from __future__ import annotations

import asyncio
import math
import json
import os
import re
import shutil
import subprocess
import uuid
from html import unescape
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request as UrlRequest
from urllib.request import urlretrieve
from urllib.request import urlopen

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field


APP_ROOT = Path(__file__).resolve().parents[1]
FRONTEND_DIR = APP_ROOT / "frontend"
DATA_DIR = Path(os.environ.get("DATA_DIR", APP_ROOT / "data")).resolve()
UPLOADS_DIR = DATA_DIR / "uploads"
RESULTS_DIR = DATA_DIR / "results"
REFERENCE_DIR = DATA_DIR / "reference"
REFERENCE_FASTA = REFERENCE_DIR / "transcriptome.fa.gz"
SALMON_INDEX_DIR = REFERENCE_DIR / "salmon_index"
REFERENCE_READY_FLAG = REFERENCE_DIR / "ready.flag"
ENSEMBL_TRANSCRIPTOME_URL = (
    "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/"
    "Homo_sapiens.GRCh38.cdna.all.fa.gz"
)

MAX_FILES_PER_RUN = 50
MAX_UPLOAD_GB = int(os.environ.get("MAX_UPLOAD_GB", "10"))
MAX_UPLOAD_BYTES = MAX_UPLOAD_GB * 1024 * 1024 * 1024
FASTQC_WORKERS = max(1, min(int(os.environ.get("FASTQC_WORKERS", "4")), 4))
ALLOWED_SUFFIXES = (".fastq", ".fastq.gz")
CHUNK_SIZE = 1024 * 1024
OPENROUTER_API_URL = "https://openrouter.ai/api/v1/chat/completions"
OPENROUTER_MODEL = os.environ.get("OPENROUTER_MODEL", "deepseek/deepseek-chat-v3-0324:free")
OPENROUTER_SITE_URL = os.environ.get("OPENROUTER_SITE_URL", "http://localhost:8000")
OPENROUTER_APP_NAME = os.environ.get("OPENROUTER_APP_NAME", "FASTQ QC App")
SUMMARY_CONTEXT_LIMIT = 60000

UPLOADS_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

fastqc_executor = ThreadPoolExecutor(max_workers=FASTQC_WORKERS)
app = FastAPI(title="FASTQ QC App", version="1.0.0")
app.mount("/results", StaticFiles(directory=RESULTS_DIR), name="results")


@dataclass
class JobState:
    job_id: str
    created_at: str
    upload_dir: Path
    result_dir: Path
    files: list[Path]
    original_names: list[str]
    status: str = "running"
    events: list[dict[str, Any]] = field(default_factory=list)
    subscribers: set[asyncio.Queue] = field(default_factory=set)
    errors: list[str] = field(default_factory=list)


@dataclass
class StreamState:
    status: str = "idle"
    events: list[dict[str, Any]] = field(default_factory=list)
    subscribers: set[asyncio.Queue] = field(default_factory=set)
    errors: list[str] = field(default_factory=list)
    started: bool = False


@dataclass
class PipelineJob:
    job_id: str
    created_at: str
    upload_dir: Path
    result_dir: Path
    files: list[Path]
    original_names: list[str]
    condition_labels: list[str]
    fastp_trim: bool
    status: str = "created"
    state: StreamState = field(default_factory=StreamState)


class ChatMessage(BaseModel):
    role: str = Field(pattern="^(user|assistant)$")
    content: str = Field(min_length=1, max_length=4000)


class ChatRequest(BaseModel):
    message: str = Field(min_length=1, max_length=4000)
    history: list[ChatMessage] = Field(default_factory=list, max_length=12)


jobs: dict[str, JobState] = {}
pipeline_jobs: dict[str, PipelineJob] = {}
reference_setup_state = StreamState()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def allowed_fastq_name(filename: str) -> bool:
    lower_name = filename.lower()
    return lower_name.endswith(ALLOWED_SUFFIXES)


def safe_filename(filename: str) -> str:
    name = Path(filename).name.strip()
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    name = name.strip("._")
    return name or f"sample_{uuid.uuid4().hex[:8]}.fastq"


def split_fastq_suffix(filename: str) -> tuple[str, str]:
    lower_name = filename.lower()
    if lower_name.endswith(".fastq.gz"):
        return filename[:-9], filename[-9:]
    return Path(filename).stem, Path(filename).suffix


def unique_destination(directory: Path, filename: str) -> Path:
    stem, suffix = split_fastq_suffix(filename)
    candidate = directory / f"{stem}{suffix}"
    index = 2
    while candidate.exists():
        candidate = directory / f"{stem}_{index}{suffix}"
        index += 1
    return candidate


def manifest_path(result_dir: Path) -> Path:
    return result_dir / "job.json"


def summary_path(result_dir: Path) -> Path:
    return result_dir / "ai_summary.json"


def report_path(job_id: str) -> Path:
    return RESULTS_DIR / job_id / "multiqc_report.html"


def report_url(job_id: str) -> str:
    return f"/results/{job_id}/multiqc_report.html"


def write_manifest(job: JobState) -> None:
    manifest = {
        "job_id": job.job_id,
        "created_at": job.created_at,
        "file_count": len(job.files),
        "files": [path.name for path in job.files],
        "original_names": job.original_names,
        "status": job.status,
        "errors": job.errors,
        "report_url": report_url(job.job_id) if report_path(job.job_id).exists() else None,
        "updated_at": utc_now(),
    }
    manifest_path(job.result_dir).write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def cleanup_uploads(job: JobState) -> None:
    shutil.rmtree(job.upload_dir, ignore_errors=True)


async def emit(job: JobState, payload: dict[str, Any]) -> None:
    event = {
        "seq": len(job.events) + 1,
        "at": utc_now(),
        **payload,
    }
    job.events.append(event)
    for subscriber in list(job.subscribers):
        await subscriber.put(event)


async def emit_stream(state: StreamState, payload: dict[str, Any]) -> None:
    event = {
        "seq": len(state.events) + 1,
        "at": utc_now(),
        **payload,
    }
    state.events.append(event)
    if payload.get("type") == "error":
        state.status = "failed"
        state.errors.append(str(payload.get("message", "Unknown error")))
    elif payload.get("type") == "done":
        state.status = "completed"
    elif payload.get("type") == "progress":
        state.status = "running"
    for subscriber in list(state.subscribers):
        await subscriber.put(event)


def encode_sse(event: dict[str, Any]) -> str:
    return f"id: {event['seq']}\ndata: {json.dumps(event)}\n\n"


def process_error_message(stderr: str, stdout: str) -> str:
    combined = "\n".join(part.strip() for part in (stderr, stdout) if part and part.strip())
    if not combined:
        return "The command exited without details."
    lines = [line.strip() for line in combined.splitlines() if line.strip()]
    return lines[-1][:500]


def run_fastqc_process(file_path: Path, result_dir: Path) -> dict[str, Any]:
    completed = subprocess.run(
        ["fastqc", "--extract", "--outdir", str(result_dir), str(file_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


async def run_fastqc(file_path: Path, result_dir: Path) -> dict[str, Any]:
    loop = asyncio.get_running_loop()
    return await loop.run_in_executor(fastqc_executor, run_fastqc_process, file_path, result_dir)


async def run_multiqc(result_dir: Path) -> dict[str, Any]:
    process = await asyncio.create_subprocess_exec(
        "multiqc",
        str(result_dir),
        "--outdir",
        str(result_dir),
        "--filename",
        "multiqc_report.html",
        "--force",
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout, stderr = await process.communicate()
    return {
        "returncode": process.returncode,
        "stdout": stdout.decode(errors="replace"),
        "stderr": stderr.decode(errors="replace"),
    }


async def run_job(job: JobState) -> None:
    try:
        await emit(
            job,
            {
                "type": "progress",
                "message": f"Starting FastQC for {len(job.files)} file(s).",
            },
        )
        completed_count = 0
        successful_count = 0

        async def run_one(file_path: Path) -> None:
            nonlocal completed_count, successful_count
            try:
                result = await run_fastqc(file_path, job.result_dir)
            except FileNotFoundError:
                result = {
                    "returncode": 127,
                    "stdout": "",
                    "stderr": "FastQC executable was not found. Run the app with Docker or install FastQC locally.",
                }
            except Exception as exc:
                result = {
                    "returncode": 1,
                    "stdout": "",
                    "stderr": str(exc),
                }
            completed_count += 1
            if result["returncode"] == 0:
                successful_count += 1
                await emit(
                    job,
                    {
                        "type": "progress",
                        "message": f"FastQC done: {file_path.name} ({completed_count}/{len(job.files)})",
                    },
                )
                return

            message = f"FastQC failed on {file_path.name}: {process_error_message(result['stderr'], result['stdout'])}"
            job.errors.append(message)
            await emit(
                job,
                {
                    "type": "error",
                    "fatal": False,
                    "message": message,
                },
            )

        await asyncio.gather(*(run_one(file_path) for file_path in job.files))

        if successful_count == 0:
            job.status = "failed"
            job.errors.append("No FastQC reports were generated, so MultiQC could not run.")
            write_manifest(job)
            await emit(
                job,
                {
                    "type": "error",
                    "fatal": True,
                    "message": "No FastQC reports were generated, so MultiQC could not run.",
                },
            )
            return

        await emit(job, {"type": "progress", "message": "Running MultiQC aggregation."})
        try:
            multiqc_result = await run_multiqc(job.result_dir)
        except FileNotFoundError:
            multiqc_result = {
                "returncode": 127,
                "stdout": "",
                "stderr": "MultiQC executable was not found. Run the app with Docker or install MultiQC locally.",
            }
        if multiqc_result["returncode"] != 0:
            job.status = "failed"
            message = f"MultiQC failed: {process_error_message(multiqc_result['stderr'], multiqc_result['stdout'])}"
            job.errors.append(message)
            write_manifest(job)
            await emit(job, {"type": "error", "fatal": True, "message": message})
            return

        if not report_path(job.job_id).exists():
            job.status = "failed"
            message = "MultiQC finished, but multiqc_report.html was not created."
            job.errors.append(message)
            write_manifest(job)
            await emit(job, {"type": "error", "fatal": True, "message": message})
            return

        job.status = "completed_with_warnings" if job.errors else "completed"
        write_manifest(job)
        await emit(
            job,
            {
                "type": "done",
                "message": "QC complete. MultiQC report is ready, and uploaded FASTQ files were cleaned up.",
                "report_url": report_url(job.job_id),
                "status": job.status,
                "errors": job.errors,
            },
        )
    except Exception as exc:  # pragma: no cover - last-resort guard for a local tool
        job.status = "failed"
        message = f"Unexpected job failure: {exc}"
        job.errors.append(message)
        write_manifest(job)
        await emit(job, {"type": "error", "fatal": True, "message": message})
    finally:
        cleanup_uploads(job)


async def save_uploads(files: list[UploadFile], upload_dir: Path) -> list[Path]:
    saved_paths: list[Path] = []
    for uploaded_file in files:
        original_name = uploaded_file.filename or ""
        safe_name = safe_filename(original_name)
        destination = unique_destination(upload_dir, safe_name)
        bytes_written = 0

        try:
            with destination.open("wb") as output_file:
                while chunk := await uploaded_file.read(CHUNK_SIZE):
                    bytes_written += len(chunk)
                    if bytes_written > MAX_UPLOAD_BYTES:
                        raise HTTPException(
                            status_code=413,
                            detail=f"{original_name} is larger than the {MAX_UPLOAD_GB} GB per-file limit.",
                        )
                    output_file.write(chunk)
        finally:
            await uploaded_file.close()

        if bytes_written == 0:
            raise HTTPException(status_code=400, detail=f"{original_name} is empty.")
        saved_paths.append(destination)

    return saved_paths


def read_manifest(result_dir: Path) -> dict[str, Any] | None:
    path = manifest_path(result_dir)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def read_json_file(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def report_text_for_summary(job_id: str) -> str:
    result_dir = RESULTS_DIR / job_id
    report = report_path(job_id)
    if not report.exists():
        raise HTTPException(status_code=404, detail="MultiQC report is not available for this job.")

    parts: list[str] = []
    manifest = read_manifest(result_dir)
    if manifest:
        parts.append("Job metadata:\n" + json.dumps(manifest, indent=2))

    for summary_file in sorted(result_dir.glob("**/summary.txt")):
        try:
            parts.append(f"FastQC summary from {summary_file.relative_to(result_dir)}:\n{summary_file.read_text(encoding='utf-8', errors='replace')}")
        except OSError:
            continue

    for data_file in sorted(result_dir.glob("**/fastqc_data.txt")):
        try:
            text = data_file.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        selected_lines = []
        keep_section = False
        for line in text.splitlines():
            if line.startswith(">>Basic Statistics"):
                keep_section = True
            elif line.startswith(">>END_MODULE"):
                if keep_section:
                    selected_lines.append(line)
                keep_section = False
            elif line.startswith(">>") and not line.startswith(">>Basic Statistics"):
                keep_section = False
            if keep_section:
                selected_lines.append(line)
        if selected_lines:
            parts.append(f"FastQC basic statistics from {data_file.relative_to(result_dir)}:\n" + "\n".join(selected_lines))

    html = report.read_text(encoding="utf-8", errors="replace")
    html = re.sub(r"(?is)<script.*?</script>", " ", html)
    html = re.sub(r"(?is)<style.*?</style>", " ", html)
    html = re.sub(r"(?s)<[^>]+>", " ", html)
    html = unescape(html)
    html = re.sub(r"\s+", " ", html).strip()
    parts.append("MultiQC report text:\n" + html[:SUMMARY_CONTEXT_LIMIT])

    return "\n\n".join(parts)[:SUMMARY_CONTEXT_LIMIT]


def request_openrouter_summary(job_id: str, report_text: str) -> str:
    api_key = os.environ.get("OPENROUTER_API_KEY", "").strip()
    if not api_key:
        raise HTTPException(
            status_code=400,
            detail="Set OPENROUTER_API_KEY in your environment, then restart Docker Compose to enable AI summaries.",
        )

    body = {
        "model": OPENROUTER_MODEL,
        "temperature": 0.2,
        "max_tokens": 800,
        "messages": [
            {
                "role": "system",
                "content": (
                    "You summarize MultiQC/FastQC reports for sequencing QC. "
                    "Be concise, practical, and do not invent metrics that are not present."
                ),
            },
            {
                "role": "user",
                "content": (
                    f"Summarize this MultiQC report for job {job_id}. "
                    "Return four short sections: Overall QC, Main warnings/failures, "
                    "Likely causes or interpretation, Recommended next checks.\n\n"
                    f"{report_text}"
                ),
            },
        ],
    }
    request = UrlRequest(
        OPENROUTER_API_URL,
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
            "HTTP-Referer": OPENROUTER_SITE_URL,
            "X-Title": OPENROUTER_APP_NAME,
        },
        method="POST",
    )

    try:
        with urlopen(request, timeout=90) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise HTTPException(status_code=502, detail=f"OpenRouter request failed: {detail[:500]}") from exc
    except (URLError, TimeoutError) as exc:
        raise HTTPException(status_code=502, detail=f"Could not reach OpenRouter: {exc}") from exc

    try:
        content = payload["choices"][0]["message"]["content"].strip()
    except (KeyError, IndexError, TypeError, AttributeError) as exc:
        raise HTTPException(status_code=502, detail="OpenRouter returned an unexpected response.") from exc

    if not content:
        raise HTTPException(status_code=502, detail="OpenRouter returned an empty summary.")
    return content


def request_openrouter_chat(job_id: str, report_text: str, chat_request: ChatRequest) -> str:
    api_key = os.environ.get("OPENROUTER_API_KEY", "").strip()
    if not api_key:
        raise HTTPException(
            status_code=400,
            detail="Set OPENROUTER_API_KEY in your environment, then restart Docker Compose to enable AI chat.",
        )

    recent_history = [
        {"role": message.role, "content": message.content}
        for message in chat_request.history[-8:]
        if message.role in {"user", "assistant"}
    ]
    messages = [
        {
            "role": "system",
            "content": (
                "You are a practical bioinformatics QC assistant. Answer questions about the provided "
                "MultiQC/FastQC report, suggest next checks, and be clear when the report does not "
                "contain enough evidence. Do not invent exact metrics that are not present."
            ),
        },
        {
            "role": "user",
            "content": (
                f"Use this report context for job {job_id} when answering the user's questions:\n\n"
                f"{report_text}"
            ),
        },
        {
            "role": "assistant",
            "content": "I have the report context. Ask me what you want to understand or investigate.",
        },
        *recent_history,
        {"role": "user", "content": chat_request.message},
    ]
    body = {
        "model": OPENROUTER_MODEL,
        "temperature": 0.3,
        "max_tokens": 900,
        "messages": messages,
    }
    request = UrlRequest(
        OPENROUTER_API_URL,
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
            "HTTP-Referer": OPENROUTER_SITE_URL,
            "X-Title": OPENROUTER_APP_NAME,
        },
        method="POST",
    )

    try:
        with urlopen(request, timeout=90) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise HTTPException(status_code=502, detail=f"OpenRouter request failed: {detail[:500]}") from exc
    except (URLError, TimeoutError) as exc:
        raise HTTPException(status_code=502, detail=f"Could not reach OpenRouter: {exc}") from exc

    try:
        content = payload["choices"][0]["message"]["content"].strip()
    except (KeyError, IndexError, TypeError, AttributeError) as exc:
        raise HTTPException(status_code=502, detail="OpenRouter returned an unexpected response.") from exc

    if not content:
        raise HTTPException(status_code=502, detail="OpenRouter returned an empty answer.")
    return content


def request_openrouter_interpretation(prompt: str) -> str:
    api_key = os.environ.get("OPENROUTER_API_KEY", "").strip()
    if not api_key:
        return "OPENROUTER_API_KEY is not configured, so AI biological interpretation was skipped."

    body = {
        "model": "mistralai/mistral-7b-instruct",
        "temperature": 0.25,
        "max_tokens": 1400,
        "messages": [
            {
                "role": "system",
                "content": "You are an expert RNA-seq bioinformatician interpreting differential expression results.",
            },
            {"role": "user", "content": prompt},
        ],
    }
    request = UrlRequest(
        OPENROUTER_API_URL,
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
            "HTTP-Referer": OPENROUTER_SITE_URL,
            "X-Title": OPENROUTER_APP_NAME,
        },
        method="POST",
    )

    try:
        with urlopen(request, timeout=120) as response:
            payload = json.loads(response.read().decode("utf-8"))
        return payload["choices"][0]["message"]["content"].strip()
    except Exception as exc:
        return f"AI interpretation failed: {exc}"


async def run_streaming_command(args: list[str], cwd: Path | None = None) -> dict[str, Any]:
    process = await asyncio.create_subprocess_exec(
        *args,
        cwd=str(cwd) if cwd else None,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout, stderr = await process.communicate()
    return {
        "returncode": process.returncode,
        "stdout": stdout.decode(errors="replace"),
        "stderr": stderr.decode(errors="replace"),
    }


async def run_reference_setup() -> None:
    state = reference_setup_state
    try:
        if REFERENCE_READY_FLAG.exists():
            await emit_stream(state, {"type": "error", "fatal": True, "message": "Reference index is already built."})
            return

        await emit_stream(state, {"type": "progress", "message": "Downloading Ensembl Homo sapiens transcriptome."})
        loop = asyncio.get_running_loop()
        await loop.run_in_executor(None, urlretrieve, ENSEMBL_TRANSCRIPTOME_URL, REFERENCE_FASTA)
        await emit_stream(state, {"type": "progress", "message": "Download complete. Building Salmon index."})

        result = await run_streaming_command(
            [
                "salmon",
                "index",
                "-t",
                str(REFERENCE_FASTA),
                "-i",
                str(SALMON_INDEX_DIR),
                "--gencode",
                "-p",
                "4",
            ]
        )
        if result["returncode"] != 0:
            raise RuntimeError(process_error_message(result["stderr"], result["stdout"]))

        REFERENCE_READY_FLAG.write_text(utc_now(), encoding="utf-8")
        await emit_stream(state, {"type": "done", "message": "Reference setup complete."})
    except Exception as exc:
        await emit_stream(state, {"type": "error", "fatal": True, "message": f"Reference setup failed: {exc}"})


def pipeline_manifest(job: PipelineJob) -> dict[str, Any]:
    return {
        "job_id": job.job_id,
        "created_at": job.created_at,
        "file_count": len(job.files),
        "files": [path.name for path in job.files],
        "original_names": job.original_names,
        "condition_labels": job.condition_labels,
        "fastp_trim": job.fastp_trim,
        "status": job.status,
        "updated_at": utc_now(),
        "pipeline": "rnaseq",
    }


def write_pipeline_manifest(job: PipelineJob) -> None:
    (job.result_dir / "pipeline_job.json").write_text(json.dumps(pipeline_manifest(job), indent=2), encoding="utf-8")


def sample_name_for(path: Path, used: set[str]) -> str:
    name = safe_filename(path.name)
    stem, _ = split_fastq_suffix(name)
    sample = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("._") or f"sample_{len(used) + 1}"
    base = sample
    index = 2
    while sample in used:
        sample = f"{base}_{index}"
        index += 1
    used.add(sample)
    return sample


def build_count_matrix(job: PipelineJob, sample_names: list[str]) -> Any:
    import pandas as pd

    series_by_sample = {}
    for sample in sample_names:
        quant_path = job.result_dir / f"salmon_{sample}" / "quant.sf"
        frame = pd.read_csv(quant_path, sep="\t")
        series_by_sample[sample] = frame.set_index("Name")["NumReads"]
    counts = pd.DataFrame(series_by_sample).fillna(0).round().astype(int)
    counts.to_csv(job.result_dir / "counts_matrix.csv")
    return counts


def run_deseq2(job: PipelineJob, counts: Any, sample_names: list[str]) -> tuple[Any, Any, Any]:
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts = counts.loc[counts.sum(axis=1) >= 10]
    metadata = pd.DataFrame({"condition": job.condition_labels}, index=sample_names)
    counts_for_deseq = counts.T
    try:
        dds = DeseqDataSet(counts=counts_for_deseq, metadata=metadata, design_factors="condition")
    except TypeError:
        dds = DeseqDataSet(counts=counts_for_deseq, metadata=metadata, design="~condition")
    dds.deseq2()
    stats = DeseqStats(dds)
    stats.summary()
    results = stats.results_df.reset_index().rename(columns={"index": "gene_id"})
    if "gene_id" not in results.columns:
        results = results.rename(columns={results.columns[0]: "gene_id"})
    keep_columns = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    for column in keep_columns:
        if column not in results.columns:
            results[column] = None
    results[keep_columns].to_csv(job.result_dir / "deseq2_results.csv", index=False)
    return dds, results[keep_columns], metadata


def significant_mask(results: Any) -> Any:
    return (results["padj"] < 0.05) & (results["log2FoldChange"].abs() > 1)


def generate_pipeline_plots(job: PipelineJob, dds: Any, results: Any, metadata: Any) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns

    plots_dir = job.result_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    plot_results = results.copy()
    plot_results["padj"] = pd.to_numeric(plot_results["padj"], errors="coerce")
    plot_results["log2FoldChange"] = pd.to_numeric(plot_results["log2FoldChange"], errors="coerce")
    plot_results["baseMean"] = pd.to_numeric(plot_results["baseMean"], errors="coerce")
    safe_padj = plot_results["padj"].clip(lower=1e-300)

    sig = significant_mask(plot_results).fillna(False)
    plt.figure(figsize=(9, 6))
    plt.scatter(plot_results["log2FoldChange"], -np.log10(safe_padj), c=np.where(sig, "#dc2626", "#94a3b8"), s=16, alpha=0.85)
    for _, row in plot_results[sig].nsmallest(10, "padj").iterrows():
        plt.text(row["log2FoldChange"], -math.log10(max(row["padj"], 1e-300)), str(row["gene_id"])[:18], fontsize=7)
    plt.axvline(-1, color="#334155", linestyle="--", linewidth=1)
    plt.axvline(1, color="#334155", linestyle="--", linewidth=1)
    plt.axhline(1.3, color="#334155", linestyle="--", linewidth=1)
    plt.xlabel("log2FoldChange")
    plt.ylabel("-log10(padj)")
    plt.title("Volcano Plot - Differential Expression")
    plt.tight_layout()
    plt.savefig(plots_dir / "volcano.png", dpi=160)
    plt.close()

    plt.figure(figsize=(9, 6))
    plt.scatter(np.log10(plot_results["baseMean"].fillna(0) + 1), plot_results["log2FoldChange"], c=np.where(plot_results["padj"] < 0.05, "#dc2626", "#94a3b8"), s=16, alpha=0.85)
    plt.axhline(0, color="#334155", linestyle="--", linewidth=1)
    plt.xlabel("log10(baseMean + 1)")
    plt.ylabel("log2FoldChange")
    plt.title("MA Plot")
    plt.tight_layout()
    plt.savefig(plots_dir / "ma.png", dpi=160)
    plt.close()

    normed = dds.layers["normed_counts"]
    normed_df = pd.DataFrame(normed, index=metadata.index, columns=dds.var_names)
    centered = np.log2(normed_df + 1)
    centered = centered - centered.mean(axis=0)
    u_matrix, singular_values, _ = np.linalg.svd(centered.to_numpy(), full_matrices=False)
    coords = u_matrix[:, :2] * singular_values[:2]
    variance = (singular_values**2) / np.sum(singular_values**2) * 100
    plt.figure(figsize=(8, 6))
    for condition in metadata["condition"].unique():
        mask = metadata["condition"] == condition
        plt.scatter(coords[mask, 0], coords[mask, 1], label=condition, s=70)
    plt.xlabel(f"PC1 ({variance[0]:.1f}% variance)")
    plt.ylabel(f"PC2 ({variance[1]:.1f}% variance)" if len(variance) > 1 else "PC2")
    plt.title("PCA Plot - Sample Clustering")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots_dir / "pca.png", dpi=160)
    plt.close()

    top_genes = plot_results.dropna(subset=["padj"]).nsmallest(50, "padj")["gene_id"].tolist()
    if top_genes:
        heatmap_data = normed_df[top_genes].T
        sns.clustermap(np.log2(heatmap_data + 1), z_score=0, cmap="vlag", figsize=(10, 10), col_colors=metadata["condition"].astype("category").cat.codes.map(lambda _: "#2563eb"))
        plt.suptitle("Heatmap - Top 50 DE Genes", y=1.02)
        plt.savefig(plots_dir / "heatmap.png", dpi=160, bbox_inches="tight")
        plt.close("all")


def build_interpretation_prompt(results: Any) -> str:
    clean = results.copy()
    clean["padj"] = clean["padj"].fillna(1)
    clean["log2FoldChange"] = clean["log2FoldChange"].fillna(0)
    up = clean[(clean["padj"] < 0.05) & (clean["log2FoldChange"] > 1)].sort_values("padj")
    down = clean[(clean["padj"] < 0.05) & (clean["log2FoldChange"] < -1)].sort_values("padj")
    return (
        "You are an expert RNA-seq bioinformatician. "
        "Analyze these differential expression results:\n"
        f"- Total genes tested: {len(clean)}\n"
        f"- Significantly upregulated (padj<0.05, log2FC>1): {len(up)}\n"
        f"- Significantly downregulated (padj<0.05, log2FC<-1): {len(down)}\n"
        f"- Top 10 upregulated genes: {', '.join(up['gene_id'].head(10).astype(str).tolist()) or 'none'}\n"
        f"- Top 10 downregulated genes: {', '.join(down['gene_id'].head(10).astype(str).tolist()) or 'none'}\n"
        "Provide a detailed biological interpretation. Discuss what these expression changes might mean biologically, "
        "what pathways could be involved, and what the data quality suggests."
    )


async def run_pipeline_job(job: PipelineJob) -> None:
    state = job.state
    job.status = "running"
    write_pipeline_manifest(job)
    used_names: set[str] = set()
    sample_names = [sample_name_for(path, used_names) for path in job.files]
    try:
        await emit_stream(state, {"type": "progress", "message": "Starting RNA-seq pipeline."})

        await emit_stream(state, {"type": "progress", "message": "STEP 1 - FastQC on raw files."})
        for file_path in job.files:
            result = await run_fastqc(file_path, job.result_dir)
            if result["returncode"] != 0:
                raise RuntimeError(f"FastQC failed on {file_path.name}: {process_error_message(result['stderr'], result['stdout'])}")
        await emit_stream(state, {"type": "progress", "message": "[OK] FastQC completed"})

        quant_inputs = job.files
        if job.fastp_trim:
            await emit_stream(state, {"type": "progress", "message": "STEP 2 - fastp trimming."})
            trimmed_dir = job.result_dir / "trimmed"
            trimmed_dir.mkdir(parents=True, exist_ok=True)
            trimmed_inputs: list[Path] = []
            for file_path, sample in zip(job.files, sample_names):
                trimmed_output = trimmed_dir / f"{sample}.trimmed.fastq.gz"
                result = await run_streaming_command(
                    [
                        "fastp",
                        "-i",
                        str(file_path),
                        "-o",
                        str(trimmed_output),
                        "--json",
                        str(job.result_dir / f"fastp_{sample}.json"),
                        "--html",
                        str(job.result_dir / f"fastp_{sample}.html"),
                        "-w",
                        "4",
                    ]
                )
                if result["returncode"] != 0:
                    raise RuntimeError(f"fastp failed on {file_path.name}: {process_error_message(result['stderr'], result['stdout'])}")
                trimmed_inputs.append(trimmed_output)
            quant_inputs = trimmed_inputs
            await emit_stream(state, {"type": "progress", "message": "[OK] fastp trimming completed"})

        await emit_stream(state, {"type": "progress", "message": "STEP 3 - Salmon quantification."})
        if not REFERENCE_READY_FLAG.exists():
            raise RuntimeError("Reference index is not ready. Run Setup Reference first.")
        for file_path, sample in zip(quant_inputs, sample_names):
            result = await run_streaming_command(
                [
                    "salmon",
                    "quant",
                    "-i",
                    str(SALMON_INDEX_DIR),
                    "-l",
                    "A",
                    "-r",
                    str(file_path),
                    "-p",
                    "4",
                    "--validateMappings",
                    "-o",
                    str(job.result_dir / f"salmon_{sample}"),
                ]
            )
            if result["returncode"] != 0:
                raise RuntimeError(f"Salmon failed on {file_path.name}: {process_error_message(result['stderr'], result['stdout'])}")
        await emit_stream(state, {"type": "progress", "message": "[OK] Salmon quantification completed"})

        await emit_stream(state, {"type": "progress", "message": "STEP 4 - Building count matrix."})
        counts = build_count_matrix(job, sample_names)
        await emit_stream(state, {"type": "progress", "message": "[OK] Count matrix built"})

        await emit_stream(state, {"type": "progress", "message": "STEP 5 - Running PyDESeq2."})
        dds, results, metadata = run_deseq2(job, counts, sample_names)
        await emit_stream(state, {"type": "progress", "message": "[OK] DESeq2 differential expression completed"})

        await emit_stream(state, {"type": "progress", "message": "STEP 6 - Generating plots."})
        generate_pipeline_plots(job, dds, results, metadata)
        await emit_stream(state, {"type": "progress", "message": "[OK] Plots generated"})

        await emit_stream(state, {"type": "progress", "message": "STEP 7 - Running MultiQC."})
        multiqc_result = await run_multiqc(job.result_dir)
        if multiqc_result["returncode"] != 0:
            raise RuntimeError(f"MultiQC failed: {process_error_message(multiqc_result['stderr'], multiqc_result['stdout'])}")
        await emit_stream(state, {"type": "progress", "message": "[OK] MultiQC report completed"})

        await emit_stream(state, {"type": "progress", "message": "STEP 8 - Generating AI biological interpretation."})
        interpretation = request_openrouter_interpretation(build_interpretation_prompt(results))
        (job.result_dir / "interpretation.txt").write_text(interpretation, encoding="utf-8")
        await emit_stream(state, {"type": "progress", "message": "[OK] AI interpretation generated"})

        (job.result_dir / "pipeline_complete.flag").write_text(utc_now(), encoding="utf-8")
        job.status = "completed"
        write_pipeline_manifest(job)
        await emit_stream(state, {"type": "done", "message": "RNA-seq pipeline complete.", "job_id": job.job_id})
    except Exception as exc:
        job.status = "failed"
        write_pipeline_manifest(job)
        await emit_stream(state, {"type": "error", "fatal": True, "message": str(exc), "job_id": job.job_id})
    finally:
        shutil.rmtree(job.upload_dir, ignore_errors=True)


@app.get("/", include_in_schema=False)
async def index() -> FileResponse:
    return FileResponse(FRONTEND_DIR / "index.html")


@app.get("/health")
async def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/api/reference-status")
async def reference_status() -> dict[str, Any]:
    return {
        "ready": REFERENCE_READY_FLAG.exists(),
        "status": reference_setup_state.status,
        "started": reference_setup_state.started,
        "errors": reference_setup_state.errors,
    }


@app.get("/api/setup-reference")
async def setup_reference() -> dict[str, Any]:
    if REFERENCE_READY_FLAG.exists():
        raise HTTPException(status_code=409, detail="Reference index is already built.")
    if reference_setup_state.started and reference_setup_state.status == "running":
        return {"started": True, "status": "running"}

    reference_setup_state.started = True
    reference_setup_state.status = "running"
    reference_setup_state.events.clear()
    reference_setup_state.errors.clear()
    asyncio.create_task(run_reference_setup())
    return {"started": True, "status": "running"}


@app.get("/api/setup-reference/stream")
async def stream_reference_setup(request: Request) -> StreamingResponse:
    async def event_generator():
        yield "retry: 2000\n\n"
        last_event_id = request.headers.get("last-event-id")
        try:
            last_seen = int(last_event_id) if last_event_id else 0
        except ValueError:
            last_seen = 0

        for event in reference_setup_state.events:
            if event["seq"] > last_seen:
                yield encode_sse(event)

        if reference_setup_state.events and reference_setup_state.events[-1]["type"] in {"done", "error"}:
            return

        queue: asyncio.Queue = asyncio.Queue()
        reference_setup_state.subscribers.add(queue)
        try:
            while True:
                if await request.is_disconnected():
                    break
                event = await queue.get()
                yield encode_sse(event)
                if event["type"] in {"done", "error"}:
                    break
        finally:
            reference_setup_state.subscribers.discard(queue)

    return StreamingResponse(event_generator(), media_type="text/event-stream")


@app.post("/api/run")
async def create_run(files: list[UploadFile] = File(...)) -> JSONResponse:
    if not files:
        raise HTTPException(status_code=400, detail="Choose at least one FASTQ or FASTQ.GZ file.")
    if len(files) > MAX_FILES_PER_RUN:
        raise HTTPException(status_code=400, detail=f"Upload at most {MAX_FILES_PER_RUN} files per run.")

    invalid_names = [file.filename or "(unnamed)" for file in files if not allowed_fastq_name(file.filename or "")]
    if invalid_names:
        invalid_list = ", ".join(invalid_names[:5])
        raise HTTPException(status_code=400, detail=f"Only .fastq and .fastq.gz files are supported: {invalid_list}")

    job_id = uuid.uuid4().hex
    created_at = utc_now()
    upload_dir = UPLOADS_DIR / job_id
    result_dir = RESULTS_DIR / job_id
    upload_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    try:
        saved_paths = await save_uploads(files, upload_dir)
    except HTTPException:
        shutil.rmtree(upload_dir, ignore_errors=True)
        shutil.rmtree(result_dir, ignore_errors=True)
        raise

    job = JobState(
        job_id=job_id,
        created_at=created_at,
        upload_dir=upload_dir,
        result_dir=result_dir,
        files=saved_paths,
        original_names=[file.filename or path.name for file, path in zip(files, saved_paths)],
    )
    jobs[job_id] = job
    write_manifest(job)
    asyncio.create_task(run_job(job))
    return JSONResponse({"job_id": job_id})


@app.post("/api/run-pipeline")
async def create_pipeline_run(
    files: list[UploadFile] = File(...),
    condition_labels: str = Form(...),
    fastp_trim: bool = Form(True),
) -> JSONResponse:
    if not files:
        raise HTTPException(status_code=400, detail="Choose at least one FASTQ or FASTQ.GZ file.")
    if len(files) > MAX_FILES_PER_RUN:
        raise HTTPException(status_code=400, detail=f"Upload at most {MAX_FILES_PER_RUN} files per run.")

    invalid_names = [file.filename or "(unnamed)" for file in files if not allowed_fastq_name(file.filename or "")]
    if invalid_names:
        invalid_list = ", ".join(invalid_names[:5])
        raise HTTPException(status_code=400, detail=f"Only .fastq and .fastq.gz files are supported: {invalid_list}")

    try:
        labels = json.loads(condition_labels)
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=400, detail="condition_labels must be a JSON array.") from exc
    if not isinstance(labels, list) or len(labels) != len(files) or not all(isinstance(label, str) and label.strip() for label in labels):
        raise HTTPException(status_code=400, detail="condition_labels must contain one non-empty label per file.")
    if len(set(labels)) < 2:
        raise HTTPException(status_code=400, detail="Provide at least two condition labels for differential expression.")

    job_id = uuid.uuid4().hex
    created_at = utc_now()
    upload_dir = UPLOADS_DIR / job_id
    result_dir = RESULTS_DIR / job_id
    upload_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    try:
        saved_paths = await save_uploads(files, upload_dir)
    except HTTPException:
        shutil.rmtree(upload_dir, ignore_errors=True)
        shutil.rmtree(result_dir, ignore_errors=True)
        raise

    job = PipelineJob(
        job_id=job_id,
        created_at=created_at,
        upload_dir=upload_dir,
        result_dir=result_dir,
        files=saved_paths,
        original_names=[file.filename or path.name for file, path in zip(files, saved_paths)],
        condition_labels=[label.strip() for label in labels],
        fastp_trim=fastp_trim,
    )
    pipeline_jobs[job_id] = job
    write_pipeline_manifest(job)
    return JSONResponse({"job_id": job_id})


@app.get("/api/run-pipeline/{job_id}/stream")
async def stream_pipeline_run(job_id: str, request: Request) -> StreamingResponse:
    async def event_generator():
        yield "retry: 2000\n\n"
        job = pipeline_jobs.get(job_id)
        if job is None:
            complete_flag = RESULTS_DIR / job_id / "pipeline_complete.flag"
            if complete_flag.exists():
                yield encode_sse({"seq": 1, "at": utc_now(), "type": "done", "message": "RNA-seq pipeline complete.", "job_id": job_id})
            else:
                yield encode_sse({"seq": 1, "at": utc_now(), "type": "error", "fatal": True, "message": "Unknown pipeline job ID.", "job_id": job_id})
            return

        if not job.state.started:
            job.state.started = True
            asyncio.create_task(run_pipeline_job(job))

        last_event_id = request.headers.get("last-event-id")
        try:
            last_seen = int(last_event_id) if last_event_id else 0
        except ValueError:
            last_seen = 0

        for event in job.state.events:
            if event["seq"] > last_seen:
                yield encode_sse(event)

        if job.state.events and job.state.events[-1]["type"] in {"done", "error"}:
            return

        queue: asyncio.Queue = asyncio.Queue()
        job.state.subscribers.add(queue)
        try:
            while True:
                if await request.is_disconnected():
                    break
                event = await queue.get()
                yield encode_sse(event)
                if event["type"] in {"done", "error"}:
                    break
        finally:
            job.state.subscribers.discard(queue)

    return StreamingResponse(event_generator(), media_type="text/event-stream")


@app.get("/api/run/{job_id}/stream")
async def stream_run(job_id: str, request: Request) -> StreamingResponse:
    async def event_generator():
        yield "retry: 2000\n\n"
        last_event_id = request.headers.get("last-event-id")
        try:
            last_seen = int(last_event_id) if last_event_id else 0
        except ValueError:
            last_seen = 0

        job = jobs.get(job_id)
        if job is None:
            report = report_path(job_id)
            if report.exists():
                event = {
                    "seq": 1,
                    "at": utc_now(),
                    "type": "done",
                    "message": "QC complete. MultiQC report is ready.",
                    "report_url": report_url(job_id),
                    "status": "completed",
                    "errors": [],
                }
                yield encode_sse(event)
                return
            event = {
                "seq": 1,
                "at": utc_now(),
                "type": "error",
                "fatal": True,
                "message": "Unknown job ID.",
            }
            yield encode_sse(event)
            return

        for event in job.events:
            if event["seq"] > last_seen:
                yield encode_sse(event)

        if job.events and job.events[-1]["type"] in {"done", "error"} and job.events[-1].get("fatal", True):
            return

        queue: asyncio.Queue = asyncio.Queue()
        job.subscribers.add(queue)
        try:
            while True:
                if await request.is_disconnected():
                    break
                event = await queue.get()
                yield encode_sse(event)
                if event["type"] == "done" or (event["type"] == "error" and event.get("fatal", True)):
                    break
        finally:
            job.subscribers.discard(queue)

    return StreamingResponse(event_generator(), media_type="text/event-stream")


@app.get("/api/run-pipeline/{job_id}/plots")
async def pipeline_plots(job_id: str) -> dict[str, str]:
    plots_dir = RESULTS_DIR / job_id / "plots"
    if not plots_dir.exists():
        raise HTTPException(status_code=404, detail="Pipeline plots are not available for this job.")
    return {
        "volcano": f"/api/plots/{job_id}/volcano.png",
        "ma": f"/api/plots/{job_id}/ma.png",
        "pca": f"/api/plots/{job_id}/pca.png",
        "heatmap": f"/api/plots/{job_id}/heatmap.png",
    }


@app.get("/api/plots/{job_id}/{filename}")
async def serve_pipeline_plot(job_id: str, filename: str) -> FileResponse:
    allowed = {"volcano.png", "ma.png", "pca.png", "heatmap.png"}
    if filename not in allowed:
        raise HTTPException(status_code=404, detail="Unknown plot.")
    plot_path = (RESULTS_DIR / job_id / "plots" / filename).resolve()
    if RESULTS_DIR not in plot_path.parents or not plot_path.exists():
        raise HTTPException(status_code=404, detail="Plot not found.")
    return FileResponse(plot_path)


@app.get("/api/run-pipeline/{job_id}/interpretation")
async def pipeline_interpretation(job_id: str) -> dict[str, str]:
    path = RESULTS_DIR / job_id / "interpretation.txt"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Interpretation is not available for this job.")
    return {"text": path.read_text(encoding="utf-8", errors="replace")}


@app.get("/api/run-pipeline/{job_id}/deseq2")
async def pipeline_deseq2(job_id: str) -> list[dict[str, Any]]:
    import pandas as pd

    path = RESULTS_DIR / job_id / "deseq2_results.csv"
    if not path.exists():
        raise HTTPException(status_code=404, detail="DESeq2 results are not available for this job.")
    frame = pd.read_csv(path)
    frame = frame.where(pd.notna(frame), None)
    return frame.to_dict(orient="records")


@app.get("/api/jobs")
async def list_jobs() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result_dir in RESULTS_DIR.iterdir():
        if not result_dir.is_dir():
            continue

        manifest = read_manifest(result_dir) or {}
        job_id = manifest.get("job_id", result_dir.name)
        report_exists = (result_dir / "multiqc_report.html").exists()
        status = manifest.get("status", "completed" if report_exists else "unknown")
        if status == "running" and job_id not in jobs and not report_exists:
            status = "interrupted"

        created_at = manifest.get("created_at")
        if not created_at:
            created_at = datetime.fromtimestamp(result_dir.stat().st_mtime, timezone.utc).isoformat()

        rows.append(
            {
                "job_id": job_id,
                "created_at": created_at,
                "file_count": manifest.get("file_count", 0),
                "status": status,
                "report_url": report_url(job_id) if report_exists else None,
                "errors": manifest.get("errors", []),
            }
        )

    return sorted(rows, key=lambda row: row["created_at"], reverse=True)


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str) -> dict[str, Any]:
    job = jobs.get(job_id)
    if job and job.status == "running":
        raise HTTPException(status_code=409, detail="This run is still running.")

    result_dir = (RESULTS_DIR / job_id).resolve()
    upload_dir = (UPLOADS_DIR / job_id).resolve()
    if RESULTS_DIR not in result_dir.parents or UPLOADS_DIR not in upload_dir.parents:
        raise HTTPException(status_code=400, detail="Invalid job ID.")

    if not result_dir.exists() and not upload_dir.exists():
        raise HTTPException(status_code=404, detail="Unknown job ID.")

    shutil.rmtree(result_dir, ignore_errors=True)
    shutil.rmtree(upload_dir, ignore_errors=True)
    jobs.pop(job_id, None)
    pipeline_jobs.pop(job_id, None)
    return {"deleted": True, "job_id": job_id}


@app.post("/api/jobs/{job_id}/summary")
async def summarize_job(job_id: str) -> dict[str, Any]:
    result_dir = RESULTS_DIR / job_id
    if not result_dir.exists() or not result_dir.is_dir():
        raise HTTPException(status_code=404, detail="Unknown job ID.")

    cached = read_json_file(summary_path(result_dir))
    if cached and cached.get("summary"):
        return cached

    report_text = report_text_for_summary(job_id)
    summary = request_openrouter_summary(job_id, report_text)
    payload = {
        "job_id": job_id,
        "model": OPENROUTER_MODEL,
        "summary": summary,
        "created_at": utc_now(),
    }
    summary_path(result_dir).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


@app.get("/api/jobs/{job_id}/summary")
async def get_job_summary(job_id: str) -> dict[str, Any]:
    result_dir = RESULTS_DIR / job_id
    if not result_dir.exists() or not result_dir.is_dir():
        raise HTTPException(status_code=404, detail="Unknown job ID.")

    cached = read_json_file(summary_path(result_dir))
    if cached and cached.get("summary"):
        return cached
    raise HTTPException(status_code=404, detail="No AI summary has been generated for this job yet.")


@app.post("/api/jobs/{job_id}/chat")
async def chat_about_job(job_id: str, chat_request: ChatRequest) -> dict[str, Any]:
    result_dir = RESULTS_DIR / job_id
    if not result_dir.exists() or not result_dir.is_dir():
        raise HTTPException(status_code=404, detail="Unknown job ID.")

    report_text = report_text_for_summary(job_id)
    answer = request_openrouter_chat(job_id, report_text, chat_request)
    return {
        "job_id": job_id,
        "model": OPENROUTER_MODEL,
        "answer": answer,
        "created_at": utc_now(),
    }

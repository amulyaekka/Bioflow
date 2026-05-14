# BIOFLOW

A local web app that runs FastQC and MultiQC on uploaded FASTQ files, streams progress in the browser, embeds the report, and can generate a simple AI summary using OpenRouter.

## Start

```bash
docker compose up
```

Then open:

```text
http://localhost:8000
```

The first build downloads the Python image, FastQC, MultiQC, and dependencies. Later starts are much faster.

## Use

1. Choose one or more `.fastq` or `.fastq.gz` files.
2. Click **Run QC**.
3. Watch the progress log while FastQC runs.
4. Open the embedded MultiQC report when the run completes.
5. Click **Generate AI Summary** if OpenRouter is configured.
6. Click **Chat About Report** to ask questions about QC warnings, likely causes, or next checks.
7. Use **Past Runs** to reload earlier reports.

## AI Summary

AI summaries are optional. Your OpenRouter API key is used only by the FastAPI backend and is never sent to the browser.

Create a local `.env` file next to `docker-compose.yml`:

```text
OPENROUTER_API_KEY=your_openrouter_key_here
OPENROUTER_MODEL=deepseek/deepseek-chat-v3-0324:free
```

Then restart the app:

```bash
docker compose down
docker compose up --build
```

Keep `.env` private. It is already listed in `.gitignore`.

## Project Layout

```text
fastq-qc-app/
|-- docker-compose.yml
|-- backend/
|   |-- Dockerfile
|   |-- main.py
|   `-- requirements.txt
|-- frontend/
|   `-- index.html
`-- data/
    |-- uploads/
    `-- results/
```

## API

| Method | Path | Purpose |
| --- | --- | --- |
| `POST` | `/api/run` | Upload FASTQ files and start a QC job |
| `GET` | `/api/run/{job_id}/stream` | Server-Sent Events progress stream |
| `GET` | `/api/jobs` | List previous runs from `data/results/` |
| `POST` | `/api/jobs/{job_id}/summary` | Generate or return cached AI summary |
| `GET` | `/api/jobs/{job_id}/summary` | Return cached AI summary |
| `POST` | `/api/jobs/{job_id}/chat` | Ask a question about the MultiQC report |
| `GET` | `/results/{job_id}/multiqc_report.html` | MultiQC report |

## Data

Uploads and results are bind-mounted into `./data`, so they survive container restarts.

```text
data/
|-- uploads/{job_id}/
`-- results/{job_id}/
    |-- *_fastqc.html
    |-- *_fastqc.zip
    |-- job.json
    |-- ai_summary.json
    `-- multiqc_report.html
```

## Stop

Press `Ctrl + C` in the terminal, or run:

```bash
docker compose down
```

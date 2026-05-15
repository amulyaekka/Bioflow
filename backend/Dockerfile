FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends fastqc default-jre-headless fastp wget ca-certificates libgomp1 libtbb12 \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz \
    && tar xzf salmon-1.10.0_linux_x86_64.tar.gz \
    && SALMON_BIN="$(find . -type f -path '*/bin/salmon' -print -quit)" \
    && test -n "$SALMON_BIN" \
    && cp "$SALMON_BIN" /usr/local/bin/salmon \
    && chmod +x /usr/local/bin/salmon \
    && rm -rf salmon-* salmon-1.10.0_linux_x86_64.tar.gz

WORKDIR /app

COPY backend/requirements.txt /app/backend/requirements.txt
RUN pip install --no-cache-dir -r /app/backend/requirements.txt

COPY backend /app/backend
COPY frontend /app/frontend

RUN mkdir -p /app/data/uploads /app/data/results

EXPOSE 8000

CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]

# Verifier image, built and run after the agent's container is torn down.
FROM python:3.13-slim-bookworm

RUN pip install --no-cache-dir pytest==9.1.1 pytest-json-ctrf==0.5.2

COPY . /tests/

# Parent directory of every path in task.toml's top-level artifacts list.
RUN mkdir -p /app

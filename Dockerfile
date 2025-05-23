FROM python:3.12-slim

WORKDIR /app

COPY pyproject.toml .
COPY documentation documentation
COPY src/ src/

RUN apt-get update && apt-get install -y build-essential

# Install dependencies using pip and pyproject.toml (via pip’s PEP 517 support)
RUN pip install --upgrade pip setuptools wheel
RUN pip install .

ENTRYPOINT ["python", "src/dementpy.py"]

FROM python:3.11.11-slim-bookworm@sha256:081075da77b2b55c23c088251026fb69a7b2bf92471e491ff5fd75c192fd38e5
ENV MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=0
RUN apt-get update && apt-get install -y --no-install-recommends procps=2:4.0.2-3 && rm -rf /var/lib/apt/lists/*
COPY requirements.lock.txt /tmp/requirements.lock.txt
RUN pip install --no-cache-dir -r /tmp/requirements.lock.txt
RUN apt-get update && apt-get install -y --no-install-recommends libxrender1=1:0.9.10-1.1 libxext6=2:1.3.4-1+b1 libexpat1=2.5.0-1+deb12u5 && rm -rf /var/lib/apt/lists/*
WORKDIR /work

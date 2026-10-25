# SINA platform image — built and verified in CI (Gate 0 evidence).
FROM continuumio/miniconda3:24.1.2-0

WORKDIR /app

COPY environment.yml environment.lock.yml* ./
RUN conda env create -f environment.yml && conda clean -afy
ENV PATH=/opt/conda/envs/sina/bin:$PATH

COPY . .

RUN pip install -e . --no-deps

# Fail the image build if the unit suite does not pass inside the container.
RUN python -m pytest tests/unit -q && ruff check src tests scripts

CMD ["python", "-c", "import sina; print('sina container OK')"]

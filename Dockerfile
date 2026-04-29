FROM mambaorg/micromamba:1.5.8

COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

WORKDIR /app
COPY . /app

RUN /opt/conda/bin/python -m pip install -e ".[structural,dev]"

CMD ["/bin/bash"]

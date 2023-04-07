FROM python:3.11.2-slim-buster

RUN apt-get update -y \
    && apt-get install -y \
    git \
    vim \
    tree \
    zip \
    unzip \
    tar \
    jq \
    curl \
    zsh \
    libxrender1

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# Install poetry
ENV PATH="/root/.local/bin:$PATH"
RUN curl -sSL https://install.python-poetry.org | python -

SHELL ["/bin/zsh", "-c"]
RUN poetry config virtualenvs.in-project true

CMD ["/bin/zsh"]

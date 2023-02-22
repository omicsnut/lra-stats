FROM rust:1.67-slim AS builder

RUN apt-get update && apt-get --yes upgrade && apt-get install --yes --no-install-recommends git cmake make gcc \
    zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libreadline-dev

RUN git clone "https://github.com/omicsnut/lra-stats.git" /lra-stats && cd /lra-stats && cargo build --release

FROM rmusunuri/samtools:v1.17
LABEL maintainer="Rajeeva Musunuri <rmusunuri@nygenome.org>"

ARG DEBIAN_FRONTEND="noninteractive"
ENV TZ="US/Eastern"

COPY --from=builder /lra-stats/target/release/lra-stats /usr/bin/lra-stats

RUN curl -sSLO "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-418.0.0-linux-x86_64.tar.gz" && \
    tar -xzvf google-cloud-cli-418.0.0-linux-x86_64.tar.gz && rm -f "google-cloud-cli-418.0.0-linux-x86_64.tar.gz" && \
    ./google-cloud-sdk/install.sh --quiet

ENV PATH "/google-cloud-sdk/bin:$PATH"

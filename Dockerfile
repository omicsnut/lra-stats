FROM rust:1.67-slim AS builder

RUN apt-get update && apt-get --yes upgrade && apt-get install --yes --no-install-recommends git cmake make gcc \
    zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libreadline-dev

RUN git clone "https://github.com/omicsnut/lra-stats.git" /lra-stats && cd /lra-stats && cargo build --release

FROM ubuntu:22.04
LABEL maintainer="Rajeeva Musunuri <rmusunuri@nygenome.org>"

ARG DEBIAN_FRONTEND="noninteractive"
ENV TZ="US/Eastern"

COPY --from=builder /lra-stats/target/release/lra-stats /usr/bin/lra-stats

RUN apt-get update && apt-get install --yes --no-install-recommends curl time tzdata ca-certificates bzip2 gcc make perl \
    zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libreadline-dev && \
    apt-get autoremove --purge --yes && apt-get clean && rm -rf "/var/lib/apt/lists/*" "/tmp/*" "/var/tmp/*"

RUN curl -sSLO "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2" && \
    tar -xjvf samtools-1.16.1.tar.bz2 && rm -f samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && ./configure && make && make install

RUN curl -sSLO "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-418.0.0-linux-x86_64.tar.gz" && \
    tar -xzvf google-cloud-cli-418.0.0-linux-x86_64.tar.gz && rm -f "google-cloud-cli-418.0.0-linux-x86_64.tar.gz" && \
    ./google-cloud-sdk/install.sh --quiet

ENV PATH "/google-cloud-sdk/bin:$PATH"

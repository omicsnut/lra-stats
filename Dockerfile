FROM rust:1-slim AS builder

RUN apt-get update && apt-get --yes upgrade && apt-get install --yes --no-install-recommends git cmake make gcc \
    zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libreadline-dev

RUN git clone "https://github.com/omicsnut/lra-stats.git" /lra-stats && cd /lra-stats && cargo build --release

FROM us-central1-docker.pkg.dev/nygc-comp-s-1856/lrpipe-deps/samtools:v1.20-gcloud479
LABEL maintainer="Rajeeva Musunuri <rmusunuri@nygenome.org>"

ARG DEBIAN_FRONTEND="noninteractive"
ENV TZ="US/Eastern"

COPY --from=builder /lra-stats/target/release/lra-stats /usr/bin/lra-stats

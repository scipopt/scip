FROM debian:bullseye-slim

ARG DEBIAN_FRONTEND=noninteractive
ARG TARGETPLATFORM
ARG TAG

RUN apt-get update && apt-get -y install \
    build-essential \
    libcliquer1 \
    gfortran \
    liblapack3 \
    libopenblas-dev \
    libgsl25 \
    libtbb2 \
    curl

WORKDIR /tmp
RUN if [ "$TARGETPLATFORM" = "linux/amd64" ]; then curl -Lo SCIPOptSuite.deb https://github.com/scipopt/scip/releases/download/$(echo "v${TAG}" | tr -d '.')/SCIPOptSuite-${TAG}-Linux-debian.deb; fi
RUN if [ "$TARGETPLATFORM" = "linux/arm64" ]; then curl -Lo SCIPOptSuite.deb https://github.com/scipopt/scip/releases/download/$(echo "v${TAG}" | tr -d '.')/SCIPOptSuite-${TAG}-Linux-arm64.deb; fi
RUN dpkg -i SCIPOptSuite.deb && rm SCIPOptSuite.deb

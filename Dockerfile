# Use a base image with a pre-installed compiler (e.g., gcc)
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        cmake \
        libtbb-dev \
        zlib1g \
        zlib1g-dev \
        libboost-iostreams-dev \
        libblosc-dev
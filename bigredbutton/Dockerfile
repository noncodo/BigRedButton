# Set the base image
FROM ubuntu

# Dockerfile author / maintainer 
MAINTAINER Name <martinalexandersmith@gmail.com>
USER root
RUN apt-get update &&\
    apt-get install -y \
        build-essential

COPY ViennaRNA-2.4.1.tar.gz /tmp/viennarna.tar.gz
RUN mkdir -p /tmp/viennarna \
    && tar -xf /tmp/viennarna.tar.gz -C /tmp/viennarna --strip-components=1 \
    && cd /tmp/viennarna \
    && chmod +x ./configure \
    && ./configure --prefix=/usr/local \
    && make -j4 \
    && make install \
    && rm -rf /tmp/*

###### SETUP LOCARNA ###############
COPY locarna-1.9.2.tar.gz /tmp/locarna.tar.gz
RUN mkdir -p /tmp/locarna \
    && tar -xf /tmp/locarna.tar.gz -C /tmp/locarna --strip-components=1 \
    && cd /tmp/locarna \
    && chmod +x ./configure \
    && ./configure --prefix=/usr/local \
    && make -j4 \
    && make install \
    && rm -rf /tmp/*

FROM r-base

RUN mkdir -p /bigredbutton
ENV WORKDIR /bigredbutton
COPY launcher.sh worker.sge fa2pp.sge /bigredbutton/


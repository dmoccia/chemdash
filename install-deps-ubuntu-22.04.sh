#!/bin/bash
set -e
set -x

export DEBIAN_FRONTEND=noninteractive

apt-get update -y

apt-get install -y \
    git \
    libgcc1 \
    libstdc++6 \
    sqlite3 \
    zlib1g \
    bzip2 \
    libcairo2 \
    libdbus-1-3 \
    libexpat1 \
    fontconfig \
    libfreetype6 \
    gettext \
    libglib2.0-0 \
    libgmp10 \
    libhdf5-dev \
    libicu-dev \
    libjpeg-dev \
    libarchive13 \
    libcurl4 \
    libedit2 \
    libffi-dev \
    libc-bin \
    libpng-dev \
    libsodium-dev \
    libssh2-1 \
    libtiff5-dev \
    libxml2 \
    libxslt1-dev \
    libomp-dev \
    liblz4-dev \
    liblzo2-2 \
    libncurses5-dev \
    libpcre3 \
    qtbase5-dev \
    qtchooser \
    qt5-qmake \
    qtbase5-dev-tools \
    libsnappy-dev \
    tk-dev \
    xz-utils \
    libyaml-dev \
    libzmq3-dev \
    libzstd-dev \
    python3-pip \
    libblas3 liblapack3 liblapack-dev libblas-dev \
    vim


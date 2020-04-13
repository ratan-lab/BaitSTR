FROM ubuntu:16.04
MAINTAINER Aakrosh Ratan <ratan@bx.psu.edu>

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update 

RUN apt-get install -y \
    autotools-dev \
    automake \
    build-essential \
    git \
    python \
    python-dev \
    python-pip \
    zlibc \
    zlib1g \
    zlib1g-dev 

RUN pip install --upgrade cython

WORKDIR /www
RUN git clone --recursive https://github.com/aakrosh/BaitSTR \
    && cd BaitSTR \
    && make \
    && make install \
    && cp bin/* /usr/local/bin/ 

CMD ["bash"]


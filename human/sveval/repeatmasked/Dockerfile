FROM ubuntu:18.04

MAINTAINER jmonlong@ucsc.edu

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        wget \
        hmmer \
        unzip \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-dev \
        gcc \ 
        bcftools \
        tabix \
        build-essential \
        && rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir requests awscli snakemake

WORKDIR /usr/local/bin

RUN wget http://tandem.bu.edu/trf/downloads/trf407b.linux64 && mv trf*.linux64 trf && chmod +x trf

WORKDIR /usr/local

RUN wget http://repeatmasker.org/RepeatMasker-open-4-0-9-p2.tar.gz \
    && tar -xzvf RepeatMasker-open*.tar.gz \
        && rm -f RepeatMasker-open*.tar.gz

WORKDIR /usr/local/RepeatMasker

RUN perl ./configure -trfbin=/usr/local/bin/trf -hmmerbin=/usr/bin/nhmmscan

RUN cpan Text::Soundex

ENV PATH=/usr/local/RepeatMasker:$PATH

WORKDIR /home

ADD test.fa /home

RUN RepeatMasker --species human test.fa


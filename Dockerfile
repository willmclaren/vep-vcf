FROM python:3.8.2-slim-buster

WORKDIR /opt/vep-vcf
RUN apt-get update && apt-get install -y \
  build-essential \
  libcurl4-openssl-dev \
  libssl-dev \
  libbz2-dev \
  liblzma-dev \
  zlib1g-dev \
  curl

RUN curl -L https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 | tar xj
WORKDIR htslib-1.10.2
RUN ./configure && make && make install
WORKDIR ..

RUN pip install virtualenv
RUN python -m virtualenv venv
ENV VIRTUAL_ENV /opt/vep-vcf/venv
ENV PATH /opt/vep-vcf/venv/bin:$PATH    
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
RUN pip install -e .
RUN touch test-report.xml

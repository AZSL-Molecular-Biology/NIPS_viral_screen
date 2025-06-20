FROM registry.access.redhat.com/ubi8/ubi-minimal:8.4-208
MAINTAINER https://www.azdelta.be

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /conda/bin:$PATH
ENV NXF_VER 21.04.3

RUN mkdir /nextflow
COPY nextflow.config /nextflow/

RUN microdnf update -y \
  && microdnf install -y python39 python39-pip python39-devel gcc gcc-c++ which wget git bzip2 ca-certificates which java-1.8.0-openjdk tar postfix \
  && alternatives --set python /usr/bin/python3.9 \
  && alternatives --set python3 /usr/bin/python3.9 \
  && pip3.9 install --no-cache-dir --user pipenv \
  && microdnf clean all \
  && mkdir -p /tmp/miniconda \
  && cd /tmp/miniconda \
  && wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh -O miniconda.sh \
  && bash miniconda.sh -b -p /conda \
  && cd /tmp \
  && rm -r miniconda \
  && cd / \
  && groupadd conda_users \
  && chgrp -R conda_users /conda \
  && chmod 770 -R /conda \
  && ln -s /conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
  && echo ". /conda/etc/profile.d/conda.sh" >> ~/.bashrc \
  && echo "conda activate base" >> ~/.bashrc \
  && echo ". /conda/etc/profile.d/conda.sh" >> /etc/skel/.bashrc \
  && echo "conda activate base" >> /etc/skel/.bashrc \
  && source ~/.bashrc \
  && cd /nextflow \
  && wget https://github.com/nextflow-io/nextflow/releases/download/v21.04.3/nextflow-21.04.3-all -O nextflow.sh \
  && chmod 755 nextflow.sh \
  && bash nextflow.sh \
  && mkdir -p ~/.nextflow/ \
  && cp /nextflow/nextflow.config ~/.nextflow/config \
  && echo "#add nextflow config" >> /etc/skel/.bashrc \
  && echo "mkdir -p ~/.nextflow" >> /etc/skel/.bashrc \
  && echo "cp /nextflow/nextflow.config ~/.nextflow/config" >> /etc/skel/.bashrc

RUN microdnf install bzip2 cairo which unzip procps

RUN mkdir /code
COPY src /code/src/
COPY pipeline.nf /code/
COPY Pipfile /code/
COPY conda_yml /code/conda_yml/
RUN cd /code \
  && chmod 777 -R /code

#create the conda environments
RUN conda env create -f /code/conda_yml/trimmomatic.yml \
  && conda env create -f /code/conda_yml/mapping.yml \
  && conda clean -a -y

RUN pip3 install --no-cache-dir --user pipenv \
  && cd /code \
  && python3.9 -m pipenv install \
  && /nextflow/nextflow.sh

WORKDIR /code
ENTRYPOINT ["bash", "/nextflow/nextflow.sh", "/code/pipeline.nf"]

FROM ubuntu:22.04

#initial setup for base image
RUN apt update
RUN apt install -y zlib1g-dev
RUN apt-get install -y wget gcc build-essential bowtie2 git default-jre

#install miniconda3
RUN mkdir -p /opt/conda && mkdir -p /opt/kraken2/ && mkdir -p /usr/local/src/QC
#RUN mkdir -p /usr/local/src/QC
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/conda/miniconda.sh
RUN bash /opt/conda/miniconda.sh -b -u -p /opt/conda
RUN rm -rf /opt/conda/miniconda.sh

#use conda to install prerequisite software
RUN /opt/conda/bin/conda install -c bioconda samtools
RUN /opt/conda/bin/conda install -c bioconda fastqc

#install troublemakers from source
RUN wget http://opengene.org/fastp/fastp -O /opt/fastp && chmod a+x /opt/fastp
WORKDIR /opt/
RUN git clone https://github.com/DerrickWood/kraken2.git
WORKDIR /opt/kraken2
RUN bash install_kraken2.sh kraken2
RUN ln -s /opt/kraken2/kraken2/kraken2 /usr/local/bin && ln -s /opt/fastp /usr/local/bin


#add the QC scripts in
ADD QC_Workflow_PE.sh /usr/local/src/QC
ADD QC_Workflow_SE.sh /usr/local/src/QC

#ENTRYPOINT ["/bin/bash"]
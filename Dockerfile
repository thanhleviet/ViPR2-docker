FROM conda/miniconda2

MAINTAINER Thanh Le Viet thanhlv@oucru.org

#Prepare for conda
RUN apt-get update \
    && apt-get install -y build-essential dh-autoreconf \
    && conda config --add channels r \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && rm /bin/sh && ln -s /bin/bash /bin/sh

#Create python3 virtual env with conda
RUN conda create -n py3k python=3.5 \
    && source activate py3k \
    && conda install iva snakemake \
    && conda clean --all

#Install common bio tools
RUN conda install bwa \
                  samtools \
                  bamtools \
                  bcftools \
                  vcftools \
                  lofreq=2.1.2 \
                  freebayes=1.1.0 \
                  bedtools \
                  mummer=3.23 \
                  matplotlib \
                  graphviz \
                  trimmomatic \
                  fastqc \
                  multiqc=1.0 \
    && conda clean --all

#Install famas
RUN apt-get install -y git ttf-dejavu && cd /opt \
    && curl -L https://github.com/andreas-wilm/famas/archive/v0.0.11.tar.gz | tar xz \
    && cd famas-0.0.11 \
    && autoreconf --install --force 2>&1 > log.txt \
    && automake --add-missing \
    && ./configure \
    && make \
    && make install \
    && cd .. && rm -r famas-0.0.11 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

#Install simple_contig_joiner and python libraries: pysam, pyvcf, biopython, scipy
RUN curl -O https://raw.githubusercontent.com/andreas-wilm/simple-contig-joiner/master/simple_contig_joiner.py \
    && chmod +x simple_contig_joiner.py \
    && mv simple_contig_joiner.py /usr/local/bin/ \
    && pip install pyvcf scipy biopython pysam \
    && mkdir src

COPY src /src

RUN chmod 777 -R /src
# Create src folder for source code of ViPR2
ENV PATH /src:$PATH


# Add user biodocker with password biodocker
# https://github.com/BioContainers/containers/blob/master/biodocker/Dockerfile
RUN groupadd fuse \
    && useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse biodocker && \
    echo `echo "biodocker\nbiodocker\n" | passwd biodocker`

# Change user
USER biodocker

# ENTRYPOINT ["vipr2t.py"]

# CMD ["-h"]

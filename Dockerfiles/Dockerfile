# getting OS
FROM conda/miniconda3

MAINTAINER piroonj <"piroonj@gmail.com">

# Install git
RUN apt-get update && \
 apt-get install -y git

# Create working directory
RUN mkdir -p /home
# Set working directory
# WORKDIR /home
# RUN git clone https://gitlab.com/piroonj/eligos2.git
WORKDIR /home/CReSIL
COPY ./CReSIL/ /home/CReSIL/

# Create the environment:
RUN conda install -y -c anaconda -c bioconda -c conda-forge python=3.6.7 biopython=1.70 mappy=2.17=py36h84994c4_0 minimap2=2.17=h8b12597_1 python-intervaltree=3.0.2 tqdm=4.7.2  flye=2.6=py36he513fc3_0 pandas=0.24.2 numpy=1.12.1 pysam=0.15.3 pybedtools=0.8.0 python-graphviz=0.13.2 matplotlib=3.1.1 networkx=2.3 samtools bioawk
# RUN conda install --yes --file cresil_linux-explicit.txt

# RUN python pip install 

# Export eligos2 directory
ENV PATH /home/CReSIL:$PATH

# Activate the environment, and make sure it's activated:
# ENTRYPOINT ["eligos2"]
CMD ["echo","Welcome to CReSIL docker container"]


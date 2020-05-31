# CReSIL_NanoCircle

Python scripts and pipeline for detecting eccDNA from Nanopore reads



* Creating an environment with commands:
    ```bash
    ## Install environment
    conda create -n cresil -y -c anaconda -c bioconda -c conda-forge python=3.6.7 biopython=1.70 mappy=2.17=py36h84994c4_0 minimap2=2.17=h8b12597_1 python-intervaltree=3.0.2 tqdm=4.7.2  flye=2.6=py36he513fc3_0 pandas=0.24.2 numpy=1.12.1 pysam=0.15.3 pybedtools=0.8.0 python-graphviz=0.13.2 matplotlib=3.1.1 networkx=2.3 samtools bioawk

    ## Activate CReSIL environment
    conda activate cresil

    ## Export CReSIL to system environment
    export PATH=$PWD:$PATH

    ## Run CReSIL
    CReSIL-trim.py -h
    CReSIL-identify.py -h
    CReSIL-verify.py -h
    ```
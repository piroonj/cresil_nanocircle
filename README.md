# To access the full version of CReSIL with extra features, please click the link below. 
# https://github.com/visanuwan/cresil

--- 

# CReSIL (developed version)

Python scripts containing bioinfomatic pipeline for detecting eccDNA from long read Oxford Nanopore.

CReSIL (developing version)

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
    
* Run CReSIL:
    ```bash
    ## Run trim 
    python CReSIL/CReSIL-trim.py -i exp_reads.fastq -r hg19.25chr.mmi -o exp_reads
    
    ## Run eccDNA identification
    python CReSIL/CReSIL-identify.py -fq exp_reads.fastq -trim exp_reads.refTrim_map.txt -g hg19.25chr.fasta.fai -b exp1 -o eccdna_result
    
    ## Run verify eccDNA
    python CReSIL/CReSIL-verify.py -d eccdna_result -i eccdna_result/exp1/exp1.refTrim_map.eccDNA.txt -r hg19.25chr.mmi -o exp1.refTrim_map.eccDNA.verified.txt
    ```
---

# NanoCircle
NanoCircle is a tool developed for identifying the coordinates of both simple and chimeric circular molecules, sequenced using long-read sequencing. Contact Rasmus Amund Henriksen, wql443@alumni.ku.dk

Please refer to https://github.com/RAHenriksen/NanoCircle for the latest version of NanoCircle

---

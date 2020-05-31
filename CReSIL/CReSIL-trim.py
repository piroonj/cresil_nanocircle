#!/usr/bin/env python
import mappy
import intervaltree as tree
# import multiprocessing as mp
import argparse, re, sys
from os import path
from tqdm import tqdm
from Bio import SeqIO

def getNonOverlap(ref, name, seq):
    refDict = {}
    mergedRead = {}
    i=0
    nonOvl = []
    fst = True
    for hit in ref.map(seq): # alignments
        if hit.is_primary and hit.mapq == 60:
#             print(name,hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand)
            if fst:
                q_st_fst = hit.q_st
                q_en_fst = hit.q_en
                r_st_fst = hit.r_st
                r_en_fst = hit.r_en
                fst = False
                refDict.setdefault(hit.ctg, tree.IntervalTree()).add(tree.Interval(hit.r_st, hit.r_en))
                nonOvl.append([name,hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand])
            else:
                if hit.ctg in refDict:
                	checkOvl = refDict[hit.ctg].overlaps(hit.r_st, hit.r_en) 
                else:
                	checkOvl = False

                if checkOvl:
#                     print(refDict[hit.ctg][hit.r_st:hit.r_en])
                    for region in refDict[hit.ctg][hit.r_st:hit.r_en]:
#                         print(region.begin, region.end)
                        st = region.begin
                        en = region.end
                        sdata = sorted([hit.r_st, hit.r_en, st, en])
                        ovl = sdata[2]-sdata[1]
                        lenR = en-st
                        lenQ = hit.r_en-hit.r_st
                        covR = ovl/float(lenR)
                        covQ = ovl/float(lenQ)
                        if covQ < 0.05:
                            refDict.setdefault(hit.ctg, tree.IntervalTree()).add(tree.Interval(hit.r_st, hit.r_en, i))
#                             print(st, en, label, ovl, covR, covQ)
                            nonOvl.append([name,hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand])
                else:
                    nonOvl.append([name,hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand])
    return nonOvl

def getMerge(nonOvl_o):
    newNonOvl = []
    mergedRead = tree.IntervalTree()
    for row in nonOvl_o: 
        mergedRead.add(tree.Interval(row[1], row[2]))
    mergedRead.merge_overlaps()
    read_mapped_pos = mergedRead[nonOvl_o[0][1]:nonOvl_o[0][2]]
    for row in nonOvl_o:
        for m_pos in read_mapped_pos:
            if m_pos.overlaps(row[1], row[2]):
                newNonOvl.append(row)
    return(read_mapped_pos, newNonOvl)
        
def get_args():
    import argparse
    parser = argparse.ArgumentParser(description='Rolling circle ONT read trimming')
    general = parser.add_argument_group(title='General options')
    # general.add_argument('--threads', '-t', help="Number of threads", type=int, default=1)
    general.add_argument('-i', "--input",  
                            help="Fastq reads",
                            type=str, default=None)
    general.add_argument('-r', "--ref",  
                            help="reference Minimap2 index .mmi or reference sequence .fasta",
                            type=str, default=None)
    general.add_argument('-o', "--output",  
                            help="output basename",
                            type=str, default=False)
    general.add_argument('-fq', "--fq-out", dest='fqout', action='store_true',
                            help="output trimmed fastq (default: false)")   
    parser.add_argument('-out','--out-datatype', dest='out', type=str, default='raw',
                        choices=['raw', 'fastq', 'fasta'],
                        help="defind out data type (default: raw)")
    parser.add_argument('-in','--in-datatype', dest='intype', type=str, default='raw',
                        choices=['raw', 'fastq', 'fasta'],
                        help="defind in data type (default: raw)")
    parser.add_argument('-v','--version', action='store_true', dest='version',
                        help="")
#     args = parser.parse_args()
    return parser

def readFaidx(fastxFile):
    if path.isfile("{}.fai".format(fastxFile)):
        count=0
        with open("{}.fai".format(fastxFile)) as inf:
            for line in inf: count += 1
        return count
    else:
        exit("please index fasta/fastq using samtools faidx/fqidx")

def main(args):
    fname = args.input
    fref = args.ref
    # threads = args.threads
    ref = args.ref
    bname, ext = path.splitext(fname)

    if ext in ['.fasta', '.fa']:
        filetype = 'fasta'
    elif ext in ['.fastq', '.fq']:
        filetype = 'fastq'
    else:
        if args.intype == 'fastq':
            indataType = "fastq"
        elif args.intype == 'fasta':
            indataType = "fasta"
        else:
            exit("no fasta/fastq file")

    if args.output:
        if args.fqout:
            oname_seq = "{}.refTrim_seq{}".format(args.output, ext)
        oname_map = "{}.refTrim_map.txt".format(args.output)
    else:
        if args.fqout:
            oname_seq = "{}.refTrim_seq{}".format(bname, ext)
        oname_map = "{}.refTrim_map.txt".format(bname)

    ref = mappy.Aligner(fref, preset="map-ont")  # load or build index
    if not ref: 
    	raise Exception("ERROR: failed to load/build index")

    numSeq = readFaidx(fname)
#     print(dir(seqFile))
#     print(numSeq)
    if args.fqout:
        oname_seq_f = open(oname_seq, "w")
    oname_map_f = open(oname_map, "w")
    oname_map_f.write("#{}\n".format("\t".join(["query","q_start","q_end","ref","r_start","r_end","match","mapBlock","mapq","strand","qlenTrimmed","freqCov","order"])))
    with tqdm(total=numSeq, desc="Progress ") as pbar:
        for record in SeqIO.parse(fname, filetype):
            name = record.id
            seq = str(record.seq)
            nonOvl_o = getNonOverlap(ref, name, seq)
            if len(nonOvl_o) > 0:
                subReads, outTab = getMerge(nonOvl_o)
                subReads = list(subReads)[0]
                if len(outTab) > 0:
                    subrecord = record[subReads.begin:subReads.end]
                    subrecord.description = subrecord.description+" subSeq={}-{}".format(subReads.begin, subReads.end)
                    if args.fqout:
                        oname_seq_f.write("{:s}\n".format(subrecord.format(filetype).strip()))
                    for idx, rTab in enumerate(sorted(outTab, key=lambda x: x[2])):
                        oTab = rTab+[len(subrecord),"{:.2f}".format((rTab[2]-rTab[1])/len(subrecord)),idx]
                        oname_map_f.write("{:s}\n".format("\t".join(map(str,oTab))))
            pbar.update()

    if args.fqout:
        oname_seq_f.close()
    oname_map_f.close()
    pbar.close()
if __name__ == '__main__':
    parser = get_args()
    args = parser.parse_args()
    version="1.0.0"
    if args.version:
    	exit("{:s} v{:s}".format(sys.argv[0],version))
    if len(sys.argv) == 1:
        exit(parser.print_help())
    # cmd = "--input test.set.fa --ref chr1.exp0.fa --output outfile"
    # args = parser.parse_args(cmd.split())
    main(args)



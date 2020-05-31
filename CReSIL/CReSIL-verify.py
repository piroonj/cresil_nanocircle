#!/usr/bin/env python
import glob
import intervaltree as tree
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
import pybedtools as bt
import mappy as mp
import os, sys
# infile="check/linear/chr19_35473074.assembly.fasta"
def reorderSeq(ref, trimSeq):
    assmOutTmp = []
    for hit in ref.map(trimSeq, cs=True):
        if hit.is_primary:
            assmOutTmp.append([hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.strand])
    tmpAssmWOmerge=[]
    for hit in sorted(assmOutTmp, key=lambda x: (x[0], x[1])):
        if hit[5]>0:
            tmpAssmWOmerge.append([[hit[0],hit[1]],[hit[2],hit[3],hit[4]]])
        else:
            tmpAssmWOmerge.append([[hit[0],hit[1]],[hit[2],hit[3],hit[4]]])
    if len(tmpAssmWOmerge)>0:
        last_seg_aln = tmpAssmWOmerge[-1][0]
        reorder = trimSeq[int(last_seg_aln[0]):] + trimSeq[:int(last_seg_aln[0])]
        return(reorder)
    else:
        return(trimSeq)
        
def trimmedCircular(infile=None, circularDict = {}, refhg19=None):
    if not circularDict:
        exit("missing cirular information")
    fname, extName = os.path.splitext(infile)
    output = {}
    with open(infile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            name = record.id
            ref = str(record.seq)
#             print(name)
            oname = ".".join([fname,name,"trimmed"+extName])
            output.setdefault(name, {"refLen": len(ref), "trimmedBase": 0,'trimmedFile':oname,
                     "seq": ref, "circular": circularDict[name], 'trimReadLen': len(ref), "trimmed":"-"})
            new_seq = False
            trimmed_reads = ref
            a = mp.Aligner(seq=ref, preset="ava-pb")
            for hit in a.map(ref, cs=True): # traverse alignments
                if ((hit.q_st<100) and (hit.r_en>(hit.ctg_len-100)) and not (hit.q_st == hit.r_st and hit.q_en ==hit.r_en)):
                    ## trim from start (a half of overlap size)
                    new_seq = ref[int(hit.q_en/2):]
            if new_seq:
                a = mp.Aligner(seq=new_seq, preset="ava-pb")
                for hit in a.map(new_seq, cs=True): # traverse alignments
                    if ((hit.q_st<100) and (hit.r_en>(hit.ctg_len-100)) and not (hit.q_st == hit.r_st and hit.q_en ==hit.r_en)):
                        ## trim from end (a half of overlap size)
                        trimmed_reads=new_seq[:(hit.r_st-1)]
                output[name] = {"refLen": len(ref), "trimReadLen": len(trimmed_reads),'trimmedFile':oname,
                         "seq": trimmed_reads, "circular": "+", 'trimmedBase': len(ref)-len(trimmed_reads), "trimmed":"+"}
            if output:
#                 print(oname)
                with open(oname, "w") as outf:
                    outf.write(f">{name} length={output[name]['trimReadLen']} circular={output[name]['circular']} trimmedBase={output[name]['trimmedBase']}\n{output[name]['seq']}\n")
            else:
#                 print(oname)
                with open(oname, "w") as outf:
                    outf.write(f">{name} length={output[name]['trimReadLen']} circular={output[name]['circular']} trimmedBase={output[name]['trimmedBase']}\n{output[name]['seq']}\n")
    return(output)


def countCloseReads(assemblyF, readsF):
    with open(assemblyF) as assemblyF_h, open(readsF) as readsF_h:
        for record1 in SeqIO.parse(assemblyF_h, "fasta"):
            readid = record1.id
            asmseq = str(record1.seq)
            ref = mp.Aligner(seq=asmseq, preset="map-ont")
            totalReads = 0
            totalbases = 0
            countCloseAssem = 0
            for record2 in SeqIO.parse(readsF_h, "fasta"):
                readid = record2.id
                seq2 = str(record2.seq)
                nonOvl_o = getNonOverlap(ref, readid, seq2)
                if len(nonOvl_o) > 1:
                    ## if there is supplementary map of read to assembly
                    startEnd = [0,nonOvl_o[0][5]] # get start_pos and end_pos
                    mapEnd = []
                    mapReads = []
                    for l in nonOvl_o:
                        totalbases+=(l[8])
                        if l[7] - l[6] > 100: # size of mapped read > 100bp
                            for eachEnd in l[6:8]: ## get assem mapped start and end
                                if eachEnd < 100 or eachEnd > (startEnd[1]-100):
                                    ## collect part of the read that map end of assem
                                    mapEnd.append(eachEnd)
                                    mapReads.append(l)
                    if len(mapEnd) == 2: # if both end of assembly has read cover
                        countCloseAssem +=1
                else:
                    for l in nonOvl_o:
                        totalbases+=(l[8])
                totalReads+=1
            expectCov = totalbases/len(asmseq)
            depthCloseRatio = countCloseAssem/expectCov
            out = {'readLen':len(asmseq), 'totalReads':totalReads, 'totalBases':totalbases, 'expectCov':expectCov, 'closeSupport':countCloseAssem, 'depthCloseRatio':depthCloseRatio}
            return(out)


def getNonOverlap(ref, name, seq):
    refDict = {}
    mergedRead = {}
    i=0
    nonOvl = []
    fst = True
    for hit in ref.map(seq, cs=True): # alignments
        queryLen = len(seq)
        if hit.is_primary and hit.mapq == 60:
#             print(name,hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand)
            if fst:
                q_st_fst = hit.q_st
                q_en_fst = hit.q_en
                r_st_fst = hit.r_st
                r_en_fst = hit.r_en
                fst = False
                refDict.setdefault(hit.ctg, tree.IntervalTree()).add(tree.Interval(hit.r_st, hit.r_en))
                nonOvl.append([name,queryLen,hit.q_st,hit.q_en, hit.ctg, hit.ctg_len, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand])
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
                            nonOvl.append([name,queryLen,hit.q_st,hit.q_en, hit.ctg, hit.ctg_len, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand])
                else:
                    nonOvl.append([name,queryLen,hit.q_st,hit.q_en, hit.ctg, hit.ctg_len, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand])
    return nonOvl


def checkAssemOnGenome(ref, mergeRegions, trimSeq, assmBed, mergeBed):
    assmOutTmp = []
    assemMapOrder = []
    assemMapOrderStr = ""
    for hit in ref.map(trimSeq, cs=True):
        if hit.is_primary:
            assmOutTmp.append([hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.strand])
    assmOut = []
    for hit in sorted(assmOutTmp, key=lambda x: (x[0], x[1])):
        if hit[5]>0:
            assmOut.append([hit[2],hit[3],hit[4],"+"])
            assemMapOrder.append(f"{hit[2]}:{int(hit[3])+1}-{hit[4]}_+")
        else:
            assmOut.append([hit[2],hit[3],hit[4],"-"])
            assemMapOrder.append(f"{hit[2]}:{int(hit[3])+1}-{hit[4]}_-")
#     print(assmOut)
    assemMapOrderStr = ",".join(assemMapOrder)

    mergebeds = []
    for reg in mergeRegions.split(","):
        merge = re.split("[:-]",reg)
        merge[1] = str(int(merge[1])-1)
        mergebeds.append(merge)
#     print(mergebeds)

    with open(assmBed, "w") as assmf:
        assmf.write("\n".join(["\t".join(map(str,x)) for x in assmOut]))
    with open(mergeBed, "w") as mergef:
        mergef.write("\n".join(["\t".join(map(str,x)) for x in mergebeds]))
    
    checkAlnCount = countMergeFinal(assmBed)

    assemUniqBase = calSizeInterval(assmBed, mergeBed, "subtract")
    assemMergeIntersect = calSizeInterval(assmBed, mergeBed, "intersect")
    mergeUniqBase = calSizeInterval(mergeBed, assmBed, "subtract")
    return {"assemMapOrderStr":assemMapOrderStr, "assemUniqBase":assemUniqBase, 
            "assemMergeIntersect":assemMergeIntersect, "mergeUniqBase":mergeUniqBase,
            "assemblyMergeAln":checkAlnCount[0],"assemblyMergeAlnCount":checkAlnCount[1]}

import re
import subprocess

def calSizeInterval(a=None, b=None, method="intersect"):
    ### method could be intersct or subtract
    with subprocess.Popen(f'''bedtools {method} -a {a} -b {b}'''+'''| awk \'BEGIN{sum=0}{sum+=$3-$2}END{print sum}\'''',
                 stdout=subprocess.PIPE, shell=True) as proc:
        return(int(proc.stdout.read().decode("utf-8").strip()))
def countMergeFinal(bed):
    cmd = f'''bedtools sort -i {bed} | bedtools merge -d 10 '''+r'''| awk '{print $1":"$2+1"-"$3"\t1"}'  | bedtools groupby -g 2 -c 1,1 -o collapse,count | cut -f2-3'''
    with subprocess.Popen(cmd,
                 stdout=subprocess.PIPE, shell=True) as proc:
        return(proc.stdout.read().decode("utf-8").strip().split("\t"))

def get_args():
    import argparse
    parser = argparse.ArgumentParser(description='verify eccDNA assemblies')
    general = parser.add_argument_group(title='General options')
    # general.add_argument('--threads', '-t', help="Number of threads", type=int, default=1)
    general.add_argument('-d', "--eccdna-dir", dest='indir',
                            help="eccDNA result directory",
                            type=str, default=None)
    general.add_argument('-i', "--input",  
                            help="eccDNA identified table",
                            type=str, default=None)
    general.add_argument('-o', "--output",  
                            help="verified eccDNA output",
                            type=str, default=None)
    general.add_argument('-r', "--ref",  
                            help="reference Minimap2 index .mmi or reference sequence .fasta",
                            type=str, default=None)
    parser.add_argument('-v','--version', action='store_true', dest='version',
                        help="")

    return parser

def main(args):
    hg19idx = args.ref
    refhg19 = mp.Aligner(fn_idx_in=hg19idx, preset="map-ont")

    assemFiles = glob.glob(f'{args.indir}/**/assembly.fasta',recursive = True)
    eccAsm = pd.read_csv(args.input, sep="\t")

    ## refine code
    header = ['sample_name','circleID','mergeRegions','mergeRegionsLength','numNode','numEdge','numDiEdge',
             'assemblyName','assemblyLength','trimmedStatus','trimmedBase','lengthAfterTrim','totalMappedReads',
             'totalMappedBases','expectedCoverage','readsSupportCloseCircle','ratioCloseReads2ExpCov','circularStatus',
             'assemblyMappedOrder','intersectBasesAssembly2mergeRegion','uniqueBasesMergeRegion','uniqueBasesAssembly',
             'assemblyMergeAln','assemblyMergeAlnCount']
    assemTrimSummary = []
    fullpath = os.path.abspath(args.input)
    odir = os.path.dirname(fullpath)
    bname, ext = os.path.splitext(os.path.basename(fullpath))
    odir = os.path.dirname(fullpath)
    if args.output:
        outName = args.output
    else:
        outName = os.path.join(odir, bname+'.verified'+ext)
    with open(outName,'w') as outf:
        outf.write("{}\n".format("\t".join(map(str,header))))
        for idx, seqF in enumerate(assemFiles):
            infoF = seqF.replace("assembly.fasta","assembly_info.txt")
            fa4asm = "/".join(seqF.split("/")[0:2]+["fasta_trim"]+["{}.fa".format(seqF.split("/")[3])])
            dirPath = os.path.dirname(seqF)

            info = pd.read_csv(infoF, sep="\t")
            infoCir = {k:v for k,v in info.loc[:,['#seq_name','circ.']].values.tolist()}
            outTrim = trimmedCircular(seqF, infoCir, refhg19)
            for ctg in outTrim:
                assmBed = os.path.join(dirPath, f"assembly.{ctg}.bed")
                mergeBed = os.path.join(dirPath, f"merged_regions.{ctg}.bed")
                varidateReads = countCloseReads(outTrim[ctg]['trimmedFile'], fa4asm)
                
                finalCircleDecision = "incomplete_circle"
                if varidateReads['depthCloseRatio'] > 0.3 or varidateReads['closeSupport'] > 10:
                    finalCircleDecision = "real_circle"
                seqid = seqF.split("/")[3]
                print(ctg)
                if len(eccAsm.loc[(eccAsm['id']==seqid)&(eccAsm['assemName']==ctg),'merge_region']) == 0:
                    continue
                mergeRegions = eccAsm.loc[(eccAsm['id']==seqid)&(eccAsm['assemName']==ctg),'merge_region'].values[0]
                if finalCircleDecision == "real_circle":
                    seq_reorder = reorderSeq(refhg19, outTrim[ctg]['seq'])
                else:
                    seq_reorder = outTrim[ctg]['seq']
    #             outTrim[ctg]['trimmedFile']
                print("-------")
                with open(outTrim[ctg]['trimmedFile'], "r") as read_handle:
                    fnameFa, ext = os.path.splitext(outTrim[ctg]['trimmedFile'])
                    seq_reorderF = "{}.reorder{}".format(fnameFa, ext)
                    with open(seq_reorderF, "w") as out_handle:
                        for record in SeqIO.parse(read_handle, "fasta"):
                            record.seq = Seq(seq_reorder)
                            SeqIO.write(record, out_handle, "fasta")

                checkAssmAln = checkAssemOnGenome(refhg19, mergeRegions, seq_reorder, assmBed, mergeBed)

                row = eccAsm.loc[(eccAsm['id']==seqid)&(eccAsm['assemName']==ctg),
                                 ['sample_name','id','merge_region','merge_len',
                                  'numNode','numEdge','numDiEdge','assemName']].values[0]
                row = list(row)
                row += [outTrim[ctg]['refLen'],outTrim[ctg]['trimmed'],outTrim[ctg]['trimmedBase'],
                        outTrim[ctg]['trimReadLen']]
                row += [varidateReads['totalReads'], varidateReads['totalBases'], varidateReads['expectCov'], 
                        varidateReads['closeSupport'],varidateReads['depthCloseRatio'],finalCircleDecision]
                row += [checkAssmAln['assemMapOrderStr'],checkAssmAln['assemMergeIntersect'],
                        checkAssmAln['mergeUniqBase'],checkAssmAln['assemUniqBase']]
                row += [checkAssmAln['assemblyMergeAln'], checkAssmAln['assemblyMergeAlnCount']]

                outf.write("{}\n".format("\t".join(map(str,row))))
        
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

#!/usr/bin/env python
import os, sys
# env_var = os.environ.get('ENV_VAR', 'refTrimOnt')
import operator
import pandas as pd
import numpy as np
import pysam, shutil
import pybedtools as bt
import subprocess
import graphviz as gv
import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import to_agraph 
#from wand.image import Image as WImage
import pathlib
# %matplotlib inline
# from subprocess import Popen, PIPE

#####################################################################
## utility fuctions #################################################

def lenLoci(loci):
    length = 0
    for r in loci.split(","):
        chrom, start, end = r.split("_")
        length += int(end)-int(start)
    return length

def draw_graph(graph, filename="test"):
    import graphviz as gv
    d = gv.Digraph(name=filename)
    for k in graph:
        node = eval(k)
        d.edge(node[0], node[1], label=str(graph[k]))
    d.render()
    return filename

import math
def splitSeq(readLen, mergeLen, ratio=0.6, minLen=1100):
    expSize = round(mergeLen*ratio)
    numLoop = math.ceil(readLen/expSize)
    mem = 0
    out = []
    for i in range(numLoop):
        if mem+expSize > readLen:
            if readLen-mem >= minLen:
                out.append([mem, readLen])
                mem=readLen
        else:
            out.append([mem, mem+expSize])
            mem+=(expSize)
    return(out)

def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)

def get_args():
    import argparse
    parser = argparse.ArgumentParser(description='Identify eccDNA')
    general = parser.add_argument_group(title='General options')
    # general.add_argument('--threads', '-t', help="Number of threads", type=int, default=1)
    general.add_argument('-fq', "--fq-input",  dest='fqinput',
                            help="input fasta/fastq",
                            type=str, default=None)
    general.add_argument('-trim', "--trim-input", dest='triminput',
                            help="eccDNA result directory",
                            type=str, default=None)
    general.add_argument('-g', "--genome-size", dest='gsize',
                            help="genome size file e.g. hg19.chrom.sizes or file.fasta.fai",
                            type=str, default=None)
    general.add_argument('-b', "--bname",  
                            help="basename prefix",
                            type=str, default=False)
    general.add_argument('-o', "--output",  
                            help="output directory",
                            type=str, default="eccdna_result")
    parser.add_argument('-v','--version', action='store_true', dest='version',
                        help="")
    return parser


def main(args):
    ## input ############################################################
    ## python eccDNA_ident.py data/hg19.chrom.sizes data/.fastq.gz data/.refTrim_map.txt bnameGlobal outDir

    # chromSizes = "hg19.chrom.sizes"
    # bnameGlobal= "lib2_8s_BC01"
    # readTrim = pd.read_csv("{}.refTrim_map.txt".format(bnameGlobal),sep="\t",header=0)
    # fastaName = "{}.fastq.gz".format(bnameGlobal)

    chromSizes = args.gsize
    fastaName = args.fqinput
    readTrim = pd.read_csv(args.triminput,sep="\t",header=0) 
    if not args.bname:
        bnameGlobal, ext = os.path.splitext(os.path.basename(args.fqinput))
    else:
        bnameGlobal = args.bname

    tmpDir="{}/{}/tmp".format(args.output,bnameGlobal)
    outDir="{}/{}".format(args.output,bnameGlobal)

    # focusListF = sys.argv[6]
    # focusList_set = set()
    # if focusListF:
    #     with open(focusListF) as inlist:
    #         for l in inlist:
    #             focusList_set.add(l.strip())


    # hg19.chrom.sizes
    # bnameGlobal
    # .refTrim_map.txt
    # .fastq.gz


    ## input ############################################################

    #####################################################################
    ## output ###########################################################

    pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(tmpDir).mkdir(parents=True, exist_ok=True)
    fastaOut = "{}/fasta_trim".format(outDir,bnameGlobal)
    fastaOutFinal = "{}/fasta_final".format(outDir,bnameGlobal)
    assemOut = "{}/assem_trim_flye".format(outDir,bnameGlobal)
    outFileName="{}/{}.refTrim_map.eccDNA.txt".format(outDir,bnameGlobal)
    ## output ###########################################################

    ## utility fuctions #################################################

    #####################################################################
    ## Identify potential eccDNA regions ################################

    ## Calculate read coverage and filter out restuion with depth <= 5x
    header = list(readTrim.columns)
    header[0] = "readid"
    readTrim.columns = header
    ord_header = ['ref', 'r_start', 'r_end', 'readid', 'q_start', 'q_end', 'match',
                  'mapBlock', 'mapq', 'strand', 'qlenTrimmed', 'freqCov', 'order']
    readTrim = readTrim.loc[:,ord_header]

    ## prepare aligned reads
    aln_reads = bt.BedTool.from_dataframe(readTrim).sort().moveto("{}/aln_reads.bed".format(tmpDir))

    ## calculate reads coverage by betools genomecov
    genome_cov = aln_reads.genome_coverage(bg=True, g=chromSizes).moveto("{}/reads_genomecov.bdg".format(tmpDir))

    ## extract merged regions based on algned reads
    aln_reads_merge = aln_reads.merge().moveto("{}/reads_merged.bed".format(tmpDir))

    ## join merged regions with read coverage and perform filtering
    merge_genomecov = aln_reads_merge.intersect(genome_cov, output='{}/merge_genomecov.txt'.format(tmpDir), wo=True)
    merge_genomecov_df = pd.read_csv("{}/merge_genomecov.txt".format(tmpDir),sep="\t",header=None,dtype=str)
    merge_genomecov_df.columns = ['m_chrom','m_start','m_end','bg_chrom','bg_start','bg_end','depth','d_ovl']
    merge_genomecov_df[['depth','bg_start','d_ovl']] = merge_genomecov_df[['depth','bg_start','d_ovl']].apply(pd.to_numeric, errors='coerce')

    ## generate mergeID
    merge_genomecov_df['mergeid'] = merge_genomecov_df['m_chrom'] + "_" + merge_genomecov_df['m_start'] + "_" + merge_genomecov_df['m_end']

    ## trim low-coverage regions (depth <=5x)
    merge_genomecov_df_filt = merge_genomecov_df[merge_genomecov_df['depth']>5]
    merge_genomecov_df = merge_genomecov_df.sort_values(['mergeid','bg_start'])
    ## group reads coverage of each mergeid and calculate mean coverage
    merge_bg_filt = merge_genomecov_df_filt.groupby(['mergeid']).agg({'bg_chrom': np.max,'bg_start' : np.min,
                                                           'bg_end' : np.max,'depth' : np.mean}).reset_index()
    ## create final potential eccDNA regions and 
    ## re-create mergeid based on trimming low-coverage regions
    merge_bg_filt[['bg_start','bg_end']] = merge_bg_filt[['bg_start','bg_end']].apply(pd.to_numeric, errors='coerce')
    merge_bg_filt['length'] = merge_bg_filt.bg_end - merge_bg_filt.bg_start
    merge_bg_filt["bg_start"]= merge_bg_filt["bg_start"].astype(str) 
    merge_bg_filt["bg_end"]= merge_bg_filt["bg_end"].astype(str) 
    merge_bg_filt['mergeid'] = merge_bg_filt.bg_chrom+"_"+merge_bg_filt.bg_start+"_"+merge_bg_filt.bg_end
    merge_bg_filt = merge_bg_filt.loc[:,['bg_chrom','bg_start','bg_end','depth','length','mergeid']]
    ## Identify potential eccDNA regions ################################

    #####################################################################
    ## Assign readid to each potential eccDNA regions (mergeid) #########

    ## create merge region
    # print("{}/aln_reads.bed".format(tmpDir))
    aln_reads = bt.BedTool.from_dataframe(readTrim).sort().moveto("{}/aln_reads.bed".format(tmpDir))
    aln_reads_merge = bt.BedTool.from_dataframe(merge_bg_filt)
    read_merged_intersect = aln_reads_merge.intersect(aln_reads, output='{}/reads_merge_intersect.bed'.format(tmpDir), wo=True)
    ## read merge region and create mergeid
    read_merged_ins_df_org = pd.read_csv("{}/reads_merge_intersect.bed".format(tmpDir),sep="\t",header=None,dtype=str)
    header = list(merge_bg_filt.columns)+list(readTrim.columns)+["ovl"]
    read_merged_ins_df_org.columns = header


    header_select = ['ref','r_start','r_end','mergeid','depth','length','readid', 
                     'q_start','q_end','match','mapBlock','mapq','strand','qlenTrimmed', 
                     'freqCov','order','ovl']
    read_merged_ins_df_org = read_merged_ins_df_org.loc[:,header_select]
    ## Assign readid to each potential eccDNA regions (mergeid) #########

    #####################################################################
    ## Annotate 5'end and 3'end region of each readid ###################

    ## collect 200bp from 5'end and 3'end of each potential eccDNA
    end5 = open("{}/end5_merge_region.bed".format(tmpDir),'w')
    end3 = open("{}/end3_merge_region.bed".format(tmpDir),'w')
    for idx_, row in merge_bg_filt.iterrows():
        chrom = row.bg_chrom
        start = row.bg_start
        end = row.bg_end
        end_size = 200
        end5.write("{}\n".format("\t".join(map(str, [chrom, start, int(start)+end_size]))))
        end3.write("{}\n".format("\t".join(map(str, [chrom, int(end)-end_size, end]))))
    end5.close()
    end3.close()
    bt_5end = bt.BedTool("{}/end5_merge_region.bed".format(tmpDir))
    bt_3end = bt.BedTool("{}/end3_merge_region.bed".format(tmpDir))
    ## Annotate each read that has 5'end and 3'end overlap
    ## add column ovl_5end and ovl_3end (0=no overlap, 1=presence of overlap) 
    read_merged_ins_df_org_bt = bt.BedTool.from_dataframe(read_merged_ins_df_org)
    read_merged_final_bt = read_merged_ins_df_org_bt.annotate(files=["{}/end5_merge_region.bed".format(tmpDir),"{}/end3_merge_region.bed".format(tmpDir)], counts=True)
    read_merged_ins_df = read_merged_final_bt.to_dataframe(names=header_select+['ovl_5end','ovl_3end'])
    ## sum 5'end and 3'end overlap (0=no overlap, 1=either 5'end or 5'end overlap, 2=both ends overlap)
    read_merged_ins_df['sum_ends'] = read_merged_ins_df.ovl_5end + read_merged_ins_df.ovl_3end
    ## Annotate 5'end and 3'end region of each mergeid ##################

    #####################################################################
    ## Identify reads containing breakpoints ############################

    ## by grouping readid and count number of mergeid
    read_merged_ins_df = read_merged_ins_df.sort_values(['readid','mergeid'], ascending = [True,True])
    groupLoci = read_merged_ins_df
    groupLoci = groupLoci.sort_values(by=['readid', 'order'], ascending=[True,True])
    groupLoci = groupLoci.groupby(['readid']).agg({'mergeid':lambda x: ','.join(x),'readid': 'count' })
    groupLoci.columns = ['mergeid','count']
    groupLoci = groupLoci.reset_index()
    numReads = len(groupLoci.readid.unique())
    ## separate reads with single mergeid and with multiple mergeid (breakpoints)
    groupLoci_single=read_merged_ins_df.loc[read_merged_ins_df.readid.isin(groupLoci.loc[groupLoci['count']==1,'readid']),
                                            ['readid','mergeid','q_start','q_end']]
    groupLoci_single['strand'] = 1
    groupLoci_single.set_index('readid', inplace=True)
    groupLoci_breakpoints = groupLoci.loc[groupLoci['count']>1,:]
    ## Identify reads containing breakpoints ############################

    #####################################################################
    ## Refine continuity of breakpoints #################################
    # print(groupLoci_breakpoints)
    ## only collect breakpoint that occur in 5'end or 3'end (200bp) 
    reads_breakpoints = []
    j = 0
    for idx, row in groupLoci_breakpoints.iterrows():
        read = read_merged_ins_df.loc[read_merged_ins_df['readid']==row['readid']].sort_values(by=['readid', 'order'], ascending=[True,True])
        collect = []
        for idx_, row in read.iterrows():
            if len(collect) == 0 and row.sum_ends > 0:
                collect.append(list(row))
            elif len(collect) > 0 and row.sum_ends == 2:
                collect.append(list(row))
            elif len(collect) > 0 and row.sum_ends == 1:
                collect.append(list(row))
                break
        if len(collect) > 1:
            reads_breakpoints = reads_breakpoints + collect
    reads_breakpoints = pd.DataFrame(reads_breakpoints)
    # print(reads_breakpoints)
    reads_breakpoints.columns = read_merged_ins_df.columns

    # create final table for reads containing breakpoints
    breakpointsFilt = reads_breakpoints.sort_values(by=['readid', 'order'], ascending=[True,True])
    breakpointsFilt = breakpointsFilt.groupby(['readid']).agg({'mergeid':lambda x: ','.join(x) , 'q_start': 'min', 'q_end':'max','strand':lambda x: x.nunique()})

    finalReadSet = pd.concat([breakpointsFilt,groupLoci_single])
    groupLociFilt = breakpointsFilt.reset_index()
    finalReadSet.reset_index().to_csv("{}/{}_final_read_set.txt".format(outDir,bnameGlobal), sep="\t", index=False)
    ## Refine continuity of breakpoints #################################

    #####################################################################
    ## create network graphs from all breakpoints #######################

    graph = {}
    for idx, row in breakpointsFilt.iterrows():
        nodes = row['mergeid'].split(",")
        for i in range(0,len(nodes)-1):
            inter = repr([nodes[i], nodes[i+1]])
            graph[inter] = graph.setdefault(inter, 0) + 1
    ## Exclude breakpoints that less than 3 evidences ****
    graphFilt = {k: v for k, v in graph.items() if v > 2}
    draw_graph(graphFilt, '{}/{}_all_graphs'.format(outDir,bnameGlobal))

    G = nx.MultiDiGraph()
    for k in graphFilt:
        nodes = eval(k)
        G.add_edge(nodes[0], nodes[1], weight=graphFilt[k])

    # nx.draw_spring(G, with_labels=False,node_size=5)
    # plt.draw()
    # plt.show()
    ## create network graphs from all breakpoints #######################

    #####################################################################
    ## Perform de novo assembly for each potential eccDNA ###############
    ## and create summary table

    ## define input and output file
    # output potential eccDNA are listed in file with .refTrim_map.eccDNA.txt
    # outFileName="lib4_bc2_21washed.refTrim_map.eccDNA.txt"
    outFile = open(outFileName, "w")
    outFile.write("{}\n".format("\t".join(["sample_name","id", "merge_region", "merge_len", "numreads", 
                             "totalbase","expectCov", "numNode", "numEdge", "numDiEdge","assemby", "assembySizeDiff",
                             "assemName","assemLength", "assemCov", "assemCirc", 
                             "assemRepeat","assemMult", "assemGraph_path","graph_detail","edge_detail"])))

    # fastaName = "lib4_bc2_21washed.fastq.gz"
    out_assem = []
    success = 0
    notsucces = 0
    # i=1
    # fastaOut = "fasta_trim_tmp"
    # assemOut = "assem_trim2_flye_tmp"
    pathlib.Path(fastaOut).mkdir(parents=True, exist_ok=True)
    pathlib.Path(assemOut).mkdir(parents=True, exist_ok=True)
    pathlib.Path(fastaOutFinal).mkdir(parents=True, exist_ok=True)

    ## go to each connected graph and perform de novo assembly using Flye software
    subgraphs = list(connected_component_subgraphs(G.to_undirected()))
    for g in subgraphs:
        nodes = list(g.nodes())
        bname = ".".join([ "_".join(x.split("_")[:2]) for x in nodes ])
        # if focusListF:
        #     if not bname in focusList_set:
        #         continue
        g.name = bname
        selectEdges = {}
        ## collect breakpoints from a given network graph
        for k in graphFilt:
            checkEdges = eval(k)
            if checkEdges[0] in nodes and checkEdges[1] in nodes:
                selectEdges[k] = graphFilt[k]
        minEdge = min(selectEdges.values())
        maxEdge = max(selectEdges.values())
        numNodes = g.number_of_nodes()
        numEdges = g.number_of_edges()
        numDiEdges = len(selectEdges)
        print(nx.info(g))
        print('Number of directional edges:',numDiEdges)
        print('nodes:',g.nodes())
        mergeRegions = []
        for node in nodes:
            (c, s, e) = node.split("_")
            mergeRegions.append(c+":"+str(int(s)+1)+"-"+e)
        mergeRegion = ",".join(mergeRegions) 
        ## draw network graph
        # imgFile = draw_graph(selectEdges, '{}/multi'.format(tmpDir))
        # img = WImage(filename="{}/multi.gv.pdf".format(tmpDir))
        # display(img)
        
        ## collect all reads from the network graph
        readInNode=read_merged_ins_df.loc[read_merged_ins_df.mergeid.isin(nodes),'readid']
        common = set(finalReadSet.index.values)&set(readInNode.values)
        if len(nodes) == 1:
            readO = finalReadSet.loc[common,].loc[finalReadSet.strand==1,]
        else:
            readO = finalReadSet.loc[common,]
        uniqReadNum = readO.shape[0]
        lengthLoci = lenLoci(",".join(nodes))
        baseCov = 0
        print("total reads:",readO.shape[0])
        readO['size'] = readO['q_end'] - readO['q_start']
        ## perform de novo assembly if the network graph has at least 5 reads 
        if uniqReadNum >= 5:
            fa = pysam.FastaFile(fastaName)
            out_filename = "{}/{}.fa".format(fastaOut,bname)
            outSeq = []
            for idx, row in readO.sort_values(by='size', ascending=True).iterrows():
                getSeq = fa.fetch(row.name)[row.q_start:row.q_end]
                for subIdx, subRead in enumerate(splitSeq(len(getSeq),lengthLoci,ratio=0.6)):
                    subS = subRead[0]
                    subE = subRead[1]
                    outSeq.append([subE-subS, ">{}-{}".format(row.name,subIdx), getSeq[subS:subE]])
                    baseCov += subE-subS
            with open(out_filename, mode='w') as fout:
                outSeq.sort(key=lambda x: x[0])
                for subOut in outSeq:
                    fout.write("{} len={}\n{}\n".format(subOut[1],subOut[0],subOut[2]))

            out_filename_final = "{}/{}_final.fa".format(fastaOutFinal,bname)
            with open(out_filename_final, mode='w') as foutF:
                numReads = 0
                for idx, row in readO.sort_values(by='size', ascending=True).iterrows():
                    getSeq = fa.fetch(row.name)[row.q_start:row.q_end]
                    foutF.write(">{}\n{}\n".format(row.name,getSeq))
                    numReads += 1
            threads=15
            ## set command line for flye assembly
            cmd = "flye --nano-raw {} --out-dir {}/{} --genome-size {} --threads {} --plasmids --meta".format(out_filename,assemOut,bname,lengthLoci,threads)
            ##cmd = "flye --nano-raw {} --out-dir {}/{} --genome-size {} --threads {} --plasmids --meta -i 2".format(out_filename,assemOut,bname,lengthLoci,threads)
            ## perform assembly
            process = subprocess.call(cmd, shell=True)
            if process:
                notsucces+=1
                print("assembly not success : {} {} {} {}".format(notsucces, bname, lengthLoci, numReads))
    #             shutil.rmtree("{}/{}".format(assemOut,bname))
    #             os.remove(out_filename)
                expectCov = "{:.2f}".format((float(baseCov)/lengthLoci))
                outList = [bnameGlobal,bname, mergeRegion, lengthLoci, numReads, baseCov, expectCov,
                      numNodes, numEdges, numDiEdges,"FALSE"]+['-']+['-', 0, 0, '-', '-', 0, 0, 0]+[str(repr(selectEdges)).replace(" ","")]
                outFile.write("{}\n".format("\t".join(map(str, outList))))
                print("{}\n".format("\t".join(map(str, outList))))
            else:
                success+=1
                print("assembly success : {} {} {} {}".format(success, bname, lengthLoci, numReads))
                info = pd.read_csv("{}/{}/assembly_info.txt".format(assemOut,bname),sep="\t")
                print(info)
                expectCov = "{:.2f}".format((float(baseCov)/lengthLoci))
                for assmIdx, assemRow in info.iterrows():
                    assem_summary = list(assemRow)
                    diffAssm = str(int(assem_summary[1])-int(lengthLoci))
                    outList = [bnameGlobal,bname, mergeRegion, lengthLoci, numReads, baseCov, expectCov,
                         numNodes, numEdges, numDiEdges, "TRUE"]+[diffAssm]+assem_summary+[str(repr(selectEdges)).replace(" ","")]
                    outFile.write("{}\n".format("\t".join(map(str, outList))))
                    print("{}\n".format("\t".join(map(str, outList))))
    #         if i==5:
    #             break
    #         i+=1
    #         break
        print("-------------------------")
        ## bname 
        ## id merge_region merge_len numreads totalbase expect_cov node edge assemby assemName assemLength assemCov assemCirc assemRepeat assemMult assemGraph_path
    outFile.close()

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

#!/usr/bin/env python3

import Utils
import pysam as ps
import numpy as np
from collections import Counter
from itertools import groupby
import pandas as pd
import numpy as np
import pybedtools as bt
import operator

class Circles:
    #Command 1
    def __init__(self,regions,bamfile,output,MapQ):
        self.regions = regions
        self.bamfile = bamfile
        self.output = output
        self.MapQ = MapQ

    def Prim_dict(self,bamfile,reg, start, end):
        """Creates a dictionary of the primary alignments for the soft-clipped reads"""
        Prim_dict = {}
        for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
            # checks if soft-clipped
            if Utils.IS_SA(read, self.MapQ) == True:
                pos = ps.AlignedSegment.get_reference_positions(read)

                # creates dict with key being read_id and value being read ref positions
                Prim_dict[read.query_name] = [pos[0], pos[-1]]
        return Prim_dict

    def Supplement_dict(self,bamfile,reg, start, end):
        """Creates a dictionary of the primary soft-clipped reads and
        their corresponding supplementary alignments"""

        Primary_dict = self.Prim_dict(bamfile,reg, start, end)
        SA_dict = {}

        for read in bamfile.fetch(reg, start, end, multiple_iterators=True):

            # only extract those supp in primary
            if read.query_name in Primary_dict.keys():
                if Utils.IS_supp(read, self.MapQ) == True:

                    # extracts positions
                    supp_pos = ps.AlignedSegment.get_reference_positions(read)
                    prim_pos = Primary_dict[read.query_name]

                    # Adding to positions based on direction to the empty dict
                    if read.query_name not in SA_dict.keys():

                        # From left to right
                        if Utils.Right_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [prim_pos[0], supp_pos[-1]]

                        # From right to left
                        if Utils.Left_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [supp_pos[0], prim_pos[-1]]

                        # From left to right once
                        if Utils.Right_circ_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]

                        # From right to left once
                        if Utils.Left_circ_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]

                    # Appends for reads with several supplementary alignments their position
                    elif read.query_name in SA_dict.keys():

                        if Utils.Right_dir(prim_pos, supp_pos) == True:
                            if supp_pos[-1] not in SA_dict[read.query_name]:
                                SA_dict[read.query_name].append(supp_pos[-1])
                        # From right to left
                        if Utils.Left_dir(prim_pos, supp_pos) == True:
                            if prim_pos[-1] not in SA_dict[read.query_name]:
                                SA_dict[read.query_name].append(prim_pos[-1])

                        # From left to right once
                        if Utils.Right_circ_dir(prim_pos, supp_pos) == True:
                            if prim_pos[-1] not in SA_dict[read.query_name]:
                                SA_dict[read.query_name].append(prim_pos[-1])

                        # From right to left once
                        if Utils.Left_circ_dir(prim_pos, supp_pos) == True:
                            if prim_pos[-1] not in SA_dict[read.query_name]:
                                SA_dict[read.query_name].append(prim_pos[-1])

        return SA_dict

    def reduce_end_coords(self,bamfile,reg, start, end):
        """ Reduce the number of end coordinates for those reads with several
        supplementary alignments by comparing the position to the other end positions"""

        Coord = self.Supplement_dict(bamfile,reg, start, end)
        # Counter of second list element for 2-element lists.
        # I dont have several primary so there is no need for counting them
        count = Counter(v[1] for v in Coord.values() if len(v) == 2)
        # Result dict
        reduce_dic = {}
        # Iterate data entries
        for k, v in Coord.items():
            # Modify lists longer than two with at least one element in the counter
            if len(v) > 2 and any(elem in count for elem in v[1:]):
                # Replace list with first element and following element with max count
                v = [v[0], max(v[1:], key=lambda elem: count.get(elem, 0))]
            # Add to result
            reduce_dic[k] = v
        return reduce_dic

    def most_frequent(self,List, pos):
        count1 = Counter(List)
        dic_val = list(count1.values())
        max_val = max(count1.values())
        occurence = dic_val.count(max_val)

        # single most occuring coordinate
        if occurence == 1:
            occ = 1
            return occ, count1.most_common(1)[0][0]

        # creates the interval immediately
        else:
            most_common_no = list(map(operator.itemgetter(1), count1.most_common()))
            most_common_idx = [i for i, x in enumerate(most_common_no) if x == max(most_common_no)]

            # if there are several equal common values
            if max(most_common_no) > 1:
                # print(max(most_common_no))
                # print(most_common_idx)
                # creating the most maximum value
                common_val = [count1.most_common()[i][0] for i in most_common_idx]
                occ = 1
                return occ, max(common_val)

            # if there isnt a most occuring coordinate, the interval is created
            else:
                occ = 2
                # start
                if pos == 0:
                    return occ, min(List)
                # end
                if pos == 1:
                    return occ, max(List)

    def Reads(self,bamfile, reg, start, end, reads_IDs):
        """ extract those reads which is found using the single coordinate dictionary"""
        Read_list = []
        for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
            if read.query_name in reads_IDs:
                if Utils.IS_SA(read, self.MapQ) == True:
                    Read_list.append(read)
        return Read_list

    def chr_coord_sort(self,chrlist, coordlist):
        """ Sort a list of chromosomes and their coordinates using the index of the numerically sorted coordinates"""
        coord_idx = np.argsort(coordlist)
        Coord_sort = [coordlist[i] for i in coord_idx]
        chr_sort = [chrlist[i] for i in coord_idx]
        return Coord_sort, chr_sort

    def Grouping_chr(self,Chr_sort, Group_size):
        """ Groups the chromosomes in to match the grouped coordinates """
        Grouped_chr = []
        Step = 0
        for i in Group_size:
            Grouped_chr.append(Chr_sort[Step:Step + i])
            Step += i
        return Grouped_chr

    def Grouping(self,Coord_list, overlap_bp):
        """ Groups a given list, for which elements within the overlap range is grouped together"""
        First = [[Coord_list[0]]]
        for i in Coord_list[1:]:
            # sets the previous and current elem as the first elem
            prev = First[-1]
            current = prev[-1]
            dist = abs(i - current)

            # if dist between is below overlapping bp, then it it belong to the same group
            if dist < overlap_bp:
                prev.append(i)
            else:
                First.append([i])
        return First

    def Single_coord(self, bamfile,reg, start, end):
        """ returns the most frequent start and end coordinate for a circle in a dictionary with the key
        being all the reads. If not having a coordinate more frequent than others, a dictionary with all reads
        and the interval for potential circle """
        Circle_dict = self.reduce_end_coords(bamfile,reg, start, end)
        start = [i[0] for i in Circle_dict.values()]
        # taking all possible end coordinates and merge into one list
        end = sum([i[1:] for i in Circle_dict.values()], [])
        # extract the most common
        if start or end != []:
            occ1, start_freq = self.most_frequent(start, 0)
            occ2, end_freq = self.most_frequent(end, 1)

            if occ1 and occ2 == 1:
                reads = []
                chr = [reg]

                if start_freq < end_freq:
                    new_val = chr + [start_freq, end_freq]
                else:
                    new_val = chr + [end_freq, start_freq]
                for k, v in Circle_dict.items():
                    if any(i == start_freq or i == end_freq for i in v):
                        reads.append(k)
                final_dict = {tuple(sorted(reads)): new_val}

                # Multiple reads
                if len(list(final_dict.keys())[0]) != 1:
                    type = 1
                    return type, final_dict

                # Single read
                else:
                    type = 2
                    return type, final_dict

            # not a more supported read
            else:
                type = 2
                chr = [reg]
                new_val = chr + [start_freq, end_freq]
                final_dict = {tuple(sorted(Circle_dict.keys())): new_val}
                return type, final_dict
        else:
            type = 0
            # these dicts are empty, and serves to continue as some regions might not create a circle
            return type, Circle_dict

    def Is_Simple(self, bamfile,reg, start, end):
        """ Check for potential complex circles by reads aligning across the genome. Afterwards it returns the simple circles
         with 1 set of coordinates (type I) several set of coordinates (type II) and the complex circles (type III)"""
        Type, Sing_dict = self.Single_coord(bamfile,reg, start, end)

        Complex_dict = {}
        Total_coord = []
        Total_chr = []
        Total_overlap = []

        if len(Sing_dict) == 1:
            reads_IDs = list(Sing_dict.keys())[0]
            Read_list = self.Reads(bamfile, reg, start, end, reads_IDs)
        else:
            reads_IDs = list(Sing_dict.keys())
            Read_list = self.Reads(bamfile, reg, start, end, reads_IDs)

        for read in Read_list:
            coord1 = list(Sing_dict.values())[0][-2]
            coord2 = list(Sing_dict.values())[0][-1]
            Total_chr.extend((reg, reg))
            Total_coord.extend((coord1, coord2))
            Tag = read.get_tag("SA").split(';')[:-1]
            Coord_list = []
            chroms = []
            cigar_len = []

            # examining for complexity
            for Tagelem in Tag:
                # splitting the individual aligment info up
                Column_list = Tagelem.split(',')
                chrom = Column_list[0]
                length = Utils.CIGAR_len(Column_list[3])
                cigar_len.append(length)
                # 1 - based positions
                pos_start = int(Column_list[1]) - 1
                pos_end = int(Column_list[1]) + length - 1
                # the overlaps between coordinates for grouping
                overlap = sum(cigar_len) * 4
                Total_overlap.append(overlap)
                # if the supp align is in between the circular input region then we already know the breakpoint from Sing_dict
                if chrom == reg and start - overlap <= pos_start <= end + overlap and start - overlap <= pos_end <= end + overlap:
                    continue

                elif int(Column_list[4]) >= self.MapQ:
                    # creates a coordinate list
                    Coord_list.append(pos_start)
                    Coord_list.append(pos_end)

                    # append chr twice to ensure same length as coord_list
                    chroms.append(chrom)
                    chroms.append(chrom)

                    Total_chr.extend((chrom, chrom))

                    Total_coord.extend((pos_start, pos_end))

                if Coord_list != []:

                    # sorts the chr and coordinates
                    Coord_sort, chr_sort = self.chr_coord_sort(chroms, Coord_list)
                    # first entry with the input region
                    Complex_dict[read.query_name] = [reg, coord1, coord2]

                    if len(Coord_sort) == 2:
                        # one coordinate pair supporting another region
                        Complex_dict[read.query_name].extend(
                            [chr_sort[0], int(min(Coord_sort)), int(max(Coord_sort))])
                    else:
                        Grouped_coord = self.Grouping(Coord_sort, max(Total_overlap))
                        Group_size = [len(i) for i in Grouped_coord]
                        Grouped_chr = self.Grouping_chr(chr_sort, Group_size)
                        for i in range(len(Grouped_coord)):
                            Complex_dict[read.query_name].extend(
                                [Grouped_chr[i][0], int(min(Grouped_coord[i])), int(max(Grouped_coord[i]))])

        # the simple circles
        if Complex_dict == {}:
            if Type == 0 or 1 or 2:
                return Sing_dict, Type, None, None
        else:
            # Sorting again to create the coordinate and chromosome list
            Coord_sort, chr_sort = self.chr_coord_sort(Total_chr, Total_coord)
            Grouped_coord = self.Grouping(Coord_sort, max(Total_overlap))
            Group_size = [len(i) for i in Grouped_coord]
            Grouped_chr = self.Grouping_chr(chr_sort, Group_size)
            # circle type 3, complex circles comprised of several chromosomal regions across the genome.
            Type = 3

            return Complex_dict, Type, Grouped_coord, Grouped_chr

    def Simple_circ_df(self,beddict, Circ_no, circ_type):
        """ returns a dataframe with circular information for the simple circles, type 1 and 2 """
        df_col = ["Chr", "Start", "End"]
        simple_df = pd.DataFrame.from_dict(beddict, orient='index', columns=df_col)
        simple_df = simple_df.sort_values(by=df_col)
        simple_df['Length'] = simple_df.apply(lambda x: x['End'] - x['Start'], axis=1)
        # add_col = ['Read_No','Read_IDs','Circle_type','Circle_ID']
        add_col = ['Read_No', 'Circle_type', 'Circle_ID']
        # convert list of ID's to 1 long string as to insert it as a single column in the df

        add_val = [len(list(beddict.keys())[0]), circ_type, "simple_circ_%d" % Circ_no]

        for i in range(len(add_col)):
            simple_df.insert(loc=len(simple_df.columns), column=add_col[i], value=add_val[i])
        return simple_df

    def Complex_full_length(self,file_name, bamfile):
        """ returns the total number of columns needed for the most complex circle"""
        total_col_len = 0
        with open(file_name) as f:
            for line in f:
                line_value = line.strip().split()
                coord = line_value[0]
                start = line_value[1]
                end = line_value[2]
                complex_dict, circ_type, circ_coord, circ_chr = self.Is_Simple(bamfile, str(coord), int(start), int(end))
                if circ_type == 3:
                    if len(circ_coord) > total_col_len:
                        total_col_len = len(circ_coord)
        return total_col_len

    def Complex_circ_df(self,dict, coord_full, chr_full, Circ_no):
        """ returns a dataframe with circular information for the complex circles, type 3 """
        d = {}
        tot_len = 0

        for i in range(len(coord_full)):
            if len(d) == 0:
                d["Circle_temp"] = [chr_full[i][0], int(min(coord_full[i])), int(max(coord_full[i])),
                                    int(max(coord_full[i])) - int(min(coord_full[i]))]
                tot_len += (max(coord_full[i]) - min(coord_full[i]))
            else:
                d["Circle_temp"].extend([chr_full[i][0], int(min(coord_full[i])), int(max(coord_full[i])),
                                         int(max(coord_full[i])) - int(min(coord_full[i]))])
                tot_len += (int(max(coord_full[i])) - int(min(coord_full[i])))
        df_col = ["Chr", "Start", "End", "Length"]

        # number of times df_col needed to be repeated
        rep = len(coord_full)

        # creates a dataframe with chromosome information
        Coord_Col = [j + "_no._" + str(i) for i in range(1, rep + 1) for j in df_col]
        complex_df = pd.DataFrame.from_dict(d, orient='index', columns=Coord_Col)

        # creates a dataframe with additional information
        add_col = ['Total_len', 'Read_No', 'Read_IDs', 'Circle_ID']
        Read_ID_str = str(list(list(dict.keys()))).replace(" ", "")
        add_val = [int(tot_len), len(list(dict.keys())), Read_ID_str,"Complex_circ_%d" % Circ_no]
        chr_col = {add_col[0]: add_val[0], add_col[1]: add_val[1], add_col[2]: add_val[2],
                   add_col[3]: add_val[3]}
        df_final5 = pd.DataFrame(columns=add_col)
        df_final5 = df_final5.append(chr_col, ignore_index=True)

        return complex_df, df_final5

    def Circle_output(self):
        Simple_count = 1
        Simple_circ = pd.DataFrame()

        Complex_count = 1
        Complex_df = pd.DataFrame()
        Final_5 = pd.DataFrame()

        Read_File = Utils.SamBamParser(self.bamfile)

        with open(self.regions) as f:
            for line in f:
                region = line.strip().split()
                coord = region[0]
                start = region[1]
                end = region[2]

                circle_dict, circ_type, circ_coord, circ_chr = self.Is_Simple(Read_File,str(coord), int(start), int(end))

                if circ_type == 1:
                    circ_bed = self.Simple_circ_df(circle_dict, Simple_count, "high_conf")
                    rows = pd.concat([Simple_circ, circ_bed])
                    Simple_circ = rows
                    Simple_count += 1

                elif circ_type == 2:
                    if len(list(circle_dict.keys())[0]) > 1:
                        circ_bed = self.Simple_circ_df(circle_dict, Simple_count, "conf")
                        rows = pd.concat([Simple_circ, circ_bed])
                        Simple_circ = rows
                        Simple_count += 1
                    else:
                        circ_bed = self.Simple_circ_df(circle_dict, Simple_count, "low_conf")
                        rows = pd.concat([Simple_circ, circ_bed])
                        Simple_circ = rows
                        Simple_count += 1

                elif circ_type == 3:
                    cird_df, Last5 = self.Complex_circ_df(circle_dict, circ_coord, circ_chr, Complex_count)

                    Complex_df = Complex_df.append(cird_df, sort=False)
                    Final_5 = Final_5.append(Last5)

                    Complex_count += 1
                else:
                    continue

            Simple_bed = bt.BedTool.from_dataframe(Simple_circ)
            Simple_bed.saveas("{0}_Simple.bed".format(str(self.output)))

            resulting_df = pd.concat([Complex_df.reset_index(drop=True), Final_5.reset_index(drop=True)], axis=1)
            pd.set_option("display.max_rows", None, "display.max_columns", None)

            cols = list(resulting_df)[0:-4]
            del cols[0::4]

            # ensuring coordinates are int
            resulting_df[cols] = resulting_df[cols].astype(float).astype('Int64')

            Bed_test = bt.BedTool.from_dataframe(resulting_df)
            Bed_test.saveas("{0}_Complex.bed".format(str(self.output)))


if __name__ == '__main__':
    print("main")


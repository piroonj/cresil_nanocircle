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

class Merge_config:
    #Command 2
    def __init__(self,circ_bed,output,dist):
        self.circ_bed = circ_bed
        self.output = output
        self.dist = dist

    def Region_seperate(self):
        """
        Splits the chimeric regions into seperate lines and merge the overlapping lines
        :return: a merged bed format file with the 4th col being the merged chimeric ID's
        """
        Merge_df = pd.read_csv(self.circ_bed, sep='\t',header=None)

        dataframe_val = Merge_df.values

        df = pd.DataFrame()
        circ_no = 0
        for line in dataframe_val:
            # the different column values from the created bed file
            chr_regions = line[:-4]
            chr_val = chr_regions[0::4]
            start_val = chr_regions[1::4]
            end_val = chr_regions[2::4]

            circ_no += 1

            chr_col = pd.DataFrame({"Chr": chr_val, "Start": start_val, "End": end_val,"ID": "Chimeric_%s" % circ_no})
            df = df.append(chr_col)

        #replace empty cells (.) with NaN and remove these rows
        df = df.replace('\.+', np.nan, regex=True).dropna()
        df['Start'] = df['Start'].astype(float).astype('Int64')
        df['End'] = df['End'].astype(float).astype('Int64')

        df = df.sort_values(by=["Chr","Start","End"],ascending=[True,True,True])
        bedformat = bt.BedTool.from_dataframe(df)

        Collapsed = bedformat.merge(d=self.dist,c=[4], o="collapse")

        #print(df.sort_values(['Chr','Start','End'], ascending=[True,True,True]))

        #chr_order = ["chr" + str(i) for i in range(1, 23)]+['chrX','chrY']
        #df['Chr'] = pd.Categorical(df['Chr'], chr_order)

        #Bed_test = bt.BedTool.from_dataframe(df)
        #Collapsed.saveas("lol.bed")

        return Collapsed

    def Collapsed_split(self):
        """
        splits the collapsed bed file up into different information
        :return:
        """

        input = self.Region_seperate()
        line_val_list = []
        id_val_list = []

        # opens and splits up the collapsed columns with regions and circle ID's
        for line_value in input:
            # line_value[0:3] => each chromosome region
            line_val_list.append(line_value[0:3])

            #list with the collapsed circle ID
            id_val_list.append(line_value[3].split(","))

        #print("line val", line_val_list)
        #print("ID val", id_val_list)

        #print("id",id_val_list)
        #print(line_val_list)
        return id_val_list, line_val_list

    def overlap_region(self,list1,list2):
        """
        :param list1: first chromosome region e.g ['chr1', 37777195, 37778735]
        :param list2: second chromosome region e.g. ['chr18', '4219198', '4222801']
        :param nt_overlap: The number of nt the regions can have between them before
        :return: list of interval [chr,min,max] if the distance between the regions are below nt_overlap. Otherwise
        the regions are concatenated e.g.  ['chr1', '37777195', '37778735','chr18', '4219198', '4222801']
        """

        #print("Chromosome", list1, list2,list1[-3],list2[0])

        # if the chromosome is the same there is potential for overlap, compare previous region list1[-3] to the next
        # region list[0]
        if list1[-3] == list2[0]:
            # if there are overlap between the regions create the interval
            if abs(int(list1[-2]) - int(list2[1])) < self.dist and abs(int(list1[-1]) - int(list2[2])) < self.dist:
                start = min(int(list1[-2]), int(list2[1]))
                end = max(int(list1[-1]), int(list2[2]))
                return (list1[:-3]+[str(list1[-3]),start,end])
            #with no overlap concatenate the list to add the next regions
            else:
                return(list1 + list2)

        else:
            return(list1 + list2)

    def modulo_length(self,list):
        """
        :param list: list with the different chromosomal fragmens chr:start:end:chr:start:end...
        :return: a list with the structure chr:start:end:length:chr:start:end:length...
        """
        final_list = []
        for i in range(len(list)):
            if i % 3 == 0:
                # chromosome
                final_list.append(list[i])

            elif i % 3 == 1:
                # start coordinate
                final_list.append(int(list[i]))

            elif i % 3 == 2:
                #last coordinate and final length
                final_list.append(int(list[i]))
                final_list.append(int(final_list[-1] - final_list[-2]))

        return final_list, sum(final_list[3::4])

    def Merge_circles(self):
        """
        Iterates through the collapsed regions and compare those regions from same chimeric ID to merge these regions
        :return:list of list of all regions,
                list of total lengths of the chimeric ecccDNA,
                list of list of chimeric ID for merged circle
        """
        id_val_list, line_val_list = self.Collapsed_split()

        #number of different regions
        idx_list = [i for i in range(len(id_val_list))]
        j_idx_list = []

        Chr_reg_list = []
        total_ID = []

        #each individual region
        for idx in idx_list:
            #print("iterate", idx, line_val_list[idx], id_val_list[idx])

            ID_list = []

            if idx in j_idx_list:
                #print("pass",idx,j_idx_list)
                # if the regions has already been merged pass
                pass
            else:
                #extracts the chr region
                line_values = line_val_list[idx]

                #compare this single region with all the others
                for j_idx in range(len(idx_list)):

                    #checks if the ID for this regions matches any of the other regions ID
                    if any(e in id_val_list[idx] for e in id_val_list[j_idx]):
                        #print("Jindex", j_idx, id_val_list[idx], id_val_list[j_idx], list(set(id_val_list[idx]+ id_val_list[j_idx])))

                        #appends the index of the matching regions based on the ID
                        j_idx_list.append(j_idx)

                        #list of all ID for the concatenated regions
                        ID_list.append(id_val_list[j_idx])

                        #if the any ID matches, check if the regions are overlapping
                        line_values = self.overlap_region(line_values, line_val_list[j_idx])
                        line_values = line_values

                        #updates the list of ID to the unique set of merged ID
                        id_val_list[idx] = list(set(id_val_list[idx] + id_val_list[j_idx]))
                    else:
                        continue

                #All the concatenated chromosome pieces
                Chr_reg_list.append(line_values)

            if ID_list != []:
                Several_IDs = list(set([Circ_ID for subreg in ID_list for Circ_ID in subreg]))
                Several_IDs = repr(Several_IDs).replace(" ","")
                #list of unique merged ID for each merged circle
                total_ID.append(Several_IDs)

        max_col = max(len(list) for list in Chr_reg_list)
        # list of full chromosome regions for each merged chimeric eccDNA
        Chr_reg_full = [self.modulo_length(i)[0] for i in Chr_reg_list]
        # list of total length of the merged circles
        Total_length = [self.modulo_length(i)[1] for i in Chr_reg_list]

        return Chr_reg_full,Total_length,max_col,total_ID

    def circle_df(self):
        """
        :return: a bed file with the merged chimeric eccDNA information
        """
        Chr_reg_list, Total_length, max_col, id_val_list = self.Merge_circles()
        print(id_val_list)
        ID = ["Merged_circ_{0}".format(i) for i in range(1, len(Chr_reg_list) + 1)]

        rep = int(max_col / 3)
        df_col = []
        for i in range(1, rep + 1):
            df_col.extend(["Chr_{0}".format(i), "Start_{0}".format(i), "End_{0}".format(i), "Length_{0}".format(i)])

        df = pd.DataFrame(Chr_reg_list, columns=df_col)

        del df_col[0::4]
        # ensuring coordinates are int
        df[df_col] = df[df_col].astype(float).astype('Int64')

        add_col = ['Total_len', 'Config_ID', 'Circle_ID']
        add_val = [Total_length, id_val_list, ID]

        for i in range(len(add_col)):
            print("II",i)
            print(add_val[i])
            df.insert(loc=len(df.columns), column=add_col[i], value=add_val[i])

        # ensure the chromosome lengths, and other numbers are integers
        Merged_complex = bt.BedTool.from_dataframe(df)
        Merged_complex.saveas(self.output)



if __name__ == '__main__':
    print("main")


import sys
import pysam
import pandas as pd
import numpy as np
from numpy import zeros
import pickle

#This script contains functions for calculating the binding intensity of ChIP-exo or ChIP-mini data (RPPM).

def load_csv(raw_csv_dir):
    raw_df = pd.read_csv(raw_csv_dir, sep = ",", header= 0, names =["sample_id","bam_1","bam_2","peak","output", "stat"])
    return raw_df

def count_reads_position(samfile):      
    total_read = 0 
    tem_count = []
    total_count = {}

    if "SQ" in samfile.header:
        chromosome_sizes = {}
        for entry in samfile.header["SQ"]:
            chromosome_sizes[entry["SN"]] = int(entry["LN"]) + 1
    else:
        for reference in samfile.references:
            chromosome_sizes[reference] = chromosome_size
    for reference in samfile.references:                                  #Make array for each chromosome
        tem_count.append(zeros(((chromosome_sizes[reference],))))
    for i, read in enumerate(samfile):
        if read.is_unmapped:
            continue        
        if not read.is_proper_pair:
            tem_count[read.tid][read.pos-1] += 1
            total_read +=1
        width = read.qlen      

    for x, reference in enumerate(samfile.references):
        total_count[reference] = {}    
        total_count[reference] = tem_count[x]
    return total_count, width, total_read

def cal_rppm(bam_filename_1, bam_filename_2, peak_dir, rppm_output):
    samfile_1 = pysam.Samfile(bam_filename_1)
    samfile_2 = pysam.Samfile(bam_filename_2)
    peak = pd.read_csv(peak_dir, sep='\t', names = ["genome_id", "source", "feature", "start", "end","score", "strand", "frame", "attr"])
    total_count_1, width, total_read_1 = count_reads_position(samfile_1)
    total_count_2, width, total_read_2 = count_reads_position(samfile_2)
    rppm = pd.DataFrame()
    rppm_mean = pd.DataFrame()

    for  genome_id in total_count_1.keys():
        read_1 = []
        read_2 = [] 
        start_lst = []
        end_lst = []
        genome_lst = []
        
        for i,idx in peak.iterrows():
            start = peak.start[i]
            end = peak.end[i]
            peak_width = int(end-start)
            if peak_width == 0:
                sys.exit("The same position information of Start and end exists in the Peak file")

            if genome_id == peak.genome_id[i]:
                read_num_1 = (sum(total_count_1[genome_id][start-width-1:end-1])*1e6)
                read_num_2 = (sum(total_count_2[genome_id][start-width-1:end-1])*1e6)
                read_1.append(float(read_num_1/(peak_width*total_read_1)))
                read_2.append(float(read_num_2/(peak_width*total_read_2)))                      
                genome_lst.append(genome_id)
                start_lst.append(start)
                end_lst.append(end)

        genome_df = pd.DataFrame(genome_lst)
        start_df = pd.DataFrame(start_lst)
        end_df = pd.DataFrame(end_lst)
        rppm_1 = pd.DataFrame(read_1)
        rppm_2 = pd.DataFrame(read_2)    
        rppm_div = pd.DataFrame(((rppm_1 + rppm_2)/2))
        
        rppm_cat = pd.concat((genome_df, start_df, end_df, rppm_1, rppm_2), axis = 1)
        rppm_mean_cat = pd.concat((genome_df, start_df, end_df, rppm_div), axis = 1)
        
        rppm = pd.concat((rppm, rppm_cat), axis = 0)
        rppm_mean = pd.concat((rppm_mean, rppm_mean_cat), axis = 0)
        
    df_rppm = rppm
    df_rppm.columns = ["Genome_id", "Start", "End", "RPPM_1", "RPPM_2"]

    df_rppm.to_csv(rppm_output + "_rppm.csv" )
    del read_1, read_2, rppm_1, rppm_2
    return rppm, rppm_mean

def initiate_rppm_cal(raw_df):
    for i,idx in raw_df.iterrows():
        sample_id = raw_df.sample_id[i]
        bam_filename_1 = raw_df.bam_1[i]
        bam_filename_2 = raw_df.bam_2[i]
        peak_dir = raw_df.peak[i]
        rppm_output = raw_df.output[i] + sample_id
        rppm_tem, rppm_mean_tem = cal_rppm(bam_filename_1, bam_filename_2, peak_dir, rppm_output)

        rppm_tem.columns = [sample_id+"_genome_id", sample_id+"_start", sample_id+"_end", sample_id + "_1",sample_id + "_2"]
        rppm_mean_tem.columns = [sample_id+"_genome_id", sample_id+"_start", sample_id+"_end", sample_id ]
        
        if  i == 0:
            rppm_df = rppm_tem.reset_index()
            rppm_mean_df = rppm_mean_tem.reset_index()
        else:
            rppm_df = pd.concat((rppm_df, rppm_tem.reset_index()), axis = 1)
            rppm_mean_df = pd.concat((rppm_mean_df, rppm_mean_tem.reset_index()), axis = 1)

        print (("Sample ID: %s RPPM calculation is done!") % (sample_id))
    rppm_df = rppm_df.drop(['index'],axis = 1)
    rppm_mean_df = rppm_mean_df.drop(['index'],axis = 1)
    
    return rppm_df, rppm_mean_df

def start_PNC(raw_df, result_name):
    rppm_df, rppm_mean_df  = initiate_rppm_cal(raw_df)    
    output_dir = raw_df.output[0] + result_name + "_rppm.pickle"
    output_mean_dir = raw_df.output[0] + result_name + "_mean_rppm.pickle"    
    with open (output_dir, "wb") as rppm:
        pickle.dump(rppm_df, rppm)   
    with open (output_mean_dir, "wb") as rppm_mean:
        pickle.dump(rppm_mean_df, rppm_mean)
    return rppm_df, rppm_mean_df
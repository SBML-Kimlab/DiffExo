import pandas as pd
import pickle

#This script contains functions for making raw files for DEseq2 in R and assembling the DEseq2 result with the RPPM result.

def make_DBPref(raw_df,rppm_dir):
    if type(rppm_dir) == pd.DataFrame:
        rppm_df = rppm_dir
    else:
        with open (rppm_dir, "rb") as dataframe:
            rppm_df = pickle.load(dataframe)       

    rppm_dic1 = {}
    for row,idx in raw_df.iterrows():
        rppm_dic2 = {}
        tem = {}
        stat = raw_df.stat[row]
        sample = raw_df.sample_id[row]
        result_dir = raw_df.output[row]

        for i,idx in rppm_df.iterrows():
            if stat == 'control':  # start: [id, end, rppm_1, rppm_2]
                if type(rppm_df[sample+"_genome_id"][i]) is not float:
                    rppm_dic1[rppm_df[sample + "_start"][i]] = [rppm_df[sample+"_genome_id"][i], rppm_df[sample+"_end"][i], 
                                                              rppm_df[sample +"_1"][i], rppm_df[sample +"_2"][i]]
            else: 
                if type(rppm_df[sample+"_genome_id"][i]) is not float:
                    rppm_dic2[rppm_df[sample + "_start"][i]] = [rppm_df[sample+"_genome_id"][i], rppm_df[sample+"_end"][i],
                                                              rppm_df[sample+"_1"][i], rppm_df[sample+"_2"][i]]   

        if stat != "control":    #Find common peaks.
            for x in rppm_dic1.keys():
                id_1 = rppm_dic1[x][0]
                start_1 = int(x)
                end_1 = int(rppm_dic1[x][1])
                for y in rppm_dic2.keys():
                    id_2 = rppm_dic2[y][0]
                    start_2 = int(y)
                    end_2 = int(rppm_dic2[y][1])
                    if id_1 == id_2:
                        if start_1 == start_2 and end_1 == end_2:          #c_start = overlap start position, c_end = overlap end position
                            c_start = start_1
                            c_end = end_1
                            tem[start_1] = [start_2, id_1, c_start, c_end]
                        elif start_1 <= start_2 and end_1 >= end_2:
                            c_start = start_1
                            c_end = end_1
                            tem[start_1] = [start_2, id_1, c_start, c_end]  
                        elif start_1 >= start_2 and end_1 <= end_2:
                            c_start = start_2
                            c_end = end_2
                            tem[start_1] = [start_2, id_1, c_start, c_end]
                        elif start_1 <= end_2 and end_1 >= end_2:
                            c_start = start_2
                            c_end = end_1  
                            tem[start_1] = [start_2, id_1, c_start, c_end]         
                        elif start_2 <= end_1 and end_2 >= end_1:
                            c_start = start_1
                            c_end = end_2
                            tem[start_1] = [start_2, id_1, c_start, c_end]                        

            num = 0            
            peak_df = pd.DataFrame()
            with open(result_dir + sample + "_DEseq2_common_peaks_ref.gff", "w") as f:
                for start in sorted(tem):
                    num += 1
                    g_id = rppm_dic1[start][0]

                    c_start = int(tem[start][2]) 
                    c_end = int(tem[start][3])

                    start_1 = int(start)
                    end_1 = int(rppm_dic1[start][1])     
                    start_2 = int(tem[start][0])
                    end_2 = int(rppm_dic2[start_2][1])
                    rppm_1 = float(rppm_dic1[start][2])  
                    rppm_2 = float(rppm_dic1[start][3])  

                    rppm_3 = float(rppm_dic2[start_2][2])  
                    rppm_4 = float(rppm_dic2[start_2][3])
                    peak = ("Peak_" +str(num).zfill(4)) 

                    exp_list = [peak, g_id, c_start, c_end,start_1, end_1, rppm_1, rppm_2, g_id, start_2, end_2, rppm_3, rppm_4]
                    exp_df = pd.DataFrame(exp_list).transpose()
                    peak_df = pd.concat((peak_df, exp_df), axis = 0)

                    peak_attr = ('gene_id "%s"; transcript_id "%s"')%("p" + str(num).zfill(4), "p" + str(num).zfill(4))      
                    f.write(("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n")%(g_id, "feature", "exon", c_start, c_end, ".", "+" ,".", peak_attr))

            peak_df = peak_df.reset_index(drop = True)
            peak_df.columns = ["peak_num", "id1", "overlap_start", "overlap_end", "start_1", "end_1", "control_rppm_1", "control_rppm_2","id2", "start_2","end_2","sample_rppm_1", "sample_rppm_2"]
            with open (result_dir + sample + "_peak_rppm.pickle", "wb") as peak_rppm:
                pickle.dump(peak_df, peak_rppm)
    print ("Making reference peak file for DESeq2 is done :)")    
        
def convert_DEseq2(peak_pickle_dir, dbp_dir, result_dir):
    dbp = pd.read_csv(dbp_dir, header= 0, sep = ",", names = ["peak_num","baseMean", "lfc", "lfcSE", "stat", "pval","qval"])
    
    with open (peak_pickle_dir, "rb") as dataframe:
        peak_df = pickle.load(dataframe)   
    peak_info = {}
    dbp_info = {}
    
    for i,idx in peak_df.iterrows():
        peak_info[peak_df.peak_num[i]] = list(peak_df.iloc[i])
    
    for j,idx in dbp.iterrows():
        dbp_peak = dbp.peak_num[j].replace("p", "Peak_")
        lfc = float(dbp.lfc[j])
        pval = float(dbp.pval[j])
        qval = float(dbp.qval[j])
        dbp_info[dbp_peak] = [lfc, pval, qval]
    
    final_dic = {}
    for x in sorted(dbp_info):
        if x in peak_info:
            peak_info[x].extend(dbp_info[x])
            final_dic[x] = peak_info[x]
  
    final = pd.DataFrame(final_dic.values())
    final.columns = ["peak_num", "id1", "overlap_start", "overlap_end","start_1", "end_1", "control_rppm_1", "control_rppm_2","id2", "start_2","end_2","sample_rppm_1", "sample_rppm_2", "LFC", "p-value", "q-value"]
    final.to_csv(result_dir, mode = 'w')
    
    print ("DEseq2 convertion is done :)")
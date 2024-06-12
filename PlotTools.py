import sys
import pickle
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from bioinfokit import analys, visuz

#This script contains functions for generating plots of ChIP-exo or ChIP-mini intensity data.

def make_boxplot(raw_df, rppm_mean_dir, result_name, log_scale = "log10", y_range = None, show_fliers= False):  #mean_dir = DataFrame or pickle files
    if type(rppm_mean_dir) == pd.DataFrame:
        rppm_mean_df = rppm_mean_dir
    else:
        with open (rppm_mean_dir, "rb") as dataframe:
            rppm_mean_df = pickle.load(dataframe)        
    
    fig_result_dir = raw_df.output[0] + "plots/" + result_name 
    del_cols = []   
    for i,idx in raw_df.iterrows():
        sample_id = raw_df.sample_id[i]
        del_cols.extend([sample_id+"_genome_id", sample_id+"_start", sample_id+"_end"])
    rppm_box = rppm_mean_df.drop(del_cols,axis = 1)   

    boxplot = rppm_box.boxplot(column=list(rppm_box.columns),showfliers=show_fliers)
    plt.xlabel("Samples")
    plt.ylabel("RPPM")
    if y_range != None:
        plt.ylim(y_range)
    plt.savefig(("%s RPPM boxplot result.png")%(fig_result_dir))
    plt.savefig(("%s RPPM boxplot result.svg")%(fig_result_dir))
    plt.show()
    
    if log_scale == "log10":
        log_rppm_box = np.log10(rppm_box.astype(float) + 1)
        boxplot = log_rppm_box.boxplot(column=list(log_rppm_box.columns),showfliers=show_fliers)
        plt.xlabel("Samples")
        plt.ylabel("Log10(RPPM)")
        plt.savefig(("%s Log10 RPPM boxplot result.png")%(fig_result_dir))
        plt.savefig(("%s Log10 RPPM boxplot result.svg")%(fig_result_dir))
        plt.show()
        print ("Making a Boxplot (Log10) is done :)")
    elif log_scale == "log2":
        log_rppm_box = np.log2(rppm_box.astype(float) + 1)
        boxplot = log_rppm_box.boxplot(column=list(log_rppm_box.columns),showfliers=show_fliers)
        plt.xlabel("Samples")
        plt.ylabel("Log2(RPPM)")
        plt.savefig(("%s Log2 RPPM boxplot result.png")%(fig_result_dir))
        plt.savefig(("%s Log2 RPPM boxplot result.svg")%(fig_result_dir))
        plt.show()
        print ("Making a Boxplot (Log2) is done :)")
    elif log_scale == "ln":
        log_rppm_box = np.log(rppm_box.astype(float) + 1)
        boxplot = log_rppm_box.boxplot(column=list(log_rppm_box.columns),showfliers=show_fliers)
        plt.xlabel("Samples")
        plt.ylabel("LN(RPPM)")
        plt.savefig(("%s LN RPPM boxplot result.png")%(fig_result_dir))
        plt.savefig(("%s LN RPPM boxplot result.svg")%(fig_result_dir))
        plt.show()
        print ("Making a Boxplot (LN) is done :)")

def make_scatterplot(raw_df, rppm_mean_dir, replicates = False, stat = "pearson"):     #If replicates = true => rppm_pickle 
    if type(rppm_mean_dir) == pd.DataFrame:
        rppm_df = rppm_mean_dir
    else:
        with open (rppm_mean_dir, "rb") as dataframe:
            rppm_df = pickle.load(dataframe) 
    output_dir = raw_df.output[0]

    if replicates == True:
        if len(rppm_df.columns) == 10:  
            for i,idx in raw_df.iterrows():
                sample_id = raw_df.sample_id[i]
                rppm_1 = list(rppm_df[sample_id + "_1"].dropna())
                rppm_2 = list(rppm_df[sample_id + "_2"].dropna())
                z = np.polyfit(rppm_1, rppm_2, 1)
                p = np.poly1d(z)
                plt.scatter(rppm_1, rppm_2, marker='o', color = "sandybrown", linewidths = 0.75, edgecolors = "black" )
                plt.plot(rppm_1, p(rppm_1), color="red", linewidth=1, linestyle="--")
                plt.title("Correlation of Replicates")
                plt.xlabel(sample_id + "_1")
                plt.ylabel(sample_id + "_2")
                plt.savefig(output_dir + "plots/"+ sample_id + "_replicates_scatter.svg")
                plt.show()    

                if stat == "pearson":
                    print ("pearson p-value " ,(pearsonr(rppm_1, rppm_2)))
                elif stat == "spearman":
                    print ("spearman p-value ", (spearmanr(rppm_1, rppm_2)))
        else:
            print ("Change input pickle file to [Test_name]_rppm.pickle")
    else:
        if len(rppm_df.columns) != 10:  
            inter_rppm = intersect_peaks(raw_df, rppm_df)  
            print ("The number of common is " + str(len(inter_rppm)))     
            con_rppm = list((inter_rppm[inter_rppm.columns[0]]))
            exp_rppm = list((inter_rppm[inter_rppm.columns[1]]))
            z = np.polyfit(con_rppm, exp_rppm, 1)
            p = np.poly1d(z)
            plt.scatter(con_rppm, exp_rppm, marker='o', color = "pink", linewidths = 0.75, edgecolors = "black" )
            plt.plot(con_rppm, p(con_rppm), color="red", linewidth=1, linestyle="--")
            plt.title("Control RPPM vs Sample RPPM")
            plt.xlabel('Control RPPM')
            plt.ylabel('Sample RPPM')
            plt.savefig(output_dir + "plots/"+ "Control_vs_Sample_scatter.svg")
            plt.show()

            if stat == "pearson":
                print ("pearson p-value " ,(pearsonr(con_rppm, exp_rppm)))
            elif stat == "spearman":
                print ("spearman p-value ", (spearmanr(con_rppm, exp_rppm)))
        else:
            print ("Change input pickle file to [Test_name]_mean_rppm.pickle")

def intersect_peaks(raw_df, rppm_mean_df):    
    rppm_dic1 = {}
    common = {}
    for row,idx in raw_df.iterrows():
        stat = raw_df.stat[row]
        sample = raw_df.sample_id[row]
        rppm_dic2 = {}
        if stat == 'control':  # start: [id, end, rppm_Mean]
            for i,idx in rppm_mean_df.iterrows():
                    if type(rppm_mean_df[sample+"_genome_id"][i]) is not float:
                        rppm_dic1[rppm_mean_df[sample + "_start"][i]] = [rppm_mean_df[sample+"_genome_id"][i], rppm_mean_df[sample+"_end"][i], 
                                                                  rppm_mean_df[sample][i]] 
            common = rppm_dic1.copy()
        else:    
            inter_dic = {}
            for j,idx in rppm_mean_df.iterrows():
                    if type(rppm_mean_df[sample+"_genome_id"][j]) is not float:

                        rppm_dic2[rppm_mean_df[sample + "_start"][j]] = [rppm_mean_df[sample+"_genome_id"][j], rppm_mean_df[sample+"_end"][j],
                                                                     rppm_mean_df[sample][j]]          
            for x in common.keys():        #Find overlappin peaks
                id_1 = common[x][0]
                start_1 = int(x)
                end_1 = int(common[x][1])
                for y in rppm_dic2.keys():
                    id_2 = rppm_dic2[y][0]
                    start_2 = int(y)
                    end_2 = int(rppm_dic2[y][1])
                    if id_1 == id_2:
                        if start_1 <= start_2 and end_1 >= end_2:
                            inter_dic[start_1] = [start_2, id_1]
                        elif start_1 <= start_2 and end_1 >= start_2:
                            inter_dic[start_1] = [start_2, id_1]
                        elif start_1 <= end_2 and end_1 >= end_2:
                            inter_dic[start_1] = [start_2, id_1]
                        elif start_1 > start_2 and end_1 < end_2:
                            inter_dic[start_1] = [start_2, id_1]
                            
            del_list = []
            for i in common.keys():
                if i not in inter_dic.keys():
                    del_list.append(i)
            for j in del_list:
                del (common[j]) 
    
    num = 0
    sample_name = [] 
    for row,idx in raw_df.iterrows():
        stat = raw_df.stat[row]
        sample = raw_df.sample_id[row]
        sample_name.append(sample)
        rppm_dic2 = {}
        if stat != 'control':  # start: [id, end, rppm_mean]
            num += 1
            for j,idx in rppm_mean_df.iterrows():
                    if type(rppm_mean_df[sample+"_genome_id"][j]) is not float:
                        rppm_dic2[rppm_mean_df[sample + "_start"][j]] = [rppm_mean_df[sample+"_genome_id"][j], rppm_mean_df[sample+"_end"][j],
                                                                     rppm_mean_df[sample][j]]     
            con_rppm = {}
            exp_rppm = {}
            for x in sorted(common):
                id_1 = common[x][0]
                start_1 = int(x)
                end_1 = int(common[x][1])
                rppm_mean1 = float(common[x][2])
                con_rppm[start_1] = rppm_mean1
                for y in rppm_dic2.keys():
                    id_2 = rppm_dic2[y][0]
                    start_2 = int(y)
                    end_2 = int(rppm_dic2[y][1])
                    rppm_mean2 = float(rppm_dic2[y][2])
                    if id_1 == id_2:
                        if start_1 <= start_2 and end_1 >= end_2:
                            exp_rppm[start_1] = rppm_mean2
                        elif start_1 <= start_2 and end_1 >= start_2:
                            exp_rppm[start_1] = rppm_mean2
                        elif start_1 <= end_2 and end_1 >= end_2:
                            exp_rppm[start_1] = rppm_mean2
                        elif start_1 > start_2 and end_1 < end_2:
                            exp_rppm[start_1] = rppm_mean2
                        
            if num == 1:
                con_df = pd.DataFrame(list(con_rppm.values()))
                exp_df = pd.DataFrame(list(exp_rppm.values()))
                inter_df = pd.concat((con_df, exp_df), axis = 1)
            else:
                exp_df = pd.DataFrame(list(exp_rppm.values()))
                inter_df = pd.concat((inter_df, exp_df), axis = 1)

    inter_df.columns = sample_name
    return inter_df

def make_heatmap(raw_df, rppm_mean_dir, result_name, color = "Reds", size = (1,50), font_size = 4, fig_type = "svg"): #rppm_mean_dir = DataFrame or pickle files
    output_dir = raw_df.output[0]
    fig_result_dir = output_dir + "plots/"+ result_name

    if type(rppm_mean_dir) == pd.DataFrame:
        rppm_mean_df = rppm_mean_dir
    else:
        with open (rppm_mean_dir, "rb") as dataframe:
            rppm_mean_df = pickle.load(dataframe) 
    
    peak_lst = []
    inter_rppm = intersect_peaks(raw_df, rppm_mean_df)  
    print ("The number of common is " + str(len(inter_rppm)))

    for i in range (len(inter_rppm)):
        num = i + 1
        peak_num = "Peak_" + str(num)
        peak_lst.append(peak_num)

    peak_df = pd.DataFrame(peak_lst)
    peak_df.columns = ["common_peak"]
    final_df = pd.concat((peak_df, inter_rppm), axis = 1)
    diff_df = final_df.set_index('common_peak', drop = True)

    visuz.gene_exp.hmap(df=diff_df, figname = fig_result_dir, rowclus=False, colclus=False, r = 720, cmap = color, dim = size, tickfont=(font_size, font_size), figtype = fig_type, show =None)
    print ("Making a Heatmap is done :)")

def make_volcano(raw_df, converted_de_dir, result_name, test = "q-value", lfc_thr = 1.0, test_thr = 0.05, peak_names = None, lengend = True, fig_type = "svg"):
    output_dir = raw_df.output[0]
    fig_result_dir = output_dir + "plots/"+ result_name

    df = pd.read_csv(converted_de_dir, sep = ",")    
    
    df = df.replace({"p-value":0}, 1e-314)
    df = df.replace({"q-value":0}, 1e-314)

    if test == "p-value":
        test_name = "pval"
        axylabel_name = "log10(p-value)"
    elif test == "q-value":
        test_name = "q-value"
        axylabel_name = "-log10(p-value)"
    
    visuz.GeneExpression.volcano(df=df, geneid = 'peak_num', figname = fig_result_dir, lfc='LFC', pv=test_name, lfc_thr = (lfc_thr,lfc_thr), pv_thr = (test_thr, test_thr), 
                                 axylabel = axylabel_name, axxlabel = "log2(FoldChange)", genenames = peak_names, plotlegend = lengend, figtype = fig_type)    
    
    print ("Making a volcano plot is done :)")


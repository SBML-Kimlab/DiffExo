import sys
import pickle
import pandas as pd
import scipy.stats

#This script contains functions for statistical analysis of ChIP-exo or ChIP-mini intensity data.

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
            for x in common.keys():        #Find common peaks
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

def parametric_test(rppm_mean_df, stat_order):
    for x in stat_order:
        stat_tem = []
        if x =="control":
            control = rppm_mean_df[stat_order[x]]
            con_id = str(stat_order[x])
        else:
            sample = rppm_mean_df[stat_order[x]]
            sam_id = str(stat_order[x])
            stat_tem.append(control)
            stat_tem.append(sample)
            stat_command = tuple(stat_tem)
            leven_pval =  scipy.stats.levene(*stat_command).pvalue
            
            if leven_pval > 0.05:
                print (("%s and %s Leven's test p-value = %s")%(con_id, sam_id,str(leven_pval)))
                equal = True
            else:
                print (("%s and %s Leven's test p-value = %s")%(con_id, sam_id,str(leven_pval)))
                print ("Change T-test to Welch's T-test")
                equal = False

            ttest_pval = (scipy.stats.ttest_ind(*stat_command, equal_var= equal)).pvalue
            print (("%s and %s T-test p-value = %s")%(con_id, sam_id, str(ttest_pval)))
                                        
def nonparametric_test(rppm_mean_df, stat_order):                    
    for x in stat_order:
        stat_tem = []
        if x =="control":
            control = rppm_mean_df[stat_order[x]]
            con_id = str(stat_order[x])
        else:
            sample = rppm_mean_df[stat_order[x]].dropna()
            sam_id = str(stat_order[x])
            stat_tem.append(control)
            stat_tem.append(sample)
            stat_command = tuple(stat_tem)
            ranksum_pval = (scipy.stats.ranksums(*stat_command)).pvalue
            print (("%s and %s Wilcoxon rank sum test p-value = %s")%(con_id, sam_id, str(ranksum_pval)))
                                  
def stat_test(raw_df, rppm_mean_dir, test = "auto", replicates = False): #mean_dir = Dataframe or pickle file, test => parametric or nonparametric.
    if type(rppm_mean_dir) == pd.DataFrame:
        rppm_mean_df = rppm_mean_dir
    else:
        with open (rppm_mean_dir, "rb") as dataframe:
            rppm_mean_df = pickle.load(dataframe)        
    
    peak_num = []
    stat_order = {}
    for i,idx in raw_df.iterrows():
        peak_dir = raw_df.peak[i]
        peak = pd.read_csv(peak_dir, sep='\t', names = ["genome_id", "source", "feature", "start", "end","score", "strand", "frame", "attr"])
        peak_num.append(len(peak))
        stat_order[raw_df.stat[i]] = raw_df.sample_id[i]

    if replicates == True:
        if len(rppm_mean_df.columns) == 10:
            for name in stat_order.values():
                rep_1 = name + "_1"
                rep_2 = name + "_2"
                auto_df = rppm_mean_df[[rep_1, rep_2]]
                
                x = len(auto_df)
                shapiro_break = False
                shapiro_switch = False
                if x <= 30:
                    for num in range(len(auto_df.columns)):
                        stat_tem = []
                        num += 1

                        shapiro_name = name + "_" + str(num)
                        shapiro_pval = scipy.stats.shapiro(auto_df[shapiro_name].dropna()).pvalue
                        print (("Shapiro-wilk test p-value of %s replicates = %s")%(shapiro_name, str(shapiro_pval)))
                        if shapiro_pval <= 0.05:
                            shapiro_switch = True

                    if shapiro_switch == True:
                        rppm_1 = list(auto_df[rep_1].dropna()) 
                        rppm_2 = list(auto_df[rep_2].dropna())
                        stat_tem.append(rppm_1)
                        stat_tem.append(rppm_2)
                        stat_command = tuple(stat_tem)
                        ranksum_pval = (scipy.stats.ranksums(*stat_command)).pvalue
                        print (("Wilcoxon rank sum test p-value of %s replicates = %s")%(name, str(ranksum_pval)))
                        shapiro_break = True

                if shapiro_break == True:
                    continue

                stat_tem = []
                rppm_1 = list(auto_df[rep_1].dropna()) 
                rppm_2 = list(auto_df[rep_2].dropna())
                stat_tem.append(rppm_1)
                stat_tem.append(rppm_2)
                stat_command = tuple(stat_tem)
                leven_pval =  scipy.stats.levene(*stat_command).pvalue

                if leven_pval > 0.05:
                    print (("Leven's test p-value of %s replicates = %s")%(name, str(leven_pval)))
                    equal = True
                else:
                    print (("Leven's test p-value of %s replicates = %s")%(name, str(leven_pval)))
                    equal = False

                ttest_pval = (scipy.stats.ttest_ind(*stat_command, equal_var = equal)).pvalue
                print (("T-test p-value of %s replicates = %s")%(name, str(ttest_pval)))
            
        else:
            sys.exit('Change input pickle file to [Test_name]_rppm.pickle')
    else:        
        if len(rppm_mean_df.columns) == 8:
            del_cols = []
            for i,idx in raw_df.iterrows():
                sample_id = raw_df.sample_id[i]
                del_cols.extend([sample_id+"_genome_id", sample_id+"_start", sample_id+"_end"])         
            auto_df = rppm_mean_df.drop(del_cols,axis = 1)  

            if test == "auto":
                shapiro_break = False
                for x in peak_num:
                    if x <= 30:
                        for sample_id in stat_order.values():
                            shapiro_pval = scipy.stats.shapiro(auto_df[sample_id]).pvalue
                            if shapiro_pval <= 0.05:
                                nonparametric_test(auto_df, stat_order)  #Nonparametric Test
                                shapiro_break = True
                                break
                        if shapiro_break == True:
                            break
                    else:
                        if len(set(peak_num)) == 1:  #Parametric Test
                            parametric_test(auto_df, stat_order)
                            break
                        else:
                            nonparametric_test(auto_df, stat_order) #Nonparametric Test
                            break
            
            elif test == "nonparametric":
                nonparametric_test(auto_df, stat_order)
            
            elif test == "parametric":
                intersect_df = intersect_peaks(raw_df,rppm_mean_df)
                print ("The number of overlapping peaks is " + str(len(intersect_df)))
                parametric_test(intersect_df, stat_order)
                
                #tem = raw_df.output[0] + raw_df.sample_id[0]
                #inter_result_dir = tem + "_intersection.pickle"  
                #with open (inter_result_dir, "wb") as inter:
                    #pickle.dump(intersect_df, inter)
        else:
            sys.exit('Change input pickle file to [Test_name]_mean_rppm.pickle')
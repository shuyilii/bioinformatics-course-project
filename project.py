#!/usr/bin/python3
import pandas as pd
from scipy import stats
from scipy.stats import levene
import sys

with open('GSE20292_series_matrix.txt', 'r') as fin:
    # To link the Sample_ID and the Sample_information with a dictionary
    Sample_information_dict = {}
    for line in fin:
        line = line.rstrip()
        if '!Sample_title' in line:
            Sample_information_list = line.replace('"','').split('\t')[1:]
        if '!Sample_geo_accession' in line:
            Sample_ID = line.replace('"','').split('\t')[1:]
Sample_information_dict = dict(zip(Sample_ID, Sample_information_list))

#To make dictionaries of control sample and disease sample
control_sample_dict = {}
disease_sample_dict = {}
for key in Sample_information_dict:
    if 'P' in Sample_information_dict[key]:
        disease_sample_dict[key] = Sample_information_dict[key]
    if 'C' in Sample_information_dict[key]:
        control_sample_dict[key] = Sample_information_dict[key]
        
#Use read_csv to turn the raw data into dataframe
data = pd.read_csv('new_file.txt', sep ='\t', header = 0, index_col = 0)
#Create the control dataframe
control_IDs = list(control_sample_dict.keys())
control_data = data[control_IDs]
#Create the disease dataframe
disease_IDs = list(disease_sample_dict.keys())
disease_data = data[disease_IDs]

#Do t-test between control sample and disease sample of each probe ID
results_dict = {}
#Get probe ID list
Geneid_list = list(control_data.index)
#Turn the expression value of each probe ID into a list
for i in range(0,len(Geneid_list)):
    control_array_for_each_ID = control_data.ix[i].values.tolist()
    disease_array_for_each_ID = disease_data.ix[i].values.tolist()
    #Do homogenetiy test
    homogeneity_test = levene(control_array_for_each_ID, disease_array_for_each_ID)
    #Do t-test and stored the results into a dictionary
    if homogeneity_test[1] > 0.05:
        tp_val = stats.ttest_ind(control_array_for_each_ID, disease_array_for_each_ID)
    else:
        tp_val = stats.ttest_ind(control_array_for_each_ID, disease_array_for_each_ID, equal_var = False)
    results_dict[Geneid_list[i]] = tp_val
#Sort the results with descending order of p_value
sorted_results_dictkey = sorted(results_dict, key=lambda k:results_dict[k][1])

#Write the results to a file
with open('Results_0.txt', 'w') as fout:
    for key in sorted_results_dictkey:
        print(key, '\t', results_dict[key], file = fout)

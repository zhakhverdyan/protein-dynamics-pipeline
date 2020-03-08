#!/usr/bin/python
import numpy as np

# this function accepts an array of peptides with corresponding MS/MS count and
# H/(H+L) values from 1 protein, filters out outliers, calculates
# the mean and standard deviation, finds the number of unique peptides and spectral counts 
# observed for each protein in each time point

def stat_analysis(peptide_array):
    peptides=[];
    msms_counts=0;
    ratios=[];
    for line in peptide_array:
        line.rstrip('\n');
        values = line.split('\t');
        peptides.append(values[0]);
        msms_counts=msms_counts+int(values[1])
        ratios.append(float(values[2]));
    num_unique_peptides=len(set(peptides))
    if num_unique_peptides>5:
        stats=find_outlier(ratios);
    elif num_unique_peptides==1:
        stats=(ratios[0], 0, 0, 0);
    else:
        stats=mean_sterror(ratios);
    return num_unique_peptides, msms_counts,  stats
        
          
# this part calculates the first and third quartile for interquartile range calculation
# to identify outliers (< Q1-1.5*IQR or >Q3+1.5*IQR). The function accepts an input of a list
# of H/(H+L) ratios. After the outlier elimination (only 1 peptide measurement per protein)
# the function calls onto mean_sterror function       
        
def find_outlier(ratios):

    p25=np.percentile(ratios,25);
    p75=np.percentile(ratios,75);
    iqr=p75-p25;
    lower_limit=p25-1.5*iqr;
    upper_limit=p75+1.5*iqr;
    if min(ratios) < lower_limit:
        ratios.remove(min(ratios));
        stats=mean_sterror(ratios);
    elif max(ratios) > upper_limit:
        ratios.remove(max(ratios));
        stats=mean_sterror(ratios);
    else:
        stats=mean_sterror(ratios);
    return stats

# this part calculates the mean and standard deviation of protein H/(H+L) ratio based on its constituent peptides 
    
def mean_sterror(ratios):
    mn=round(np.average(ratios), 5);
    std=round(np.std(ratios), 5);
    return mn, std
    
    
    
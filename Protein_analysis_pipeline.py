# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 11:34:11 2017

@author: zhannahakhverdyan
"""

#!/usr/bin/python
"""This script accepts the MaxQuant evidence.txt file(s). Make sure that sequential RAW files going 
into MaxQuant analysis are sequentially (alphanumerically) named. After initial filtering for 
contaminants and reverse matches, the ratio of H/(H+L) is calulated for each peptide identification 
based on the H/L entry. Next, the data is grouped by the RAW file, where each peptide was observed. 
Next, in the data corresponding to each RAW file, peptides are grouped into the 'leading razor protein'. 
In addition, Stat_analysis_final.py is used to calculate various statistics (average and standard 
deviation of the H/(H+L) value) and pull out the number of unique peptides, which signifies how 
reliable each identification was (for each protein in each time point). Next, there is a log-transformation
of average H/(H+L) values. A linear regression is performed on the time course data wherever there are >=3 
data points available. If there are 4 or 5 data points available and the initial fit does not have a good r 
cotrrelation coefficient (i.e. r<0.9), the data points are sequenttialy removed (only one at a time) and the 
regression performed on a subset of the original dataset to improve r. The slope and intercept of the linear 
fit are recorded, as well as the r correlation coefficient, P value (two-sided p-value for a hypothesis test 
whose null hypothesis is that the slope is zero) and standard error of the estimated gradient. Finally, the 
systematic names of proteins are matched to the common names in the SGD database and, based on gene descriptions, 
the proteins are split into 3 goups: nucleoporin, transport facotrs and other. The three groups are plotted in 
different colors for visual inspection. All the files are organized in new folders: MaxQuant, Processed, 
Processed/Assigned, Processed/Excluded, Processed/Graphs, Processed/Intermediate."""
from datetime import datetime
startTime = datetime.now()
import sys
import Stat_analysis_final_log
import os
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
if not os.path.exists('Processed'):
    os.makedirs('Processed')

f = open(sys.argv[1], 'r');
# open the evidence.txt file in read only mode

file_name = sys.argv[1];
file_name = file_name[:-4];
# create a variable name from the file name

filt_name = file_name + '_filtered_out';
analyz_name = file_name + '_analyzed';
col_data = file_name+'_collated';
# create names for various files were the calculated, analyzed or excluded data will
# be stored, which is based on the original file name

filt = open('Processed/'+filt_name+'.txt', 'w');
analyzed = open('Processed/'+analyz_name+'.txt', 'w');
collated = open('Processed/'+col_data+'.txt', 'w');
# now create writable text files with all of the names

first_line = f.readline();
# gets the first line, column headers

first_line=first_line.rstrip('\r\n');

filt.write(first_line);
filt.write('\n');
analyzed.write(first_line);
analyzed.write('\tRatio H/(H+L)\n')
# the data that cannot be used (e.g. contaminant, reverse or no K peptides)
# will be simply copied to _filtered_out.txt file to keep track, while
# all the data that will go into analysis will be written into _analyzed.txt file
# with an aditional column Ratio H/(H+L) that will be calculated from Ratio H/L

# make all lower letters
first_line=first_line.lower();
headers = first_line.split('\t');
#this is for case insensitive matching

# get the column indices of interest
i0=headers.index("sequence");
i1=headers.index("k count");
i2=headers.index("leading razor protein");
i3=headers.index("raw file");
i4=headers.index("ms/ms count");
i5=headers.index("ratio h/l");
i6=headers.index("reverse");
i7=headers.index("potential contaminant");
#i7=headers.index("contaminant"); depending on the version of MaxQuant, this header may be used

"""
this part of the program writes all the unusable data to a _filtered_out.txt 
file, all the usable data into _analysed.txt file, and keeps track of all the different raw files.
Also we create a dictionary, where the keys are raw files and the values are arrays containing all 
the peptide characteristics to be analyzed from the corresponding raw file.
the dictionary is printed out based on keys into separate txt files which have 4 column entries - "Sequence", 
"Leading Razor Protein", "MS/MS Count", "log(Ratio H/total)"""

raw_files = {};
for lines in f:
    entries=lines.split('\t');
    if entries[i1] == '0':
        filt.write('\t'.join(entries));
    elif entries[i4] == '0':
        filt.write('\t'.join(entries));
    elif entries[i5] == '' or entries[i5] == "NaN":
        filt.write('\t'.join(entries));
    elif entries[i6] == "+":
        filt.write('\t'.join(entries));
    elif entries[i7] == "+":
        filt.write('\t'.join(entries));
    else:
        line = '\t'.join(entries);
        analyzed.write(line.rstrip('\r\n'));
        analyzed.write('\t');
        ratio = float(entries[i5])/(float(entries[i5])+float(1));
        if ratio > 0:
            ratio = np.log(ratio)
            ratio=round(ratio, 5)
            analyzed.write(str(ratio));
            analyzed.write('\n');
            line = (entries[i0]+'\t'+entries[i2]+'\t'+entries[i4]+'\t'+str(ratio)+'\n');
            if entries[i3] not in raw_files:
                raw_files[entries[i3]] = [line];
            else:
                raw_files[entries[i3]].append(line);
        else:
            filt.write('\t'.join(entries));
f.close();
filt.close();
analyzed.close();

# this part writes the contents from the corresponding RAW files into separate text files

raw_file_header = ("Sequence"+'\t'+"Leading Razor Protein"+'\t'+"MS/MS Count"+'\t'+"ln(H/(H+L))"+'\n')

for files in sorted(raw_files):
    f = open('Processed/'+file_name+files+'.txt', 'w');
    f.write(raw_file_header)
    for line in raw_files[files]:
        f.write(line);
    f.close();
# Next, the text files are read and organised by proteins. A dictionary is created
# for each file, the keys are the leading razor proteins,
# the values are arrays with corresponding peptides and characteristics

protein_analysis_headers = ("Protein systematic name"+'\t'+"Number of unique peptides"+'\t'+"Total MS/MS counts"+'\t'+"Average H/(H+L)"+'\t'+"Standard deviation"+'\n')

for files in sorted(raw_files):
    f = open('Processed/'+file_name+files+'.txt', 'r');
    f1 = open('Processed/'+file_name+files+'_processed.txt', 'w');
    f1.write(protein_analysis_headers);
    first_line = f.readline();
    headers = first_line.split('\t');
    i0=headers.index("Sequence");
    i2=headers.index("Leading Razor Protein");
    i4=headers.index("MS/MS Count");
    i8=headers.index("ln(H/(H+L))\n");
    proteins={};
    for lines in f:
        entries=lines.split('\t');
        line = (entries[i0]+'\t'+entries[i4]+'\t'+entries[i8]+'\n');
        if entries[i2] not in proteins:
            proteins[entries[i2]]=[line];
        else:
            proteins[entries[i2]].append(line);
    for protein in proteins:
        l=Stat_analysis_final_log.stat_analysis(proteins[protein]);
        f1.write(protein+'\t');
        f1.write(str(l[0])+'\t'+str(l[1])+'\t'+str(l[2][0])+'\t'+str(l[2][1])+'\n')
    f1.close()
    
# the _final_analysis text file will contain the time course data for each protein that will be fit to an exponential curve
f2=open('Processed/'+file_name+'_final_analysis.txt', 'w');
# create an empty list of protein names to be populated with protein names from all files
protein_names=[];    

# create an empty list of dictioneries, that are going to contain protein=key, 
# the rest of descriptors=value pairs
file_names=[];

for files in sorted(raw_files):
    f3 = open('Processed/'+file_name+files+'_processed.txt', 'r');
    f3.readline(); # get rid of the header row
    # make an empty dictionary named "filename" and add it to the file_name list
    file_name1={};
    # get rid of the new line character, split the line on tabs, add the protein name 
    #to the list of proteins, create a key-value pair from the protein name and the rest 
    #of the characteristics
    for lines in f3:
        lines=lines.rstrip('\n');
        entries=lines.split('\t');
        protein_names.append(entries[0])
        file_name1[entries[0]]=('\t'.join(entries[1:]));
    f3.close();
    file_names.append(file_name1);
# go through the list of protein names and write the protein names and values from each file
# next to each other
num_files=len(file_names); # this counts the number of time points in the experiment
protein_analysis_header = ("Protein systematic name"+('\t'+"Number of unique peptides"+'\t'+"Total MS/MS counts"+'\t'+"Average ln(H/(H+L))"+'\t'+"Standard deviation ln(H/H+L)")*len(raw_files)+'\n')
f2.write(protein_analysis_header);

# this part of the code pulls out the protein measuremnt characteristics from each file
# and prints them next to each other in the _final_analysis.txt file for linear regression.
# the order of values is based on string sorted file_name, that is why sequential time point samples 
# need to be named in an alphanumerical ascending sequeantial order
unique_protein_names=set(protein_names);
four_tabs=('\t0'*4);
for proteins in unique_protein_names:
    f2.write(proteins);
    sorted(file_names, key=str);
    for file_name1 in file_names:
        if proteins not in file_name1:
            f2.write(four_tabs);
            continue;
        values=file_name1[proteins];
        val_entries = values.split('\t');
        f2.write('\t'+values);
    f2.write('\n');
f2.close();

# this part deletes, unused text files and and can be optionally removed in case the individual samples need to be separated
for files in sorted(raw_files):
    os.remove('Processed/'+file_name+files+'.txt')
    os.remove('Processed/'+file_name+files+'_processed.txt')
    
"""
This part performs a linear regression on the log-transformed average H/(H+L) measurements over the timecourse. 
I do not enforce a (0,log(1)) intercept, since there is incomplete labeling and excess tagged protein in the  
first hour of the time course. 
"""

f = open('Processed/'+file_name+'_final_analysis.txt', 'r');
f1 = open('Processed/'+file_name+'_fitted.txt', 'w');
f2 = open('Processed/'+file_name+'_notfitted.txt', 'w');
f3 = open('Processed/'+file_name+'_poorfit.txt', 'w');
f4 = open('Processed/'+file_name+'_lowconfidence.txt','w');    
average_timepoint='';
stdev_timepoint='';
for i in range(num_files):
    m=str(i+1);
    average_timepoint+='\tAverage '+m;
    stdev_timepoint+='\tStandard deviation '+m;
f1.write("Protein systematic name"+average_timepoint+stdev_timepoint+"\tSlope\tIntercept\tr value\tp value\tstd err\tAve uniqiue peptides\tAve MS/MS count\tTime point filtered\n");
f2.write("Protein systematic name\tTime points\tAve unique peptide\tAve MS/MS count\n");
f3.write("Protein systematic name"+average_timepoint+stdev_timepoint+"\tSlope\tIntercept\tr value\tp value\tstd err\tAve uniqiue peptides\tAve MS/MS count\tTime point filtered\n");
f4.write("Protein systematic name\tCommon name\tDescription"+average_timepoint+stdev_timepoint+"\tSlope\tIntercept\tr value\tp value\tstd err\tAve uniqiue peptides\tAve MS/MS count\tTime point filtered\n"); 
t1=np.array([1,2,3,4,5]);
t2=np.array([1.5,3,4.5,6,7.5]); # Nup120-null strain grew slower and hence the time cours is longer for this strain

f.readline();
for lines in f:
    lines=lines.rstrip('\r\n');
    entries=lines.split('\t');
    stats1 = np.array(entries[1:len(entries)], dtype=float);
    unique=stats1[0:len(stats1):4];
    counts=stats1[1:len(stats1):4];
    ave_unique = np.average(unique, axis=None);
    ave_counts = np.average(counts, axis=None);
    average1=stats1[2:len(stats1):4];
    stdev1=stats1[3:len(stats1):4];
    average1_filt=stats1[2:len(stats1):4];
    average1_filt=average1_filt.ravel()[np.flatnonzero(average1_filt)];
    average1=average1.tolist();
    for n,i in enumerate(average1):
        if i==0.:
            average1[n]='nan';
            stdev1[n]='nan';
    if re.search("Nup120-null",file_name):   
        t=t2;
    else:
        t=t1;
    if average1_filt.shape[0]>=3: 
        unique=unique.ravel()[np.flatnonzero(average1_filt)]
        counts=counts.ravel()[np.flatnonzero(average1_filt)]
        tmod=t.ravel()[np.flatnonzero(average1_filt)]
        slope, intercept, r_value, p_value, std_err = stats.linregress(tmod,average1_filt)
        """at least 3 points are needed for the regression, proteins with less data points are added to the not-fitted list.
        Fits with r correlation coefficient<0.9 and p value >0.05 are either discarded (added to the poor fit list) or
        further manipulated if there are 4 or 5 points. Linear regression is performed on the subset of the data and if r and p
        values do not improve, then the protein is added to the poor fit list otherwise it is recorded in the assigned list. 
        Proteins with >=3 time points, r correlation coefficient >=0.9 and p value <=0.05 are added to the assigned list."""
        if -r_value >= 0.9 and p_value<=0.05:
            f1.write(entries[0]+'\t'+'\t'.join(map(str, average1))+'\t'+'\t'.join(map(str, stdev1))+'\t'+str(round(slope,5))+'\t'+str(round(intercept,5))+'\t'+str(-round(r_value,5))+'\t'+str(round(p_value,5))+'\t'+str(round(std_err,5))+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\t'+'\n');
        else:
            if (average1_filt.shape[0]==3):
                f3.write(entries[0]+'\t'+'\t'.join(map(str, average1))+'\t'+'\t'.join(map(str, stdev1))+'\t'+str(round(slope,5))+'\t'+str(round(intercept,5))+'\t'+str(-round(r_value,5))+'\t'+str(round(p_value,5))+'\t'+str(round(std_err,5))+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\n');
            elif (average1_filt.shape[0]>3):
                j=0;
                average3=[];
                stdev3=[];
                for i in range(len(average1_filt)):
                    average2=np.delete(average1_filt, i)
                    tmod_new=np.delete(tmod,i)
                    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(tmod_new,average2)
                    if -r_value1 >= 0.9 and p_value1 <= 0.05:
                        j=i;
                        slope, intercept, r_value, p_value, std_err= slope1, intercept1, r_value1, p_value1, std_err1;
                if -r_value >= 0.9 and p_value <= 0.05: 
                    average3=np.delete(average1, j); # if the exclusion of the value improved the fit, the value is discarded
                    average3=average3.tolist();
                    average3.insert(j, 'nan');
                    stdev2=np.delete(stdev1,j);
                    stdev3=stdev2.tolist();
                    stdev3.insert(j,'nan');
                    f1.write(entries[0]+'\t'+'\t'.join(map(str,average3))+'\t'+'\t'.join(map(str,stdev3))+'\t'+str(round(slope,5))+'\t'+str(round(intercept,5))+'\t'+str(-round(r_value,5))+'\t'+str(round(p_value,5))+'\t'+str(round(std_err,5))+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\tYes'+'\n');
                else:
                    f3.write(entries[0]+'\t'+'\t'.join(map(str, average1))+'\t'+'\t'.join(map(str, stdev1))+'\t'+str(round(slope,5))+'\t'+str(round(intercept,5))+'\t'+str(-round(r_value,5))+'\t'+str(round(p_value,5))+'\t'+str(round(std_err,5))+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\t'+'\n');
    else:
        f2.write(entries[0]+'\t'+str(average1_filt.shape[0])+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\n')
f.close();
f1.close();
f2.close();  
f3.close();      
"""
This part searches the SGD database for the systematic names of proteins and groups 
the final set of proteins with fitted values into 3 groups: nups, transport factors and other.
This is done  for assigned, not fitted and poor fit files, separately.
"""

f = open('Processed/'+file_name+'_fitted.txt', 'r');
f1 = open('Processed/'+file_name+"_nups_grouped.txt", 'w');
f2 = open('Processed/'+file_name+"_kaps_grouped.txt", 'w');
f3 = open('Processed/'+file_name+"_other_grouped.txt", 'w');
fp = open('Processed/'+file_name+'_poorfit.txt', 'r')
fp1 = open('Processed/'+file_name+"_nups_poorfit.txt", 'w');
fp2 = open('Processed/'+file_name+"_kaps_poorfit.txt", 'w');
fp3 = open('Processed/'+file_name+"_other_poorfit.txt", 'w');
fn = open('Processed/'+file_name+'_notfitted.txt', 'r');
fn1 = open('Processed/'+file_name+"_nups_notfitted.txt", 'w');
fn2 = open('Processed/'+file_name+"_kaps_notfitted.txt", 'w');
fn3 = open('Processed/'+file_name+"_other_notfitted.txt", 'w');

protein_desc={};
prot=open("../../../orf_trans_all_2015_Nup145.txt", 'r');
for line in prot:
    if re.match('^>', line):
        line=line.rstrip();
        line=line[1:];
        entries=line.split(' ');
        protein_description=(' '.join(entries[2:]));
        protein_desc[entries[0]]=(entries[1]+'\t'+protein_description);
    else:
        next;
              
header=f.readline();
names= header.split('\t');
i9=names.index("Ave uniqiue peptides");
new_header=(names[0]+'\t'+"Common name"+'\t'+"Description"+'\t'+('\t'.join(names[1:])));
f1.write(new_header);
f2.write(new_header);
f3.write(new_header);

for line in f:
    entries=line.split('\t');
    protein_stats=('\t'.join(entries[1:]));
    first_name=entries[0].split(';');
    if re.search('nuclear pore complex', protein_desc[first_name[0]], flags=re.IGNORECASE):
        if first_name[0]=='YLR347C' or first_name[0]=='YMR308C' or first_name[0]=='YHR170W' or first_name[0]=='YDR101C':
            if float(entries[i9])>=2:
                f2.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
            elif float(entries[i9])<2:
                f4.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        elif first_name[0]=='YMR290C' or first_name[0]=='YML107C':
            if float(entries[i9])>=4:
                f3.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
            elif float(entries[i9])<4:
                f4.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        else:
            f1.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
    elif re.search('karyopherin|transport factor|exportin|importin|mRNA export|nuclear import|nucleocytoplasmic transport', protein_desc[first_name[0]], flags=re.IGNORECASE):
        if first_name[0]=='YDL075W' or first_name[0]=='YLR406C' or first_name[0]=='YGL086W' or first_name[0]=='YDL084W' or first_name[0]=='YDL160C':
            if float(entries[i9])>=4:
                f3.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
            elif float(entries[i9])<4:
                f4.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        else:
            if float(entries[i9])>=2:
                f2.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
            elif float(entries[i9])<2:
                f4.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
    elif first_name[0]=='YDR424C' or first_name[0]=='HA_GFP':
        f1.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);     
    else:
        if float(entries[i9])>=4:
                f3.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        elif float(entries[i9])<4:
                f4.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);


f.close();
f1.close();
f2.close();
f3.close();
f4.close();

headerp=fp.readline();
namesp= headerp.split('\t');
new_headerp=(namesp[0]+'\t'+"Common name"+'\t'+"Description"+'\t'+('\t'.join(namesp[1:])));
fp1.write(new_headerp);
fp2.write(new_headerp);
fp3.write(new_headerp);

for line in fp:
    entries=line.split('\t');
    protein_stats=('\t'.join(entries[1:]));
    first_name=entries[0].split(';');
    if re.search('nuclear pore complex', protein_desc[first_name[0]], flags=re.IGNORECASE):
        if first_name[0]=='YLR347C' or first_name[0]=='YMR308C' or first_name[0]=='YHR170W' or first_name[0]=='YDR101C':
            fp2.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        elif first_name[0]=='YMR290C' or first_name[0]=='YML107C':
            fp3.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        else:
            fp1.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
    elif re.search('karyopherin|transport factor|exportin|importin|mRNA export|nuclear import|nucleocytoplasmic transport', protein_desc[first_name[0]], flags=re.IGNORECASE):
        if first_name[0]=='YDL075W' or first_name[0]=='YLR406C'  or first_name[0]=='YGL086W' or first_name[0]=='YDL084W' or first_name[0]=='YDL160C':
            fp3.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        else:
            fp2.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
    elif first_name[0]=='YDR424C' or first_name[0]=='HA_GFP':
        fp1.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);     
    else:
        fp3.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        
fp.close();
fp1.close();
fp2.close();
fp3.close();

headern=fn.readline();
namesn= headern.split('\t');
new_headern=(namesp[0]+'\t'+"Common name"+'\t'+"Description"+'\t'+('\t'.join(namesn[1:])));
fn1.write(new_headern);
fn2.write(new_headern);
fn3.write(new_headern);

for line in fn:
    entries=line.split('\t');
    protein_stats=('\t'.join(entries[1:]));
    first_name=entries[0].split(';');
    if re.search('nuclear pore complex', protein_desc[first_name[0]], flags=re.IGNORECASE):
        if first_name[0]=='YLR347C' or first_name[0]=='YMR308C' or first_name[0]=='YHR170W' or first_name[0]=='YDR101C':
            fn2.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        elif first_name[0]=='YMR290C' or first_name[0]=='YML107C'  or first_name[0]=='YGL086W' or first_name[0]=='YDL084W' or first_name[0]=='YDL160C':
            fn3.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        else:
            fn1.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
    elif re.search('karyopherin|transport factor|exportin|importin|mRNA export|nuclear import|nucleocytoplasmic transport', protein_desc[first_name[0]], flags=re.IGNORECASE):
        if first_name[0]=='YDL075W' or first_name[0]=='YLR406C':
            fn3.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        else:
            fn2.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
    elif first_name[0]=='YDR424C' or first_name[0]=='HA_GFP':
        fn1.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);     
    else:
        fn3.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats);
        

fn.close();
fn1.close();
fn2.close();
fn3.close();
prot.close();

os.remove('Processed/'+file_name+'_fitted.txt');
os.remove('Processed/'+file_name+'_poorfit.txt');
os.remove('Processed/'+file_name+'_notfitted.txt');

"""
finally the data and fits in the assigned (fitted) file are plotted for visual inspection
"""

f1=open('Processed/'+file_name+"_nups_grouped.txt", 'r');
f2=open('Processed/'+file_name+"_kaps_grouped.txt", 'r');
f3=open('Processed/'+file_name+"_other_grouped.txt", 'r');



t3=np.arange(0, t[-1]+0.5, 0.0025);

first_line=f1.readline();
entries=first_line.split('\t');
i10=entries.index("Slope");
i11=entries.index("Intercept");
f2.readline();
f3.readline();
lw=0.75;
timepoints_kap=[];
if num_files==4:
    t=t[:-1];
error_kap = [];
for line in f2:
    line=line.rstrip('\r\n');
    entries=line.split('\t');
    param2=float(entries[i10]);
    y=param2*t3;
    plt.plot(t3,y,'mediumpurple', label='Transport factors', lw=lw);
    t1=[];
    if num_files==4:
        timepoints_kap=entries[3:7];
        for n,i in enumerate(timepoints_kap):
            timepoints_kap[n]=float(timepoints_kap[n]);
        timepoints_kap[:] = [x - float(entries[i11]) for x in timepoints_kap];
        error_kap = entries[7:11];
        for n,i in enumerate(error_kap):
            error_kap[n]=float(error_kap[n]);
        plt.errorbar(t, timepoints_kap, yerr=error_kap, fmt='o', color='lightgrey');
    elif num_files==5:
        timepoints_kap=entries[3:8];
        for n,i in enumerate(timepoints_kap):
            timepoints_kap[n]=float(timepoints_kap[n]);
        timepoints_kap[:] = [x - float(entries[i11]) for x in timepoints_kap];
        error_kap = entries[8:13];
        for n,i in enumerate(error_kap):
            error_kap[n]=float(error_kap[n]);
        plt.errorbar(t, timepoints_kap, yerr=error_kap, fmt='o', color='lightgrey');
        
timepoints_nup=[];
error_nup = [];     
for line in f1:
    line=line.rstrip('\r\n');
    entries=line.split('\t');
    param2=float(entries[i10]);
    y=param2*t3;
    plt.plot(t3,y,'skyblue', label='Nucleoporins', lw=lw);
    if num_files==4:
        timepoints_nup=entries[3:7];
        for n,i in enumerate(timepoints_nup):
            timepoints_nup[n]=float(timepoints_nup[n]);
        timepoints_nup[:] = [x - float(entries[i11]) for x in timepoints_nup];
        error_nup = entries[7:11];
        for n,i in enumerate(error_nup):
            error_nup[n]=float(error_nup[n]);
        plt.errorbar(t, timepoints_nup, yerr=error_nup, fmt='o', color='lightgrey');
    elif num_files==5:
        timepoints_nup=entries[3:8];
        for n,i in enumerate(timepoints_nup):
            timepoints_nup[n]=float(timepoints_nup[n]);
        timepoints_nup[:] = [x - float(entries[i11]) for x in timepoints_nup];
        error_nup = entries[8:13];
        for n,i in enumerate(error_nup):
            error_nup[n]=float(error_nup[n]);
        plt.errorbar(t, timepoints_nup, yerr=error_nup, fmt='o', color='lightgrey');

# to have a basis for comparison with "other" proteins in the pull out, largely abundant housekeeping
#proteins, the average +/-2 standard deviation of the observed population is plotted as well
       
params=np.array([]);    
for line in f3:
    line=line.rstrip('\r\n');
    entries=line.split('\t');
    param2=float(entries[i10]);
    params=np.append(params, param2);
ave_param=np.mean(params);
std_param=np.std(params);
y=ave_param*t3;
plt.plot(t3,y,'black', label='Average', lw=1);
ave_minus=ave_param-2*std_param;
y1=ave_minus*t3;
plt.plot(t3,y1,'black', label='Average +/- 2 st dev', lw=1, linestyle='dashed');
ave_plus=ave_param+2*std_param;
y2=ave_plus*t3;
plt.plot(t3,y2,'black', lw=1, linestyle='dashed');

# following, fast turnover proteins from Christiano et al, 2014 are plotted for comparison
fast_turn=[-5.07237, -4.14418, -4.01516, -3.57477, -3.47162, -3.44578, -3.17846, -2.99037, -2.91265, -2.84860];
y_fast=[];
for item in fast_turn:
    y_fast=item*t3;
    plt.plot(t3,y_fast, color='tomato', linestyle=':', label='Fast turnover proteins', lw=1.5)

plt.xlabel('Time (h)', fontsize=20);
plt.ylabel('Log(H/(H+L))',fontsize=20);
plt.title(file_name+' turnover\n', fontsize=24); # replace ' turnover\n' with ' exchange\n' for exchange experiment
plt.xticks(np.arange(0, max(t)+1, 1), fontsize=18);
plt.xlim(0,t[-1]+0.5);
plt.yticks(np.arange(-5,0,1), fontsize=18);
plt.ylim(-5,0);

#ax=plt.gcf();
#ax.axis["right"].set_visible(False)
#ax.axis["top"].set_visible(False)

plt.box('off');
#plt.axhline(y=0);
#plt.axvline(x=0);
plt.tick_params(
    axis='x',          
    which='both',      
    bottom='off',     
    top='off') ;
plt.tick_params(
    axis='y',          
    which='both',      
    left='off',     
    right='off') ; 
#ax = plt.axes()       
#ax=pylab.axis()
#print(ax)
#ax.axis["right"].set_visible(False)
#ax.axis["top"].set_visible(False)

# Hide the right and top spines
#ax=plt.figure(); 
#plt.spines['right'].set_visible(False)
#plt.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
#plt.yaxis.set_ticks_position('left')
#plt.xaxis.set_ticks_position('bottom')


# include 1 legend for all graphs.
handles, labels = plt.gca().get_legend_handles_labels()
newLabels, newHandles = [], []
for handle, label in zip(handles, labels):
  if label not in newLabels:
    newLabels.append(label)
    newHandles.append(handle)
plt.legend(newHandles, newLabels, loc=0, frameon=False)

plt.savefig('Processed/'+file_name+'graph.eps' , bbox_inches='tight', format='eps', dpi=1000)
f1.close();
f2.close();
f3.close();  

# make some new folders to organize all the files
if not os.path.exists('MaxQuant'):
    os.makedirs('MaxQuant')
if not os.path.exists('Processed/Graphs'):
    os.makedirs('Processed/Graphs')
if not os.path.exists('Processed/Excluded'):
    os.makedirs('Processed/Excluded')
if not os.path.exists('Processed/Intermediate'):
    os.makedirs('Processed/Intermediate')
if not os.path.exists('Processed/Assigned'):    
    os.makedirs('Processed/Assigned')
if not os.path.exists('Processed/Excluded/Notfitted'):
    os.makedirs('Processed/Excluded/Notfitted')
if not os.path.exists('Processed/Excluded/Poorfit'):
    os.makedirs('Processed/Excluded/Poorfit')
if not os.path.exists('Processed/Excluded/Lowconfidence'):
    os.makedirs('Processed/Excluded/Lowconfidence')

src = "Processed"
dst1 = "Processed/Excluded"
dst2 = "Processed/Graphs"
dst3 = "Processed/Intermediate"
dst4 = "Processed/Assigned"
dst5 = "Processed/Excluded/Notfitted"
dst6 = "Processed/Excluded/Poorfit"
dst7 = "MaxQuant" 
dst8 = "Processed/Excluded/Lowconfidence"

listOfFiles = os.listdir(src)
for f in listOfFiles:
    if re.search("_filtered",f):
        fullPath = src + "/" + f
        os.system ("mv"+ " " + fullPath + " " + dst1)
    elif re.search("graph",f):
        fullPath = src + "/" + f
        os.system ("mv"+ " " + fullPath + " " + dst2)
    elif re.search("_analyzed|_analysis",f):
        fullPath = src + "/" + f
        os.system ("mv"+ " " + fullPath + " " + dst3)    
    elif re.search("grouped",f):
        fullPath = src + "/" + f
        os.system ("mv"+ " " + fullPath + " " + dst4)  
    elif re.search("notfitted",f):
        fullPath = src + "/" + f
        os.system ("mv"+ " " + fullPath + " " + dst5)
    elif re.search("poorfit",f):
        fullPath = src + "/" + f
        os.system ("mv"+ " " + fullPath + " " + dst6)
    elif re.search("lowconfidence",f):
        fullPath = src + "/" + f
        os.system ("mv"+ " " + fullPath + " " + dst8)

fullPath = "./" + file_name + ".txt"
os.system ("mv"+ " " + fullPath + " " + dst7)

print(datetime.now() - startTime)  
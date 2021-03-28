#!/Users/zhannahakhverdyan/opt/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun March 14 11:34:11 2021

@author: zhannahakhverdyan
"""

"""
This script accepts two command line arguments: 1. MaxQuant evidence.txt file 2. a json file, a dictionary with 2 keys: the 1st key is 'time' 
which corresponds to a list with time course values in hours, e.g.[1, 2, 3, 4, 5] for a five-hour time course where samples were collected hourly; 
2nd key is 'keywords' - a list of keywords to pull the proteins of interest, e.g. ['nup', 'nucleoporin', nuclear pore complex']. 
Make sure that sequential RAW files going into MaxQuant analysis are sequentially (alphanumerically) named, e.g. timepoint 1, 2, 3.. are labeled 
a1, a2, a3.., etc. Analyze one experiment at a time, e.g. all time points for turnover or exchange but not both.

All outputs of this script are saved into data_out folder. Note, the analysis is intended for heavy lysine labeling, if you
have an alternate label, e.g. heavy arginine, you have to edit the script to include the corresponding columns. Additionally, the final
grouping is based on yeast gene descriptions from SGD. You have to provide an alternative FASTA file for other organisms, the one used
for MaxQuant search will work. Throughout the script points for substitution are highlighted in comments.

Analysis workflow:

1. Filter out contaminants, reverse matches and peptides with no lysine (the latter provides no quantitative information). Write the
filtered out peptides to data_out/filename_filtered_out.txt file.

2. Compute H/(H+L) from H/L column estimated by MaxQuant for all remaining peptides.

3. Group the data by RAW file, i.e. by each time point.

4. Within each group (step 3), group peptides into the 'leading razor protein' identified by Maxquant, which is usually 1 per peptide.

5. Stat_analysis_final.py script is applied to every group of peptides constituting a protein to eliminate outliers, compute the
average and standard deviation of H/(H+L) fractions, and count the number of unique peptides (the higher the number, the more reliable
the protien identification and H/(H+L) fraction measurement.

6. Log-transform the average H/(H+L) values and perform linear regression with the time as independent variable and the heavy fraction
as dependednt variable:
    * at least 3 timepoints are required for fitting
    * if there are 4 or 5 data points available and the initial fit does not have a good r 
      cotrrelation coefficient (i.e. r<0.9), the data points are sequenttialy removed
      (only one at a time) and the regression performed on the resulting subset of the original dataset.
    * the fit with best correlation coefficient is saved
    * record the slope, intercept, correlation coefficient and P value (two-sided null hypothesis test that the slope=0
7. the systematic names of proteins are matched to the common names in the SGD database and, based on gene descriptions,
   the proteins are split into 2 goups: 1. proteins of interest, gene descriptions contain keywords from the provided list
   2. all other protsins - usually composed of anundant contaminants that can be used as a baseline/ null distribution.
8. a summary plot is constructed for the 2 groups. This is intended for a quick visual assessment of data and fit and not as
   a final analysis figure.


 The three groups are plotted in 
different colors for visual inspection. All the files are organized in new folders: MaxQuant, data_out, 
data_out/Assigned, data_out/Excluded, data_out/Graphs, data_out/Intermediate."""

import sys
import stat_analysis_final_log
import os
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
import json
if not os.path.exists('data_out'):
    os.makedirs('data_out')

# open the evidence.txt file resulting from MaxQuant run
f = open(sys.argv[1], 'r')

# create a variable name from the file name
file_name = sys.argv[1]
file_name = file_name.split('data_in/')[1]
file_name = file_name.split('.txt')[0]

# create names for various files were the calculated, analyzed or excluded data will
# be stored, which is based on the original file name
filt_name = file_name + '_filtered_out'
analyz_name = file_name + '_analyzed'

# now create writable text files with all of the names
filt = open('data_out/'+filt_name+'.txt', 'w') # the data that cannot be used (e.g. contaminant, reverse or
# no K peptides) will be stored here
analyzed = open('data_out/'+analyz_name+'.txt', 'w') #data that will go into analysis will be written here
# with an aditional column Ratio H/(H+L) that will be calculated from Ratio H/L

# gets the first line, column headers
first_line = f.readline()
first_line=first_line.rstrip('\r\n')

# coppy the header line to filtered and analyzed files
filt.write(first_line)
filt.write('\n')
analyzed.write(first_line)
analyzed.write('\tRatio H/(H+L)\n')


# make all lower letter for case insensitive search
first_line=first_line.lower()
headers = first_line.split('\t')
#this is for case insensitive matching

# get the column indices of interest
i0=headers.index("sequence")
i1=headers.index("k count") # note, if you are using arginine instead of lysine, edit this column to "r count"
i2=headers.index("leading razor protein")
i3=headers.index("raw file")
i4=headers.index("ms/ms count")
i5=headers.index("ratio h/l")
i6=headers.index("reverse")
i7=headers.index("potential contaminant")

"""
this part of the script writes all the unusable data to a _filtered_out.txt 
file, all the usable data into _analysed.txt file, and keeps track of all the different raw files.
Also we create a dictionary, where the keys are raw files and the values are arrays containing all 
the peptide characteristics to be analyzed from the corresponding raw file.
the dictionary is printed out based on keys into separate txt files which have 4 column entries - "Sequence", 
"Leading Razor Protein", "MS/MS Count", "ln(Ratio H/(H+L)"""

raw_files = {}
for lines in f:
    entries=lines.split('\t')
    if entries[i1] == '0':
        filt.write('\t'.join(entries))
    elif entries[i4] == '0':
        filt.write('\t'.join(entries))
    elif entries[i5] == '' or entries[i5] == "NaN":
        filt.write('\t'.join(entries))
    elif entries[i6] == "+":
        filt.write('\t'.join(entries))
    elif entries[i7] == "+":
        filt.write('\t'.join(entries))
    else:
        line = '\t'.join(entries)
        analyzed.write(line.rstrip('\r\n'))
        analyzed.write('\t')
        ratio = float(entries[i5])/(float(entries[i5])+float(1))
        if ratio > 0:
            ratio = np.log(ratio)
            ratio=round(ratio, 5)
            analyzed.write(str(ratio))
            analyzed.write('\n')
            line = (entries[i0]+'\t'+entries[i2]+'\t'+entries[i4]+'\t'+str(ratio)+'\n')
            if entries[i3] not in raw_files:
                raw_files[entries[i3]] = [line]
            else:
                raw_files[entries[i3]].append(line)
        else:
            filt.write('\t'.join(entries))
f.close()
filt.close()
analyzed.close()

# this part writes the contents from the corresponding RAW files into separate text files

raw_file_header = ("Sequence"+'\t'+"Leading Razor Protein"+'\t'+"MS/MS Count"+'\t'+"ln(H/(H+L))"+'\n')

for files in sorted(raw_files):
    f = open('data_out/'+file_name+files+'.txt', 'w')
    f.write(raw_file_header)
    for line in raw_files[files]:
        f.write(line)
    f.close()

# Next, the text files are read and organised by proteins. A dictionary is created
# for each file, the keys are the leading razor proteins,
# the values are arrays with corresponding peptides and characteristics

protein_analysis_headers = "Protein systematic name"+'\t'+"Number of unique peptides"+'\t'+"Total MS/MS counts"+'\t'+"Average H/(H+L)"+'\t'+"Standard deviation"+'\n'

for files in sorted(raw_files):
    f = open('data_out/'+file_name+files+'.txt', 'r')
    f1 = open('data_out/'+file_name+files+'_processed.txt', 'w') # average, standard deviation of heavy fractions and unique peptide count written here
    f1.write(protein_analysis_headers)
    first_line = f.readline()
    headers = first_line.split('\t')
    i0=headers.index("Sequence")
    i2=headers.index("Leading Razor Protein")
    i4=headers.index("MS/MS Count")
    i8=headers.index("ln(H/(H+L))\n")
    proteins={}
    for lines in f:
        entries=lines.split('\t')
        line = (entries[i0]+'\t'+entries[i4]+'\t'+entries[i8]+'\n')
        if entries[i2] not in proteins:
            proteins[entries[i2]]=[line]
        else:
            proteins[entries[i2]].append(line)
    for protein in proteins:
        l=stat_analysis_final_log.stat_analysis(proteins[protein])
        f1.write(protein+'\t')
        f1.write(str(l[0])+'\t'+str(l[1])+'\t'+str(l[2][0])+'\t'+str(l[2][1])+'\n')
    f1.close()
    f.close()
    
# the _final_analysis text file will contain the time course data for each protein that will be fit with linear regression
f2=open('data_out/'+file_name+'_final_analysis.txt', 'w')

# create an empty list of protein names to be populated with protein names from all files
protein_names=[]    

# create an empty list of dictioneries, that are going to contain protein-key, descriptor-value pairs
prot_descs=[]

for files in sorted(raw_files):
    f3 = open('data_out/'+file_name+files+'_processed.txt', 'r')
    f3.readline() # get rid of the header row
    # make an empty dictionary named "prot_desc" and add it to the file_name list
    prot_desc={}
    # get rid of the new line character, split the line on tabs, add the protein name 
    #to the list of proteins, create a key-value pair from the protein name and the rest 
    #of the characteristics
    for lines in f3:
        lines=lines.rstrip('\n')
        entries=lines.split('\t')
        protein_names.append(entries[0])
        prot_desc[entries[0]]=('\t'.join(entries[1:]))
    f3.close()
    prot_descs.append(prot_desc)
# go through the list of protein names and write the protein names and values from each file
# next to each other
num_files=len(prot_descs) # this counts the number of time points in the experiment
protein_analysis_header = ("Protein systematic name"+('\t'+"Number of unique peptides"+'\t'+"Total MS/MS counts"+'\t'+"Average ln(H/(H+L))"+'\t'+"Standard deviation ln(H/H+L)")*len(raw_files)+'\n')
f2.write(protein_analysis_header)

# this part of the code pulls out the protein measurement characteristics from each file
# and prints them next to each other in the _final_analysis.txt file for linear regression.
# the order of values is based on string sorted file_name, that is why sequential time point samples 
# need to be named in an alphanumerical ascending sequeantial order
unique_protein_names=set(protein_names)
four_tabs=('\t0'*4)
for proteins in unique_protein_names:
    f2.write(proteins)
    sorted(prot_descs, key=str)
    for prot_desc in prot_descs:
        if proteins not in prot_desc:
            f2.write(four_tabs)
            continue
        values=prot_desc[proteins]
        val_entries = values.split('\t')
        f2.write('\t'+values)
    f2.write('\n')
f2.close()

# this part deletes, unused text files and can be optionally removed in case the individual samples need to be separated
for files in sorted(raw_files):
    os.remove('data_out/'+file_name+files+'.txt')
    os.remove('data_out/'+file_name+files+'_processed.txt')
    
"""
This part performs a linear regression on the log-transformed average H/(H+L) measurements over the time course. 
Note, (0,log(1)) intercept is not enforced, since there is incomplete labeling and excess tagged protein in the  
first hour of the time course. 
"""

f = open('data_out/'+file_name+'_final_analysis.txt', 'r')
f1 = open('data_out/'+file_name+'_fitted.txt', 'w')
f2 = open('data_out/'+file_name+'_notfitted.txt', 'w')
f3 = open('data_out/'+file_name+'_poorfit.txt', 'w')
f4 = open('data_out/'+file_name+'_lowconfidence.txt','w')    
average_timepoint=''
stdev_timepoint=''
for i in range(num_files):
    m=str(i+1)
    average_timepoint+='\tAverage '+m
    stdev_timepoint+='\tStandard deviation '+m
f1.write("Protein systematic name"+average_timepoint+stdev_timepoint+"\tSlope\tIntercept\tr value\tp value\tstd err\tAve uniqiue peptides\tAve MS/MS count\tTime point filtered\n")
f2.write("Protein systematic name\tTime points\tAve unique peptide\tAve MS/MS count\n")
f3.write("Protein systematic name"+average_timepoint+stdev_timepoint+"\tSlope\tIntercept\tr value\tp value\tstd err\tAve uniqiue peptides\tAve MS/MS count\tTime point filtered\n")
f4.write("Protein systematic name\tCommon name\tDescription"+average_timepoint+stdev_timepoint+"\tSlope\tIntercept\tr value\tp value\tstd err\tAve uniqiue peptides\tAve MS/MS count\tTime point filtered\n") 

with open(sys.argv[2], 'r') as arg_file:
    arg_dict = json.load(arg_file)
#print(arguments)
#arg_dict = json.loads(arguments)
#t = sys.argv[2]
t = np.array(arg_dict['time'])

f.readline()
for lines in f:
    lines=lines.rstrip('\r\n')
    entries=lines.split('\t')
    vals = np.array(entries[1:len(entries)], dtype=float)
    unique=vals[0:len(vals):4]
    counts=vals[1:len(vals):4]
    ave_unique = np.average(unique, axis=None)
    ave_counts = np.average(counts, axis=None)
    average=vals[2:len(vals):4]
    stdev=vals[3:len(vals):4]
    average_filt=vals[2:len(vals):4]
    average_filt=average_filt.ravel()[np.flatnonzero(average_filt)]
    average=average.tolist()
    for n,i in enumerate(average):
        if i==0.:
            average[n]='nan'
            stdev[n]='nan'
    if average_filt.shape[0]>=3: 
        unique=unique.ravel()[np.flatnonzero(average_filt)]
        counts=counts.ravel()[np.flatnonzero(average_filt)]
        tmod=t.ravel()[np.flatnonzero(average_filt)]
        slope, intercept, r_value, p_value, std_err = stats.linregress(tmod,average_filt)
        """at least 3 points are needed for the regression, proteins with less data points are added to the not-fitted list.
        Fits with r correlation coefficient<0.9 and p value >0.05 are either discarded (added to the poor fit list) or
        further manipulated if there are 4 or 5 points. Linear regression is performed on the subset of the data and if r and p
        values do not improve, then the protein is added to the poor fit list otherwise it is recorded in the assigned list. 
        Proteins with >=3 time points, r correlation coefficient >=0.9 and p value <=0.05 are added to the assigned list."""
        if -r_value >= 0.9 and p_value<=0.05:
            f1.write(entries[0]+'\t'+'\t'.join(map(str, average))+'\t'+'\t'.join(map(str, stdev))+'\t'+str(round(slope,5))+'\t'+str(round(intercept,5))+'\t'+str(-round(r_value,5))+'\t'+str(round(p_value,5))+'\t'+str(round(std_err,5))+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\t'+'\n')
        else:
            if (average_filt.shape[0]==3):
                f3.write(entries[0]+'\t'+'\t'.join(map(str, average))+'\t'+'\t'.join(map(str, stdev))+'\t'+str(round(slope,5))+'\t'+str(round(intercept,5))+'\t'+str(-round(r_value,5))+'\t'+str(round(p_value,5))+'\t'+str(round(std_err,5))+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\n')
            elif (average_filt.shape[0]>3):
                j=0
                average3=[]
                stdev3=[]
                for i in range(len(average_filt)):
                    average2=np.delete(average_filt, i)
                    tmod_new=np.delete(tmod,i)
                    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(tmod_new,average2)
                    if -r_value1 >= 0.9 and p_value1 <= 0.05:
                        j=i
                        slope, intercept, r_value, p_value, std_err= slope1, intercept1, r_value1, p_value1, std_err1
                if -r_value >= 0.9 and p_value <= 0.05: 
                    average3=np.delete(average, j) # if the exclusion of the value improved the fit, the value is discarded
                    average3=average3.tolist()
                    average3.insert(j, 'nan')
                    stdev2=np.delete(stdev,j)
                    stdev3=stdev2.tolist()
                    stdev3.insert(j,'nan')
                    f1.write(entries[0]+'\t'+'\t'.join(map(str,average3))+'\t'+'\t'.join(map(str,stdev3))+'\t'+str(round(slope,5))+'\t'+str(round(intercept,5))+'\t'+str(-round(r_value,5))+'\t'+str(round(p_value,5))+'\t'+str(round(std_err,5))+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\tYes'+'\n')
                else:
                    f3.write(entries[0]+'\t'+'\t'.join(map(str, average))+'\t'+'\t'.join(map(str, stdev))+'\t'+str(round(slope,5))+'\t'+str(round(intercept,5))+'\t'+str(-round(r_value,5))+'\t'+str(round(p_value,5))+'\t'+str(round(std_err,5))+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\t'+'\n')
    else:
        f2.write(entries[0]+'\t'+str(average_filt.shape[0])+'\t'+str(round(ave_unique,5))+'\t'+str(round(ave_counts,5))+'\n')
f.close()
f1.close()
f2.close()  
f3.close()      
"""
This part searches the SGD database for the systematic names of proteins and groups 
the final set of proteins with fitted values into 2 groups: proteins of interest and "other".
This is done  for assigned, not fitted and poor fit files, separately.
"""

f = open('data_out/'+file_name+'_fitted.txt', 'r')
f1 = open('data_out/'+file_name+"_poi_grouped.txt", 'w')
f2 = open('data_out/'+file_name+"_other_grouped.txt", 'w')
fp = open('data_out/'+file_name+'_poorfit.txt', 'r')
fp1 = open('data_out/'+file_name+"_poi_poorfit.txt", 'w')
fp2 = open('data_out/'+file_name+"_other_poorfit.txt", 'w')
fn = open('data_out/'+file_name+'_notfitted.txt', 'r')
fn1 = open('data_out/'+file_name+"_poi_notfitted.txt", 'w')
fn2 = open('data_out/'+file_name+"_other_notfitted.txt", 'w')

# Extract gene descriptions from SGD
# Note, if using another organism, download the proteome and extract gene descriptions. Need to match the systematic name of the 'leading rasor protein' with gene description.
protein_desc={}
prot=open("data_in/orf_trans.fasta", 'r')
for line in prot:
    if re.match('^>', line):
        line=line.rstrip()
        line=line[1:]
        entries=line.split(' ')
        protein_description=(' '.join(entries[2:]))
        protein_desc[entries[0]]=(entries[1]+'\t'+protein_description)
    else:
        next
              
header=f.readline()
names= header.split('\t')
i9=names.index("Ave uniqiue peptides")
new_header=(names[0]+'\t'+"Common name"+'\t'+"Description"+'\t'+('\t'.join(names[1:])))
f1.write(new_header)
f2.write(new_header)

keyword_list = arg_dict['keywords']
#keyword_list = sys.argv[3]
for line in f:
    entries=line.split('\t')
    protein_stats=('\t'.join(entries[1:]))
    first_name=entries[0].split(';')
    match = False
    for keyword in keyword_list:
        if re.search(keyword, protein_desc[first_name[0]], flags=re.IGNORECASE):
            match = True
            if float(entries[i9])>=4:
                f1.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats)
            elif float(entries[i9])<4:
                f4.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats)
    if match == False:
        if float(entries[i9])>=4:
            f2.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats)
        elif float(entries[i9])<4:
            f4.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats)


f.close()
f1.close()
f2.close()

f4.close()

headerp=fp.readline()
namesp= headerp.split('\t')
new_headerp=(namesp[0]+'\t'+"Common name"+'\t'+"Description"+'\t'+('\t'.join(namesp[1:])))
fp1.write(new_headerp)
fp2.write(new_headerp)

for line in fp:
    entries=line.split('\t')
    protein_stats=('\t'.join(entries[1:]))
    first_name=entries[0].split(';')
    match = False
    for keyword in keyword_list:
        if re.search(keyword, protein_desc[first_name[0]], flags=re.IGNORECASE):
            fp1.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats)
            match = True
    if match == False:
        fp2.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats)
        
fp.close()
fp1.close()
fp2.close()

headern=fn.readline()
namesn= headern.split('\t')
new_headern=(namesp[0]+'\t'+"Common name"+'\t'+"Description"+'\t'+('\t'.join(namesn[1:])))
fn1.write(new_headern)
fn2.write(new_headern)

for line in fn:
    entries=line.split('\t')
    protein_stats=('\t'.join(entries[1:]))
    first_name=entries[0].split(';')
    match = False
    for keyword in keyword_list:
        if re.search(keyword, protein_desc[first_name[0]], flags=re.IGNORECASE):
            fn1.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats)
            match = True
        if match == False:
            fn2.write(first_name[0]+'\t'+protein_desc[first_name[0]]+'\t'+protein_stats)
        

fn.close()
fn1.close()
fn2.close()
prot.close()

os.remove('data_out/'+file_name+'_fitted.txt')
os.remove('data_out/'+file_name+'_poorfit.txt')
os.remove('data_out/'+file_name+'_notfitted.txt')

"""
finally the data and fits in the assigned (fitted) file are plotted for visual inspection
"""

f1=open('data_out/'+file_name+"_poi_grouped.txt", 'r')
f2=open('data_out/'+file_name+"_other_grouped.txt", 'r')

x=np.arange(0, t[-1]+0.5, 0.0025)

first_line=f1.readline()
entries=first_line.split('\t')
i10=entries.index("Slope")
i11=entries.index("Intercept")
f2.readline()
lw=0.75
        
timepoints=[]
errors = []     
for line in f1:
    line=line.rstrip('\r\n')
    entries=line.split('\t')
    param2=float(entries[i10])
    y=param2*x
    plt.plot(x,y,'skyblue', label='proteins of interest', lw=lw)
    timepoints=entries[3:3+num_files]
    for n,i in enumerate(timepoints):
        timepoints[n]=float(timepoints[n])
    timepoints[:] = [x - float(entries[i11]) for x in timepoints]
    errors = entries[7:7+num_files]
    for n,i in enumerate(errors):
        errors[n]=float(errors[n])
    plt.errorbar(t, timepoints, yerr=errors, fmt='o', color='lightgrey')
    

# to have a basis for comparison with "other" proteins in the pullout, largely abundant housekeeping
#proteins, the average +/-2 standard deviation of the observed population is plotted as well
       
params=np.array([])    
for line in f2:
    line=line.rstrip('\r\n')
    entries=line.split('\t')
    param2=float(entries[i10])
    params=np.append(params, param2)
ave_param=np.mean(params)
std_param=np.std(params)
y=ave_param*x
plt.plot(x,y,'black', label='Average', lw=1)

ave_minus=ave_param-2*std_param
y1=ave_minus*x
plt.plot(x,y1,'black', label='Average +/- 2 st dev', lw=1, linestyle='dashed')
ave_plus=ave_param+2*std_param
y2=ave_plus*x
plt.plot(x,y2,'black', lw=1, linestyle='dashed')

plt.xlabel('Time (h)', fontsize=20)
plt.ylabel('Log(H/(H+L))',fontsize=20)
plt.title(file_name+' time course\n', fontsize=24) # this is either turnover or exchange time-course
plt.xticks(np.arange(0, max(t)+1, 1), fontsize=18)
plt.xlim(0,t[-1]+0.5)
plt.box('off')
plt.tick_params(
    axis='x',          
    which='both',      
    bottom='off',     
    top='off') 
plt.tick_params(
    axis='y',          
    which='both',      
    left='off',     
    right='off')  

# include 1 legend for all graphs.
handles, labels = plt.gca().get_legend_handles_labels()
newLabels, newHandles = [], []
for handle, label in zip(handles, labels):
  if label not in newLabels:
    newLabels.append(label)
    newHandles.append(handle)
plt.legend(newHandles, newLabels, loc=0, frameon=False)

plt.savefig('data_out/'+file_name+'graph.eps' , bbox_inches='tight', format='eps', dpi=1000)
f1.close()
f2.close()

# make some new folders to organize all the files

if not os.path.exists('data_out/Graphs'):
    os.makedirs('data_out/Graphs')
if not os.path.exists('data_out/Excluded'):
    os.makedirs('data_out/Excluded')
if not os.path.exists('data_out/Intermediate'):
    os.makedirs('data_out/Intermediate')
if not os.path.exists('data_out/Assigned'):    
    os.makedirs('data_out/Assigned')
if not os.path.exists('data_out/Excluded/Notfitted'):
    os.makedirs('data_out/Excluded/Notfitted')
if not os.path.exists('data_out/Excluded/Poorfit'):
    os.makedirs('data_out/Excluded/Poorfit')
if not os.path.exists('data_out/Excluded/Lowconfidence'):
    os.makedirs('data_out/Excluded/Lowconfidence')

src = "data_out"
dst1 = "data_out/Excluded"
dst2 = "data_out/Graphs"
dst3 = "data_out/Intermediate"
dst4 = "data_out/Assigned"
dst5 = "data_out/Excluded/Notfitted"
dst6 = "data_out/Excluded/Poorfit"
dst7 = "data_out/Excluded/Lowconfidence"

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
        os.system ("mv"+ " " + fullPath + " " + dst7)

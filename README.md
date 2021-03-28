# Procedure
This repo accompanies the publication "Protein dynamics quantification" in cell protocols. The two objectives are to organize the heavy label quantification values produced by Maxquant and fit the time course data to estimate the apparent exchange/turnover rates of proteins of interest.
Protein_analysis_pipeline.py file provides step by step description of the procedure with inline comments.

1. clone the protein-dynamics-pipeline repository<br>
```git clone https://github.com/zhakhverdyan/protein-dynamics-pipeline.git```<br>

2. Create and activate an environment with conda<br>
```conda create conda env create -f environment.yml```<br>
```conda activate protein-dynamics-pipeline```<br>

* Note alternatively, pip/conda install scipy>=1.6 and matplotlib>=3.3
3. Prepare data, sample files are provided for testing:<br>
	a. Place the evidence.txt file resulting from Maxquant run into data_in folder<br>
	b. Place the orf_trans.fasta (or alternative pasta file with the proteome of the organism into data_in folder<br>
	c. Prepare a json file with two input arguments: 1. key = "time", value = list of timepoints in hours, 2. key = "keywords", value = list of keywords used to pull proteins of interest.<br>
4. Run the protein_analysis_pipeline.py script with two input arguments: the evidence.txt file and arguments.json file<br>
```./protein_analysis_pipeline.py data_in/evidence.txt data_in/arguments.json```<br>
5. Analyze the output in data_out folder:
* Excluded - contains data excluded from analysis at various stages, e.g, no lysine, not enough peptides, poor fit, etc. see protein_analysis_pipeline.py for detailed description
* Intermediate - contains two files: evidence_analyzed.txt - all peptides that made it into analysis, evidence_final_analysis.txt - all protein timepoints data collated
* Assigned - contains two files - evidence_poi_grouped.txt - summary of time point measurements and linear regression fit values for proteins of interest, evidence_other_grouped.txt - summary of time point measurements and linear regression fit values for all other proteins
* Graphs - graph of timepoints measurements and fitted curves for protein of interest (from evidence_poi_grouped.txt); the 3 diagonals are average +/- 2 standard deviations of the baseline protein fitted values (from evidence_other_grouped.txt). This graph can be used to quickly assess the data, e.g. the further away fitted curves are from the baseline, the slower the exchange of the protein of interest

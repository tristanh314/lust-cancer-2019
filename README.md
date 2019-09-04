#### lust-cancer-2019 ####
This repository contains the tools and results of a data analysis of gene expression and clinical outcomes available from the Broad institute at https://gdac.broadinstitute.org/ using the LUST algorithm (see Dbasinapp.pdf). We have aligned clinical data with genetic expression data for patients in this database. The first part of the algorithm then analyzes the genetic expression data for each patient and discretizes the data by evaluating whether a gene in a given sample is expressed at a "low," "normal," or "high" rate. At this stage our programs (see MATLAB_functions and data_cleaning) then look for genetic signiatures, combinations of gense and expression values, that are common in the genetic expression data (see metagenes). We have ranked these genetic signiatures for each type of tumor by how often each signiature is present among the samples (see comparisons_and_overlaps/Master Metagene Comparison.xlsx). Additionally we have created tables to show the overlap between metagenes (see comparisons_and_overlaps/Master Metagene Overlaps.xlsx and comparisons_and_overlaps/Significant Overlaps.xlsx).

The second part of the algorithm takes clinical data and attempts to find genetic signiatures that are good predictors for long term patient survival. So far this has been done for liver and stomach cancers and the signiatures found are of various effectiveness. The files for this stage of the analysis will be uploaded soon.


#### File Descriptions ####

clean_clinical_data/

	*_clin_clean.tar.xz: Compressed csv file of cleaned clinical data for the given type of 	tumor (note that the various kinds of sarcoma are not separated by organ type, but that 	this information is available in the file)

clean_gene_data/

	*_gene_clean.tar.xz: Compressed csv file of cleaned genetic expression data for the given 		type of tumor

comparisons_and_overlaps/

	Master Metagene Comparison.xlsx: Table showing the rank of the prominence of each metagene 		for each type of tumor

	Master Metagene Overlaps.xlsx: Shows the intersection of genetic factors for all 		significant metagenes

	Significant Overlaps.xlsx: Shows those metagenese where significant overlaps occur

data_cleaning/

	Firehose index.txt: Summary of the available data available from the Firehose repository 		from the Broad institute

	master_data_cleaner.py: Script use to clean the data from Firehose to align clinical and 		genetic data

MATLAB_functions/

	LUST_manual.pdf: The entire process of using the LUST algorithm, beginning to end,
	is a long one. This document is intended to help the user get started.

	MATLAB_function_density_calc.txt: Subroutine for density_fix. Saved as a text file that 	may be copy and pasted for use in MATLAB or Octave

	MATLAB_function_density_fix.txt: This program adjusts the 'factor' paramter so that in 		descretizing, expression values > 1/factor become +1, those < 1/factor become -1. The 		secant method is used to achieve the target density 'target'. It is called by the main 		LUST program, though could be used separately. Saved as a text file that may be copy and 		pasted for use in MATLAB or Octave

	MATLAB_function_eigsurvplot_topbot_noage.txt: Input is a signature, gene expression file, 		clinical data nodays is mean survival estimate from ecdf. Output is prediction accuracy 	for various thresholds on the escore

	MATLAB_function_logrank2.txt: Log-rank test(Approximate). The logrank test is a hypothesis 		test to compare the survival distributions of two samples. It is a nonparametric test and 		appropriate to use when the data are right censored (technically, the censoring must be 	non-informative).
	
	MATLAB_function_lust_find_metagenes.txt: This version finds metagenes from the large 		expression file. The groups are ranked using a graph density objective function, and 		overlaps listed. Input:  rids of genes, expression file, density, conftol; dens is the 		desired density of +1,-1; conftol is the confidence parameter. There are also some 		internal parameters For merging the groups which are hard-wired below.

	MATLAB_function_lust_find_metagenes_july19.txt: This version finds metagenes from the 		large expression file. The groups are ranked using a graph density objective function, and 		overlaps listed. Input:  rids of genes, expression file, density, conftol; dens is the 		desired density of +1,-1; conftol is the confidence parameter. There are also some 		internal parameters For merging the groups which are hard-wired below.

	MATLAB_function_lust_find_meta_to_sigs_july19.txt: This takes a metagene, identified by 	the rids_file, extracts its expression file, and uses the LUST algorithm to find 		signatures. Signatures are ranked by their Fisher score. Kaplan-Meier curves are plotted 		for the top ones. Input: rids_file is rids for the metagene; exp_file is expression for 	the cancer, clinical file is survival, censoring, etc.; metax is the label for the 		metagene, e.g., 'stomach meta R2'; dens is density; algorithm is run at several conftol 	values, starting with conftolmin. Output: Top 10 signatures; Their Fisher scores and 		Kaplan-Meir survival curves; Differential expression on immune system genes for groups 		identified by the signature.

	MATLAB_function_preproc.txt: Takes expression data, adds 1, log2, quantilenorm, center, 	frobenius norm. Saved as a text file that may be copy and pasted for use in MATLAB or 		Octave.

	MATLAB_function_symb2indx.txt: 

	TCGA_trimrids.mat:

metagenes/

	Metagene_*_list.xlsx: Lists the genes present in the given metagene with their expressions 		levels

papers/

	Dbasinapp.pdf: We introduce the parameter of relevance of an attribute of a binary 		table to another attribute of the same table, computed with respect to an implicational 	basis of a closure system associated with the table. This enables a ranking of all 		attributes, by relevance parameter to the same fixed attribute, and, as a consequence, 		reveals the implications of the basis most relevant to this attribute. As an application 		of this new metric, we test the algorithm for D-basis extraction presented in Adaricheva 		and Nation [1] on biomedical data related to the survival groups of patients with 		particular types of cancer. Each test case requires a specialized approach in converting 		the real-valued data into binary data and careful analysis of the transformed data in a 	multi-disciplinary environment of cross-field collaboration.

	lust_2019_part_1.pdf: The Lattice Up-Stream Targeting (LUST) algorithm is a discrete 		mathematical method for analyzing expression data. It uses a variation of association 		rules to find groups of genes whose expression is correlated. These sets of genes are 		called metagenes, as they are associated with a common biological process or func-
	tion. Metagenes can be refined to smaller subsets called signatures that represent the 		entire metagene. This study uses the LUST algorithm to find metagenes and predictive 		signatures for the 33 different types of cancer in the TCGA database, based on mRNA 		expression data. This allows us to identify the metagenes that are common to multiple 		cancers, and those that occur with only one or two types. Knowing which metagenes are 		associated with a particular cancer enables us to determine factors that are significant 		in determining the progress of the disease, which are then candidates for further study.
	This first part concentrates on the mathematical background for the LUST algorithm. A 		second part will present the biological results, analyzing the different metagenes that 	arise.

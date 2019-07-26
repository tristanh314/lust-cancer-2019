#### lust-cancer-2019 ####
This repository contains the tools and results of a data analysis of gene expression and clinical outcomes available from the Broad institute at https://gdac.broadinstitute.org/ using the LUST algorithm (see Dbasinapp.pdf). We have aligned clinical data with genetic expression data for patients in this database. The first part of the algorithm then analyzes the genetic expression data for each patient and discretizes the data by evaluating whether a gene in a given sample is expressed at a "low," "normal," or "high" rate. At this stage our programs a then look for genetic signiatures, combinations of gense and expression values, that are common in the genetic expression data (see Metagene_*_list.xlsx). We have ranked these genetic signiatures for each type of tumor by how often each signiature is present among the samples (see Master Metagene Comparison.xlsx). Additionally we have created tables to show the overlap between metagenes (see Master Metagene Overlaps.xlsx and Significant Overlaps.xlsx).

The second part of the algorithm takes clinical data and attempts to find genetic signiatures that are good predictors for long term patient survival. So far this has been done for liver and stomach cancers and the signiatures found are of various effectiveness. The files for this stage of the analysis will be uploaded soon.


#### File Descriptions ####

*_clin_clean.tar.xz: Compressed csv file of cleaned clinical data for the given type of tumor (note that the various kinds of sarcoma are not separated by organ type, but that this information is available in the file)

*_gene_clean.tar.xz: Compressed csv file of cleaned genetic expression data for the given type of tumor

Dbasinapp.pdf Article: regarding using the D-basis in data analysis in biomedical studies

Firehose index.txt: Summary of the available data available from the Firehose repository from the Broad institute

master_data_cleaner.py: Script use to clean the data from Firehose to align clinical and genetic data

Master Metagene Comparison.xlsx: Table showing the rank of the prominence of each metagene for each type of tumor

Master Metagene Overlaps.xlsx: Shows the intersection of genetic factors for all significant metagenes

Metagene_*_list.xlsx: Lists the genes present in the given metagene with their expressions levels

Significant Overlaps.xlsx: Shows those metagenese where significant overlaps occur

Scripts used to identify metagenes coming soon!

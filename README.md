This repository contains the necessary codes to reproduce the figures presented in the PHACTboost manuscript. Each figure's code and related data are provided in their respective named folders.

## PHACTboost_Model

This folder contains the PHACTboost model and variants in the training set and TS1-5, which are alternative test sets used in the manuscript. The definition of these sets are given below:

![Picture1](https://github.com/CompGenomeLab/PHACTboost_manuscript/assets/68369488/3750f6e3-e1a7-499c-86a3-a5d883eb748d)


## Figure 1

Figure 1C illustrates the evaluation of the PHACTboost algorithm's performance across all alternative test sets, revealing the enhancements achieved over the baseline PHACT model, and two alternative models MSAboost and TREEboost. Both MSAboost and TREEboost employs the same machine learning algorithm but differing in their feature sets.

Models of MSAboost and TREEboost is also provided in folders, [MSAboost](Figure1/Data/MSAboost) and [TREEboost](Figure1/Data/TREEboost), respectively.

Figure 1D is dedicated to exploring Feature Importance – Shapley Values for the PHACTboost Algorithm.

Figure 1C are generated using the script "Figure_1C.R". The script for plotting Figure 1D is named "Figure_1D.R".

The data provided in the [Figure1/Data](Figure1/Data) folder was generated using the script "Figure1_ProduceData.R".

## Figure 2

Figure 2A displays AUROC comparisons between PHACTboost and tools reported in dbNSFP across TS1, TS2, TS3, TS4, and TS5, respectively.

To generate Part A and B, you can use the script "Figure2A.R," and "Figure2B.R," which are located in the same directory, respectively.

The necessary files are provided in the Figure2/Data folder.

Figure 2C displays examples of variants misclassified by VEST4 and REVEL but correctly predicted by PHACTboost. The script "Figure2C.R" is used to produce the subplots of the phylogenetic tree and MSA.

The data found in the [Figure2/Data](Figure2/Data) folder was generated using the script "Figure2_ProduceData.R.

## Figure 3


Figures 3A, 3B, and 3C illustrate the AUROC difference between PHACTboost and AlphaMissense, EVE, and CPT-1 across all alternative test sets, respectively. You can find the code to produce these results, along with AUPR comparison, in 'Figure3A_C.R.'

Figure 3D compares the average Spearman correlation between PHACTboost and DMS results with other tools and DMS data, encompassing all variants and hard cases (cases with conflicting pathogenicity predictions). You can access the code to reproduce these results in 'Figure3D.R.'

In Figure 3E, a Chi-square test examines whether there is a significant difference between pathogenic variants mislabeled by PHACTboost and others in terms of their submission to ClinVar, following ACMG guidelines. The code for testing this significance is available as 'Figure3E.R.'

Figure 3F highlights a pathogenic variation mislabeled by PHACTboost due to a lack of consideration for coevolution. The code for generating the MSA and phylogenetic tree plots can be found in 'Figure3F.R.'

The data in the [Figure3/Data](Figure3/Data) folder was generated using the script 'Figure3_ProduceData.R.'

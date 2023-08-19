# PHACTboost_manuscript

This repository contains the necessary codes to reproduce the figures presented in the PHACTboost manuscript. Each figure's code and related data are provided in their respective named folders.

## Figure 2

Figure 2A illustrates the evaluation of the PHACTboost algorithm's performance across all alternative test sets, revealing the enhancements achieved over the baseline PHACT model.

Figure 2B is dedicated to exploring Feature Importance – Shapley Values for the PHACTboost Algorithm.

In Figure 2C, we present a comparative analysis between PHACTboost and MSAboost, both employing the same machine learning algorithm but differing in their feature sets. This comparison was conducted across the entire dataset and encompassed all alternative test sets.

Both Figure 2A and 2C are generated using the script "Figure_2AC.R". The script for plotting Figure 2B is named "Figure_2B.R".

The data provided in the [Figure2/Data](Figure2/Data) folder was generated using the script "Figure2_ProduceData.R".

## Figure 3

Figure 3A, 3B, 3C, 3D, and 3E display AUROC comparisons between PHACTboost and tools reported in dbNSFP across TS1, TS2, TS3, TS4, and TS5, respectively.

To generate Part A, B, C, D, and E, you can use the script "Figure3_Comparison_TS.R," which is located in the same directory.

The necessary files are provided in the Figure2/Data folder.

The data found in the [Figure3/Data](Figure3/Data) folder was generated using the script "Figure3_ProduceData.R.

## Figure 4

A. AUROC and AUPR difference between PHACTboost and EVE, CPT-1, VARITY, gMVP and ENVISION over the whole test set. B. AUROC and AUPR difference between PHACTboost and EVE, CPT-1, VARITY, gMVP and ENVISION over hard cases. C. Protein-by-protein comparison of PHACTboost with EVE over 60 proteins with at least 3 neutral and 3 pathogenic variations. D. Comparison of average Spearman correlation between PHACTboost and DMS results with the average Spearman correlation between other tools and DMS data.




# PHACTboost_manuscript

This repository includes the required codes to reproduce the figures at PHACTboost manuscript.

## Figure 2

Figure 2A illustrates the evaluation of the PHACTboost algorithm's performance across all alternative test sets, revealing the enhancements achieved over the baseline PHACT model.

Figure 2B is dedicated to exploring Feature Importance – Shapley Values for the PHACTboost Algorithm.

In Figure 2C, we present a comparative analysis between PHACTboost and MSAboost, both employing the same machine learning algorithm but differing in their feature sets. This comparison was conducted across the entire dataset and included all alternative test sets.

Figure 2A and 2C are plotted by using "Figure_2AC.R".
The code to plot 2C is names as "Figure_2B.R".

The data given in Figure2/Data folder is produced by using "Figure2_ProduceData.R".

## Figure 3

Figure 3A, 3B, 3C, 3D, 3E show AUROC comparisons of PHACTboost and tools reported in dbNSFP across TS1, TS2, TS3, TS4 and TS5, respectively. 

Part A, B, C, D and E can be plotted by using "Figure3_Comparison_TS.R".
The required files are given in the same folder.

The data given in Figure3/ folder is produced by using "Figure3_ProduceData.R".

## Figure 4

A. AUROC and AUPR difference between PHACTboost and EVE, CPT-1, VARITY, gMVP and ENVISION over the whole test set. B. AUROC and AUPR difference between PHACTboost and EVE, CPT-1, VARITY, gMVP and ENVISION over hard cases. C. Protein-by-protein comparison of PHACTboost with EVE over 60 proteins with at least 3 neutral and 3 pathogenic variations. D. Comparison of average Spearman correlation between PHACTboost and DMS results with the average Spearman correlation between other tools and DMS data.




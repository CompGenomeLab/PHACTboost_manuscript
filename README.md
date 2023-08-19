# PHACTboost_manuscript

This repository includes the required codes to reproduce the figures at ...

## Figure 2

This figure shows the performance difference between PHACTboost, PHACT and MSAboost.

A. The performance of the PHACTboost algorithm was evaluated over the complete set of variations in the test set. B. Comparison of PHACTboost results over full test set, test set excluding same positions with the training set, and test set excluding the same proteins as training set. C. Feature Importance – Shapley Values for PHACTboost algorithm. D. Comparison of PHACTboost with MSAboost which corresponds to the same machine learning algorithm with only MSA-based features over the whole dataset. E. Comparison of PHACTboost with MSAboost over the hard cases.

## Figure 3

Figure 3 shows AUPR comparisons of PHACTboost and tools reported in dbNSFP across TS1, TS2, TS3, TS4 and TS5, respectively. 

Part A, B, C, D and E can be plotted by using "Figure3_Comparison_TS.R".
The required files are given in the same folder.

The data given in Figure3/ folder is produced by using "Figure3_ProduceData.R".

## Figure 4

A. AUROC and AUPR difference between PHACTboost and EVE, CPT-1, VARITY, gMVP and ENVISION over the whole test set. B. AUROC and AUPR difference between PHACTboost and EVE, CPT-1, VARITY, gMVP and ENVISION over hard cases. C. Protein-by-protein comparison of PHACTboost with EVE over 60 proteins with at least 3 neutral and 3 pathogenic variations. D. Comparison of average Spearman correlation between PHACTboost and DMS results with the average Spearman correlation between other tools and DMS data.




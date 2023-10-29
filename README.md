# PHACTboost_manuscript

This repository contains the necessary codes to reproduce the figures presented in the PHACTboost manuscript. Each figure's code and related data are provided in their respective named folders.

## Figure 1

Figure 1A illustrates the evaluation of the PHACTboost algorithm's performance across all alternative test sets, revealing the enhancements achieved over the baseline PHACT model.

Figure 1B is dedicated to exploring Feature Importance – Shapley Values for the PHACTboost Algorithm.

In Figure 1C, we present a comparative analysis between PHACTboost and MSAboost, both employing the same machine learning algorithm but differing in their feature sets. This comparison was conducted across the entire dataset and encompassed all alternative test sets.

Both Figure 1A and 1C are generated using the script "Figure_1AC.R". The script for plotting Figure 1B is named "Figure_1B.R".

The data provided in the [Figure1/Data](Figure1/Data) folder was generated using the script "Figure1_ProduceData.R".

## Figure 2

Figure 2A displays AUROC comparisons between PHACTboost and tools reported in dbNSFP across TS1, TS2, TS3, TS4, and TS5, respectively.

To generate Part A and B, you can use the script "Figure2A.R," and "Figure2B.R," which are located in the same directory, respectively.

The necessary files are provided in the Figure2/Data folder.

The data found in the [Figure2/Data](Figure2/Data) folder was generated using the script "Figure2_ProduceData.R.

## Figure 3

Figure 3A demonstrates the AUROC difference between PHACTboost and EVE across all alternative test sets. Similarly, Figure 3B displays the AUROC levels of PHACTboost and CPT-1 across all alternative test sets.

In Part C, we present a protein-by-protein comparison of PHACTboost with EVE over 60 proteins, each containing at least 3 neutral and 3 pathogenic variations.

Finally, Part D is dedicated to comparing the average Spearman correlation between PHACTboost and DMS results with the average Spearman correlation between other tools and DMS data.

The data found in the [Figure3/Data](Figure3/Data) folder was generated using the script "Figure3_ProduceData.R.



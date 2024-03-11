##### FIND THE COORDINATES OF CONFLICTING INTERP. VARIANTS

library(data.table)
library(stringr)
library(readxl)
library(dplyr)

clinvar_bef <- "/Users/nurdankuru/Desktop/PHACTboost_Revision1/variant_summary_2023-01.txt"

clinvar_bef <- fread(file = clinvar_bef, sep = "\t", header = T, data.table = F)
clinvar_bef_38 <- clinvar_bef[clinvar_bef$Assembly == "GRCh38", ]
clinvar_bef_38 <- clinvar_bef_38[clinvar_bef_38$Type == "single nucleotide variant", ]
clinvar_bef_38$chr_vars <- paste0(clinvar_bef_38$Chromosome, "-", clinvar_bef_38$Start, clinvar_bef_38$ReferenceAlleleVCF, ">", clinvar_bef_38$AlternateAlleleVCF)

clinvar_aft <- "/Users/nurdankuru/Desktop/PHACTboost_Revision1/variant_summary.txt"

clinvar_aft <- fread(file = clinvar_aft, sep = "\t", header = T, data.table = F)
clinvar_aft_38 <- clinvar_aft[clinvar_aft$Assembly == "GRCh38", ]
clinvar_aft_38 <- clinvar_aft_38[clinvar_aft_38$Type == "single nucleotide variant", ]
clinvar_aft_38$chr_vars <- paste0(clinvar_aft_38$Chromosome, "-", clinvar_aft_38$Start, clinvar_aft_38$ReferenceAlleleVCF, ">", clinvar_aft_38$AlternateAlleleVCF)

rm(clinvar_aft, clinvar_bef)

label1 <- "Pathogenic"
label2 <- "Benign"
label3 <- "Likely pathogenic"
label4 <- "Pathogenic/Likely pathogenic"
label5 <- "Benign/Likely benign"
label6 <- "Likely benign"

ll <- c(which(clinvar_bef_38$ClinicalSignificance==label1), which(clinvar_bef_38$ClinicalSignificance==label2),
        which(clinvar_bef_38$ClinicalSignificance==label3), which(clinvar_bef_38$ClinicalSignificance==label4),
        which(clinvar_bef_38$ClinicalSignificance==label5), which(clinvar_bef_38$ClinicalSignificance==label6))

clinvar_bef_38 <- clinvar_bef_38[ll,]

loc <- grep("onflicting", clinvar_aft_38$ClinicalSignificance)
clinvar_aft_38 <- clinvar_aft_38[loc,]

com <- intersect(clinvar_aft_38$chr_vars, clinvar_bef_38$chr_vars)
clinvar_aft_38 <- clinvar_aft_38[match(com, clinvar_aft_38$chr_vars),]
clinvar_bef_38 <- clinvar_bef_38[match(com, clinvar_bef_38$chr_vars),]

chosen_vars <- clinvar_aft_38$chr_vars

vect <- c()
for (i in 1:length(clinvar_aft_38$chr_vars)){
  var <- clinvar_aft_38$chr_vars[i]
  var <- unlist(strsplit(var, "-"))
  chr <- var[1]
  var2 <- unlist(strsplit(var[2], ">"))
  alt <- var2[2]
  loc <- substr(var2[1], 1, (nchar(var2[1])-1))
  ref <- substr(var2[1], nchar(var2[1]), nchar(var2[1]))
  
  vect <- rbind(vect, c(chr, loc, ref, alt))
}

write.table(vect, quote = F, col.names = F, row.names = F, "ClinVar_VUS_Coordinates.txt")


##### Our Code used 


















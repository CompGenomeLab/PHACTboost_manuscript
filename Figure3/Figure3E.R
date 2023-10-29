library(data.table)
library(stringr)
library(readxl)
library(dplyr)

clinvar_data <- "variant_summary.txt"
clinvar_data <- fread(file = clinvar_data, sep = "\t", header = T, data.table = F)
clinvar_data_38 <- clinvar_data[clinvar_data$Assembly == "GRCh38", ]
clinvar_data_38 <- clinvar_data_38[clinvar_data_38$Type == "single nucleotide variant", ]
clinvar_data_38$chr_vars <- paste0(clinvar_data_38$Chromosome, "-", clinvar_data_38$Start, clinvar_data_38$ReferenceAlleleVCF, ">", clinvar_data_38$AlternateAlleleVCF)

load("TestSet.RData")
load("lightgbm_replication_1_prediction.RData")
test <- cbind(prediction$test_prediction, test)

test <- test[which(test$variant_info==1),]

db <- grep(";", test$chr_vars)
test2 <- test[db,]
test <- test[-db,]

guideline_vec1 <- clinvar_data_38$Guidelines[match(test$chr_vars, clinvar_data_38$chr_vars)]

guideline_vec2 <- c()
for (var in test2$chr_vars){
  chr_var <- unlist(strsplit(var, ";"))
  gui <- c()
  for (i in chr_var){
    gui <- c(gui, clinvar_data_38$Guidelines[which(clinvar_data_38$chr_vars==i)])
  }
  
  if (length(unique(gui))>1){
    gui <- paste(gui, collapse = "-")
  } else {
    gui <- unique(gui)
  }
  guideline_vec2 <- c(guideline_vec2, gui)
}

test$Guideline <- guideline_vec1
test2$Guideline <- guideline_vec2

test <- rbind(test, test2)

notfound_pat <- intersect(which(test$variant_info==1), which(test$`prediction$test_prediction`<0.5))
test_nf <- test[notfound_pat,]

data1 <- as.numeric(table(test$Guideline))
data2 <-  as.numeric(table(test_nf$Guideline))

n1 <- length(which(test$Guideline!="-"))
n2 <- length(which(test$Guideline=="-"))
m1 <- length(which(test_nf$Guideline!="-"))
m2 <- length(which(test_nf$Guideline=="-"))

# Create a contingency table
data <- matrix(c(n1, n2, m1, m2), nrow = 2, byrow = TRUE)
colnames(data) <- c("With_Guidelines", "Without_Guidelines")
rownames(data) <- c("Main_Dataset", "Subset")

# Chi-squared test
chi_sq_result <- chisq.test(data)
cat("Chi-Squared Test:\n")
print(chi_sq_result)






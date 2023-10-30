library(data.table)
library(ggplot2)
library(dplyr)
library(AUC)
library(PRROC)
library(wesanderson)
library(ggpubr)
library(readxl)
library(grid)
save_path <- "./"
data_path <- "../PHACTboost_Model/"
my_theme1 =  theme_bw() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # no gridlines
  theme(panel.border = element_blank(), axis.line = element_line(color = "black"))+ # no border, just major axis
  theme(axis.text = element_text(color = "black", size = 11), axis.title = element_text(size = 11)) + # axis text size and color
  theme(plot.title = element_text(size = 11, face = "bold")) + 
  theme(text = element_text(size = 11))+
  theme(plot.subtitle = element_text(size=10)) + # subtitle size
  theme(plot.caption=element_text(size=8)) + # caption size, if present
  theme(legend.text=element_text(size=10), legend.position = "right", legend.title = element_blank()) 

my_theme2 <- theme_void() +
  theme(axis.text.x = element_text(color = "black", size = 11, margin = margin(r = 6)),
        axis.text.y = element_text(color = "black", size = 11, hjust = 1, margin = margin(r = 6)),#
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(color = "black", size = 11),
        plot.subtitle = element_text(color = "black", size = 11),
        panel.grid.major.y = element_line(color = "grey90", size = 0.6),
        plot.background = element_rect(fill = "white", color = "white"), plot.margin = margin(c(0,0,5,0)))

load(sprintf("%s/TrainingSet.RData", data_path))
load(sprintf("%s/TestSet.RData", data_path))
dbnsfp_testset <- fread(sprintf("%s/dbNSFP_44a_TestSet.csv", data_path))

train_set <- train
test_set <- test
train_set$variant_info <- as.numeric(train_set$variant_info)
test_set$variant_info <- as.numeric(test_set$variant_info)

dbnsfp_data <- dbnsfp_testset

prediction <- fread(sprintf("%s/PHACTboost_TestPrediction.csv", data_path))
test_set$PHACTboost <- prediction$PHACTboost_Score
test_set$chr_var_with_aa <- paste0(test_set$chr_vars, "-", test_set$Ref_AA, test_set$Alt_AA)

test_set <- merge(test_set, dbnsfp_data[, c("chr_var_with_aa", colnames(dbnsfp_data)[grepl("rankscore", colnames(dbnsfp_data))])], by = "chr_var_with_aa", all.x = T)

selected_methods <- colnames(dbnsfp_data)[grep("rankscore", colnames(dbnsfp_data))]

pmv_train <- train_set %>% group_by(UNIPROTKB) %>% 
  summarise(Pathogenic = sum(variant_info == 1), 
            Neutral = sum(variant_info == -1)) %>% as.data.frame()

pmv_train$MV <- pmv_train$Pathogenic / rowSums(pmv_train[,-1])
pmv_test <- test_set[, c("UNIPROTKB","chr_var_with_aa", "prot_vars", "Positions", "Ref_AA", "Alt_AA", "variant_info", "PHACTboost", selected_methods)]
pmv_test$MV <- pmv_train[match(pmv_test$UNIPROTKB, pmv_train$UNIPROTKB), "MV"]
pmv_test[is.na(pmv_test$MV), "MV"] <- 0.5

test_set <- merge(test_set, pmv_test[, c("chr_var_with_aa", "MV")], by = "chr_var_with_aa", all.x = T)

pure_pathogenic_proteins <- unique(pmv_train[pmv_train$MV == 1, "UNIPROTKB"])
variants_in_pure_pathogenic_proteins <- test_set[test_set$UNIPROTKB %in% pure_pathogenic_proteins, ]
dim(variants_in_pure_pathogenic_proteins)
pure_neutral_proteins <- unique(pmv_train[pmv_train$MV == 0, "UNIPROTKB"])
variants_in_pure_neutral_proteins <- test_set[test_set$UNIPROTKB %in% pure_neutral_proteins, ]
dim(variants_in_pure_neutral_proteins)
mixed_proteins <- unique(pmv_train[pmv_train$MV != 1 & pmv_train$MV != 0, "UNIPROTKB"])

pure_test_set <- pmv_test[pmv_test$UNIPROTKB %in% c(pure_neutral_proteins, pure_pathogenic_proteins), ]
mixed_test_set <- pmv_test[(pmv_test$UNIPROTKB %in% c(mixed_proteins)), ]

pmv_0.1_0.9 <- pmv_test[(pmv_test$MV >= 0.1 & pmv_test$MV <= 0.9), ]
pmv_0.2_0.8 <- pmv_test[(pmv_test$MV >= 0.2 & pmv_test$MV <= 0.8), ]
pmv_0.3_0.7 <- pmv_test[(pmv_test$MV >= 0.3 & pmv_test$MV <= 0.7), ]
pmv_0.4_0.6 <- pmv_test[(pmv_test$MV >= 0.4 & pmv_test$MV <= 0.6), ]

method_colors <- c(wes_palettes$GrandBudapest1, wes_palettes$GrandBudapest2, wes_palettes$Moonrise1, wes_palettes$Moonrise2, wes_palettes$Darjeeling1, wes_palettes$Darjeeling2,
                   wes_palettes$BottleRocket1, wes_palettes$BottleRocket2, wes_palettes$Royal1, wes_palettes$Royal2, wes_palettes$FantasticFox1, wes_palettes$Rushmore)

methods <- selected_methods
methods <- methods[!grepl("BayesDel|ClinPred|MetaRNN|MetaSVM|MetaLR|Eigen|LINSIGHT|M-CAP|MutPred", methods)]
methods <- c(methods, "MV")

results_MV <- data.frame(Data_MV = NULL, Method = NULL, Delta_AUC = NULL, Delta_PR = NULL)
data_size <- data.frame(Data_MV = NULL, Method = NULL, N_Protein = NULL, Neutral = NULL, Pathogenic = NULL)
supplementary_table_pmv <- data.frame(Dataset = NULL, Tool = NULL, N_Protein = NULL, N_Neutral = NULL, N_Pathogenic = NULL, AUROC_Tool = NULL, 	AUROC_PHACTboost = NULL,	AUPR_Tool = NULL,	AUPR_PHACTboost = NULL)
for (method in methods) {
  method_scores <- as.numeric(test_set[[method]])
  method_vars <- test_set[, "variant_info"]
  method_proteins <- test_set$UNIPROTKB
  pb_scores <- as.numeric(test_set[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Data_MV = "All", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size <- rbind(data_size, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  dat0 <- data.frame(Data_MV = "All", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "All", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_pmv <- rbind(supplementary_table_pmv, df_res)
  
  method_scores <- as.numeric(pure_test_set[[method]])
  method_vars <- pure_test_set[, "variant_info"]
  method_proteins <- pure_test_set[, "UNIPROTKB"]
  
  pb_scores <- as.numeric(pure_test_set[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Data_MV = "Pure", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size <- rbind(data_size, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat1 <- data.frame(Data_MV = "Pure", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "Pure", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_pmv <- rbind(supplementary_table_pmv, df_res)
  
  method_scores <- as.numeric(mixed_test_set[[method]])
  method_vars <- mixed_test_set[, "variant_info"]
  method_proteins <- mixed_test_set[, "UNIPROTKB"]
  pb_scores <- as.numeric(mixed_test_set[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Data_MV = "Mixed", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size <- rbind(data_size, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat2 <- data.frame(Data_MV = "Mixed", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "Mixed", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_pmv <- rbind(supplementary_table_pmv, df_res)
  
  method_scores <- as.numeric(pmv_0.1_0.9[[method]])
  method_vars <- pmv_0.1_0.9[, "variant_info"]
  method_proteins <- pmv_0.1_0.9[, "UNIPROTKB"]
  pb_scores <- as.numeric(pmv_0.1_0.9[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Data_MV = "0.1_0.9", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size <- rbind(data_size, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat3 <- data.frame(Data_MV = "0.1_0.9", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "0.1_0.9", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_pmv <- rbind(supplementary_table_pmv, df_res)
  
  method_scores <- as.numeric(pmv_0.2_0.8[[method]])
  method_vars <- pmv_0.2_0.8[, "variant_info"]
  method_proteins <- pmv_0.2_0.8[, "UNIPROTKB"]
  pb_scores <- as.numeric(pmv_0.2_0.8[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Data_MV = "0.2_0.8", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size <- rbind(data_size, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat4 <- data.frame(Data_MV = "0.2_0.8", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "0.2_0.8", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_pmv <- rbind(supplementary_table_pmv, df_res)
  
  method_scores <- as.numeric(pmv_0.3_0.7[[method]])
  method_vars <- pmv_0.3_0.7[, "variant_info"]
  method_proteins <- pmv_0.3_0.7[, "UNIPROTKB"]
  pb_scores <- as.numeric(pmv_0.3_0.7[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Data_MV = "0.3_0.7", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size <- rbind(data_size, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat5 <- data.frame(Data_MV = "0.3_0.7", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "0.3_0.7", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_pmv <- rbind(supplementary_table_pmv, df_res)
  
  method_scores <- as.numeric(pmv_0.4_0.6[[method]])
  method_vars <- pmv_0.4_0.6[, "variant_info"]
  method_proteins <- pmv_0.4_0.6[, "UNIPROTKB"]
  pb_scores <- as.numeric(pmv_0.4_0.6[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Data_MV = "0.4_0.6", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size <- rbind(data_size, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat6 <- data.frame(Data_MV = "0.4_0.6", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "0.4_0.6", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_pmv <- rbind(supplementary_table_pmv, df_res)
  
  results_MV <- rbind(results_MV, dat0, dat1, dat2, dat3, dat4, dat5, dat6)
}

results_MV$Data_MV <- factor(results_MV$Data_MV, levels = c("All", "Pure", "Mixed", "0.1_0.9", "0.2_0.8", "0.3_0.7", "0.4_0.6"))

sel_methods_MV <- results_MV[results_MV$Data_MV == "All", "Method"]
sel_methods_MV <- sel_methods_MV[order(results_MV[results_MV$Data_MV == "All", "Delta_AUC"], decreasing = T)]

results_MV2 <- results_MV[results_MV$Method %in% c(sel_methods_MV), ]

sel_methods_MV <- gsub("_rankscore|_raw_rankscore|_converted_rankscore", "", sel_methods_MV)
results_MV2$Method <- gsub("_rankscore|_raw_rankscore|_converted_rankscore", "", results_MV2$Method)
results_MV2$Method <- factor(results_MV2$Method, levels = c(sel_methods_MV))

TS1 <- fread(sprintf("%s/TS1.csv", data_path))
TS2 <- fread(sprintf("%s/TS2.csv", data_path))
TS3 <- fread(sprintf("%s/TS3.csv", data_path))
TS4 <- fread(sprintf("%s/TS4.csv", data_path))
TS5 <- fread(sprintf("%s/TS5.csv", data_path))

training_proteins <- unique(train_set$UNIPROTKB)
training_positions <- paste0(train_set$UNIPROTKB, "-", train_set$Positions)
ts1 <- test_set[test_set$prot_vars %in% TS1$prot_vars, ]
ts2 <- ts1[(ts1$prot_vars %in% TS2$prot_vars), ]
ts3 <- ts1[(ts1$prot_vars %in% TS3$prot_vars), ]
ts4 <- ts1[(ts1$prot_vars %in% TS4$prot_vars), ]
ts5 <- ts1[(ts1$prot_vars %in% TS5$prot_vars), ]

hard_cases <- ts4
hard_cases_mixed <- ts5

pure_test_set <- pure_test_set[pure_test_set$prot_vars %in% hard_cases$prot_vars, ]
mixed_test_set <- mixed_test_set[mixed_test_set$prot_vars %in% hard_cases$prot_vars, ]
pmv_0.1_0.9 <- pmv_0.1_0.9[pmv_0.1_0.9$prot_vars %in% hard_cases$prot_vars, ]
pmv_0.2_0.8 <- pmv_0.2_0.8[pmv_0.2_0.8$prot_vars %in% hard_cases$prot_vars, ]
pmv_0.3_0.7 <- pmv_0.3_0.7[pmv_0.3_0.7$prot_vars %in% hard_cases$prot_vars, ]
pmv_0.4_0.6 <- pmv_0.4_0.6[pmv_0.4_0.6$prot_vars %in% hard_cases$prot_vars, ]

methods <- selected_methods
methods <- methods[!grepl("BayesDel|ClinPred|MetaRNN|MetaSVM|MetaLR|Eigen|LINSIGHT", methods)]
methods <- methods[!grepl("M-CAP_rankscore|MutPred_rankscore", methods)]
methods <- c(methods, "MV")

results_dbnsfp <- data.frame(Dataset = NULL, Method = NULL, Delta_AUC = NULL, Delta_PR = NULL)
data_size_dbnsfp <- data.frame(Dataset = NULL, Method = NULL, N_Protein = NULL, Neutral = NULL, Pathogenic = NULL)
supplementary_table_dbnsfp <- data.frame(Dataset = NULL, Tool = NULL, N_Protein = NULL, N_Neutral = NULL, N_Pathogenic = NULL, AUROC_Tool = NULL, 	AUROC_PHACTboost = NULL,	AUPR_Tool = NULL,	AUPR_PHACTboost = NULL)

for(method in methods){
  
  method_scores <- as.numeric(ts1[[method]])
  method_vars <- ts1[, "variant_info"]
  method_proteins <- gsub("-.*", "", ts1$prot_vars)
  pb_scores <- as.numeric(ts1[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Dataset = "TS1", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size_dbnsfp<- rbind(data_size_dbnsfp, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  dat0 <- data.frame(Dataset = "TS1", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "TS1", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_dbnsfp <- rbind(supplementary_table_dbnsfp, df_res)
  
  method_scores <- as.numeric(ts2[[method]])
  method_vars <- ts2[, "variant_info"]
  method_proteins <- ts2[, "UNIPROTKB"]
  
  pb_scores <- as.numeric(ts2[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Dataset = "TS2", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size_dbnsfp <- rbind(data_size_dbnsfp, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat1 <- data.frame(Dataset = "TS2", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "TS2", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_dbnsfp <- rbind(supplementary_table_dbnsfp, df_res)
  
  method_scores <- as.numeric(ts3[[method]])
  method_vars <- ts3[, "variant_info"]
  method_proteins <- ts3[, "UNIPROTKB"]
  
  pb_scores <- as.numeric(ts3[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Dataset = "TS3", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size_dbnsfp <- rbind(data_size_dbnsfp, df)
  
  if(method == "MV") {
    method_auc <- 0.5
    v <- 1 * (method_vars[elims] == 1)
    method_pr <- sum(v == 1) / length(v)
  }else{
    method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
    method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  }
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat2 <- data.frame(Dataset = "TS3", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "TS3", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_dbnsfp <- rbind(supplementary_table_dbnsfp, df_res)
  
  method_scores <- as.numeric(ts4[[method]])
  method_vars <- ts4[, "variant_info"]
  method_proteins <- ts4[, "UNIPROTKB"]
  pb_scores <- as.numeric(ts4[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Dataset = "TS4", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size_dbnsfp <- rbind(data_size_dbnsfp, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  if(length(unique(method_scores[elims])) != 1) {
    method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  }else{
    method_pr <- var_table[2] / sum(var_table)
  }
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat3 <- data.frame(Dataset = "TS4", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "TS4", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_dbnsfp <- rbind(supplementary_table_dbnsfp, df_res)
  
  method_scores <- as.numeric(ts5[[method]])
  method_vars <- ts5[, "variant_info"]
  method_proteins <- ts5[, "UNIPROTKB"]
  pb_scores <- as.numeric(ts5[, "PHACTboost"])
  
  elims <- !(is.na(as.numeric(method_scores)))
  var_table <- table(method_vars[elims])
  prot_table <- method_proteins[elims]
  df <- data.frame(Dataset = "TS5", Method = method, N_protein = length(unique(prot_table)), Neutral = var_table[1], Pathogenic = var_table[2])
  data_size_dbnsfp <- rbind(data_size_dbnsfp, df)
  
  method_auc <- auc(roc(method_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  method_pr <- pr.curve(scores.class0 = method_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  pb_auc <- auc(roc(pb_scores[elims], as.factor(1 * (method_vars[elims] == 1))))
  pb_pr <- pr.curve(scores.class0 = pb_scores[elims], weights.class0 = (1 * (method_vars[elims] == 1)), curve =T)$auc.integral
  
  dat4 <- data.frame(Dataset = "TS5", Method = method, Delta_AUC = method_auc - pb_auc, Delta_PR = method_pr - pb_pr)
  df_res <- data.frame(Dataset = "TS5", Tool = method, N_Protein = length(unique(prot_table)), 
                       N_Neutral = var_table[1], N_Pathogenic = var_table[2], AUROC_Tool = method_auc, 	
                       AUROC_PHACTboost = pb_auc,	AUPR_Tool = method_pr,	AUPR_PHACTboost = pb_pr)
  supplementary_table_dbnsfp <- rbind(supplementary_table_dbnsfp, df_res)
  
  results_dbnsfp <- rbind(results_dbnsfp, dat0, dat1, dat2, dat3, dat4)
}

results_dbnsfp$Dataset <- factor(results_dbnsfp$Dataset, levels = c("TS1", "TS2", "TS3", "TS4", "TS5"))

sel_methods_dbnsfp <- results_dbnsfp[results_dbnsfp$Dataset == "TS1", "Method"]
sel_methods_dbnsfp <- sel_methods_dbnsfp[order(results_dbnsfp[results_dbnsfp$Dataset == "TS1", "Delta_AUC"], decreasing = T)]

results_dbnsfp2 <- results_dbnsfp[results_dbnsfp$Method %in% c(sel_methods_dbnsfp), ]

sel_methods_dbnsfp <- gsub("_rankscore|_raw_rankscore|_converted_rankscore", "", sel_methods_dbnsfp)
results_dbnsfp2$Method <- gsub("_rankscore|_raw_rankscore|_converted_rankscore", "", results_dbnsfp2$Method)
results_dbnsfp2$Method <- factor(results_dbnsfp2$Method, levels = c(sel_methods_dbnsfp))

sel_methods_pmv <- sel_methods_MV[!grepl("M-CAP|MutPred", sel_methods_MV)]

sel_methods_dbnsfp <- sel_methods_dbnsfp[!grepl("M-CAP|MutPred", sel_methods_dbnsfp)]

datasets <- c("TS1", "TS2", "TS3", "TS4", "TS5")
results_MV2 <- results_MV2[results_MV2$Method %in% sel_methods_pmv, ]
results_dbnsfp2 <- results_dbnsfp2[results_dbnsfp2$Method %in% sel_methods_dbnsfp, ]

names(method_colors) <- sel_methods_pmv
method_colors <- method_colors[sel_methods_pmv]
results_MV2$Method <- factor(results_MV2$Method, levels = sel_methods_pmv)
results_dbnsfp2$Method <- factor(results_dbnsfp2$Method, levels = sel_methods_dbnsfp)

ensembl_info <- read_xlsx(path = sprintf("%s/methods_info.xlsx", data_path))
ensembl_info <- as.data.frame(ensembl_info)
ensembl_info[is.na(ensembl_info)] <- 0

results_dbnsfp3 <- merge(results_dbnsfp2, ensembl_info, by = "Method", all.x = T)

unique_methods <- c("VEST4", "REVEL", "gMVP", "DEOGEN2", "MVP", "VARITY_R", "CADD", "LIST-S2", "fathmm-XF_coding",
                    "MutationTaster", "PrimateAI", "phyloP100way_vertebrate", "Polyphen2_HVAR", "FATHMM", "SIFT", 
                    "MutationAssessor", "PROVEAN", "MV", "LRT", "SiPhy_29way_logOdds", "MPC", "GERP++_RS", "DANN",
                    "phastCons100way_vertebrate", "GenoCanyon", "integrated_fitCons", "H1-hESC_fitCons", "bStatistic")
################
####plots#######
################
results_MV3 <- results_MV2[results_MV2$Data_MV != "All", ]
results_MV3 <- results_MV3[results_MV3$Method %in% unique_methods, ]
results_MV3$Method <- gsub("100way|470way|17way","", results_MV3$Method)
results_MV3$Method <- gsub("SiPhy_29way_logOdds","SiPhy", results_MV3$Method)
results_MV3$Method <- gsub("_coding","", results_MV3$Method)
mets <- as.vector(results_MV3[results_MV3$Data_MV == "0.4_0.6", "Method"])
names(method_colors) <- mets[order(results_MV3[results_MV3$Data_MV == "0.4_0.6", "Delta_AUC"], decreasing = T)]
method_colors[-c(1:10)] <- "#5A5A5A"

results_MV3$Data_MV <- gsub("_", "-",results_MV3$Data_MV)
results_MV3$Data_MV <- factor(results_MV3$Data_MV, levels = c("Pure", "Mixed", "0.1-0.9", "0.2-0.8", "0.3-0.7", "0.4-0.6"))

legend_methods1 <- results_MV3[results_MV3$Data_MV == "Pure", "Method"] 
legend_methods1 <- as.vector(legend_methods1[order(results_MV3[results_MV3$Data_MV == "Pure", "Delta_AUC"], decreasing = T)])

legend_methods2 <- results_MV3[results_MV3$Data_MV == "0.4-0.6", "Method"] 
legend_methods2 <- as.vector(legend_methods2[order(results_MV3[results_MV3$Data_MV == "0.4-0.6", "Delta_AUC"], decreasing = T)])
lim1 <- 0.125
lim2 <- -0.575
seq_legend <- seq(lim1,lim2, length.out = length(unique_methods))

lines1 <- data.frame(x1 = rep(0.60, length(unique_methods)), x2 = rep(0.95, length(unique_methods)), y1 = seq_legend, y2 = sort(results_MV3[results_MV3$Data_MV == "Pure", "Delta_AUC"], decreasing = T))
lines2 <- data.frame(x1 = rep(6.05, length(unique_methods)), x2 = rep(6.40, length(unique_methods)), y1 = sort(results_MV3[results_MV3$Data_MV == "0.4-0.6", "Delta_AUC"], decreasing = T), y2 = seq_legend)

pmv_auc <- ggplot(data=results_MV3, aes(x = Data_MV, y = Delta_AUC, group = Method, color = Method)) + my_theme2 +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") + geom_line(lwd = 0.5, position = position_identity())  + geom_point(size = 1, position = position_identity()) + 
  scale_y_continuous(limits = c(lim2,lim1)) +
  scale_x_discrete(breaks = c("Pure", "Mixed", "0.1-0.9", "0.2-0.8", "0.3-0.7", "0.4-0.6"), limits = c("Pure", "Mixed", "0.1-0.9", "0.2-0.8", "0.3-0.7", "0.4-0.6"),
                   expand = expansion(mult = c(0, 0.7))) +
  scale_color_manual(values = method_colors) + labs(y = paste0("\u0394", "AUROC\n"), x = "") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(angle = 90), plot.title = element_text(hjust = 0.5),# axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        legend.text = element_text(size = 11), axis.text.x = element_text(angle = 45), plot.margin = margin(1,0,0,0, "cm")) + 
  annotate("text", x = 6.5, y = seq_legend, color = method_colors[legend_methods2], size = 3,
           label =  as.vector(legend_methods2),
           hjust = 0) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), linetype=2, lwd = 0.1, color = method_colors[legend_methods2], data = lines2, inherit.aes = F) + 
  labs(title = "TS1 Subsets")


legend_methods1 <- results_MV3[results_MV3$Data_MV == "Pure", "Method"] 
legend_methods1 <- as.vector(legend_methods1[order(results_MV3[results_MV3$Data_MV == "Pure", "Delta_PR"], decreasing = T)])

legend_methods2 <- results_MV3[results_MV3$Data_MV == "0.4-0.6", "Method"] 
legend_methods2 <- as.vector(legend_methods2[order(results_MV3[results_MV3$Data_MV == "0.4-0.6", "Delta_PR"], decreasing = T)])
lim1 <- 0.125
lim2 <- -0.7
seq_legend <- seq(lim1,lim2, length.out = length(unique_methods))

lines1 <- data.frame(x1 = rep(0.60, length(unique_methods)), x2 = rep(0.95, length(unique_methods)), y1 = seq_legend, y2 = sort(results_MV3[results_MV3$Data_MV == "Pure", "Delta_PR"], decreasing = T))
lines2 <- data.frame(x1 = rep(6.05, length(unique_methods)), x2 = rep(6.40, length(unique_methods)), y1 = sort(results_MV3[results_MV3$Data_MV == "0.4-0.6", "Delta_PR"], decreasing = T), y2 = seq_legend)

pmv_aupr <- ggplot(data=results_MV3, aes(x = Data_MV, y = Delta_PR, group = Method, color = Method)) + my_theme2 +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") + geom_line(lwd = 0.5, position = position_identity())  + geom_point(size = 1, position = position_identity()) + 
  scale_y_continuous(limits = c(lim2,lim1)) +
  scale_x_discrete(breaks = c("Pure", "Mixed", "0.1-0.9", "0.2-0.8", "0.3-0.7", "0.4-0.6"), limits = c("Pure", "Mixed", "0.1-0.9", "0.2-0.8", "0.3-0.7", "0.4-0.6"),
                   expand = expansion(mult = c(0, 0.7))) +
  scale_color_manual(values = method_colors) + labs(y = paste0("\u0394", "AUPR\n"), x = "") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(angle = 90), plot.title = element_text(hjust = 0.5),# axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        legend.text = element_text(size = 11), axis.text.x = element_text(angle = 45), plot.margin = margin(1,0,0,0, "cm")) +
  annotate("text", x = 6.5, y = seq_legend, color = method_colors[legend_methods2], 
           label =  as.vector(legend_methods2), size = 3,
           hjust = 0) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), linetype=2, lwd = 0.1, colour = method_colors[legend_methods2], data = lines2, inherit.aes = F) + 
  labs(title = "TS1 Subsets")

data1 <- results_dbnsfp3[results_dbnsfp3$Dataset == "TS1", ]
data1 <- data1[data1$Method %in% unique_methods, ]
data1$Method <- gsub("100way|470way|17way","", data1$Method)
data1$Method <- gsub("SiPhy_29way_logOdds","SiPhy", data1$Method)
data1$Method <- gsub("_coding","", data1$Method)
data1$Method <- factor(data1$Method, levels = as.vector(data1[order(data1[data1$Dataset == "TS1", "Delta_AUC"], decreasing = T), "Method"]))

gray_dark <- rgb(72/255,73/255,77/255)

g1 <- ggplot(data1, aes(x = Method, y = Delta_AUC)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data1$Method), width = 0.5) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 1)), breaks = c(-0.4, -0.2, 0), limits = c(min(data1$Delta_AUC)*1.05, max(data1$Delta_AUC, 0.0))) +
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none",  plot.margin = margin(1,0,0.2,0, "cm"), axis.title.y = element_text(angle = 90),
        axis.title.x = element_blank(),  axis.ticks.x = element_blank(), axis.line.x = element_blank(),  plot.title = element_text(hjust = 0.5),  
        axis.text.x.top = element_text(size = 8, angle = 90, vjust=0.5, hjust = 0, colour = method_colors[as.vector(data1$Method)])) +
  geom_text(aes(y = 0, group = Method, label = ifelse(`Meta-Predictor` == 1, "\u25CF", ""), vjust = -1.75), colour = method_colors[as.vector(data1$Method)], size = 3) +
  geom_text(aes(y = 0, group = Method, label = ifelse(`ML Method` == 1, "\u25CF", ""), vjust = -0.5), colour = method_colors[as.vector(data1$Method)], show.legend = T, size = 3) +
  annotate("text", x = 0.5, y = 0,
           label =  "Meta-Pred.", size = 3,
           hjust = 1, vjust = -1.75) +
  annotate("text", x = 0.5, y = 0,
           label =  "ML", size = 3,
           hjust = 1, vjust = -0.5) +
  labs(title = "TS1", y = "", x = "") 
g1

data2 <- results_dbnsfp2[results_dbnsfp2$Dataset == "TS2", ]
data2 <- data2[data2$Method %in% unique_methods, ]

g2 <- ggplot(data2, aes(x = Method, y = Delta_AUC)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data2$Method), width = 0.5) + 
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0), limits = c(min(data2$Delta_AUC)*1.05, max(data2$Delta_AUC, 0))) + 
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none",  plot.title = element_text(hjust = 0.5), axis.title.y = element_text(angle = 90), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(title = "TS2", y = "", x = "")


data3 <- results_dbnsfp2[results_dbnsfp2$Dataset == "TS3", ]
data3 <- data3[data3$Method %in% unique_methods, ]

g3 <- ggplot(data3, aes(x = Method, y = Delta_AUC)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data3$Method), width = 0.5) + 
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0),limits = c(min(data3$Delta_AUC)*1.05, max(data3$Delta_AUC, 0)))  +
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),  axis.title.y = element_text(angle = 90), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(title = "TS3", y = "", x = "")

data4 <- results_dbnsfp2[results_dbnsfp2$Dataset == "TS4", ]
data4 <- data4[data4$Method %in% unique_methods, ]

g4 <- ggplot(data4, aes(x = Method, y = Delta_AUC)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data4$Method), width = 0.5) + 
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0), limits = c(min(data4$Delta_AUC)*1.05, max(data4$Delta_AUC, 0))) +
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.y = element_text(angle = 90),  
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(title = "TS4", y = "", x = "")

data5 <- results_dbnsfp2[results_dbnsfp2$Dataset == "TS5", ]
data5 <- data5[data5$Method %in% unique_methods, ]

g5 <- ggplot(data5, aes(x = Method, y = Delta_AUC)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data5$Method), width = 0.5) + 
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0), limits = c(min(data5$Delta_AUC)*1.05, max(data5$Delta_AUC, 0))) +
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none",  plot.title = element_text(hjust = 0.5), axis.title.y = element_text(angle = 90),  
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())+ 
  labs(title = "TS5", y = "", x = "")

r1 <- diff(c(min(data1$Delta_AUC)*1.05, max(data1$Delta_AUC, 0)))
r2 <- diff(c(min(data2$Delta_AUC)*1.05, max(data2$Delta_AUC, 0)))
r3 <- diff(c(min(data3$Delta_AUC)*1, max(data3$Delta_AUC, 0)))
r4 <- diff(c(min(data4$Delta_AUC)*0.9, max(data4$Delta_AUC, 0)))
r5 <- diff(c(min(data5$Delta_AUC)*0.9, max(data5$Delta_AUC, 0)))
dbnsfp_auc <- ggarrange(g1, g2, g3, g4, g5, labels = c("", "", "", "", ""), nrow = 5, ncol = 1, heights = c(r1*3.6, r2, r3, r4, r5)) 
annotate_figure(dbnsfp_auc, left = textGrob(paste0("\u0394", "AUROC"), rot = 90, vjust = 1))

space_plot <- ggplot() + theme_void()
figure2 <- ggarrange(dbnsfp_auc, space_plot, pmv_auc, labels = c("A", "", "B"), ncol = 3, widths = c(4.4, 0.1, 4.1))
figure2 <- annotate_figure(figure2, left = textGrob(paste0("\u0394", "AUROC"), rot = 90, vjust = 1, hjust = 0.62, gp = gpar(cex = 0.9)))

ggsave(sprintf("%s/Figure2_AB.pdf", save_path), figure2, width=8.5, height=5.5, device = cairo_pdf)

g6 <- ggplot(data1, aes(x = Method, y = Delta_PR)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data1$Method), width = 0.5) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 1)), breaks = c(-0.4, -0.2, 0), limits = c(min(data1$Delta_PR)*1.05, max(data1$Delta_PR, 0.0))) +
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none",  plot.margin = margin(1,0,0.2,0, "cm"), axis.title.y = element_text(angle = 90),
        axis.title.x = element_blank(),  axis.ticks.x = element_blank(), axis.line.x = element_blank(),  plot.title = element_text(hjust = 0.5),  
        axis.text.x.top = element_text(size = 8, angle = 90, vjust=0.5, hjust = 0, colour = method_colors[as.vector(data1$Method)])) +
  geom_text(aes(y = 0, group = Method, label = ifelse(`Meta-Predictor` == 1, "\u25CF", ""), vjust = -1.75), colour = method_colors[as.vector(data1$Method)], size = 3) +
  geom_text(aes(y = 0, group = Method, label = ifelse(`ML Method` == 1, "\u25CF", ""), vjust = -0.5), colour = method_colors[as.vector(data1$Method)], show.legend = T, size = 3) +
  annotate("text", x = 0.5, y = 0,
           label =  "Meta-Pred.", size = 3,
           hjust = 1, vjust = -1.75) +
  annotate("text", x = 0.5, y = 0,
           label =  "ML", size = 3,
           hjust = 1, vjust = -0.5) +
  labs(title = "TS1", y = "", x = "")


g7 <- ggplot(data2, aes(x = Method, y = Delta_PR)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data2$Method), width = 0.5) + 
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0), limits = c(min(data2$Delta_PR)*1.05, max(data2$Delta_PR, 0))) +
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none",  plot.title = element_text(hjust = 0.5), axis.title.y = element_text(angle = 90), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(title = "TS2", y = "", x = "")


g8 <- ggplot(data3, aes(x = Method, y = Delta_PR)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data3$Method), width = 0.5) + 
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0),limits = c(min(data3$Delta_PR)*1.05, max(data3$Delta_PR, 0))) +
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),  axis.title.y = element_text(angle = 90),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(title = "TS3", y = "", x = "")

g9 <- ggplot(data4, aes(x = Method, y = Delta_PR)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data4$Method), width = 0.5) + 
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0), limits = c(min(data4$Delta_PR)*1.05, max(data4$Delta_PR, 0))) +
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.y = element_text(angle = 90), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(title = "TS4", y = "", x = "")

g10 <- ggplot(data5, aes(x = Method, y = Delta_PR)) + my_theme2 +
  geom_bar(stat = "identity", aes(fill = data5$Method), width = 0.5) + 
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0), limits = c(min(data5$Delta_PR)*1.05, max(data5$Delta_PR, 0))) + 
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.25, 0))) + scale_fill_manual(values = method_colors) +
  theme(legend.position = "none",  plot.title = element_text(hjust = 0.5), axis.title.y = element_text(angle = 90), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())+ 
  labs(title = "TS5", y = "", x = "")

r6 <- diff(c(min(data1$Delta_PR)*1.05, max(data1$Delta_PR, 0)))
r7 <- diff(c(min(data2$Delta_PR)*1.05, max(data2$Delta_PR, 0)))
r8 <- diff(c(min(data3$Delta_PR)*0.9, max(data3$Delta_PR, 0)))
r9 <- diff(c(min(data4$Delta_PR)*0.9, max(data4$Delta_PR, 0)))
r10 <- diff(c(min(data5$Delta_PR)*0.9, max(data5$Delta_PR, 0)))
dbnsfp_pr <- ggarrange(g6, g7, g8, g9, g10, labels = c("A", "", "", "", ""), nrow = 5, ncol = 1, heights = c(r6*3.6, r7, r8, r9, r10))


space_plot <- ggplot() + theme_void()
supp_fig1 <- ggarrange(dbnsfp_pr, space_plot, pmv_aupr, labels = c("A", "", "B"), ncol = 3, widths = c(4.4, 0.1, 4.1))
supp_fig1 <- annotate_figure(supp_fig1, left = textGrob(paste0("\u0394", "AUPR"), rot = 90, vjust = 1, hjust = 0.62, gp = gpar(cex = 0.9)))

ggsave(sprintf("%s/Supplementary_Figure_1.pdf", save_path), supp_fig1, width=8.5, height=5.5, device = cairo_pdf)

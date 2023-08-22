library(PRROC)
library(AUC)

load("/Users/nurdankuru/Desktop/PHACTboost_August/PHACTboost_GnomADSHARED_TRAIN_with_phact_vars_with_0_CountNodes_3/TrainingSet.RData")
load("/Users/nurdankuru/Desktop/PHACTboost_August/PHACTboost_GnomADSHARED_TRAIN_with_phact_vars_with_0_CountNodes_3/TestSet.RData")
load("/Users/nurdankuru/Desktop/PHACTboost_August/PHACTboost_GnomADSHARED_TRAIN_with_phact_vars_with_0_CountNodes_3/lightgbm_replication_1_prediction.RData")
test <- cbind(prediction$test_prediction, test)
phactboost <- test
rm(prediction)

load("/Users/nurdankuru/Downloads/EVE_Scores.RData")
eve <- all
rm(all)

load("/Users/nurdankuru/Downloads/CPT1_EVEData.RData")
cpt1 <- full
load("/Users/nurdankuru/Downloads/CPT1_XGMData.RData")
cpt1 <- rbind(cpt1, full)
rm(full)

# TS3
var_train <- paste(train$UNIPROTKB, train$Positions, sep = "-")
var_test <- paste(test$UNIPROTKB, test$Positions, sep = "-")
inds_ts3 <- setdiff(var_test, var_train)
ind1 <- c()
for (id in inds_ts3){
  ind1 <- c(ind1, which(var_test==id))
}
test_unseenvar <- test[ind1, ]

# TS4
ids_ts4 <- setdiff(unique(test$UNIPROTKB), unique(train$UNIPROTKB))
ind2 <- c()
for (id in ids_ts4){
  ind2 <- c(ind2, which(test$UNIPROTKB==id))
}
test_unseenprot <- test[ind2, ]
rm(test, train)

load("/Users/nurdankuru/Hard_Cases.RData")
load("/Users/nurdankuru/Hard_Cases_Mixed.RData")

vars_all <- phactboost$prot_vars
vars_diffprots <- test_unseenprot$prot_vars
vars_diffvars <- unlist(test_unseenvar$prot_vars)
vars_hardcases1 <- hard_cases$prot_vars
vars_hardcases2 <- hard_cases2$prot_vars
rm(id, ids_ts4, ind1, ind2, inds_ts3, var_test, var_train, vars_all)
rm(hard_cases, hard_cases2, test_unseenprot, test_unseenvar)

data_main <- cbind(phactboost$variant_info, phactboost$`prediction$test_prediction`,
                   phactboost$Phact_wl_param_CountNodes_3_alt,
                   eve$EVE_scores_ASM[match(phactboost$prot_vars, eve$prot_vars)],
                   cpt1$CPT[match(phactboost$prot_vars, cpt1$prot_vars)])

data_main <- as.data.frame(data_main)
colnames(data_main) <- c("variant_info", "PHACTboost", "PHACT", "EVE", "CPT-1")

ts1 <- data_main
ts2 <- data_main[match(vars_hardcases1, data_main$prot_vars),]
ts3 <- data_main[match(vars_hardcases2, data_main$prot_vars),]
ts4 <- data_main[match(vars_diffprots, data_main$prot_vars),]
ts5 <- data_main[match(vars_diffvars, data_main$prot_vars),]

save(ts1, file = "TS1_EVE_CPT1.RData")
save(ts2, file = "TS2_EVE_CPT1.RData")
save(ts3, file = "TS3_EVE_CPT1.RData")
save(ts4, file = "TS4_EVE_CPT1.RData")
save(ts5, file = "TS5_EVE_CPT1.RData")

write.csv(ts1, quote = F, row.names = F, "Figure2_TS1_EVE_CPT1.csv")
write.csv(ts2, quote = F, row.names = F, "Figure2_TS2_EVE_CPT1.csv")
write.csv(ts3, quote = F, row.names = F, "Figure2_TS3_EVE_CPT1.csv")
write.csv(ts4, quote = F, row.names = F, "Figure2_TS4_EVE_CPT1.csv")
write.csv(ts5, quote = F, row.names = F, "Figure2_TS5_EVE_CPT1.csv")





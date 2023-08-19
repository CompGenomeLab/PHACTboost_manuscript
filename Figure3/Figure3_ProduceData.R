library(PRROC)
library(AUC)

load("/Users/nurdankuru/Desktop/PHACTboost/TrainingSet.RData")
load("/Users/nurdankuru/Desktop/PHACTboost/TestSet.RData")
load("/Users/nurdankuru/Desktop/PHACTboost/lightgbm_replication_1_prediction.RData")
test <- cbind(prediction$test_prediction, test)
phactboost <- test
load("/Users/nurdankuru/Desktop/PHACTboost/dbNSFP_44a_TestSet.RData")

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

load("/Users/nurdankuru/PHACTboost/Hard_Cases.RData")
load("/Users/nurdankuru/PHACTboost/Hard_Cases_Mixed.RData")

vars_all <- phactboost$prot_vars
vars_diffprots <- test_unseenprot$prot_vars
vars_diffvars <- unlist(test_unseenvar$prot_vars)
vars_hardcases1 <- hard_cases$prot_vars
vars_hardcases2 <- hard_cases2$prot_vars
rm(id, ids_ts4, ind1, ind2, inds_ts3, var_test, var_train, vars_all)
rm(hard_cases, hard_cases2, test_unseenprot, test_unseenvar)

dbnsfp_testset <- cbind(dbnsfp_testset$prot_vars, dbnsfp_testset[, c(grep("ank", colnames(dbnsfp_testset)))])
colnames(dbnsfp_testset)[1] <- "prot_vars"
tools_considering_AF <- c("MetaSVM_rankscore", "MetaLR_rankscore", "MetaRNN_rankscore",
                          "BayesDel_addAF_rankscore", "BayesDel_noAF_rankscore", "M-CAP_rankscore",
                          "MutPred_rankscore", "ClinPred_rankscore")

dbnsfp_testset <- dbnsfp_testset[,-match(tools_considering_AF, colnames(dbnsfp_testset))]

com <- intersect(phactboost$prot_vars, dbnsfp_testset$prot_vars)
data_main <- cbind(phactboost$variant_info[match(com, phactboost$prot_vars)], 
                   phactboost$`prediction$test_prediction`[match(com, phactboost$prot_vars)],
                   dbnsfp_testset[match(com, dbnsfp_testset$prot_vars),],
                   phactboost$prot_vars[match(com, phactboost$prot_vars)])

data_main <- as.data.frame(data_main)
colnames(data_main)[c(1,2,47)] <- c("variant_info", "PHACTboost", "prot_vars")
data_main <- data_main[,1:46]

vars_hardcases1 <- intersect(vars_hardcases1, dbnsfp_testset$prot_vars)
vars_hardcases2 <- intersect(vars_hardcases2, dbnsfp_testset$prot_vars)
vars_diffprots <- intersect(vars_diffprots, dbnsfp_testset$prot_vars)
vars_diffvars <- intersect(vars_diffvars, dbnsfp_testset$prot_vars)

ts1 <- data_main
ts2 <- data_main[match(vars_hardcases1, data_main$prot_vars),]
ts3 <- data_main[match(vars_hardcases2, data_main$prot_vars),]
ts4 <- data_main[match(vars_diffprots, data_main$prot_vars),]
ts5 <- data_main[match(vars_diffvars, data_main$prot_vars),]

save(ts1, file = "TS1_dbNSFP.RData")
save(ts2, file = "TS2_dbNSFP.RData")
save(ts3, file = "TS3_dbNSFP.RData")
save(ts4, file = "TS4_dbNSFP.RData")
save(ts5, file = "TS5_dbNSFP.RData")

write.csv(ts1, quote = F, row.names = F, "Figure3_TS1_dbNSFP.csv")
write.csv(ts2, quote = F, row.names = F, "Figure3_TS2_dbNSFP.csv")
write.csv(ts3, quote = F, row.names = F, "Figure3_TS3_dbNSFP.csv")
write.csv(ts4, quote = F, row.names = F, "Figure3_TS4_dbNSFP.csv")
write.csv(ts5, quote = F, row.names = F, "Figure3_TS5_dbNSFP.csv")


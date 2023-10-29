library(PRROC)
library(AUC)

phactboost <- read.csv("../PHACTboost_Model/PHACTboost_TestPrediction.csv")
msaboost <- read.csv("Data/MSAboost/MSAboost_TestPrediction.csv")
treeboost <- read.csv("Data/TREEboost/TREEboost_TestPrediction.csv")

ts1 <- read.csv("../PHACTboost_Model/TS1.csv")
ts2 <- read.csv("../PHACTboost_Model/TS2.csv")
ts3 <- read.csv("../PHACTboost_Model/TS3.csv")
ts4 <- read.csv("../PHACTboost_Model/TS4.csv")
ts5 <- read.csv("../PHACTboost_Model/TS5.csv")

data_main <- cbind(phactboost$variant_info, phactboost$PHACTboost_Score,
                   phactboost$PHACT_Score, treeboost$TREEboost_Score,
                   msaboost$MSAboost_Score, phactboost$prot_vars)

data_main <- as.data.frame(data_main)
colnames(data_main) <- c("variant_info", "PHACTboost", "PHACT", "TREEboost", "MSAboost", "prot_vars")

ts1 <- data_main
ts2 <- data_main[match(ts2$prot_vars, data_main$prot_vars),]
ts3 <- data_main[match(ts3$prot_vars, data_main$prot_vars),]
ts4 <- data_main[match(ts4$prot_vars, data_main$prot_vars),]
ts5 <- data_main[match(ts5$prot_vars, data_main$prot_vars),]

write.csv(ts1, quote = F, row.names = F, "Figure1_TS1.csv")
write.csv(ts2, quote = F, row.names = F, "Figure1_TS2.csv")
write.csv(ts3, quote = F, row.names = F, "Figure1_TS3.csv")
write.csv(ts4, quote = F, row.names = F, "Figure1_TS4.csv")
write.csv(ts5, quote = F, row.names = F, "Figure1_TS5.csv")

library(PRROC)
library(AUC)

phactboost <- read.csv("../PHACTboost_Model/PHACTboost_TestPrediction.csv")
alphamissense <- read.csv("Data/AlphaMissense_TestPrediction.csv")
eve <- read.csv("Data/EVE_TestPrediction.csv")
cpt1 <- read.csv("Data/CPT1_TestPrediction.csv")

ts1 <- read.csv("../PHACTboost_Model/TS1.csv")
ts2 <- read.csv("../PHACTboost_Model/TS2.csv")
ts3 <- read.csv("../PHACTboost_Model/TS3.csv")
ts4 <- read.csv("../PHACTboost_Model/TS4.csv")
ts5 <- read.csv("../PHACTboost_Model/TS5.csv")

data_main <- cbind(phactboost$variant_info, phactboost$PHACTboost_Score,
                   phactboost$PHACT_Score, 
                   alphamissense$AlphaMissense[match(phactboost$prot_vars, alphamissense$prot_vars)], 
                   eve$EVE[match(phactboost$prot_vars, eve$prot_vars)],
                   cpt1$CPT1[match(phactboost$prot_vars, cpt1$prot_vars)], phactboost$prot_vars)

data_main <- as.data.frame(data_main)
colnames(data_main) <- c("variant_info", "PHACTboost", "PHACT", "AlphaMissense", "EVE", 
                         "CPT-1", "prot_vars")

ts1 <- data_main
ts2 <- data_main[match(ts2$prot_vars, data_main$prot_vars),]
ts3 <- data_main[match(ts3$prot_vars, data_main$prot_vars),]
ts4 <- data_main[match(ts4$prot_vars, data_main$prot_vars),]
ts5 <- data_main[match(ts5$prot_vars, data_main$prot_vars),]

write.csv(ts1, quote = F, row.names = F, "Figure3_TS1.csv")
write.csv(ts2, quote = F, row.names = F, "Figure3_TS2.csv")
write.csv(ts3, quote = F, row.names = F, "Figure3_TS3.csv")
write.csv(ts4, quote = F, row.names = F, "Figure3_TS4.csv")
write.csv(ts5, quote = F, row.names = F, "Figure3_TS5.csv")

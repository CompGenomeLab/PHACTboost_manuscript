library(AUC)
library(PRROC)

load("/Users/nurdankuru/Desktop/PHACTboost/Test_gnomAD_Shared.RData")
load("/Users/nurdankuru/Desktop/PHACTboost/lightgbm_replication_1_prediction.RData")
test <- cbind(prediction$test_prediction, test)
load("/Users/nurdankuru/Desktop/PHACTboost/Train_gnomAD_Shared.RData")
rm(prediction)

# TS1
PHACTboost <- as.numeric(test$`prediction$test_prediction`)
y_pb <- as.numeric(test$variant_info)

# TS3
var_train <- paste(train$UNIPROTKB, train$Positions, sep = "-")
var_test <- paste(test$UNIPROTKB, test$Positions, sep = "-")
inds_ts3 <- setdiff(var_test, var_train)
ind1 <- c()
for (id in inds_ts3){
  ind1 <- c(ind1, which(var_test==id))
}
test_ts3 <- test[ind1, ]
PHACTboost_ts3 <- as.numeric(test_ts3$`prediction$test_prediction`)
y_ts3 <- as.numeric(test_ts3$variant_info)

# TS4
ids_ts4 <- setdiff(unique(test$UNIPROTKB), unique(train$UNIPROTKB))
ind2 <- c()
for (id in ids_ts4){
  ind2 <- c(ind2, which(test$UNIPROTKB==id))
}
test_ts4 <- test[ind2, ]
PHACTboost_ts4 <- as.numeric(test_ts4$`prediction$test_prediction`)
y_ts4 <- as.numeric(test_ts4$variant_info)


output_roc <- "PHACTboost_DiffVariantSet_ROCCurves.pdf"
output_pr <- "PHACTboost_DiffVariantSet_PRCurves.pdf"

cols <- c("#ca0020", "#4393c3", "#2166ac", "#ffa056")

auc_val <- list()
pdf(file = sprintf(output_roc), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)
auc_values <- c()
ind <- 1
ind <-  ind + 1
phact_roc <- roc(as.numeric(PHACTboost_ts3), as.factor(1 * (y_ts3 == 1)))
plot(phact_roc, col = cols[ind], lwd = 2)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACTboost - TS3 = %.3f", auc(phact_roc)), col = cols[ind], cex = 0.9)

ind <- ind + 1
phactboost_roc <- roc(as.numeric(PHACTboost_ts4), as.factor(1 * (y_ts4 == 1)))
lines(phactboost_roc$fpr, phactboost_roc$tpr, col = cols[ind], lwd = 2)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACTboost - TS4 = %.3f", auc(phactboost_roc)), col = cols[ind], cex = 0.9)

ind <- 1
phactboost_roc <- roc(as.numeric(PHACTboost), as.factor(1 * (y_pb == 1)))
lines(phactboost_roc$fpr, phactboost_roc$tpr, col = cols[ind], lwd = 2)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACTboost - TS1 = %.3f", auc(phactboost_roc)), col = cols[ind], cex = 1)

mtext(side = 2, at = 0.5 , text = "TPR (sensitivity)", line = 2.5, cex = 1)
mtext(side = 1, at = 0.5 , text = "FPR (1 - specificity)", line = 2.5, cex = 1)

dev.off()  


auc_val <- list()
pdf(file = sprintf(output_pr), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)

auc_values <- c()
ind <- 1
ind <-  ind + 1
phact_pr <- pr.curve(scores.class0 = PHACTboost_ts3, weights.class0 = (1 * (y_ts3 == 1)), curve = TRUE)
plot(phact_pr, col = cols[ind], lwd = 2, main = "", auc.main = FALSE)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACTboost - TS3 = %.3f", (phact_pr$auc.integral)), col = cols[ind], cex = 0.9)

ind <- ind + 1
phactboost_pr <- pr.curve(scores.class0 = as.numeric(PHACTboost_ts4), weights.class0 = (1 * (y_ts4 == 1)), curve = TRUE)
lines(phactboost_pr$curve[,1], phactboost_pr$curve[,2], col = cols[ind], lwd = 2)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACTboost - TS4 = %.3f", phactboost_pr$auc.integral), col = cols[ind], cex = 0.9)

ind <- 1
phactboost_pr <- pr.curve(scores.class0 = as.numeric(PHACTboost), weights.class0 = (1 * (y_pb == 1)), curve = TRUE)
lines(phactboost_pr$curve[,1], phactboost_pr$curve[,2], col = cols[ind], lwd = 2)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACTboost - TS1 = %.3f", phactboost_pr$auc.integral), col = cols[ind], cex = 1)

mtext(side = 2, at = 0.5 , text = "Precision", line = 2.5, cex = 1)
mtext(side = 1, at = 0.5 , text = "Recall", line = 2.5, cex = 1)
dev.off() 

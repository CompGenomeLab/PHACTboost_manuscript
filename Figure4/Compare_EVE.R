library(bio3d)
library(stringr)
library(AUC)
library(PRROC)
library(bio3d)
library(RColorBrewer)
library(caret)

load("/Users/nurdankuru/Desktop/PHACTboost/Test_gnomAD_Shared.RData")
load("/Users/nurdankuru/Desktop/PHACTboost/lightgbm_replication_1_prediction.RData")
test <- cbind(prediction$test_prediction, test)

load("/Users/nurdankuru/PHACTboost/EVE_Scores.RData")
all <- all[-which(is.na(all$EVE_scores_ASM)),]
com <- intersect(all$prot_vars, test$prot_vars)

phactboost_full <- as.numeric(test[match(com, test$prot_vars),1])
alt_full <- as.numeric(all$EVE_scores_ASM[match(com, all$prot_vars)])
y <- as.numeric(test[match(com, test$prot_vars),2])

cols <- c("mediumblue", "#5392A0")

    auc_val <- list()
    pdf(file = sprintf("%s_ROC.pdf", "EVE"), width = 4, height = 4)
    par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)
    auc_values <- c()
    ind <- 1
    ind <-  ind + 1
    pp2 <- roc(alt_full, as.factor(1 * (y == 1)))
    plot(pp2, col = cols[ind], lwd = 2)
    text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("%s = %.3f", "EVE", auc(pp2)), col = cols[ind], cex = 1)
    ind <- 1
    phylas_roc <- roc(phactboost_full, as.factor(1 * (y == 1)))
    lines(phylas_roc$fpr, phylas_roc$tpr, col = cols[ind], lwd = 2)
    text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("PHACTboost = %.3f", auc(phylas_roc)), col = cols[ind], cex = 1)
    mtext(side = 2, at = 0.5 , text = "TPR (sensitivity)", line = 2.5, cex = 1)
    mtext(side = 1, at = 0.5 , text = "FPR (1 - specificity)", line = 2.5, cex = 1)
    dev.off()  
    
    
    auc_val <- list()
    pdf(file = sprintf("%s_PR.pdf", "EVE"), width = 4, height = 4)
    par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)
    auc_values <- c()
    ind <- 1
    ind <-  ind + 1
    pp2 <- pr.curve(scores.class0 = alt_full, weights.class0 = (1 * (y == 1)), curve = TRUE)
    plot(pp2, col = cols[ind], lwd = 2, main = "", auc.main = FALSE)
    text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("%s = %.3f", "EVE", (pp2$auc.integral)), col = cols[ind], cex = 1)
    ind <- 1
    phylas_prauc <- pr.curve(scores.class0 = phactboost_full, weights.class0 = (1 * (y == 1)), curve = TRUE)
    lines(phylas_prauc$curve[,1], phylas_prauc$curve[,2], col = cols[ind], lwd = 2)
    text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("PHACTboost = %.3f", phylas_prauc$auc.integral), col = cols[ind], cex = 1)
    
    mtext(side = 2, at = 0.5 , text = "Precision", line = 2.5, cex = 1)
    mtext(side = 1, at = 0.5 , text = "Recall", line = 2.5, cex = 1)
    dev.off() 




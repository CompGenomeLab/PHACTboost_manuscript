library(AUC)
library(PRROC)

data <- read.csv("PHACT_Original_vs_SingleChar.csv")
y <- as.numeric(data$variant_info)
val1 <- as.numeric(data$Phact_wl_param_CountNodes_3_alt)
val2 <- as.numeric(data$parsimony)

cols <- c("blue", "red")
cols <- c("#3498db", "#e74c3c")

auc_val <- list()
pdf(file = sprintf("ROC.pdf"), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)
auc_values <- c()
ind <- 1
ind <-  ind + 1
pp2 <- roc(1 - val2, as.factor(1 * (y == 1)))
plot(pp2, col = cols[ind], lwd = 2)
text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("Maximum Parsimony = %.3f", auc(pp2)), col = cols[ind], cex = 1)
ind <- 1
phylas_roc <- roc(1 - val1, as.factor(1 * (y == 1)))
lines(phylas_roc$fpr, phylas_roc$tpr, col = cols[ind], lwd = 2)
text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("Maximum Likelihood = %.3f", auc(phylas_roc)), col = cols[1], cex = 1)
mtext(side = 2, at = 0.5 , text = "TPR (sensitivity)", line = 2.5, cex = 1)
mtext(side = 1, at = 0.5 , text = "FPR (1 - specificity)", line = 2.5, cex = 1)

dev.off()  


auc_val <- list()
pdf(file = sprintf("PR.pdf"), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)
auc_values <- c()
ind <- 1
ind <-  ind + 1
pp2 <- pr.curve(scores.class0 = 1 - val2, weights.class0 = (1 * (y == 1)), curve = TRUE)
plot(pp2, col = cols[ind], lwd = 2, main = "", auc.main = FALSE)
text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("Maximum Parsimony = %.3f", (pp2$auc.integral)), col = cols[ind], cex = 1)
ind <- 1
phylas_prauc <- pr.curve(scores.class0 = 1 - val1, weights.class0 = (1 * (y == 1)), curve = TRUE)
lines(phylas_prauc$curve[,1], phylas_prauc$curve[,2], col = cols[ind], lwd = 2)
text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("Maximum Likelihood = %.3f", phylas_prauc$auc.integral), col = cols[1], cex = 1)
mtext(side = 2, at = 0.5 , text = "Precision", line = 2.5, cex = 1)
mtext(side = 1, at = 0.5 , text = "Recall", line = 2.5, cex = 1)
dev.off() 




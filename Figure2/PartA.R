library(AUC)
library(PRROC)


# Part 2A

output_roc <- "TS1_PHACT_ROCCurves.pdf"
output_pr <- "TS1_PHACT_PRCurves.pdf"

data <- read.csv("PHACTboost_vs_PHACT_TS1.csv")
y <- as.numeric(data$Type)

cols <- c("#ca0020", "#4393c3", "#2166ac", "#ffa056")

auc_val <- list()
pdf(file = sprintf(output_roc), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)
auc_values <- c()
ind <- 1
ind <-  ind + 1
msaboost_roc <- roc(1 - as.numeric(data$PHACT), as.factor(1 * (y == 1)))
plot(msaboost_roc, col = cols[ind], lwd = 2)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACT = %.3f", auc(msaboost_roc)), col = cols[ind], cex = 0.9)

ind <- 1
phactboost_roc <- roc(as.numeric(data$PHACTboost), as.factor(1 * (y == 1)))
lines(phactboost_roc$fpr, phactboost_roc$tpr, col = cols[ind], lwd = 2)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACTboost = %.3f", auc(phactboost_roc)), col = cols[1], cex = 1)

mtext(side = 2, at = 0.5 , text = "TPR (sensitivity)", line = 2.5, cex = 1)
mtext(side = 1, at = 0.5 , text = "FPR (1 - specificity)", line = 2.5, cex = 1)

dev.off()  


auc_val <- list()
pdf(file = sprintf(output_pr), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)

auc_values <- c()
ind <- 1
ind <-  ind + 1
msaboost_pr <- pr.curve(scores.class0 = 1 - as.numeric(data$PHACT), weights.class0 = (1 * (y == 1)), curve = TRUE)
plot(msaboost_pr, col = cols[ind], lwd = 2, main = "", auc.main = FALSE)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACT = %.3f", (msaboost_pr$auc.integral)), col = cols[ind], cex = 1)

ind <- 1
phactboost_pr <- pr.curve(scores.class0 = as.numeric(data$PHACTboost), weights.class0 = (1 * (y == 1)), curve = TRUE)
lines(phactboost_pr$curve[,1], phactboost_pr$curve[,2], col = cols[ind], lwd = 2)
text(0.5, 0.4 - (ind - 1) * 0.07, label = sprintf("PHACTboost = %.3f", phactboost_pr$auc.integral), col = cols[1], cex = 1)

mtext(side = 2, at = 0.5 , text = "Precision", line = 2.5, cex = 1)
mtext(side = 1, at = 0.5 , text = "Recall", line = 2.5, cex = 1)
dev.off() 

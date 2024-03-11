library(PRROC)
library(AUC)

set.seed(123)

test_data <- read.csv("ConfInt_PHACTboostResults.csv")

TS1 <- read.csv("https://github.com/CompGenomeLab/PHACTboost_manuscript/blob/main/Figure1/Data/Figure1_TS1.csv")
TS1 <- TS1[,c(4:14, 2)]

pat <- TS1[which(TS1$variant_info==1),]
net <- TS1[which(TS1$variant_info==-1),]

i1 <- sample(1:length(pat$variant_info), length(which(test_data$Variant_Info==1)))
i2 <- sample(1:length(net$variant_info), length(which(test_data$Variant_Info==-1)))

test_data_original <- rbind(pat[i1, ], net[i2, ])


y1 <- as.numeric(test_data_original$variant_info)
y2 <- as.numeric(test_data$Variant_Info)
val1 <- as.numeric(test_data_original$PHACTboost)
val2 <- as.numeric(test_data$PHACTboost)

cols <- c("#3498db", "#e74c3c")

auc_val <- list()
pdf(file = sprintf("ROC.pdf"), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)
auc_values <- c()
ind <- 1
ind <-  ind + 1
pp2 <- roc(val2, as.factor(1 * (y2 == 1)))
plot(pp2, col = cols[ind], lwd = 2)
text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("ConflictingInt = %.3f", auc(pp2)), col = cols[ind], cex = 1)
ind <- 1
phylas_roc <- roc(val1, as.factor(1 * (y1 == 1)))
lines(phylas_roc$fpr, phylas_roc$tpr, col = cols[ind], lwd = 2)
text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("TS1_Subset = %.3f", auc(phylas_roc)), col = cols[1], cex = 1)
mtext(side = 2, at = 0.5 , text = "TPR (sensitivity)", line = 2.5, cex = 1)
mtext(side = 1, at = 0.5 , text = "FPR (1 - specificity)", line = 2.5, cex = 1)

dev.off()  

auc_val <- list()
pdf(file = sprintf("PR.pdf"), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)
auc_values <- c()
ind <- 1
ind <-  ind + 1
pp2 <- pr.curve(scores.class0 = val2, weights.class0 = (1 * (y2 == 1)), curve = TRUE)
plot(pp2, col = cols[ind], lwd = 2, main = "", auc.main = FALSE)
text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("ConflictingInt = %.3f", (pp2$auc.integral)), col = cols[ind], cex = 1)
ind <- 1
phylas_prauc <- pr.curve(scores.class0 = val1, weights.class0 = (1 * (y1 == 1)), curve = TRUE)
lines(phylas_prauc$curve[,1], phylas_prauc$curve[,2], col = cols[ind], lwd = 2)
text(0.6, 0.35 - (ind - 1) * 0.07, label = sprintf("TS1_Subset = %.3f", phylas_prauc$auc.integral), col = cols[1], cex = 1)
mtext(side = 2, at = 0.5 , text = "Precision", line = 2.5, cex = 1)
mtext(side = 1, at = 0.5 , text = "Recall", line = 2.5, cex = 1)
dev.off() 


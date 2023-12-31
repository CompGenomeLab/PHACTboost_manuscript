ts1 <- read.csv("Data/Figure3_TS1.csv")
ts2 <- read.csv("Data/Figure3_TS2.csv")
ts3 <- read.csv("Data/Figure3_TS3.csv")
ts4 <- read.csv("Data/Figure3_TS4.csv")
ts5 <- read.csv("Data/Figure3_TS5.csv")

all_sets <- list()
all_sets$ts1 <- ts1
all_sets$ts2 <- ts2
all_sets$ts3 <- ts3
all_sets$ts4 <- ts4
all_sets$ts5 <- ts5

values_roc_alpmis <- c()
values_roc_eve <- c()
values_roc_cpt1 <- c()
values_pr_alpmis <- c()
values_pr_eve <- c()
values_pr_cpt1 <- c()

for (i in 1:5){
  chosen <- all_sets[[paste("ts", i, sep = "")]]
  
  y <- as.numeric(chosen$variant_info)
  pb <- as.numeric(chosen$PHACTboost)
  
  alpmis <- as.numeric(chosen$AlphaMissense)
  el <- which(is.na(alpmis))
  if (length(el)>0){
    pb_roc <- auc(roc(pb[-el], as.factor(1 * (y[-el] == 1))))
    alpmis_roc <- auc(roc(alpmis[-el], as.factor(1 * (y[-el] == 1))))
    pb_pr <- pr.curve(scores.class0 = pb[-el], weights.class0 = (1 * (y[-el] == 1)), curve = TRUE)
    alpmis_pr <- pr.curve(scores.class0 = alpmis[-el], weights.class0 = (1 * (y[-el] == 1)), curve = TRUE)
    nvar <- table(y[-el])
  } else {
    pb_roc <- auc(roc(pb, as.factor(1 * (y == 1))))
    alpmis_roc <- auc(roc(alpmis, as.factor(1 * (y == 1))))
    pb_pr <- pr.curve(scores.class0 = pb, weights.class0 = (1 * (y == 1)), curve = TRUE)
    alpmis_pr <- pr.curve(scores.class0 = alpmis, weights.class0 = (1 * (y == 1)), curve = TRUE)
    nvar <- table(y)
  }
  values_roc_alpmis <- cbind(values_roc_alpmis, c(nvar, pb_roc, alpmis_roc))
  values_pr_alpmis <- cbind(values_pr_alpmis, c(nvar, pb_pr$auc.integral, alpmis_pr$auc.integral))
  
  eve <- as.numeric(chosen$EVE)
  el <- which(is.na(eve))
  if (length(el)>0){
    pb_roc <- auc(roc(pb[-el], as.factor(1 * (y[-el] == 1))))
    eve_roc <- auc(roc(eve[-el], as.factor(1 * (y[-el] == 1))))
    pb_pr <- pr.curve(scores.class0 = pb[-el], weights.class0 = (1 * (y[-el] == 1)), curve = TRUE)
    eve_pr <- pr.curve(scores.class0 = eve[-el], weights.class0 = (1 * (y[-el] == 1)), curve = TRUE)
    nvar <- table(y[-el])
  } else {
    pb_roc <- auc(roc(pb, as.factor(1 * (y == 1))))
    eve_roc <- auc(roc(eve, as.factor(1 * (y == 1))))
    pb_pr <- pr.curve(scores.class0 = pb, weights.class0 = (1 * (y == 1)), curve = TRUE)
    eve_pr <- pr.curve(scores.class0 = eve, weights.class0 = (1 * (y == 1)), curve = TRUE)
    nvar <- table(y)
  }
  values_roc_eve <- cbind(values_roc_eve, c(nvar, pb_roc, eve_roc))
  values_pr_eve <- cbind(values_pr_eve, c(nvar, pb_pr$auc.integral, eve_pr$auc.integral))
  
  cpt1 <- as.numeric(chosen$`CPT-1`)
  el <- which(is.na(cpt1))
  if (length(el)>0){
    pb_roc <- auc(roc(pb[-el], as.factor(1 * (y[-el] == 1))))
    cpt1_roc <- auc(roc(cpt1[-el], as.factor(1 * (y[-el] == 1))))
    pb_pr <- pr.curve(scores.class0 = pb[-el], weights.class0 = (1 * (y[-el] == 1)), curve = TRUE)
    cpt1_pr <- pr.curve(scores.class0 = cpt1[-el], weights.class0 = (1 * (y[-el] == 1)), curve = TRUE)
    nvar <- table(y[-el])
  } else {
    pb_roc <- auc(roc(pb, as.factor(1 * (y == 1))))
    cpt1_roc <- auc(roc(cpt1, as.factor(1 * (y == 1))))
    pb_pr <- pr.curve(scores.class0 = pb, weights.class0 = (1 * (y == 1)), curve = TRUE)
    cpt1_pr <- pr.curve(scores.class0 = cpt1, weights.class0 = (1 * (y == 1)), curve = TRUE)
    nvar <- table(y)
  }
  values_roc_cpt1 <- cbind(values_roc_cpt1, c(nvar, pb_roc, cpt1_roc))
  values_pr_cpt1 <- cbind(values_pr_cpt1, c(nvar, pb_pr$auc.integral, cpt1_pr$auc.integral))
}

values_roc_alpmis <- as.data.frame(values_roc_alpmis)
values_pr_alpmis <- as.data.frame(values_pr_alpmis)
values_roc_eve <- as.data.frame(values_roc_eve)
values_pr_eve <- as.data.frame(values_pr_eve)
values_roc_cpt1 <- as.data.frame(values_roc_cpt1)
values_pr_cpt1 <- as.data.frame(values_pr_cpt1)

rownames(values_roc_alpmis) <- c("Neutral", "Pathogenic", "PHACTboost", "AlphaMissense")
rownames(values_pr_alpmis) <- c("Neutral", "Pathogenic", "PHACTboost", "AlphaMissense")
rownames(values_roc_eve) <- c("Neutral", "Pathogenic", "PHACTboost", "EVE")
rownames(values_pr_eve) <- c("Neutral", "Pathogenic", "PHACTboost", "EVE")
rownames(values_roc_cpt1) <- c("Neutral", "Pathogenic", "PHACTboost", "CPT1")
rownames(values_pr_cpt1) <- c("Neutral", "Pathogenic", "PHACTboost", "CPT1")


dataset_names <- c("TS1", "TS2", "TS3", "TS4", "TS5")
colnames(values_roc_alpmis) <- dataset_names
colnames(values_pr_alpmis) <- dataset_names
colnames(values_roc_eve) <- dataset_names
colnames(values_pr_eve) <- dataset_names
colnames(values_roc_cpt1) <- dataset_names
colnames(values_pr_cpt1) <- dataset_names

cols <- c("#18401E", "#A9D9A7")
pdf(file = "AlphaMissense_ROC.pdf", width = 9, height = 5)
barplot(as.matrix(values_roc_alpmis[3:4,]), col = cols, beside = TRUE, ylim = c(0,1.3), yaxt = "n")
legend(12, 1.25, legend = rownames(values_roc_alpmis[3:4,]), fill = cols, bty = "n", y.intersp = 0.7)
axis(side = 2, at = seq(0,1,0.2), labels = sprintf("%.2f",seq(0,1,0.2)), las = 2)
mtext(side = 2, text = "AUROC", line = 3, cex = 1)
dev.off()

pdf(file = "AlphaMissense_PR.pdf", width = 9, height = 5)
barplot(as.matrix(values_pr_alpmis[3:4,]), col = cols, beside = TRUE, ylim = c(0,1.3), yaxt = "n")
legend(12, 1.25, legend = rownames(values_pr_alpmis[3:4,]), fill = cols, bty = "n", y.intersp = 0.7)
axis(side = 2, at = seq(0,1,0.2), labels = sprintf("%.2f",seq(0,1,0.2)), las = 2)
mtext(side = 2, text = "AUPR", line = 3, cex = 1)
dev.off()

cols <- c("#18401E", "#9A99BF")
pdf(file = "EVE_ROC.pdf", width = 9, height = 5)
barplot(as.matrix(values_roc_eve[3:4,]), col = cols, beside = TRUE, ylim = c(0,1.3), yaxt = "n")
legend(12, 1.25, legend = rownames(values_roc_eve[3:4,]), fill = cols, bty = "n", y.intersp = 0.7)
axis(side = 2, at = seq(0,1,0.2), labels = sprintf("%.2f",seq(0,1,0.2)), las = 2)
mtext(side = 2, text = "AUROC", line = 3, cex = 1)
dev.off()

pdf(file = "EVE_PR.pdf", width = 9, height = 5)
barplot(as.matrix(values_pr_eve[3:4,]), col = cols, beside = TRUE, ylim = c(0,1.3), yaxt = "n")
legend(12, 1.25, legend = rownames(values_pr_eve[3:4,]), fill = cols, bty = "n", y.intersp = 0.7)
axis(side = 2, at = seq(0,1,0.2), labels = sprintf("%.2f",seq(0,1,0.2)), las = 2)
mtext(side = 2, text = "AUPR", line = 3, cex = 1)
dev.off()

cols <- c("#18401E", "#F2B694")
pdf(file = "CPT1_ROC.pdf", width = 9, height = 5)
barplot(as.matrix(values_roc_cpt1[3:4,]), col = cols, beside = TRUE, ylim = c(0,1.3), yaxt = "n")
legend(12, 1.25, legend = rownames(values_roc_cpt1[3:4,]), fill = cols, bty = "n", y.intersp = 0.7)
axis(side = 2, at = seq(0,1,0.2), labels = sprintf("%.2f",seq(0,1,0.2)), las = 2)
mtext(side = 2, text = "AUROC", line = 3, cex = 1)
dev.off()

pdf(file = "CPT1_PR.pdf", width = 9, height = 5)
barplot(as.matrix(values_pr_cpt1[3:4,]), col = cols, beside = TRUE, ylim = c(0,1.3), yaxt = "n")
legend(12, 1.25, legend = rownames(values_pr_cpt1[3:4,]), fill = cols, bty = "n", y.intersp = 0.7)
axis(side = 2, at = seq(0,1,0.2), labels = sprintf("%.2f",seq(0,1,0.2)), las = 2)
mtext(side = 2, text = "AUPR", line = 3, cex = 1)
dev.off()


ts1 <- read.csv("Data/Figure1_TS1.csv")
ts2 <- read.csv("Data/Figure1_TS2.csv")
ts3 <- read.csv("Data//Figure1_TS3.csv")
ts4 <- read.csv("Data//Figure1_TS4.csv")
ts5 <- read.csv("Data//Figure1_TS5.csv")

all_sets <- list()
all_sets$ts1 <- ts1
all_sets$ts2 <- ts2
all_sets$ts3 <- ts3
all_sets$ts4 <- ts4
all_sets$ts5 <- ts5

values_roc <- c()
values_pr <- c()
for (i in 1:5){
  chosen <- all_sets[[paste("ts", i, sep = "")]]
  
  y <- as.numeric(chosen$variant_info)
  pb <- as.numeric(chosen$PHACTboost)
  tb <- as.numeric(chosen$TREEboost)
  mb <- as.numeric(chosen$MSAboost)
  p <- 1-as.numeric(chosen$PHACT)
  
  pb_roc <- auc(roc(pb, as.factor(1 * (y == 1))))
  tb_roc <- auc(roc(tb, as.factor(1 * (y == 1))))
  mb_roc <- auc(roc(mb, as.factor(1 * (y == 1))))
  p_roc <- auc(roc(p, as.factor(1 * (y == 1))))
  values_roc <- cbind(values_roc, c(pb_roc, tb_roc, mb_roc, p_roc))
  
  pb_pr <- pr.curve(scores.class0 = pb, weights.class0 = (1 * (y == 1)), curve = TRUE)
  tb_pr <- pr.curve(scores.class0 = tb, weights.class0 = (1 * (y == 1)), curve = TRUE)
  mb_pr <- pr.curve(scores.class0 = mb, weights.class0 = (1 * (y == 1)), curve = TRUE)
  p_pr <- pr.curve(scores.class0 = p, weights.class0 = (1 * (y == 1)), curve = TRUE)
  values_pr <- cbind(values_pr, c(pb_pr$auc.integral, tb_pr$auc.integral, mb_pr$auc.integral, p_pr$auc.integral))
}

values_roc <- as.data.frame(values_roc)
values_pr <- as.data.frame(values_pr)
rownames(values_roc) <- c("PHACTboost", "TREEboost", "MSAboost", "PHACT")
rownames(values_pr) <- c("PHACTboost", "TREEboost", "MSAboost", "PHACT")
colnames(values_roc) <- c("TS1", "TS2", "TS3", "TS4", "TS5")
colnames(values_pr) <-  c("TS1", "TS2", "TS3", "TS4", "TS5")

cols <- c("#18401E", "#9A99BF", "#F2B694", "#A9D9A7")

pdf(file = "Figure1_ROC.pdf", width = 9, height = 5)
barplot(as.matrix(values_roc), col = cols, beside = TRUE, ylim = c(0,1.3), yaxt = "n")
legend(17, 1.25, legend = rownames(values_roc), fill = cols, bty = "n", y.intersp = 0.7)
axis(side = 2, at = seq(0,1,0.2), labels = sprintf("%.2f",seq(0,1,0.2)), las = 2)
mtext(side = 2, text = "AUROC", line = 3, cex = 1)
dev.off()

pdf(file = "Figure1_PR.pdf", width = 9, height = 5)
barplot(as.matrix(values_pr), col = cols, beside = TRUE, ylim = c(0,1.3), yaxt = "n")
legend(21, 1.3, legend = rownames(values_pr), fill = cols, bty = "n", y.intersp = 0.7)
axis(side = 2, at = seq(0,1,0.2), labels = sprintf("%.2f",seq(0,1,0.2)), las = 2)
mtext(side = 2, text = "AUPR", line = 3, cex = 1)
dev.off()

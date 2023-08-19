library(bio3d)
library(stringr)
library(AUC)
library(PRROC)
library(bio3d)
library(RColorBrewer)
library(caret)

# Update this part to reproduce the results for test set TSX
name <- "TS1"

scores <- read.csv(sprintf("Figure3_AUC_AUPR_%s_dbNSFP.csv", name))

names <- c("SIFT","SIFT4G","Polyphen2_HDIV","Polyphen2_HVAR","LRT","MutationTaster",
           "MutationAssessor","FATHMM","PROVEAN","VEST4", "REVEL",
           "MVP", "gMVP", "MPC","PrimateAI","DEOGEN2","LIST-S2", "VARITY_R", "VARITY_ER", 
           "VARITY_R_LOO", "VARITY_ER_LOO", "CADD_raw","CADD_raw_hg19","DANN",
           "fathmm-MKL_coding","fathmm-XF_coding","Eigen-raw_coding","Eigen-PC-raw_coding",
           "GenoCanyon","integrated_fitCons","GM12878_fitCons","H1-hESC_fitCons","HUVEC_fitCons",
           "LINSIGHT", "GERP++","phyloP100way","phyloP470way","phyloP17way","phastCons100way",
           "phastCons470way","phastCons17way","SiPhy","bStatistic")


cc <- c()
k <- 0
keep <- c()
values <- c()
for (i in 4:length(scores[1,])){
  k <- k + 1
  scores_pair <- scores
  
  y_full <- as.numeric(scores_pair$variant_info)
  phactboost_full <- as.numeric(scores_pair$PHACTboost)
  alt_full <- scores_pair[,i]
  
  el <- which(alt_full==".")
  
  y_full <- y_full[-el]
  phactboost_full <- phactboost_full[-el]
  alt_full <- as.numeric(alt_full[-el])
  
  yy <- unique(y_full)
  
  if (length(yy)>1){
    phylas_roc <- roc(alt_full, as.factor(1 * (y_full == 1)))
    auc1 = auc(phylas_roc)
    phylas_roc <- roc(phactboost_full, as.factor(1 * (y_full == 1)))
    auc2 = auc(phylas_roc)
    
    phylas_prauc <- pr.curve(scores.class0 = alt_full, weights.class0 = (1 * (y_full == 1)), curve = TRUE)
    aupr1 = phylas_prauc$auc.integral
    phylas_prauc <- pr.curve(scores.class0 = phactboost_full, weights.class0 = (1 * (y_full == 1)), curve = TRUE)
    aupr2 = phylas_prauc$auc.integral
    
    keep <- rbind(keep, c(names[k], length(which(y_full==1)), length(which(y_full==-1)),
                          auc1, auc2, aupr1, aupr2))
    values <- rbind(values, c(names[k], auc1, auc2, (auc1-auc2), aupr1, aupr2, (aupr1-aupr2)))
  }
  
}
keep <- as.data.frame(keep)
colnames(keep) <- c("ToolName", "Num_Pat", "Num_Net", "AUROC_Tool", "AUROC_PHACTboost", "AUPR_Tool", "AUPR_PHACTboost")

ts1_result <- keep
write.csv(ts1_result, quote = F, row.names = F, sprintf("Figure3_AUC_AUPR_%s_dbNSFP.csv", name))

f <- sort(as.numeric(values[,4]), decreasing = T, index.return = T)
values <- values[f$ix,]
dd <- sort(as.numeric(values[,3]), decreasing =T, index.return = T)

cols <- c("#21457A", "#E8F3FD")

thr <- as.numeric(max(values[,7]))-0.1
# thr <- -0.5
pdf(file = sprintf("PR_%s.pdf", name), width = 9, height = 5)
par(oma = c(2, 1, 1, 0),  mar = c(1,4,1,0), cex = 11/12, lwd = 1, cex.axis = 1.0, cex.main = 1)
x1 <- barplot(as.numeric(values[,4]), col = cols[1], ylim = c(thr,0.6), horiz = F, yaxt = "n",
              las = 2, xaxt="n")
text(x1-0.4, 0.01*matrix(1,43,1), labels = values[,1], srt = 90, cex = 1, col = "black", pos = 4)
axis(side = 2, at = seq(-1,0,0.1), labels = sprintf("%.2f",seq(-1,0,0.1)), las = 2)
mtext(side = 2, at = thr/2, text = expression(Delta ~ "AUROC"), line = 3.2, cex = 1)
dev.off()

f <- sort(as.numeric(values[,7]), decreasing = T, index.return = T)
values <- values[f$ix,]
dd <- sort(as.numeric(values[,6]), decreasing =T, index.return = T)

thr <- as.numeric(max(values[,7]))-0.1
# thr <- -0.75
pdf(file = sprintf("ROC_%s.pdf", name), width = 9, height = 5)
par(oma = c(2, 1, 1, 0),  mar = c(1,4,1,0), cex = 11/12, lwd = 1, cex.axis = 1.0, cex.main = 1)
x1 <- barplot(as.numeric(values[,7]), col = cols[1], ylim = c(thr,0.6), horiz = F, yaxt = "n",
              las = 2, xaxt="n")
text(x1-0.4, 0.01*matrix(1,43,1), labels = values[,1], srt = 90, cex = 1, col = "black", pos = 4)
axis(side = 2, at = seq(-1,0,0.1), labels = sprintf("%.2f",seq(-1,0,0.1)), las = 2)
mtext(side = 2, at = thr/2, text = expression(Delta ~ "AUPR"), line = 3.2, cex = 1)
dev.off()


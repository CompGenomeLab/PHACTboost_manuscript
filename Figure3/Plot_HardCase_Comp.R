library(bio3d)
library(stringr)
library(AUC)
library(PRROC)
library(bio3d)
library(RColorBrewer)
library(caret)

#load("/Users/nurdankuru/Hard_Cases.RData")
load("/Users/nurdankuru/Hard_Cases_Mixed.RData")

output_name_roc <- "ROC_HardCases2.pdf"
output_name_pr <- "PR_HardCases2.pdf"

names <- c("SIFT","SIFT4G","Polyphen2_HDIV","Polyphen2_HVAR","LRT","MutationTaster",
           "MutationAssessor","FATHMM","PROVEAN","VEST4", "REVEL",
           "MVP", "gMVP", "MPC","PrimateAI","DEOGEN2","LIST-S2", "VARITY_R", "VARITY_ER", 
           "VARITY_R_LOO", "VARITY_ER_LOO", "CADD_raw","CADD_raw_hg19","DANN",
           "fathmm-MKL_coding","fathmm-XF_coding","Eigen-raw_coding","Eigen-PC-raw_coding",
           "GenoCanyon","integrated_fitCons","GM12878_fitCons","H1-hESC_fitCons","HUVEC_fitCons",
           "LINSIGHT", "GERP++","phyloP100way","phyloP470way","phyloP17way","phastCons100way",
           "phastCons470way","phastCons17way","SiPhy","bStatistic")

scores <- hard_cases2[,2:46]
cc <- c()
k <- 0
keep <- c()
values <- c()
for (i in 3:length(scores[1,])){
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
    
    values <- rbind(values, c(names[k], auc1, auc2, (auc1-auc2), aupr1, aupr2, (aupr1-aupr2)))
  }
  
}

f <- sort(as.numeric(values[,4]), decreasing = T, index.return = T)
values <- values[f$ix,]
dd <- sort(as.numeric(values[,3]), decreasing =T, index.return = T)

cols <- c("#21457A", "#E8F3FD")

pdf(file = output_name_roc, width = 9, height = 5)
par(oma = c(2, 1, 1, 0),  mar = c(1,4,1,0), cex = 11/12, lwd = 1, cex.axis = 1.0, cex.main = 1)
x1 <- barplot(as.numeric(values[,4]), col = col_upd[1], ylim = c(-0.8,0.6), horiz = F, yaxt = "n",
              las = 2, xaxt="n")
text(x1-0.4, 0.01*matrix(1,43,1), labels = values[,1], srt = 90, cex = 1, col = "black", pos = 4)
axis(side = 2, at = seq(-1,0,0.1), labels = sprintf("%.2f",seq(-1,0,0.1)), las = 2)
mtext(side = 2, at = -0.35, text = "AUROC Difference", line = 3.2, cex = 1)
dev.off()

f <- sort(as.numeric(values[,7]), decreasing = T, index.return = T)
values <- values[f$ix,]
dd <- sort(as.numeric(values[,6]), decreasing =T, index.return = T)

pdf(file = output_name_pr, width = 9, height = 5)
par(oma = c(2, 1, 1, 0),  mar = c(1,4,1,0), cex = 11/12, lwd = 1, cex.axis = 1.0, cex.main = 1)
x1 <- barplot(as.numeric(values[,7]), col = col_upd[1], ylim = c(-0.7,0.6), horiz = F, yaxt = "n",
              las = 2, xaxt="n")
text(x1-0.4, 0.01*matrix(1,43,1), labels = values[,1], srt = 90, cex = 1, col = "black", pos = 4)
axis(side = 2, at = seq(-1,0,0.1), labels = sprintf("%.2f",seq(-1,0,0.1)), las = 2)
mtext(side = 2, at = -0.3, text = "AUPR Difference", line = 3.2, cex = 1)
dev.off()

folder_name <- "Data/DMS_AllVariants"
# folder_name <- "Data/DMS_HardCases"
files <- list.files(folder_name)

cor_mat <- c()
for (file in files){
  load(sprintf("%s/%s", folder_name, file))
  print(dim(comb))
  
  id <- unlist(strsplit(unlist(strsplit(file, ".RData"))[1], "_"))[4]
  
  el2 <- which(is.na(comb$PHACTboost))
  if (length(el2)>0){
    comb <- comb[-el2,]
  }
  el3 <- which(is.na(comb$EVE))
  if (length(el3)>0){
    comb <- comb[-el3,]
  }
  el4 <- which(is.na(comb$AlphaMissense))
  if (length(el4)>0){
    comb <- comb[-el4,]
  }
  el5 <- which(is.na(comb$REVEL))
  if (length(el5)>0){
    comb <- comb[-el5,]
  }
  el6 <- which(is.na(comb$CPT1))
  if (length(el6)>0){
    comb <- comb[-el6,]
  }
  el7 <- which(is.na(comb$DMS))
  if (length(el7)>0){
    comb <- comb[-el7,]
  }
  
  cor1 <- cor(as.numeric(comb$DMS), as.numeric(comb$PHACTboost), method = "spearman")
  cor2 <- cor(as.numeric(comb$DMS), as.numeric(comb$EVE), method = "spearman")
  cor3 <- cor(as.numeric(comb$DMS), as.numeric(comb$REVEL), method = "spearman")
  cor4 <- cor(as.numeric(comb$DMS), as.numeric(comb$CPT1), method = "spearman")
  cor5 <- cor(as.numeric(comb$DMS), as.numeric(comb$AlphaMissense), method = "spearman")
  
  cor_mat <- rbind(cor_mat, c(file, abs(cor1), abs(cor2), abs(cor3), abs(cor4), abs(cor5), length(as.numeric(comb$DMS))))
  cor_mat <- as.data.frame(cor_mat)
  colnames(cor_mat) <- c("File", "PHACTboost", "EVE", "REVEL", "CPT1", "AlphaMissense",  "#")
  
  rm(comb, el2, el3, el4, el5, el6)
}

pb <- sum(as.numeric(cor_mat[,2]))/length(cor_mat[,1])
ev <- sum(as.numeric(cor_mat[,3]))/length(cor_mat[,1])
re <- sum(as.numeric(cor_mat[,4]))/length(cor_mat[,1])
am <- sum(as.numeric(cor_mat[,6]))/length(cor_mat[,1])
cp <- sum(as.numeric(cor_mat[,5]))/length(cor_mat[,1])

vect <- c(pb, cp, am, re, ev)
vect <- rbind(vect,  c("PHACTboost", "CPT-1", "AlphaMissense", "REVEL", "EVE"))

max.temp <- sort(as.numeric(vect[1,]), decreasing = T, index.return=T)
# barchart with added parameters
pdf(file = sprintf("BarPlot_AllCases.pdf"), width = 12, height = 8)
par(mar=c(4,5,4,1))
barplot(max.temp$x,
        main = "Spearman Correlation with DMS Data",
        ylab = expression("DMS assays correlation (mean)"),
        xlab = "",
        names.arg = vect[2,max.temp$ix],
        col = c("darkblue"),
        horiz = FALSE,
        cex.axis = 1.5, cex.names = 1.5, cex.main = 1.5, ylim = c(0,0.6),
        cex.lab=1.5)
dev.off()


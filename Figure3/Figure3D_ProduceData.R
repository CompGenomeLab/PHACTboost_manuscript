type <- 1
# 0 if Hard Cases

library(lightgbm)
library(AUC)

genes <- c("ADRB2", "A4", "SYUA", "BRCA1", "P53", "PTEN", "MSH2", "VKOR1")
ids <- c("P07550", "P05067", "P37840", "P38398", "P04637", "P60484", "P43246", "Q9BQB6")

aa_to_num <- function(aa) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  # num <- ifelse(sum(amino_acids %in% aa) == 1, as.numeric(which(amino_acids %in% aa)), 21)
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  aa <- ifelse(num == 21, 21, amino_acids[num])
  return(aa)
}

mat <- c()
for (ind in 1:length(genes)){
  gene <- genes[ind]
  id <- ids[ind]
  
  load(sprintf("Data/PHACTboost_Results/PHACTboost_%s.RData", id))
  phactboost <- data
  rm(data)
  load(sprintf("Data/EVE_Results/%s.RData", id))
  eve <- data
  rm(data)
  load(sprintf("Data/REVEL_Results/%s.RData", id))
  revel <- sub
  revel <- as.data.frame(revel)
  revel$label <- matrix(0, length(revel$V1), 1)
  revel$label[which(revel$V3>0.5)] <- "pathogenic"
  revel$label[which(revel$V3<0.5)] <- "neutral"
  rm(sub)
  load(sprintf("Data/AlphaMissense_Results/%s.RData", id))
  alphamis <- sub
  rm(sub)
  cpt1 <- read.csv(sprintf("Data/CPT1_Results/%s_HUMAN.csv", gene))
  
  revel_var <- revel[,1]
  vect_r <- c()
  for (i in 1:length(revel_var)){
    vec <- unlist(strsplit(revel_var[i], "-"))[2]
    vec <- unlist(strsplit(vec, "/"))
    ref <- substr(vec[1], nchar(vec[1]), nchar(vec[1]))
    alt <- vec[2]
    pos <- as.numeric(substr(vec[1], 1, (nchar(vec[1])-1)))
    
    vect_r <- rbind(vect_r, c(pos, ref, alt))
  }
  
  alphamis_var <- alphamis[,2]
  vect_am <- c()
  for (i in 1:length(alphamis_var)){
    vec <- alphamis_var[i]
    ref <- substr(vec, 1, 1)
    alt <- substr(vec, nchar(vec), nchar(vec))
    pos <- as.numeric(substr(vec, 2, (nchar(vec)-1)))
    
    vect_am <- rbind(vect_am, c(pos, ref, alt))
  }
  
  phactboost$vars <- paste(phactboost$Ref_AA, phactboost$Positions, phactboost$Alt_AA, sep = "")
  cpt1$vars <- cpt1$mutant
  
  cpt1_vals <- cpt1$CPT1_score
  # Calculate the 40th and 60th percentiles
  threshold_40 <- quantile(cpt1_vals, 0.4)
  threshold_60 <- quantile(cpt1_vals, 0.6)
  cpt1$label <- matrix(0, length(cpt1$mutant), 1)
  cpt1$label[which(cpt1_vals<=threshold_40)] <- "pathogenic"
  cpt1$label[which(cpt1_vals>=threshold_60)] <- "neutral"
  cpt1$label[which(cpt1$label==0)] <- "uncertain"
  
  eve$vars <- paste(eve$wt_aa, eve$position, eve$mt_aa, sep = "")
  revel$vars <- paste(vect_r[,2], vect_r[,1], vect_r[,3], sep = "")
  alphamis$vars <- paste(vect_am[,2], vect_am[,1], vect_am[,3], sep = "")
  rm(vect_am, vect_r)
  
  
  if (gene == "MSH2"){
    dms <- readxl::read_xlsx("Data/DMS/MSH2_IDown.xlsx")
    dms <- dms[-which(is.na(dms$`LOF score`)),]
    dms$vars <- dms$Variant
    dms$Score <- dms$`LOF score`
  } else if (gene == "A4") {
    dms <- readxl::read_xlsx("Data/DMS/A4_IDown.xlsx")
    dms$Pos <- dms$Pos + 671
    dms$vars <- paste(dms$WT_AA, dms$Pos, dms$Mut, sep = "")
    dms <- dms[-which(dms$Mut=="*"),]
    dms$Score <- dms$nscore
  } else if (gene == "SYUA") {
    dms <- readxl::read_xlsx("Data/DMS/SYUA_Data.xlsx")
    vect <- c()
    for (i in 2:141){
      sub <- dms[,i]
      val <- as.numeric(unlist(sub))
      names <- unlist(dms[,1])
      ref <- names[which(sub==0)]
      vv <- paste(ref, (i-1), names, sep = "")
      vect <- rbind(vect, cbind(vv, val))
    }
    vect <- as.data.frame(vect)
    dms <- vect
    colnames(dms) <- c("vars", "Score")
    dms <- dms[-which(dms$Score=="NaN"),]
  } else if (gene == "VKOR1"){
    dms <- read.csv("Data/DMS/VKOR1_IDown.csv")
    dms$vars <- dms$variant
    dms$Score <- dms$abundance_score
  } else if  (gene == "P53") {
    dms <- readxl::read_xlsx("Data/DMS/P53_Giacomelli_IDown.xlsx")
    dms$vars <- dms$Allele
    dms$Score <- dms$`A549_p53WT_Nutlin-3_Z-score`
  } else if (gene == "PTEN") {
    dms <- readxl::read_xlsx("/Users/nurdankuru/Desktop/PHACTboost_Files/PHACTboost_Datasets/DMS/PTEN_Data_Mighell.xlsx")
    dms$vars <- dms$`Variant (one letter)`
    dms$Score <- dms$Imputed_Score
    el <- which(dms$Score=="NA")
    dms <- dms[-el, ]
  } else if (gene == "ADRB2") {
    dms <- readxl::read_xls("Data/DMS/ADRB2_IDown.xls")
    cond <- 0
    dms <- dms[which(dms$Condition==cond),]
    dms$vars <- paste("A", dms$Pos, dms$AA, sep = "")
    dms$Score <- dms$Norm
  } else  if (gene == "BRCA1") {
    dms <- readxl::read_xlsx("Data/DMS/BRCA1_Findlay_IDown.xlsx")
    dms <- dms[which(dms$consequence=="Missense"),]
    dms$vars <- paste(dms$aa_ref, dms$aa_pos, dms$aa_alt, sep = "")
    dms$Score <- dms$function.score.mean
  }
  
  com <- dms$vars
  if (gene=="ADRB2"){
    substr(phactboost$vars, 1,1) <- "A"
    substr(dms$vars, 1,1) <- "A"
    substr(alphamis$vars, 1,1) <- "A"
    substr(cpt1$vars, 1,1) <- "A"
    substr(revel$vars, 1,1) <- "A"
    substr(eve$vars, 1,1) <- "A"
  } 
  comb <- cbind(gene, id, com, phactboost$PhactBoost_Scores[match(com, phactboost$vars)],dms$Score[match(com, dms$vars)],
                eve$EVE_scores_ASM[match(com, eve$vars)], eve$EVE_classes_75_pct_retained_ASM[match(com, eve$vars)],
                revel$V3[match(com, revel$vars)], revel$label[match(com, revel$vars)],
                alphamis$V3[match(com, alphamis$vars)], alphamis$V4[match(com, alphamis$vars)],
                cpt1$CPT1_score[match(com, cpt1$vars)], cpt1$label[match(com, cpt1$vars)])
  
  comb <- as.data.frame(comb)
  colnames(comb) <- c("gene", "id", "vars", "PHACTboost", "DMS", "EVE", "EVE_Label", "REVEL", "REVEL_Label",
                      "AlphaMissense", "AlphaMissense_Label", "CPT1", "CPT1_Label")
  i1 <- which(is.na(comb$EVE))
  i2 <- which(is.na(comb$REVEL))
  i3 <- which(is.na(comb$AlphaMissense))
  i4 <- which(is.na(comb$CPT1))
  ii <- unique(c(i1, i2, i3, i4))
  if (length(ii)>0){
    comb <- comb[-ii,]
  }

  if (type==0){
    keep <- c()
    labels <- cbind(comb$EVE_Label, comb$REVEL_Label, comb$AlphaMissense_Label, comb$CPT1_Label)
    for (i in 1:length(comb$gene)){
       labs <- labels[i, ]
       nets <- length(unique(c(grep("eutral", labs), grep("enign", labs))))
       pats <- length(unique(c(grep("athogen", labs))))
    
       if (nets==2 && pats==2){
          keep <- c(keep, i)
      }
     }
     comb <- comb[keep,]
  }
  save(comb, file = sprintf("DMS_%s_%s.RData", gene, id))
  
  rm(dms, com, comb, eve, phactboost, revel, cpt1)
}

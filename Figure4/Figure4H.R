library(ape)
library(bio3d)
library(ggtree)
library(paletteer)

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

vars <- c("Q9UQ88-93R/W", "Q6ZTR7-22E/K", "Q07666-421P/L", "O60500-1M/T")

chosen_pos <- c()
for (var in vars){
  k <- k+1
  var <- vars[k]
  id <- unlist(strsplit(var, "-"))[1]
  rem <- unlist(strsplit(var, "-"))[2]
  alt <- unlist(strsplit(rem, "/"))[2]
  rem <- unlist(strsplit(rem, "/"))[1]
  ref <- substr(rem, nchar(rem), nchar(rem))
  pos <- as.numeric(substr(rem, 1, (nchar(rem)-1)))
  
  file_nwk <- sprintf("/Users/nurdankuru/%s.treefile", id)
  file_fasta <- sprintf("/Users/nurdankuru/%s_MaskedMSA.fasta", id)
  
  chosen_pos <- c(pos, pos)
  
  fasta <- read.fasta(file = file_fasta)
  msa <- fasta$ali
  
  # Read Tree
  tree <- read.tree(file_nwk)
  
  # Location of the leaf of human
  human_plc <- grep(id, tree$tip.label)
  my_label <- matrix(, (length(tree$tip.label)+tree$Nnode), 1)
  my_label[human_plc] <- ".....Homo sapiens"
  
  # A new MSA file including the chosen positions
  write.fasta(alignment = NULL, ids = fasta$id, seqs = fasta$ali[,chosen_pos], file = sprintf("%s.fasta", id), append = FALSE)
  
  aas <- unique(msa[,chosen_pos[1]])
  aas <- sort(aas)
  
  i1 <- which(aas==ref)
  i2 <- which(aas==alt)
  
  cols <- c("white", paletteer_c("ggthemes::Orange Light", 50))
  cols[i1] <- "#8D2526"
  cols[i2] <- "#18401E"
  
  pdf(file = sprintf("%s_2.pdf", id), width = 4, height = 8)
  msaplot(p = ggtree(tree) + geom_tiplab(aes(label=my_label)), fasta = sprintf("%s.fasta", id),
          offset = 1, col = cols)
  dev.off()
}





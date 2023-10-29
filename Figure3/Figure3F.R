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

id <- "P54803"
pos <- c(62, 65)

file_nwk <- sprintf("Data/%s.treefile", id)
file_fasta <- sprintf("Data/%s_MaskedMSA.fasta", id)

fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali

# Read Tree
tree <- read.tree(file_nwk)

# Location of the leaf of human
human_plc <- grep(id, tree$tip.label)
my_label <- matrix(, (length(tree$tip.label)+tree$Nnode), 1)
my_label[human_plc] <- ".....Homo sapiens"

# A new MSA file including the chosen positions
write.fasta(alignment = NULL, ids = fasta$id, seqs = fasta$ali[,pos], file = sprintf("%s.fasta", id), append = FALSE)

cols <- c("white", "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3",
          "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3",
          "#fff9e3", "#fff9e3", "#fff9e3", "#fff9e3")
cols[5] <- "#F2B694"
cols[15] <- "#A9D9A7"
cols[2] <- "#D7AECC"
cols[16] <- "#9A99BF"
cols[17] <- "#18401E"

pdf(file = sprintf("Figure3E.pdf"), width = 4, height = 8)
msaplot(p = ggtree(tree, size = 0.1) + geom_tiplab(aes(label=my_label)), fasta = sprintf("%s.fasta", id), 
        offset = 1, color = cols, bg_line = F, width = 0.05)
dev.off()


#!/usr/bin/env Rscript
library(ape)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
library(Peptides)
args = commandArgs(trailingOnly=TRUE)

uniprot_id <- args[4]

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

compute_score <- function(file_nwk, file_rst, file_fasta, output_name, human_id, pos_chosen, parameters) {
  
  # Read tree file
  tr_org <- read.tree(file_nwk)
  x <- read.table(file = file_rst, sep = '\t', header = TRUE, fill = TRUE)
  x[,1] <- str_remove(x[,1], "Node")
  colnames(x)[4:ncol(x)] <- gsub("p_", replacement = "", x = colnames(x)[4:ncol(x)], fixed = TRUE )
  x[,4:ncol(x)] <- matrix(0, nrow(x), (ncol(x)-3))
  el <- which(x$State=="-")
  oth <- which(x$State!="-")
  x[cbind(oth, match(x$State[oth], colnames(x)))] <- 1
  

  # Tree_info: node-node, node-leaf connections
  tree_info <- as.data.frame(as_tibble(tr_org))
  
  # Read fasta file, MSA
  fasta <- read.fasta(file = file_fasta)
  msa <- fasta$ali
  
  # connections_1: Parent node, connections_2: connected node/leaf
  connections_1 <- tree_info$parent
  connections_2 <- tree_info$node
  
  # Names of leaves
  names_all <- tr_org[["tip.label"]]
  msa <- msa[names_all, ]
  # Number of total leaves&nodes
  num_leaves <- length(tr_org[["tip.label"]])
  num_nodes <- tr_org[["Nnode"]]
  
  # Distance between leaves & nodes
  dd_node <- dist.nodes(tr_org)
  dist_leaf <- dd_node[1:num_leaves, 1:num_leaves]
  dist_node <- dd_node[(num_leaves+1):(num_leaves + num_nodes), (num_leaves+1):(num_leaves + num_nodes)]
  
  # Human position (leaf & node)
  h_name <- human_id
  human_codeml <- names_all[grep(pattern = h_name, x = names_all, fixed = TRUE)]
  leaf_human <- tree_info[which(tree_info$label == human_codeml), "node"]
  human_plc <- leaf_human
  node_human <- tree_info[which(tree_info$label == human_codeml), "parent"]
  nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
  names(nodes_raxml) <- tree_info[num_leaves+1:num_nodes, "node"]
  
  # Total number of positions from ancestralProbs file
  total_pos <- max(x$Site)
  
  # Chosen positions (all or some)
  if (pos_chosen[1] == "all"){
    positions <- 1:total_pos
    score_all <- matrix(0, total_pos, 21)
  } else {
    positions <- pos_chosen
    score_all <- matrix(0, length(positions), 21)
  }
  
  human_gaps <- which(msa[match(human_codeml, row.names(msa))]=="-")
  positions <- setdiff(positions, human_gaps)
  
    
  ####################################################
  ####################################################
  
  # Connections between leaves & nodes
  chosen_leaves <- tree_info[1:num_leaves,c("parent", "node")]
  # Connections between nodes & nodes
  chosen_nodes <- tree_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
  leaf_names <- tree_info$label
  
  human_leaf_len <- as.double(tree_info[human_plc, "branch.length"])
  
  if (num_nodes == 1) {
    d_n <- dist_node + human_leaf_len
  } else {
    d_n <- dist_node[as.character(node_human),] + human_leaf_len
  }
  
  d_l <- dist_leaf[leaf_human,]
  
  # chosen_nodes2: ordered connections (for probability differences)
  chosen_nodes2 <- matrix(0, num_nodes-1, 2)
  
  n1 <- as.numeric(chosen_nodes$parent)
  n2 <- as.numeric(chosen_nodes$node)
  dist_f <- d_n[as.character(n1)]
  dist_s <- d_n[as.character(n2)]
  
  # chosen_nodes2: ordered connections (for probability differences)
  chosen_nodes2[which(dist_f < dist_s), 1] <- n2[which(dist_f < dist_s)]
  chosen_nodes2[which(dist_f < dist_s), 2] <- n1[which(dist_f < dist_s)]
  
  chosen_nodes2[which(dist_f >= dist_s), 1] <- n1[which(dist_f >= dist_s)]
  chosen_nodes2[which(dist_f >= dist_s), 2] <- n2[which(dist_f >= dist_s)]
  
  # Number of nodes between nodes & leaf of human
  nodes_conn <- numeric(num_nodes)
  nodes_conn[node_human-num_leaves] <- 1
  names(nodes_conn) <- names(d_n)
  chs <- c()
  chs2 <- c()
  ##########################
  inds <- chosen_nodes2[chosen_nodes2[,2]==(node_human),1]
  nodes_conn[as.character(inds)] <- 2
  chs <- inds
  
  s0 <- sapply(3:num_leaves, function(i){
    for (j in chs){
      inds <- chosen_nodes2[chosen_nodes2[,2]==j,1]
      if (length(inds)!=0){
        nodes_conn[as.character(inds)] <<- i
        chs2 <- c(chs2, inds)
      }
    }
    chs <<- chs2
    chs2 <- c()
  })
  
  # Number of nodes between leaves & leaf of human
  if (num_nodes == 1) {
    leaves_conn <- nodes_conn*matrix(1, 1, num_leaves)
  } else {
    leaves_conn <- nodes_conn[as.character(chosen_leaves[,1])]
  }  
  #########################################
  
  parameters <- unlist(str_split(parameters, pattern = ","))
  
  score_norm <- t(mapply(function(ps){position_score(ps, x, msa, trim_final, human_id, names_all, tr_org, num_nodes, num_leaves, tree_info, num_nodes, nodes_raxml, num_leaves, total_pos, human_plc, node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf, parameters)}, rep(positions)))
  
  pl <- 1
  scores <- list()
  for (p in 1:length(parameters)) {
    score_norm_with_leaf <- matrix(unlist(score_norm[ ,pl]), nrow = length(positions), ncol = 20, byrow = TRUE)
    score_norm_without_leaf <- matrix(unlist(score_norm[ ,(pl+1)]), nrow = length(positions), ncol = 20, byrow = TRUE)
    score_norm_diversity <- matrix(unlist(score_norm[ ,(pl+2)]), nrow = length(positions), ncol = 1, byrow = TRUE)
    score_norm_nogap <- matrix(unlist(score_norm[ ,(pl+3)]), nrow = length(positions), ncol = 20, byrow = TRUE)
    
    pl <- pl + 4
    
    score_norm_with_leaf <- cbind(1:length(positions), score_norm_with_leaf)
    score_norm_without_leaf <- cbind(1:length(positions), score_norm_without_leaf)
    score_norm_nogap <- cbind(1:length(positions), score_norm_nogap)
    
    colnames(score_norm_with_leaf) <- c("Pos/AA", num_to_aa(1:20))
    colnames(score_norm_without_leaf) <- c("Pos/AA", num_to_aa(1:20))
    colnames(score_norm_nogap) <- c("Pos/AA", num_to_aa(1:20))
    
    filename <- ifelse(parameters[p] == "0", "max05", ifelse(parameters[p] == "X", "max05_Gauss", parameters[p]))
    filename <- ifelse(parameters[p] == "X", "max05_Gauss", parameters[p])
    scores[[sprintf("%s_wl_param_%s", output_name, filename)]] <- score_norm_with_leaf
    scores[[sprintf("%s_wol_param_%s", output_name, filename)]] <- score_norm_without_leaf
    scores[[sprintf("%s_diversity_%s", output_name, filename)]] <- score_norm_diversity
    scores[[sprintf("%s_nogap_param_%s", output_name, filename)]] <- score_norm_nogap
  }

  print("saving scores") 
  save("scores", file = sprintf("Scores/%s_scores.RData", output_name))

}

position_score <- function(ps, x, msa, trim_final, human_id, names_all, tr_org, num_nodes, num_leaves, tree_info, num_nodes_prev, nodes_raxml_prev, num_leaves_prev, total_pos, human_plc, node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf, parameters) {
  position <- ps
  b1 <- position + total_pos*(0:(num_nodes-1))
  TT <- x[b1,]
  
  node_info <- as.numeric(TT[,1])
  sort_node_info <- sort(node_info, decreasing = F, index.return=T)
  TT <- TT[sort_node_info$ix,]
  
  matrix_prob <- matrix(0, num_nodes, 20)
  
  probs <- data.matrix((TT[, (4:ncol(TT))]))
  rownames(probs) <- NULL
  rr <- aa_to_num(colnames(x)[4:ncol(TT)])
  matrix_prob[,rr] <- probs
  matrix_prob <- matrix_prob[nodes_raxml,]
  if (length(matrix_prob)==20){
    names(matrix_prob) <- names(sort(rr))
  } else {
    colnames(matrix_prob) <- names(sort(rr))
  }
  
  msa_upd <- msa
  
  names_msa <- rownames(msa)
  trims <- names_msa[which(msa[,ps]=="-")]
  
  parameters <- unlist(str_split(parameters, pattern = ","))
  
  trims <- unique(trims)
  if (length(grep(human_id, trims))>0){
    ALL_SCORES <- list()
    pl <- 1
    for (parameter in parameters) {
      scores <- list()
      scores$score_with_leaf <- matrix(0,1,20)
      scores$score_without_leaf <- matrix(0,1,20)
      scores$diversity <- 0
      scores$score_nogap <- matrix(0,1,20)
      
      ALL_SCORES <- c(ALL_SCORES, scores)
      names(ALL_SCORES)[pl:(pl+3)] <- paste(parameter, "_", names(ALL_SCORES)[pl:(pl+3)], sep = "")
      pl <- pl + 4
    }
    
    ALL_SCORES$trims <- num_leaves
    ALL_SCORES$div_pos <- 0
    ALL_SCORES$div_pos2 <- 0
  } else {
    if (length(trims)>=1) {
      tree_new <- drop.tip(tr_org, trims)
      tree_new_info <- as.data.frame(as_tibble(tree_new))
      
      num_leaves <- length(tree_new[["tip.label"]])
      num_nodes <- tree_new[["Nnode"]]
      
      nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_new_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
      names(nodes_raxml) <- tree_info[num_leaves+1:num_nodes, "node"]
      
      elim_nodes <- setdiff(nodes_raxml_prev, nodes_raxml)
      if (length(elim_nodes)>0){
        els <- t(mapply(function(jj){which(as.numeric(nodes_raxml_prev) == elim_nodes[jj])}, rep(1:length(elim_nodes))))
        matrix_prob <- matrix_prob[-els,]
      }
      
      if (length(trims)>=1){
        msa_upd <- msa_upd[-t(mapply(function(jj){which(rownames(msa_upd)==trims[jj])}, rep(1:length(trims)))),]
      }
    } else {
      tree_new <- tr_org
      tree_new_info <-  tree_info
    }
    
    ALL_SCORES <- list()
    if (num_leaves == 1){
      score <- matrix(0,1,20)
      aa_f <- aa_to_num(msa_upd[ps])
      if (aa_f!=21){
        score[aa_f] <- 1
      }
      diversity <- 0
      position_vec <- msa_upd[ps]
      position_num <- aa_to_num(position_vec)
      
      pl <- 1
      for (parameter in parameters){
        scores <- list()
        scores$score_with_leaf <- score/(num_nodes+num_leaves)
        scores$score_without_leaf <- score/num_nodes
        
        scores$diversity <- diversity
        num_gaps <- sum(position_num==21)
        scores$score_nogap <- score/(num_nodes+num_leaves-num_gaps)
        
        ALL_SCORES <- c(ALL_SCORES, scores)
        names(ALL_SCORES)[pl:(pl+3)] <- paste(parameter, "_", names(ALL_SCORES)[pl:(pl+3)], sep = "")
        pl <- pl + 4
      }
      
    } else {
      h_name <- human_id
      human_codeml <- names_all[grep(pattern = h_name, x = names_all, fixed = TRUE)]
      leaf_human <- tree_new_info[which(tree_new_info$label == human_codeml), "node"]
      human_plc <- leaf_human
      node_human <- tree_new_info[which(tree_new_info$label == human_codeml), "parent"]
      
      
      dd_node <- dist.nodes(tree_new)
      dist_leaf <- dd_node[1:num_leaves, 1:num_leaves]
      dist_node <- dd_node[(num_leaves+1):(num_leaves + num_nodes), (num_leaves+1):(num_leaves + num_nodes)]
      
      # Connections between leaves & nodes
      chosen_leaves <- tree_new_info[1:num_leaves,c("parent", "node")]
      # Connections between nodes & nodes
      chosen_nodes <- tree_new_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
      leaf_names <- tree_new_info$label
      
      human_leaf_len <- as.double(tree_new_info[human_plc, "branch.length"])
      
      if (num_nodes == 1) {
        d_n <- dist_node + human_leaf_len
      } else {
        d_n <- dist_node[as.character(node_human),] + human_leaf_len
      }
      
      d_l <- dist_leaf[leaf_human,]
      
      # chosen_nodes2: ordered connections (for probability differences)
      chosen_nodes2 <- matrix(0, num_nodes-1, 2)
      
      n1 <- as.numeric(chosen_nodes$parent)
      n2 <- as.numeric(chosen_nodes$node)
      dist_f <- d_n[as.character(n1)]
      dist_s <- d_n[as.character(n2)]
      
      # chosen_nodes2: ordered connections (for probability differences)
      chosen_nodes2[which(dist_f < dist_s), 1] <- n2[which(dist_f < dist_s)]
      chosen_nodes2[which(dist_f < dist_s), 2] <- n1[which(dist_f < dist_s)]
      
      chosen_nodes2[which(dist_f >= dist_s), 1] <- n1[which(dist_f >= dist_s)]
      chosen_nodes2[which(dist_f >= dist_s), 2] <- n2[which(dist_f >= dist_s)]
      
      # Number of nodes between nodes & leaf of human
      nodes_conn <- numeric(num_nodes)
      nodes_conn[node_human-num_leaves] <- 1
      names(nodes_conn) <- names(d_n)
      chs <- c()
      chs2 <- c()
      
      inds <- chosen_nodes2[chosen_nodes2[,2]==(node_human),1]
      nodes_conn[as.character(inds)] <- 2
      chs <- inds
      
      s0 <- sapply(3:num_leaves, function(i){
        for (j in chs){
          inds <- chosen_nodes2[chosen_nodes2[,2]==j,1]
          if (length(inds)!=0){
            nodes_conn[as.character(inds)] <<- i
            chs2 <- c(chs2, inds)
          }
        }
        chs <<- chs2
        chs2 <- c()
      })
      
      # Number of nodes between leaves & leaf of human
      if (num_nodes == 1) {
        leaves_conn <- nodes_conn*matrix(1, 1, num_leaves)
      } else {
        leaves_conn <- nodes_conn[as.character(chosen_leaves[,1])]
      }
      
      
      position_vec <- msa_upd[, ps]
      
      position_num <- aa_to_num(position_vec)
      prob_leaves <- matrix(0, num_leaves, 20)
      prob_leaves[cbind(which(position_num <= 20), position_num[which(position_num <= 20)])] <- 1
      
      gaps <- which(position_num == 21)
      
      diff_leaves <- matrix(0, num_leaves, 20)
      if (num_nodes == 1) {
        diff_leaves <- prob_leaves - do.call(rbind, replicate(num_leaves, matrix_prob, simplify=FALSE))
        diff_nodes <- matrix(0, 1, 20)
        vect_human <- matrix_prob
      } else {
        diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
        diff_nodes <- matrix(0, num_nodes-1, 20)
        diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]
        vect_human <- matrix_prob[node_human - num_leaves,]
      }
      diff_leaves[human_plc,]<- -diff_leaves[human_plc,]
      
      pl <- 1
      ################## weights
      for (parameter in parameters) {
        weights <- weight_fnc(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn, mxx) # Updated with related parameters
        weight_leaf <- weights[1:num_leaves]
        weight_node <- tail(weights,num_nodes)
        
        score <- matrix(0,1,20)
        
        if (num_nodes != 1) {
          s1 <- sapply(1:20, function(ii){
            if (num_nodes == 2) {
              dif_pr <- diff_nodes[ii]
            } else {
              dif_pr <- diff_nodes[1:(num_nodes-1),ii]
            }
            dif_pr[dif_pr<0] <- 0
            sel_node <- chosen_nodes2[1:length(dif_pr), 1] - num_leaves
            score[ii] <<- score[ii] + sum(weight_node[sel_node] * dif_pr)
          })
        }
        
        ### NOVEL 29.03
        aa_f <- position_num[human_plc]
        if (aa_f != 21) {
          vect_human[aa_f]<-0
        }
        
        score_without_leaf <- score
        score_without_leaf <- score_without_leaf + weight_node[(node_human-num_leaves)]*vect_human
        score <- score + weight_node[(node_human-num_leaves)]*vect_human
        
        s2 <- sapply(1:20, function(ii){
          diff_lf <- diff_leaves[1:num_leaves,ii]
          diff_lf[gaps] <-  0
          diff_lf[diff_lf<0] <- 0
          
          s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
          score[ii] <<- score[ii] + s1
        })
        
        aa_f <- position_num[human_plc]
        if (aa_f != 21){
          score[aa_f] <- score[aa_f] + weight_leaf[human_plc]*1
          score_without_leaf[aa_f] <- score_without_leaf[aa_f] + weight_leaf[human_plc]*1
        }
        
        sum_exc_max <- sum(score_without_leaf)-max(score_without_leaf)
        diversity <- (-(length(which(score_without_leaf<0.0001))*0.1)/20+0.1)*(sum_exc_max)
        
        scores <- list()
        scores$score_with_leaf <- (score*0.9 + diversity)/(num_nodes+num_leaves)
        scores$score_without_leaf <- (score_without_leaf*0.9 + diversity)/num_nodes
        
        scores$diversity <- diversity
        num_gaps <- sum(position_num==21)
        scores$score_nogap <- (score*0.9 + diversity)/(num_nodes+num_leaves-num_gaps)
        
        ALL_SCORES <- c(ALL_SCORES, scores)
        names(ALL_SCORES)[pl:(pl+3)] <- paste(parameter, "_", names(ALL_SCORES)[pl:(pl+3)], sep = "")
        pl <- pl + 4
      }
      
    }
    
    ALL_SCORES$trims <- length(trims)/num_leaves
    ALL_SCORES$div_pos <- length(unique(position_vec))
    if (ps==1) {
      un1 <- length(unique(msa[which(msa[,(ps+1)]!="-"),(ps+1)])) + length(unique(position_vec))
    } else if (ps == length(msa[1,])) {
      un1 <- length(unique(msa[which(msa[,(ps-1)]!="-"),(ps-1)])) + length(unique(position_vec))
    } else {
      un1 <- length(unique(msa[which(msa[,(ps-1)]!="-"),(ps-1)])) + length(unique(position_vec)) +
        length(unique(msa[which(msa[,(ps+1)]!="-"),(ps+1)]))
    }
    ALL_SCORES$div_pos2 <- un1
    
  }
  
  
  return(ALL_SCORES)
}


weight_fnc <- function(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn, mxx) {
  # print(parameter)
  if (parameter=="0"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_l_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
    
  } else if (parameter=="0_MinNode"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
    
  } else if (parameter=="0_MinNode_Mix"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node1 <- 1/d_n2
    weight_leaf1 <- 1/d_l2
    
    param <- mean(c(d_n, d_l))
    weight_leaf2 <- exp(-d_l^2/param^2)
    weight_node2 <- exp(-d_n^2/param^2)
    
    weight_node<-sqrt(weight_node1*weight_node2)
    weight_leaf<-sqrt(weight_leaf1*weight_leaf2)
    
  } else if (parameter=="0_MinNode_Mix2"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node1 <- 1/d_n2
    weight_leaf1 <- 1/d_l2
    
    param <- mean(c(d_n, d_l))
    weight_leaf2 <- exp(-d_l^2/param^2)/2
    weight_node2 <- exp(-d_n^2/param^2)/2
    
    weight_node<-sqrt(weight_node1*weight_node2)
    weight_leaf<-sqrt(weight_leaf1*weight_leaf2)
    
  } else if (parameter == "mean"){
    param <- mean(c(d_n, d_l))
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
  } else if (parameter == "median"){
    param <- median(c(d_n, d_l))
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
  } else if (parameter == "X"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_l_n)
    min_n <- min(d_n)
    min_ch <- min(min_l, min_n)
    d_l2 <- d_l-min_ch
    d_n2 <- d_n-min_ch
    weight_leaf <- exp(-d_l^2)/2
    weight_node <- exp(-d_n^2)/2
    weight_leaf[human_plc] <- 1
  } else if (parameter == "CountNodes_1"){      ### NewFunction 1 (15Nov)
    weight_node = (exp(-d_n^2) + 1/nodes_conn)/2
    weight_leaf = (exp(-d_l^2) + 1/leaves_conn)/2
  } else if (parameter == "CountNodes_2"){      ### NewFunction 2 (15Nov)
    weight_node = (exp(-d_n^2) + exp(-nodes_conn^2))/2
    weight_leaf = (exp(-d_l^2) + exp(-leaves_conn^2))/2
  } else if (parameter == "CountNodes_3"){      ### NewFunction 3 (15Nov)
    weight_node = sqrt(exp(-d_n^2)*1/nodes_conn)
    weight_leaf = sqrt(exp(-d_l^2)*1/leaves_conn)
  } else if (parameter == "CountNodes_4"){      ### NewFunction 4 (15Nov)
    weight_node = exp(-(sqrt(d_n*nodes_conn))^2)
    weight_leaf = exp(-(sqrt(d_l*leaves_conn))^2)
  } else if (parameter == "CountNodes_5"){
    param <- mean(c(d_n, d_l))
    weight_node = sqrt(exp(-d_n^2/param^2)*1/nodes_conn)
    weight_leaf = sqrt(exp(-d_l^2/param^2)*1/leaves_conn)
  } else if (parameter == "Equal"){
    
    weight_leaf <- matrix(1,1,length(d_l))
    weight_node <- matrix(1,1,length(d_n))
    
  } else if (parameter == "MinThreshold"){
    max_dis <- max(c(d_l, d_n))
    param_min <- as.double(param_min)    
    weight_node <- (-1+param_min)/max_dis*d_n + 1
    weight_leaf <- (-1+param_min)/max_dis*d_l + 1
  } else if (parameter == "MinThreshold_Gauss"){
    max_dis <- max(c(d_l, d_n))
    param_min <- as.double(param_min)
    param <- sqrt(-max_dis^2/log(param_min))
    
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
    
  } else {
    param <- as.double(parameter)
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
  }
  weights = c(weight_leaf, weight_node)
  
  
  return(weights)
}

csv_file <- compute_score(file_nwk=args[1],file_rst=args[2],file_fasta=args[3], output_name=args[4],human_id=args[5],'all', parameters = args[6])




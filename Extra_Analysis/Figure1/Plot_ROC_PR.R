library(PRROC)
library(AUC)

set.seed(123)

TS1 <- read.csv("https://github.com/CompGenomeLab/PHACTboost_manuscript/blob/main/Figure3/Data/Figure3_TS1.csv")

TS1 <- TS1[,1:4]
TS1 <- TS1[-which(is.na(TS1$AlphaMissense)==1),]

pats <- TS1[which(TS1$variant_info==1),]
nets <- TS1[which(TS1$variant_info==-1),]

ratio <- c(0.0001, 0.001, seq(0.01, 1, 0.01))

keep <- c()
keep2 <- c()
inds <- c()
for (i in 1:length(ratio)){
  num <- round(length(nets$variant_info)*ratio[i])
  num <- min(num, length(pats$variant_info))
  
  ind <- sample(1:length(pats$variant_info), num)
  pats_round <- pats[ind,]
  
  inds <- c(inds, paste(ind, collapse = "-"))
  
  all <- rbind(nets, pats_round)
  
  y <- as.numeric(all$variant_info)
  val1 <- as.numeric(all$PHACTboost)
  val2 <- as.numeric(all$AlphaMissense)
  
  pr2 <- pr.curve(scores.class0 = val1, weights.class0 = (1 * (y == 1)), curve = TRUE)
  pr3 <- pr.curve(scores.class0 = val2, weights.class0 = (1 * (y == 1)), curve = TRUE)
  keep <- rbind(keep, c(i, ratio[i], pr2$auc.integral, pr3$auc.integral))
  
  pp2 <- roc(val1, as.factor(1 * (y == 1)))
  pp3 <- roc(val2, as.factor(1 * (y == 1)))
  keep2 <- rbind(keep2, c(i, ratio[i], auc(pp2), auc(pp3)))
}


library(ggplot2)
library(ggsignif)

PHACTboost_color <- "#771B16"
AlphaMissense_color <- "#505218"

val1 <- as.numeric(keep2[,3])
val2 <- as.numeric(keep2[,4])

df <- data.frame(
  Algorithm = rep(c("PHACTboost", "AlphaMissense"), each = length(val1)),
  Value = c(val1, val2)
)

t_test_result <- t.test(val1, val2)
p_value <- t_test_result$p.value

df$Algorithm <- factor(df$Algorithm, levels = c("PHACTboost", "AlphaMissense"))

p1 <- ggplot(df, aes(x = Algorithm, y = Value, color = Algorithm)) +
      geom_jitter(position = position_jitter(width = 0.3), size = 2) +
      geom_boxplot(width = 0.8, outlier.shape = NA) +
      theme_minimal() +
      labs(title = "",
           x = "Algorithm", y = "AUROC") +
      theme(legend.position = "top",
            axis.text.x = element_text(size = 12),   # Adjust x-axis label size
            axis.text.y = element_text(size = 12),   # Adjust y-axis label size
            axis.title.x = element_text(size = 14),  # Adjust x-axis title size
            axis.title.y = element_text(size = 14)) +  # Adjust y-axis title size
      geom_signif(comparisons = list(c("PHACTboost", "AlphaMissense")), textsize = 10, vjust = -0.5, step_increase = 0.05, map_signif_level = TRUE, color = "black") +
      scale_color_manual(values = c("PHACTboost" = PHACTboost_color, "AlphaMissense" = AlphaMissense_color)) +
      theme(legend.position = "topright") +
      coord_cartesian(ylim = c(0, 1.05))
    
    
val1 <- as.numeric(keep[,3])
val2 <- as.numeric(keep[,4])
    
df <- data.frame(
  Algorithm = rep(c("PHACTboost", "AlphaMissense"), each = length(val1)),
  Value = c(val1, val2)
)
    
t_test_result <- t.test(val1, val2)
p_value <- t_test_result$p.value
    
df$Algorithm <- factor(df$Algorithm, levels = c("PHACTboost", "AlphaMissense"))
    
p2 <- ggplot(df, aes(x = Algorithm, y = Value, color = Algorithm)) +
      geom_jitter(position = position_jitter(width = 0.3), size = 2) +
      geom_boxplot(width = 0.8, outlier.shape = NA) +
      theme_minimal() +
      labs(title = "",
           x = "Algorithm", y = "AUPR") +
      theme(legend.position = "top",
            axis.text.x = element_text(size = 12),   # Adjust x-axis label size
            axis.text.y = element_text(size = 12),   # Adjust y-axis label size
            axis.title.x = element_text(size = 14),  # Adjust x-axis title size
            axis.title.y = element_text(size = 14)) +  # Adjust y-axis title size
      geom_signif(comparisons = list(c("PHACTboost", "AlphaMissense")), textsize = 10, vjust = -0.5, step_increase = 0.05, map_signif_level = TRUE, color = "black") +
      scale_color_manual(values = c("PHACTboost" = PHACTboost_color, "AlphaMissense" = AlphaMissense_color)) +
      theme(legend.position = "topright") +
      coord_cartesian(ylim = c(0, 1.05))


library(cowplot)

combined_plot <- plot_grid(p1, p2, ncol = 2)

ggsave(sprintf("AUROC_AUPR.pdf"), combined_plot, width = 14, height = 10)






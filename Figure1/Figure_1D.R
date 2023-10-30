
library(lightgbm)
library(SHAPforxgboost)
library(data.table)
library(ggplot2)

my_theme1 =  theme_bw() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # no gridlines
  theme(panel.border = element_blank(), axis.line = element_line(color = "black"))+ # no border, just major axis
  theme(axis.text = element_text(color = "black", size = 11), axis.title = element_text(size = 11)) + # axis text size and color
  theme(plot.title = element_text(size = 11, face = "bold")) + 
  theme(text = element_text(size = 11))+
  theme(plot.subtitle = element_text(size=10)) + # subtitle size
  theme(plot.caption=element_text(size=7)) + # caption size, if present
  theme(legend.text=element_text(size=11), legend.position = "right", legend.title = element_blank()) 

my_theme2 <- theme_void() +
  theme(axis.text.x = element_text(color = "black", size = 11, margin = margin(r = 6)),
        axis.text.y = element_text(color = "black", size = 11, hjust = 1, margin = margin(r = 6)),#
        axis.line.x = element_line(color = "black"),
        plot.subtitle = element_text(color = "black", size = 11),
        panel.grid.major.y = element_line(color = "grey90", size = 0.6),
        plot.background = element_rect(fill = "white", color = "white"), plot.margin = margin(c(0,0,5,0)))


load("./TrainingSet.RData")
load("./TestSet.RData")
load("./PHACTboost_FinalModel_TrainTest/lightgbm_replication_1_prediction.RData")
state <- lgb.load("./PHACTboost_FinalModel_TrainTest/lightgbm_replication_1_model.txt")

X_train <- train[, prediction$selected_features]
X_test <- test[,prediction$selected_features]

X_train <- scale(X_train)
X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
colnames(X_test) <- colnames(X_train)

shap_values <- shap.values(xgb_model = state, X_train = as.matrix(X_train))
shap_values_pb <- shap_values$shap_score

mygreen <- rgb(100/256, 185/256, 80/256)
mypurple <- rgb(75/256, 0, 130/256)


pdf(file = "./Figure_1_D.pdf", width = 9.5, height = 5.5)
sh <- shap.plot.summary.wrap2(shap_score = shap_values_pb, X = X_train, top_n = 20, dilute = 100) + my_theme1 +
  scale_color_continuous(low = mygreen, high = mypurple, labels = c("Low", "", "", "", "High")) +
  theme(legend.position = "right", legend.key.height = unit(1, "cm"), legend.key.width = unit(0.2, "cm"), legend.title = element_text(size = 11)) #+
sh$layers[[2]]$aes_params$size <- 2
sh$layers[[3]]$aes_params$size <- 3
sh
dev.off()

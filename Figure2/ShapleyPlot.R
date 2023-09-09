library(shapviz)
library(ggplot2)
library(xgboost)
library(lightgbm)

load("/Users/nurdankuru/Desktop/PHACTboost/Train_gnomAD_Shared.RData")
load("/Users/nurdankuru/Desktop/PHACTboost/Test_gnomAD_Shared.RData")
state <- lgb.load("/Users/nurdankuru/Desktop/PHACTboost/lightgbm_replication_1_model.txt")

all_train <- train
all_test <- test

y_train <- all_train$variant_info
X_train <- all_train[, -c(1:11, 498, grep("SIFT", colnames(all_train)))]
X_train <- X_train[, (!grepl("wol", colnames(X_train)))]

y_test <- all_test$variant_info
X_test <- all_test[, -c(1:11, 498, grep("SIFT", colnames(all_test)))]
X_test <- X_test[, (!grepl("wol", colnames(X_test)))]

parameters <- c("0.1", "0.5", "1", "2", "3",  "5", "mean", "median", "0", paste0("CountNodes_", c(1:5)), "0_MinNode_Mix", "0_MinNode_Mix2", "Equal", "max05_Gauss")
param_choice <- "CountNodes_3"
if (param_choice  %in% parameters) {
  left_params <- parameters[parameters != param_choice]
  param_cols <- which(grepl("param", colnames(X_train)))
  param_selected <- param_cols[which(!grepl(paste0("param_", left_params, collapse = "|"), colnames(X_train)[param_cols]))]
  
  selected_cols <- c(param_selected, which(!grepl("param", colnames(X_train))))
  X_train <- X_train[, selected_cols]
  X_test <- X_test[, selected_cols]
}

elims <- which(apply(X_train, 2, sd) == 0)
X_train <- X_train[,-elims]
X <- data.matrix(X_train)

shp <- shapviz(state, X_pred = X, X = X_train)

pdf(file = "Figure_ShapleyVal.pdf", width = 8, height = 8)
sv_importance(shp, kind = "beeswarm", max_display = 20L,
              viridis_args = list(begin = 0.8, end = 0, option = "viridis"))
dev.off()



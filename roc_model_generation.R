roc_model <- function(mtx, pData, method = c("caret", "h2o")) {
  
  method <- match.arg(method)
  
  if(!is.matrix(mtx)) {
    mtx <- as.matrix(mtx)
  }
  
  if(!is.data.frame(pData)) {
    pData <- as.data.frame(pData)
  }

  if(!"condition" %in% colnames(pData)) {
    stop("pData needs to have a column named 'condition'")
  }
  
  if(!is.factor(pData$condition)) {
    n1 <- unique(pData$condition)[1]
    n2 <- unique(pData$condition)[2]
    pData <- factor(pData$condition, levels = c(n1, n2))
  }
  
  df <- data.frame(pheno = pData$condition,
                   scale(t(mtx)))
  
  df$pheno <- as.numeric(ifelse(df$pheno == levels(df$pheno)[1], 0, 1))
  df$pheno <- factor(df$pheno, levels = c(0,1))
  
  message("Make sure ref group is set as the control group before proceeding")
  
  #### Using caret package for prediction ####
  if (method == "caret") {
    library(caret)
  
  ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 10, allowParallel = FALSE)
  ### Include other algorithms
  
  set.seed(1)
  split <- createDataPartition(df$pheno,
                               p = 0.7,
                               list = FALSE,
                               times = 1)
  
  train <- df[split, ]
  test <- df[-split, ]
  
  set.seed(2)
  model <- train(pheno ~.,
                 method = "rf",
                 data = train,
                 trControl = ctrl)
  
  cm <- confusionMatrix(predict(model, test), test$pheno, positive = levels(test$pheno)[2])
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  
  predict <- predict(model, test, type = "prob")[, 1]
  prediction <- prediction(predict, test$pheno)
  
  auc <- performance(prediction, "auc")@y.values[[1]]
  performance <- performance(prediction, "tpr", "fpr")
  
  roc_obj <- roc(test$pheno, predict)
  ci_auc <- ci.auc(roc_obj)
  
  results_df <- data.frame(
    fpr = unlist(performance@x.values),
    tpr = unlist(performance@y.values)
  )
  
  
  plot <- ggplot(results_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "blue", linewidth = 2) +
    geom_abline(linetype = "dashed", color = "black") +
    labs(
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)"
    ) +
    annotate("text", x = 0.6, y = 0.1, 
             label = sprintf("AUC = %.3f (95%% CI: %.2f-%.2f)\n Sensitivity = %.2f \n Specificity = %.2f", 
                             auc, ci_auc[1], ci_auc[3],
                             sensitivity, specificity), 
             size = 5, color = "black", hjust = .5) +
    theme(legend.position = "none",
          plot.margin = unit(c(.75,.75,.75,.75), "inches"),
          axis.text = element_text(size = 12),  
          axis.text.x = element_text(size = 14),  
          axis.text.y = element_text(size = 14))
  
  
  file_path <- getwd()
  ggsave(paste0(file_path, "/", "roc_plot.pdf"), plot = plot, width = 8, height = 8, dpi = 300)
  
  #### Using h2o package for prediction ####
  } else if(method == "h2o") {
    library(h2o)
    h2o.init()
    
    df_h2o <- as.h2o(df)
    split <- h2o.splitFrame(df_h2o, ratios = 0.7)
    
    train <- split[[1]]
    test <- split[[2]]
    
    predictors <- df[, -1]
    response <- "pheno"
    
    model <- h2o.automl(x = predictors,
                        y = response,
                        training_frame = train,
                        max_models = 10,
                        seed = 123)
    
    best_model <- aml$leader
    predict <- h2o.predict(best_model, test)
    performance <- h2o.performance(best_model, test)
    
    auc <- h2o.auc(performance)
    roc <- h2o.roc(performance)
    
    cm <- as.data.frame(h2o.confusionMatrix(performance))
    sensitivity <- cm[1, 1] / (cm[1, 1] + cm[1, 2])
    specificity <- cm[2, 2] / (cm[2, 2] + cm[2, 1])
    
    fpr <- performance@metrics$thresholds_and_metric_scores$fpr
    tpr <- performance@metrics$thresholds_and_metric_scores$tpr
    
    results_df <- data.frame(
      fpr = fpr,
      tpr = tpr
    )
    
    plot <- ggplot(results_df, aes(x = fpr, y = tpr)) +
      geom_line(color = "blue", linewidth = 2) +
      geom_abline(linetype = "dashed", color = "black") +
      labs(
        x = "1 - Specificity (False Positive Rate)",
        y = "Sensitivity (True Positive Rate)"
      ) +
      annotate("text", x = 0.6, y = 0.1, 
               label = sprintf("AUC = %.3f (95%% CI: %.2f-%.2f)\n Sensitivity = %.2f \n Specificity = %.2f", 
                               auc, ci_auc[1], ci_auc[3],
                               sensitivity, specificity), 
               size = 5, color = "black", hjust = .5) +
      theme(legend.position = "none",
            plot.margin = unit(c(.75,.75,.75,.75), "inches"),
            axis.text = element_text(size = 12),  
            axis.text.x = element_text(size = 14),  
            axis.text.y = element_text(size = 14))
    
    
    file_path <- getwd()
    ggsave(paste0(file_path, "/", "roc_plot.pdf"), plot = plot, width = 8, height = 8, dpi = 300)
    h2o.shutdown(prompt = FALSE)
  }
  
  return(model)
  
}

#### Example usage ####
roc_model(exprs_gse8586, pData_gse8586, method = "caret")

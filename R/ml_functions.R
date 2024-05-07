
calculate_simple_roc <- function(data,
                                 features,
                                 target,
                                 rand_seed = 111,
                                 train_prop = 0.7,
                                 plot = FALSE,
                                 smooth = TRUE){
  ## Data preparation
  set.seed(rand_seed)
  data$Label <- data[[target]]
  train_samples <- sample(rownames(data),round(nrow(data)*train_prop), replace = FALSE)
  test_samples <- setdiff(rownames(data), train_samples)
  train_data <- data[train_samples,]
  test_data <- data[test_samples,]
  if (n_distinct(train_data$Label) == 1 | n_distinct(test_data$Label) == 1){
    stop("Train or test data include only one type of label, choose another random seed")}
  
  # Train and predict
  model <- glm(as.formula(paste("Label~", paste(features, collapse = "+"))) , train_data ,family="binomial")
  preds <- predict(model, test_data, type="response") %>% as.numeric()
  test_data$prediction <- preds
  roc <- pROC::roc(test_data, "Label", "prediction", plot = plot, ci = TRUE, auc = TRUE)
  if (smooth) roc <- pROC::smooth(roc, method = "density")
  return(roc)
}

calculate_mean_roc <- function(data,
                                  features,
                                  target,
                                  rand_seeds = 111:116,
                                  train_prop = 0.7,
                                  plot = TRUE,
                                  smooth = FALSE) {
  roc_curves <- data.frame()
  for (rand_seed in rand_seeds){
    roc <- calculate_simple_roc(data, c(features), "response",
                                rand_seed = rand_seed,
                                train_prop = train_prop,
                                plot = plot,
                                smooth = smooth)
    res <- data.frame(specificity = roc$specificities,
                      sensitivity = roc$sensitivities,
                      auc= roc$auc[[1]], rand_seed)
    
    res$ind <- 1:nrow(res)
    roc_curves <- bind_rows(roc_curves, res)
  }
  
  roc_curves_summarised <- roc_curves %>% group_by(ind) %>% 
    summarise_at(vars(c("specificity", "sensitivity")), funs(mean, sd)) %>%
    ungroup %>% arrange(desc(ind))
  mean_auc <- roc_curves %>% select(rand_seed, auc) %>% 
    unique() %>% pull(auc) %>% mean() %>% round(.,3)
  
  return(list("roc_curves_summarised" = roc_curves_summarised,
              "mean_auc" = mean_auc))
  }
  
plot_mean_roc <- function(roc_curves_summarised, mean_auc = NULL){
  ggp <- ggplot(roc_curves_summarised, aes(specificity_mean, sensitivity_mean)) + 
    geom_line() + scale_x_reverse() +     
    geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "black") + 
    labs(x = 'Specificity', y = 'Sensitivity') + 
    geom_ribbon(aes(ymax = sensitivity_mean + sensitivity_sd,
                    ymin = sensitivity_mean - sensitivity_sd), alpha = 0.4, fill = 'steelblue')
  
  if (!is.null(mean_auc)) ggp <- ggp + labs(subtitle = paste("Mean AUC:", mean_auc))
  return(ggp)
}

.generate_cv_folds <- function(n_samples, n_ml_folds, n_ml_repeats){
  folds <- rsample::vfold_cv(
    data.frame(sample_num = 1:n_samples),
    v = n_ml_folds,
    repeats = n_ml_repeats
  )
  return(folds)
}
## Cross validation ROC curves (based on the function get_binary_predictability)
calculate_cv_roc <- function(data,
                             features,
                             target,
                             n_ml_folds = 5,
                             n_ml_repeats = 10,
                             n_ml_shuffle = 5,
                             rand_seed = 100,
                             smooth_roc = FALSE,
                             show_plot = FALSE){
  
  set.seed(rand_seed)
  data$Label <- data[[target]]
  predictions_for_ROC <- data.frame()
  cv_aucs <- data.frame()
  
  # Generate folds for cross-validation
  folds <- .generate_cv_folds(nrow(data), n_ml_folds, n_ml_repeats)

  # Also train models on shuffled labels, as a null model (0 = no shuffling)
  for (shuffle_index in 0:n_ml_shuffle) {
    tmp <- data  # Copy before shuffling label (if relevant)
    if (shuffle_index > 0) tmp <- tmp %>% mutate(Label = sample(Label)) # Shuffle if needed

    # Cross validation loop
    for (i in 1:nrow(folds)) {
      fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
      train_samples <- folds$splits[[i]]$in_id
      train_data <- tmp[train_samples,]
      test_data <- tmp[-train_samples,]
      
    # When train/test data do not include examples of both categories, no model is trained
      if (n_distinct(train_data$Label) > 1 & n_distinct(test_data$Label) > 1) {
        formula <- paste("Label~", paste(features, collapse = "+")) %>% as.formula()
        model <- glm(formula , train_data ,family="binomial")
        preds <- predict(model, test_data, type="response") %>% as.numeric()  
        predictions_for_ROC <- bind_rows(predictions_for_ROC,
                                         data.frame(shuffle_index =  shuffle_index,
                                                    fold_id = fold_id,
                                                    preds = preds,
                                                    label = test_data$Label))
      } 
    }
  }
  
  # Prepare an ROC plot for the current analyzed variable
  # First, get ROC curves for true + shuffled (null) models
  if (nrow(predictions_for_ROC) > 0) {
    
    rocs <- lapply(0:n_ml_shuffle, function(i) {
      pROC::roc(response = predictions_for_ROC %>% filter(shuffle_index == i) %>% pull(label), 
                predictor	= predictions_for_ROC %>% filter(shuffle_index == i) %>% pull(preds), 
                levels = levels(data$Label),
                # direction = '<',
                ci = TRUE, plot = show_plot)}
    )
    
    names(rocs) <- 0:n_ml_shuffle
    
    # Organize ROC curves with confidence intervals in a table
    df_roc_plot <- data.frame()
    for (i in 0:n_ml_shuffle) {
      if (smooth_roc){
      roc <- pROC::smooth(rocs[[as.character(i)]], method = "binormal")
      roc_auc <- as.numeric(roc$auc)
      } else {
      roc <- rocs[[as.character(i)]]
      roc_auc <-as.numeric(pROC::auc(roc))
      }
      df_roc_plot <- bind_rows(df_roc_plot,
                               data.frame(shuffle_index = i,
                                          shuffled_flag = ifelse(i > 0, 'Null models', 'True model'),
                                          Specificity = roc$specificities,
                                          Sensitivity = roc$sensitivities))
      cv_aucs <- bind_rows(cv_aucs, data.frame(shuffle_index =  i, auc = roc_auc))
    }
  }
  
  return(list("roc" = df_roc_plot,
              "auc" = cv_aucs,
              "rocs" = rocs))
}

plot_cv_roc <- function(df_roc_plot){       
  # Patch to make true model plotted last
  df_roc_plot <- df_roc_plot %>% mutate(shuffle_index = 100 - shuffle_index)
  # Build plot
  p <- ggplot(df_roc_plot, aes(Specificity,Sensitivity, group = shuffle_index)) +
    geom_path(aes(color = shuffled_flag)) +
    scale_color_manual(values = c('Null models' = 'darkgrey', 'True model' = 'steelblue')) +
    theme_bw() +
    coord_equal() +
    scale_x_reverse() +
    geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "black") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank())
  
  return(p)
}

summarise_cv_auc <- function(cv_aucs){
  summary_aucs <- cv_aucs %>%
    mutate(is_shuffled = ifelse(shuffle_index > 0, 'shuffled_model', 'true_model')) %>%
    group_by(is_shuffled) %>%
    summarise(n_folds = n(),
              mean_auc = mean(auc, na.rm = TRUE),
              sd_auc = sd(auc, na.rm = TRUE),
              .groups = 'drop')
  return(summary_aucs)
}
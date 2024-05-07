custom_theme <- function(is_figure_small = FALSE, remove_legend = T){
  if (is_figure_small){
    axis_text = 12
    axis_title = 14
  } else {
    axis_text = 10
    axis_title = 12
  }  
    t <- theme(
          strip.background = element_rect(fill = "white"),
          axis.text = element_text(size = axis_text, color = "black"),
          axis.title = element_text(size = axis_title, color = "black"),
          title = element_text(size = axis_text, color = "black")
    )
    if (remove_legend){
      t <- t + theme(legend.position = "none")
    }
  return(t)
}


plot_feature <- function(data, x_var, y_var, palette = NULL, paired = FALSE){
  is_numeric <- is.numeric(data[[feature]])
  if (is_numeric) ggp <- geom_point_with_corr(data, x_var, y_var)
  else ggp <- geom_boxplot_general(data, x_var, y_var, palette = palette) + 
      ggpubr::stat_compare_means(label.x.npc = 0.5, label = "p.signif", paired = paired) +
      ggpubr::stat_compare_means(label.x.npc = 0.6,label.y.npc = 0, paired = paired)
  return(ggp)
}

geom_point_with_corr <- function(data, x_var, y_var, add_rmse = FALSE, color_var = NULL){
  x_var_label <- gsub("_", " ", x_var) %>% stringr::str_to_title()
  y_var_label <- gsub("_", " ", y_var) %>% stringr::str_to_title()
  ggp <- ggplot(data, aes_string(x_var, y_var)) + 
    geom_smooth(method = "lm", color = "black", linetype = "dashed") +
    ggpubr::stat_cor() + 
    labs(x = x_var_label, y = y_var_label)
  
  if (!is.null(color_var)) {
    ggp <- ggp + geom_point(aes(color = !!sym(color_var)), size = 2, alpha = 0.8)
  } else {
    ggp <- ggp + geom_point(size = 2, alpha = 0.8)
  }
  
  if (add_rmse){
    model <- lm(as.formula(paste(y_var, "~", x_var)), data)
    rmse <- sqrt(mean(model$residuals^2))
    ggp <- ggp + labs(subtitle = paste("RMSE:",signif(rmse,3)))
  }
  return(ggp)
}

geom_boxplot_general <- function(data, x_var, y_var, color_var = NULL, palette = NULL, group_var = "subject"){
  if(is.null(color_var)) color_var <- x_var
  
  x_var_label <- gsub("_", " ", x_var) %>% stringr::str_to_title()
  y_var_label <- gsub("_", " ", y_var) %>% stringr::str_to_title()
  
  ggp <- ggplot(data, aes_string(x_var, y_var, color = color_var, fill =color_var)) + 
    geom_boxplot(alpha = 0.1, lwd = 1) + 
    geom_point(size = 3, alpha= 0.8) + 
    labs(x = x_var_label, y = y_var_label)
  
  if(!is.null(group_var)) ggp <- ggp + geom_line(aes(group = !!sym(group_var)), size = 1, alpha = 0.8)
  
  if (!is.null(palette)){
    ggp <- ggp + 
      scale_color_manual(values = palette) +
      scale_fill_manual(values = palette)}
  
  return(ggp)
}

geom_longitudinal_line <- function(data, y_var,
                                   x_var = "timepoint",
                                   group_var = "subject",
                                   color_var = "subject", 
                                   hline_intercept = 0) {
  ggp <- ggplot(data, 
                aes_string(x_var, y_var, group = group_var, color = color_var)) +
         geom_point() + geom_line()
  
  if (is.numeric(hline_intercept)) {
    ggp <- ggp + geom_hline(yintercept = hline_intercept, linetype = "dashed")}
  
  return(ggp)
  
}
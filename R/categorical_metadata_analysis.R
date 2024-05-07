get_permanova <- function(disct_metadata, 
                          dist_matrix, 
                          dist_metric_name,
                          sample_id_column = 'sample_id', 
                          metadata_cols_to_analyze = NULL, 
                          metadata_cols_fix_effects = NULL,
                          metadata_cols_rand_effects = NULL,
                          fdr_threshold = 0.1,
                          permutations = 9999,
                          plot_labels = NULL,
                          pcoa_fill = NULL,
                          ord = NULL,
                          add_centroid = TRUE,
                          seed = 111) {

    # Required libraries loading
    require(dplyr)
    require(tibble)
    require(vegan)
    
    set.seed(seed)
  
    # Define metadata_cols_to_analyze as all columns if not provided, and check if exist if provided.
    if (!is.null(metadata_cols_to_analyze)) {
      if (!all(metadata_cols_to_analyze %in% colnames(disct_metadata)))
        stop("One of metadata_cols_to_analyze not found in metadata")
    } else {
      metadata_cols_to_analyze <- colnames(disct_metadata) 
      metadata_cols_to_analyze <- metadata_cols_to_analyze[!metadata_cols_to_analyze %in% c(sample_id_column)]
      warning("metadata_cols_to_analyze not provided, all metadata columns will be taken except sample_id_column.\n")
    }
  
    # More data verifications
    if (!sample_id_column %in% colnames(dist_matrix))
      stop("sample_id_column not found in dist_matrix")
    if (!sample_id_column %in% colnames(disct_metadata))
      stop("sample_id_column not found in metadata")
    if (!all(metadata_cols_to_analyze %in% colnames(disct_metadata)))
      stop("One of metadata_cols_to_analyze not found in metadata")
    if (nrow(disct_metadata) != nrow(dist_matrix))
      warning("# of samples is different between the distance matrix and the metadata. Only samples with both features and metadata will be used for the following analysis.\n")
  
    # Make sure all random/fixed effects are indeed in data
    if (!is.null(metadata_cols_fix_effects)) {
      if (!all(metadata_cols_fix_effects %in% colnames(disct_metadata)))
        stop("One of metadata_cols_fix_effects not found in metadata")
    }
    if (!is.null(metadata_cols_rand_effects)) {
      if (!all(metadata_cols_rand_effects %in% colnames(disct_metadata)))
        stop("One of metadata_cols_rand_effects not found in metadata")
    }
  
    # Merge the dataframes
    df <- merge(x = dist_matrix, y = disct_metadata, by = sample_id_column, all.x=FALSE)  %>%
      remove_rownames %>% 
      column_to_rownames(var=sample_id_column)
    
    # Add fixed effects to a string that later will be added to the formula and check that they are not one of the metadata_cols_to_analyze.
    mixed_string = ""
    if (!is.null(metadata_cols_fix_effects)){
      for (fixed_col in metadata_cols_fix_effects){
        # Validate that all metadata_cols_fix_effects not in metadata_cols_to_analyze.
        if (fixed_col %in% metadata_cols_to_analyze)
          stop(paste0("Shared column between metadata_cols_fix_effects and metadata_cols_to_analyze: ", 
                      fixed_col, "\nTry again with different metadata_cols_fix_effects!\n"))
      mixed_string =  paste0(mixed_string, " + tmp_df[,'", fixed_col, "'] ")
      }
    }

    # Run adonis2
    pvals <- c()
    for (col in metadata_cols_to_analyze) {
      message('Working on variable: ', col)
      tmp_df <- df[!is.na(df[,col]), ]
      
      formula_string <- paste("tmp_df[, rownames(tmp_df)] ~ tmp_df[[col]]", mixed_string)
      permanova <- adonis2(as.formula(formula_string),
                           permutations = permutations)
      
      pvals[[col]] <- permanova$`Pr(>F)`[1]
    }
    
    # Call plot_colored_pcoa_discrete for significant features.
    all_fdrs <- p.adjust(pvals, 'fdr')
    plots <- list()
    if(!all(all_fdrs > fdr_threshold)) {
      for (metadata_feature in names(all_fdrs[all_fdrs <= fdr_threshold])){
        p12 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 2, plot_labels, pcoa_fill, ord, add_centroid)
        plots[[paste0(metadata_feature,"_pc1&pc2")]] <- p12
        p13 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 3, plot_labels, pcoa_fill, ord, add_centroid)
        plots[[paste0(metadata_feature,"_pc1&pc3")]] <- p13
        p23 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 2, 3, plot_labels, pcoa_fill, ord, add_centroid)
        plots[[paste0(metadata_feature,"_pc2&pc3")]] <- p23
      }
    }
    
    # Return a list of a data frame of the p-values and the FDRs, and the PCoA plots.
    results = data.frame(unlist(pvals), unlist(all_fdrs)) %>%
      rename(pvals = 1, FDR = 2) %>%
      tibble::rownames_to_column("feature")
      
    return(list(results = results, plots = plots))
}


plot_colored_pcoa_discrete <- function(dist_matrix,
                                       disct_metadata,
                                       sample_id_column,
                                       var_to_color_by,
                                       dist_metric_name,
                                       axis_x = 1,
                                       axis_y = 2,
                                       plot_labels = NULL,
                                       pcoa_fill = NULL,
                                       ord = NULL,
                                       add_centroid = TRUE){
  
    # Required libraries loading
    require(dplyr)
    require(ggplot2)
    require(ape)
    
    # Validation
    if(! var_to_color_by %in% names(disct_metadata))
      stop(var_to_color_by, " column is missing from metadata")
    if (!sample_id_column %in% colnames(dist_matrix))
      stop("sample_id_column not found in dist_matrix")
    if (!sample_id_column %in% colnames(disct_metadata))
      stop("sample_id_column not found in disct_metadata")
    if (nrow(disct_metadata) != nrow(dist_matrix))
      warning("# of samples is different between the distance matrix and the disct_metadata. Only samples with features and disct_metadata will be used for the following analysis.\n")
  
    df <- merge(x = dist_matrix, y = disct_metadata, by = sample_id_column, all.x=FALSE)
    df <- df %>% drop_na(all_of(var_to_color_by)) # TODO: instead of dropping missing values, give them a default color?
    
    # Compute PCoA using the ape package
    dist_matrix <- as.dist(df[df[,sample_id_column]])
    PCOA <- pcoa(dist_matrix)
    pcoa_values <- PCOA$vectors[,1:10] %>%
      as.data.frame() %>%
      tibble::rownames_to_column(sample_id_column)
    
    # Format titles for the axes
    ax_x_title <- paste0('Axis ', axis_x, ' (', round(PCOA$values$Relative_eig[axis_x]*100,1) , '% var. explained)')
    ax_y_title <- paste0('Axis ', axis_y, ' (', round(PCOA$values$Relative_eig[axis_y]*100,1) , '% var. explained)')
    
    # Add metadata to pcoa values, for plotting
    pcoa_values <- pcoa_values %>%
      left_join(disct_metadata, by = sample_id_column)
    
    # Convert discrete values to given plot labels by plot_labels.
    if (! is.null(plot_labels)){
      for (val in names(plot_labels))
        pcoa_values[, var_to_color_by][pcoa_values[, var_to_color_by] == strtoi(val)] <- plot_labels[[val]]
    }
    
    if (!is.null(ord))
      pcoa_values[, var_to_color_by] <- factor(pcoa_values[, var_to_color_by] , levels=ord)
    
    pcoa_values[, "var_to_color_by"] <- pcoa_values[, var_to_color_by]
    
    # Plot
    p <- ggplot(pcoa_values,
                aes_string(x = paste0('Axis.', axis_x),
                           y = paste0('Axis.', axis_y))) +
      geom_point(aes(color = var_to_color_by), alpha = 0.6,  size = 3) +
      ggtitle(paste0(dist_metric_name, ' ~ ', var_to_color_by)) +
      xlab(ax_x_title) +
      ylab(ax_y_title) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.subtitle = element_text(hjust = 0.5)) +
      theme(legend.title = element_text(size = 8)) 
    
    if (add_centroid){
      p <- p + stat_ellipse(aes_string(x = paste0('Axis.', axis_x),
                                       y = paste0('Axis.', axis_y),
                                       fill = var_to_color_by),
                            geom = "polygon",
                            alpha = 0.08) 
    
    means <- pcoa_values %>% 
      group_by(var_to_color_by) %>% 
      summarise(mean_x = mean(!!sym(paste0('Axis.', axis_x))), mean_y = mean(!!sym(paste0('Axis.', axis_y))))
    means[, paste0('Axis.', axis_x)] <- means[, "mean_x"]
    means[, paste0('Axis.', axis_y)] <- means[, "mean_y"]
    means[, var_to_color_by] <- means[, "var_to_color_by"] 
    }
    
    if (!is.null(pcoa_fill))
      p <- p + scale_fill_manual(values = pcoa_fill) + scale_color_manual(values = pcoa_fill)
  return(p)
}
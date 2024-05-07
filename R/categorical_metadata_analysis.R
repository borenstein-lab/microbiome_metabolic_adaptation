###############################################################################
#' The following script contains various functions for the analysis of 
#' microbiome-related feature tables against categorical metadata features.
#'  
#' Last updated: 01/01/2023
#' 
#' Implemented functions:
#' - get_differential_abundance
#' - get_permanova
#' - plot_colored_pcoa_discrete
#' - within_groups_beta_div_boxplot
#' - get_binary_predictability
#' - get_multiclass_predictability
###############################################################################


#' Run MaAsLin2 (differential abundance) in case there are more than 2 discrete 
#' values in the col in metadata_cols_to_analyze.
#' This function will get the differential abundance of each pair and plot all 
#' the results in one boxplot.
#' In case you want to do "one vs. all" comparison you should doe one-hote encoding
#' to the specific column and run the get_differential_abundance function.
#'    
#' @param feat_table Any table with continuous features.
#' @param disct_metadata A table with the metadata features to analyze. A 
#'   'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'   be discrete.
#' @param sample_id_column The name of the sample id column.
#' @param metadata_cols_to_analyze The discrete columns in the disct_metadata 
#'   to run MaAaLin2.
#' @param metadata_cols_fix_effects A vector of columns to add as fixed 
#'   effects (e.g. age, sex etc.). 
#'   If numeric with less than 5 or categorical, will perform one-hot encoding.
#'   NULL to ignore.
#' @param metadata_cols_rand_effects A vector of columns to add as random 
#'   effects (e.g. subject idnetifier in the case of multiple samples per 
#'   subject, or batches in case of a multi-batch dataset). NULL to ignore.
#' @param fdr_threshold FDR threshold used to define significant findings. 
#'   Plots will only be generated for significant findings.
#' @param output_dir the directory to keep MaAsLin2 output.
#' @param return_maaslin_plots Boolean. FALSE (default) will return our own 
#'   plots instead of MaAsLin2 plots.
#' @param keep_maaslin_plots Boolean. FALSE (default) will not generate 
#'   MaAsLin2's default output files.
#' @param normalization_method A normalization method supported by MaAsLin2. 
#'   Defaults to 'TSS'. Other options are: CLR, CSS, NONE, TMM
#' @param transformation_method A transformation method supported by MaAsLin2. 
#'   Defaults to 'LOG'. Other options are: LOGIT, AST, NONE
#' @param analysis_method The model MaAsLin2 will run. Defaults to 'LM'.
#'   Other options are: CPLM, NEGBIN, ZINB
#' @param category_name_map raw category values will be mapped to names 
#'   according to this named vector, for plotting purposes. Also defines the 
#'   order of categories. NULL to skip.
#' @param logscale_y_for_plot log scale the Y axis in the plots (default: FALSE).
#' 
#' @return A list of two elements:
#'   `sig_results` contains statistical test results. Data frame with rows that 
#'    represent significant results. 
#'   `all_results` contains statistical test results. Data frame with rows that 
#'    represent all the results (provide the option to perform additional FDR 
#'    correction if needed,)
#'   `plots` is a list of relevant plots: boxplots that represent the 
#'    differential abundance of feature regarding the discrete metadata feature.
#'    Plots generated for adjusted p-value < fdr_threshold.
#'    The names of the plots list are the unique values (strings of the rows in 
#'    the results dataframe.)
get_differential_abundance <- function(feat_table, 
                                       disct_metadata, 
                                       sample_id_column = 'sample_id', 
                                       metadata_cols_to_analyze = NULL,
                                       metadata_cols_fix_effects = NULL,
                                       metadata_cols_rand_effects = NULL,
                                       fdr_threshold = 0.1,
                                       output_dir = "MaAsLin2_output",
                                       return_maaslin_plots = FALSE,
                                       keep_maaslin_plots = FALSE,
                                       normalization_method = "TSS",
                                       transformation_method = "LOG",
                                       analysis_method = "LM",
                                       category_name_map = NULL,
                                       skip_fdr_correction = FALSE,
                                       logscale_y_for_plot = FALSE,
                                       quiet = FALSE) {
  # Required libraries loading
  require(Maaslin2)
  require(ggplot2)
  require(ggsignif)
  require(grid)
  require(png)
  require(dplyr)
  require(ggpubr)
  require(purrr)
  
  # Data verifications
  
  # In case maaslin plots required, change keep_maaslin_plots to TRUE.
  if (return_maaslin_plots & !(keep_maaslin_plots)){
    message("Can not return maaslin plots if keep_maaslin_plots = FALSE. keep_maaslin_plots was changed to TRUE.\n")
    keep_maaslin_plots = TRUE
  }
  
  # If no metadata_cols_to_analyze provided, run the MaAsLin2 on all the metadata columns.
  if (is.null(metadata_cols_to_analyze))
    metadata_cols_to_analyze = colnames(disct_metadata)
  
  # Check if metadata columns are discrete (1, 0 and NA).
  for (col in metadata_cols_to_analyze){ 
    if (length(unique(disct_metadata[,col])) > 3)
      warning(paste0("Too many values in disct_metadata in the ", col, " column"))
  }
  
  # Create output directory (will be deleted if not required).
  if (!dir.exists(output_dir))
    dir.create(output_dir)
  
  # Convert sample_id column to rownames for MaAsLin2 run.
  disct_metadata <- disct_metadata %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  feat_table <- feat_table %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  
  # Remove samples not in both tables, and order samples identically
  samples_to_include <- intersect(row.names(disct_metadata), row.names(feat_table))
  disct_metadata <- disct_metadata[samples_to_include,]
  feat_table <- feat_table[samples_to_include,,drop=FALSE]
  
  # TODO: warn if some samples were not in both tables
  
  # Rum MaAsLin2
  maaslin2_signif_results <- data.frame() # Initialize
  maaslin2_all_results <- data.frame() # Initialize
  
  for (col in metadata_cols_to_analyze) {
    if (!quiet) message('Working on column: ', col)
    
    # For multi-category variables, we run differential abundance separately for each pair of categories
    for (pair in combn(na.omit(unique(disct_metadata[[col]])), 2, simplify = FALSE)) {
      
      # Keep only the two selected values (also drop missing values)
      tmp_metadata <- disct_metadata[(!is.na(disct_metadata[[col]])) & (disct_metadata[[col]] %in% pair), ]
      
      # Get the final list of samples for which DA will be calculated
      samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table))
      tmp_metadata <- tmp_metadata[samples_to_include,,drop=F]
      tmp_feat_table <- feat_table[samples_to_include,,drop=F]
      tmp_metadata_rownames <- rownames(tmp_metadata)
      
      # Add the metadata_cols_to_analyze as fixed effects for the MaAsLin2 run.
      fix <- c(col)
      
      for (fixed_col in metadata_cols_fix_effects){
        
        # Validate that metadata_cols_fix_effects are not in metadata_cols_to_analyze.
        if (fixed_col %in% metadata_cols_to_analyze)
          stop(paste0("Shared column between metadata_cols_fix_effects and metadata_cols_to_analyze: ", 
                      fixed_col, "\nTry again with different metadata_cols_fix_effects!\n"))

        # Check if fixed column is numeric with less than 5 unique values or categorical. Perform one-hot encoding in these cases.
        if(((is.numeric(tmp_metadata[, fixed_col])) & (n_distinct(tmp_metadata[, fixed_col]) < 5)) | 
           (!is.numeric(tmp_metadata[, fixed_col]))){
          
          # Message if there are more than 5 categories in the metadata_cols_fix_effects because it might affect the FDR correction.
          if (n_distinct(tmp_metadata[ ,fixed_col]) > 5)
            message("Too many unique values in ", fixed_col, " column. It might effect the FDR correction and omit some significant results!! Consider other method for differential abundance analysis.")
          
          # Perform one-hot encoding for columns if metadata_cols_fix_effects if they have more then 2 categories.
          if (n_distinct(tmp_metadata[ ,fixed_col]) > 2){
            message("Performing one-hot encoding on ", fixed_col)
            fix <- c(fix, (tmp_metadata[ ,fixed_col]))
            tmp_metadata <- tmp_metadata %>% mutate(value = 1)  %>% spread(fixed_col, value,  fill = 0 )
            rownames(tmp_metadata) <- tmp_metadata_rownames
            
          } else { # Categorical with 2 values or less. One-hot is not needed 
            fix <- c(fix, fixed_col)
          }
    
        } else { # Numeric case, one-hot is not needed  
          fix <- c(fix, fixed_col)
        }
      }
      
      # Check that all the metadata_cols_rand_effects are not in metadata_cols_to_analyze.
      for (rand_col in metadata_cols_rand_effects){
        if (rand_col %in% metadata_cols_to_analyze)
          stop(paste0("Shared column between metadata_cols_rand_effects and metadata_cols_to_analyze: ", 
                      fixed_col, "\nTry again with different metadata_cols_rand_effects!\n"))
      }
      
      # Fix column names if needed
      colnames(tmp_metadata) <- make.names(colnames(tmp_metadata))
      if (!is.null(metadata_cols_rand_effects))
        make.names(metadata_cols_rand_effects)
      
      invisible(capture.output(tmp_results <- Maaslin2(
        input_data = tmp_feat_table, 
        input_metadata = tmp_metadata, 
        output = paste(output_dir, "/", col, sep=""),
        fixed_effects = make.names(fix),
        random_effects = metadata_cols_rand_effects,
        normalization = normalization_method,
        transform = transformation_method,
        analysis_method = analysis_method,
        plot_scatter = keep_maaslin_plots,
        max_significance = fdr_threshold,
        plot_heatmap = FALSE
      )))
      
      tmp_results$results$pair_a <- as.character(pair[[1]])
      tmp_results$results$pair_b <- as.character(pair[[2]])
      
      # Save significant results only
      maaslin2_signif_results <- bind_rows(
        maaslin2_signif_results,
        tmp_results$results %>% filter(qval <= fdr_threshold)
      )
      
      # Save all results
      maaslin2_all_results <- bind_rows(
        maaslin2_all_results,
        tmp_results$results
      )
    }
  }
  
  # New FDR 
  if (skip_fdr_correction){
    maaslin2_all_results$FDR <- maaslin2_all_results$pval
  } else {
    maaslin2_all_results$FDR <- p.adjust(maaslin2_all_results$pval, 'fdr')
  }
  maaslin2_signif_results <- maaslin2_all_results[maaslin2_all_results$FDR < fdr_threshold,]
  
  # Return significant results only for column in metadata_cols_to_analyze (and not in metadata_cols_fix_effects).
  maaslin2_signif_results <- maaslin2_signif_results[maaslin2_signif_results$metadata %in% metadata_cols_to_analyze, ]
  
  # Generate plots to plots list
  plots <- list()
  
  if (nrow(maaslin2_signif_results) > 0) {
    if (return_maaslin_plots) {
      png_list <- fs::dir_ls(path = output_dir, recurse = TRUE, type = "file", glob = "*.png")
      for (png_dir in png_list){
        print(png_dir)
        plots[[png_dir]] <- rasterGrob(readPNG(png_dir))
      }
    } else {
      
      # Fix feature names to match how maaslin outputs them
      feat_table2 <- feat_table
      names(feat_table2) <- make.names(names(feat_table2))
      
      
      for (col in unique(maaslin2_signif_results$metadata)){
        for (feat in unique(maaslin2_signif_results$feature)){
          
          unique_signif_results <- maaslin2_signif_results[(maaslin2_signif_results$metadata == col) & (maaslin2_signif_results$feature == feat), ]
          if (dim(unique_signif_results)[1] == 0)
            next
            

          tmp_metadata <- disct_metadata %>% drop_na(col)
          samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table2))
          tmp_metadata <- tmp_metadata[samples_to_include,]
          feat_table2 <- feat_table[samples_to_include,]
          
          tmp <- data.frame(
            col_vals = tmp_metadata[[col]],
            feat_vals = feat_table2[[feat]]
          )
          
          if (!is.null(category_name_map)) 
            tmp$col_vals <- factor(category_name_map[as.character(tmp$col_vals)],
                                   levels = unname(category_name_map))
          
          y_m <- max(tmp$feat_vals) 
          tmp$col_vals<- as.character(tmp$col_vals)
          
          # Add number of samples in each group to the x-axis ticks.
          xlabs <- list()
          for (val in as.character(unique(tmp$col_vals)))
            xlabs[[val]] = paste0(val,"\n(", dim(tmp[as.character(tmp$col_vals) == val,])[1], ")")
          
          p <- ggplot(tmp, aes(x = col_vals, group = col_vals, y = feat_vals)) +
            geom_boxplot(color = 'black', fill = 'lightblue', alpha = 0.7, outlier.shape = NA) +
            geom_jitter(height = 0, width = 0.1, alpha = 0.5, color = 'black', size = 2) +
            theme_classic() +
            xlab(col) +
            ylab(feat) +
            #ylim(min(tmp$feat_vals), max(tmp$feat_vals)*1.15) + 
            labs(subtitle = 'MaAsLin2') +
            theme(plot.subtitle = element_text(hjust = 0.5)) +
            scale_x_discrete(labels=xlabs)
          
          # Log scale y axis if needed
          if (logscale_y_for_plot) p <- p + yscale("log10", .format = TRUE)
          
          # Add annotations of significant differences (note patch for case when y axis is log scaled)
          for (i in c(1:dim(unique_signif_results)[1])) {
            y_pos <- y_m*(1-(0.1*(i-1)))
            if (logscale_y_for_plot) y_pos <- y_m+(ifelse(y_m > 0.1, 10, 50)*i*y_m)
            p <- p + ggsignif::geom_signif(annotations=round(unique_signif_results$FDR[i],3), y_position =y_pos, 
                                           xmin = unique_signif_results$pair_a[i], xmax = unique_signif_results$pair_b[i])
          }
          
          plots[[paste(col, feat, collapse = ';')]] <- p
        }
      }
    } 
  }
  
  if (! keep_maaslin_plots){
    unlink(output_dir, recursive=TRUE)
  }
    

  return(list(sig_results = maaslin2_signif_results, all_results = maaslin2_all_results, plots = plots))
}

#' Run PERMANOVA (significance test for group-level differences) with adonis2 
#' (vegan library). 
#' 
#'    
#' @param disct_metadata A table with the metadata features to analyze. A 
#'    'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'    be categorical (or NA).
#' @param dist_matrix Microbiome distance matrix (column names correspond to 
#'    sample id's in the first column and are in the same order).
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `dist_matrix` and `disct_metadata`
#' @param metadata_cols_to_analyze The discrete columns in the disct_metadata 
#'    to run permanova.
#' @param fdr_threshold FDR threshold used to define significant findings. PCoA 
#'    plots will only be generated for significant findings.
#' @param plot_labels PCoA plot labels in case the discrete metadata are numeric values.
#' @param pcoa_fill optional fill color to the PCoA (should be as the same size
#'   of n_distinct(disct_metadata$metadata_cols_to_analyze)). If not provided, 
#'   the PCoA colors will be the defaults.
#' @param ord New order to the x-axis in the PCoA plot (should be the same values
#'    as in the disct_metadata$metadata_col_to_analyze).
#' @add_centroid TRUE (default) if you wand to enter ellipse and centroid center
#' to PCoA plot.
#' @return A list of two elements:
#'    `results` Table of p-values and FDRs for all `sample_id_column`.
#'    `plots` is a list of relevant plots: PCoA colored by discrete feature for all
#'    significant results.

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



#' Creates a PCoA plot, where samples are colored by some discrete variable.
#'   Take only samples that exist in the distance matrix and in the metadata.
#'   
#' @param dist_matrix Microbiome distance matrix (column 
#'    names correspond to sample id's in the first column and are in the same 
#'    order).
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected. 
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `dist_matrix` and `disct_metadata`
#' @param vars_to_color_by discrete feature to test (column in metadata).
#' @param dist_metric_name Name of distance metric (for plot title). 
#' @param axis_x Which PC should be used for the x axis? Default: 1. 
#' @param axis_y Which PC should be used for the y axis? Default: 2.
#' @param plot_labels PCoA plot labels in case the discrete metadata are numeric values.
#'
#' @return PCoA plot

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
    
    # p <- p + geom_point(data = means,
    #                    aes_string(x = paste0('Axis.', axis_x),
    #                               y = paste0('Axis.', axis_y),
    #                               color = var_to_color_by, fill = var_to_color_by),
    #                    size = 6, alpha = 0.9, shape = 24)
    }
    
    if (!is.null(pcoa_fill))
      p <- p + scale_fill_manual(values = pcoa_fill) + scale_color_manual(values = pcoa_fill)
  return(p)
}


#' Create boxplot of distances within all pairs of samples in each group by 
#' metadata_col_to_analyze. The distances compared with one of
#' the two methods:
#' 1. Permanova
#' 2. Wilcoxon test of (two-sided) 
#' Significant annotations was added to the significant 
#' comparisons (after FDR).
#' The significant annotations: *, <= fdr_threshold
#'                              **, <= fdr_threshold/10
#'                              ***, <= fdr_threshold/100
#'                              ****, <= fdr_threshold/100
#'
#' @param disct_metadata A table with the metadata features to analyze. A 
#'   'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'   be discrete.
#' @param dist_matrix Microbiome distance matrix (column 
#'    names correspond to sample id's in the first column and are in the same 
#'    order).
#' @param dist_metric_name the distance metric name (string used in the plot)
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `disct_metadata` and `dist_matrix`
#' @param metadata_col_to_analyze The metadata column to group by.
#' @param colors Vector of colors to the boxplot (default: lightblue)
#' @param ord New order to the x-axis in the boxplot (should be the same values
#'    as in the disct_metadata$metadata_col_to_analyze).
#' @param fdr_threshold FDR threshold used to define significant findings (used for
#'    the significant annotation above the boxplot).
#' @param method one of `permanova` or `wilcoxon` to compare the distances between
#'    the two groups (default: wilcoxon).
#'
#' @return list of significant results and boxplot
within_groups_beta_div_boxplot <- function(disct_metadata, 
                                           dist_matrix, 
                                           dist_metric_name,
                                           sample_id_column = 'sample_id', 
                                           metadata_col_to_analyze,
                                           colors = "lightblue",
                                           ord = NULL,
                                           fdr_threshold = 0.1,
                                           method = "wilcoxon"){
  
  require("reshape2")
  require("dplyr")
  require("ggplot2")
  
  # Validations
  if(! metadata_col_to_analyze %in% names(disct_metadata))
    stop(metadata_col_to_analyze, " column is missing from metadata")
  if (!sample_id_column %in% colnames(dist_matrix))
    stop("sample_id_column not found in dist_matrix")
  if (!sample_id_column %in% colnames(disct_metadata))
    stop("sample_id_column not found in disct_metadata")
  if (nrow(disct_metadata) != nrow(dist_matrix))
    warning("Note # of samples is different between the distance matrix and the disct_metadata
                Only samples with features and disct_metadata will be used for the following analysis.")
  if (!method %in% c("permanova", "wilcoxon"))
    stop("Incorrect method. Try again with method='wilcoxon' or method='permanova'.")
  
  # Keep only samples with values in metadata_col_to_analyze
  orig_nrows <- nrow(disct_metadata)
  disct_metadata <- disct_metadata %>% drop_na(all_of(metadata_col_to_analyze))
  if (nrow(disct_metadata) < orig_nrows)
    print(paste0("Drop ", orig_nrows - nrow(disct_metadata), " samples with missing ", metadata_col_to_analyze, " values."))
  df <- merge(x = dist_matrix, y = disct_metadata, by = sample_id_column, all.x=FALSE) %>%
    remove_rownames %>% 
    column_to_rownames(var=sample_id_column)
  
  # Order by default of ord not provided
  if (is.null(ord))
    ord <- unique(df[,metadata_col_to_analyze])
  
  # Reshape the data frame
  tmp_df <- df[rownames(df)]
  
  melt_df <- melt(as.matrix(tmp_df)) %>% 
    filter(as.character(Var1) != as.character(Var2)) %>% 
    mutate_if(is.factor,as.character) 
  
  # Add the needed metadata to the melted data frame
  disct_metadata[, "Var1"] <- disct_metadata[, sample_id_column]  
  disct_metadata[, "Var2"] <- disct_metadata[, sample_id_column]  
  disct_metadata[, "metadata_col_to_analyze_1"] <- disct_metadata[, metadata_col_to_analyze] 
  disct_metadata[, "metadata_col_to_analyze_2"] <- disct_metadata[, metadata_col_to_analyze] 
  melt_df <- left_join(melt_df, disct_metadata[c("Var1", "metadata_col_to_analyze_1")], by=c("Var1"))
  melt_df <- left_join(melt_df, disct_metadata[c("Var2", "metadata_col_to_analyze_2")], by=c("Var2"))
  melt_df <- melt_df[melt_df$metadata_col_to_analyze_1 == melt_df$metadata_col_to_analyze_2,]
  
  # Order the x-axis 
  melt_df$metadata_col_to_analyze_1 <- factor(melt_df$metadata_col_to_analyze_1 , levels=ord)
  
  # Add number of sample in each group to the x-axis ticks.
  xlabs <- list()
  for (val in ord)
    xlabs[[val]] = paste0(val,"\n(", n_distinct(melt_df[melt_df$metadata_col_to_analyze_1 == val,]$Var1), ")")
  
  melt_df <- melt_df %>%
    rowwise() %>%      # for each row
    mutate(sorted_pairs = paste(sort(c(Var1, Var2)), collapse = "_@_")) %>%  # sort the teams alphabetically and then combine them separating with _@_
    ungroup()  %>%  # forget the row grouping
    distinct(sorted_pairs, .keep_all= TRUE) # drop duplicate to hace each distance only once
  
  # Perform two-sided wilcoxon test between the withon group distances, of each pair of groups.
  results_df = data.frame()
  for (pair in combn(na.omit(unique(melt_df$metadata_col_to_analyze_1)), 2, simplify = FALSE)){
    
    if (method == "permanova"){
      pair_samples <- unique((subset(melt_df, metadata_col_to_analyze_1 %in% c(pair[1], pair[2]))$Var1))
      permanova <- adonis2(df[pair_samples, pair_samples] ~ df[pair_samples, metadata_col_to_analyze],
                           permutations=9999)
      results_df <- rbind(results_df, data.frame(a=c(pair[1]), b=c(pair[2]),pval=c(permanova$`Pr(>F)`[1])))
    }
    
    if (method == "wilcoxon"){
      pair_a <- melt_df[melt_df$metadata_col_to_analyze_1 == pair[1],]$value
      pair_b <- melt_df[melt_df$metadata_col_to_analyze_1 == pair[2],]$value
      wilcox_res <- wilcox.test(pair_a, pair_b, alternative = "two.sided")
      results_df <- rbind(results_df, data.frame(a=c(pair[1]), b=c(pair[2]),pval=c(wilcox_res$p.value)))
    }
  }
  
  # Perform FDR correction and keep only the significant values.
  results_df$FDR <- p.adjust(results_df$pval, 'fdr')
  results_df <- results_df[results_df$FDR < fdr_threshold, ]
  
  # Get the boxplots
  p <- ggplot(melt_df, aes(x = metadata_col_to_analyze_1 , 
                           group = metadata_col_to_analyze_1, 
                           y = value)) +
    geom_boxplot(color = 'black', fill = colors, alpha = 0.7, outlier.shape = NA) +
    ylab(dist_metric_name) +
    xlab(metadata_col_to_analyze) +
    theme_classic() + 
    labs(title=dist_metric_name) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels=xlabs)
  
  # Add significant annotations.
  if (dim(results_df)[1] > 0){
    for (i in c(1:dim(results_df)[1])){
      an = "*"
      if (results_df$FDR[i] < (fdr_threshold / 10))
        an <- "**"
      if (results_df$FDR[i] < (fdr_threshold / 100))
        an <- "***"
      if (results_df$FDR[i] < (fdr_threshold / 1000))
        an <- "****"
      p <- p + ggsignif::geom_signif(annotations=an, y_position =1+(0.1*i), 
                                     xmin = results_df$a[i], xmax =results_df$b[i])
    }
  }
  return(list(results = results_df, plot = p))
}


#' Plot boxplot of values in feat_name (feat_table) grouped by 
#' metadata_col_to_analyze (disct_metadata).
#' Check if values feat_name in feat_table have different distribution in 
#' different groups of metadata_col_to_analyze in disct_metadata.
#' The different distribution is measured by a two-sided wilcoxon test. 
#' Significant annotation were were added according to fdr_threshold.
#' Significance annotations: *, <= fdr_threshold
#'                          **, <= fdr_threshold/10
#'                          ***, <= fdr_threshold/100
#'                          ****, <= fdr_threshold/100
#'
#' @param disct_metadata A table with the metadata features to analyze. A 
#'   'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'   be discrete.
#' @param feat_table A feature table (taxa abundances, pathways, alpha-diversity 
#'    metrics, etc.) with rows as samples and columns as features. A 'sample id' 
#'    column is expected in order to merge the table with the metadata.
#' @param feat_name The column in `feat_table` to analyse.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `disct_metadata` and `feat_table`
#' @param metadata_col_to_analyze The metadata column to group by.
#' @param ord New order to the x-axis in the boxplot (should be the same values
#'    as in the disct_metadata$metadata_col_to_analyze).
#' @param colors Vector of colors to the boxplot (default: lightblue)
#' @param fdr_threshold FDR threshold used to define significant findings (used for
#'    the significant annotation above the boxplot).
#'
#' @return List of significant results and boxplots
 
feat_boxplot_wilcoxon <- function(disct_metadata, 
                                feat_table, 
                                feat_name,
                                sample_id_column = 'sample_id', 
                                metadata_col_to_analyze = NULL,
                                ord = NULL,
                                colors = "lightblue",
                                fdr_threshold = 0.1){
  
  require("reshape2")
  require("dplyr")
  require("ggplot2")
  
  # Validations
  if(! metadata_col_to_analyze %in% names(disct_metadata))
    stop(metadata_col_to_analyze, " column is missing from metadata")
  if(! feat_name %in% names(feat_table))
    stop(feat_name, " column is missing from feat_table")
  if (!sample_id_column %in% colnames(feat_table))
    stop("sample_id_column not found in feat_table")
  if (!sample_id_column %in% colnames(disct_metadata))
    stop("sample_id_column not found in disct_metadata")
  if (nrow(disct_metadata) != nrow(feat_table))
    warning("# of samples is different between the distance matrix and the disct_metadata. Only samples with features and disct_metadata will be used for the following analysis.\n")
  
  # Keep only samples with values in metadata_col_to_analyze
  orig_nrows <- nrow(disct_metadata)
  disct_metadata <- disct_metadata %>% drop_na(all_of(metadata_col_to_analyze))
  if (nrow(disct_metadata) < orig_nrows)
    message("Dropped ", orig_nrows - nrow(disct_metadata), " samples with missing ", metadata_col_to_analyze, " values.")
  
  df <- merge(x = feat_table, y = disct_metadata, by = sample_id_column, all.x=FALSE)
  df[,"metadata_col_to_analyze"] <- df[,metadata_col_to_analyze]
  df["feat_name"] <- df[, feat_name]
  
  # Perform thw wilcoxon test
  wilcoxon_df = data.frame()
  for (pair in combn(na.omit(unique(df$metadata_col_to_analyze)), 2, simplify = FALSE)){
    pair_a <- df[df$metadata_col_to_analyze == pair[1],]$feat_name
    pair_b <- df[df$metadata_col_to_analyze == pair[2],]$feat_name
    wilcox_res <- wilcox.test(pair_a, pair_b, alternative = "two.sided")
    wilcoxon_df <- rbind(wilcoxon_df, data.frame(a=c(pair[1]), b=c(pair[2]),pval=c(wilcox_res$p.value)))
  }
  
  # Perform FDR correction and keep only the significant values.
  wilcoxon_df$FDR <- p.adjust(wilcoxon_df$pval, 'fdr')
  wilcoxon_df <- wilcoxon_df[wilcoxon_df$FDR < fdr_threshold, ]
  
  # Order by default of `ord` if provided
  if (!is.null(ord))
    df$metadata_col_to_analyze <- factor(df$metadata_col_to_analyze , levels=ord)
  
  # Add number of samples in each group to the x-axis ticks.
  xlabs <- list()
  for (val in as.character(unique(df$metadata_col_to_analyze)))
    xlabs[[val]] = paste0(val,"\n(", dim(df[as.character(df$metadata_col_to_analyze) == val,])[1], ")")
  
  # Get the boxplots
  p <- ggplot(df, aes(x = metadata_col_to_analyze , 
                      group = metadata_col_to_analyze, 
                      y = feat_name)) +
    geom_boxplot(color = 'black', fill = colors, alpha = 0.7, outlier.shape = NA) +
    ylab(feat_name) +
    xlab(metadata_col_to_analyze) +
    theme_classic() + 
    labs(title=feat_name) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels=xlabs)
  
  # Add significance annotations.
  ym <- max(df[, feat_name])
  if (dim(wilcoxon_df)[1] > 0){
    
    mult_fact <- 0.1
    if (dim(wilcoxon_df)[1] > 4)
      mult_fact <- 0.05
    
    for (i in c(1:dim(wilcoxon_df)[1])){
      an = "*"
      if (wilcoxon_df$FDR[i] < (fdr_threshold / 10))
        an <- "**"
      if (wilcoxon_df$FDR[i] < (fdr_threshold / 100))
        an <- "***"
      if (wilcoxon_df$FDR[i] < (fdr_threshold / 1000))
        an <- "****"
      p <- p + ggsignif::geom_signif(annotations=an, y_position =ym+(mult_fact*i*ym), 
                                     xmin = wilcoxon_df$a[i], xmax =wilcoxon_df$b[i])
    }
  }
  return(list(results = wilcoxon_df, plot = p))
}

#' Plot heatmap of mean distance between all pairs of samples within group
#' (on the diagonal) and between groups.
#'
#' @param disct_metadata A table with the metadata features to analyze. A 
#'   'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'   be discrete.
#' @param dist_matrix Microbiome distance matrix (column 
#'    names correspond to sample id's in the first column and are in the same 
#'    order).
#' @param dist_metric_name the distance metric name (string used in the plot)
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `disct_metadata` and `dist_matrix`
#' @param metadata_col_to_analyze The metadata column to group by.
#' @param ord New order to the x-axis in the boxplot (should be the same values
#'    as in the disct_metadata$metadata_col_to_analyze).
#' @param low_col the lower color to the heatmap scale (default = lightcyan).
#' @param high_col the higher color to the heatmap scale (default = midnightblue).
#'
#' @return
#' @export
#'
#' @examples
mean_beta_div_heatmap <- function(disct_metadata, 
                                  dist_matrix, 
                                  dist_metric_name,
                                  sample_id_column = 'sample_id', 
                                  metadata_col_to_analyze = NULL,
                                  ord = NULL,
                                  low_col = 'lightcyan',
                                  high_col = 'midnightblue'){
  
  require("reshape2")
  require("dplyr")
  require("ggplot2")
  
  if(! metadata_col_to_analyze %in% names(disct_metadata))
    stop(metadata_col_to_analyze, " column is missing from metadata")
  if (!sample_id_column %in% colnames(dist_matrix))
    stop("sample_id_column not found in dist_matrix")
  if (!sample_id_column %in% colnames(disct_metadata))
    stop("sample_id_column not found in disct_metadata")
  if (nrow(disct_metadata) != nrow(dist_matrix))
    warning("Note # of samples is different between the distance matrix and the disct_metadata
                Only samples with features and disct_metadata will be used for the following analysis.")
  
  # Keep only samples with values in metadata_col_to_analyze
  orig_nrows <- nrow(disct_metadata)
  disct_metadata <- disct_metadata %>% drop_na(all_of(metadata_col_to_analyze))
  if (nrow(disct_metadata) < orig_nrows)
    print(paste0("Dropped ", orig_nrows - nrow(disct_metadata), " samples with missing ", metadata_col_to_analyze, " values."))
  
  
  df <- merge(x = dist_matrix, y = disct_metadata, by = sample_id_column, all.x=FALSE) %>%
    remove_rownames %>% 
    column_to_rownames(var=sample_id_column)
  df <- df[rownames(df)]
  
  melt_df <- melt(as.matrix(df)) %>% 
    filter(as.character(Var1) != as.character(Var2)) %>% 
    mutate_if(is.factor,as.character) 
  
  disct_metadata[, "Var1"] <- disct_metadata[, sample_id_column]  
  disct_metadata[, "Var2"] <- disct_metadata[, sample_id_column]  
  
  disct_metadata[, "metadata_col_to_analyze_1"] <- disct_metadata[, metadata_col_to_analyze] 
  disct_metadata[, "metadata_col_to_analyze_2"] <- disct_metadata[, metadata_col_to_analyze] 
  
  melt_df <- left_join(melt_df, disct_metadata[c("Var1", "metadata_col_to_analyze_1")], by=c("Var1"))
  melt_df <- left_join(melt_df, disct_metadata[c("Var2", "metadata_col_to_analyze_2")], by=c("Var2"))
  
  plot_df <- melt_df %>% group_by(metadata_col_to_analyze_1,metadata_col_to_analyze_2) %>% 
    summarise(mean_beta_diversity=mean(value),
              .groups = 'drop') %>%
    as.data.frame()    
  
  if (! is.null(ord)){
    plot_df$metadata_col_to_analyze_1 <- factor(plot_df$metadata_col_to_analyze_1 , levels=ord)
    plot_df$metadata_col_to_analyze_2 <- factor(plot_df$metadata_col_to_analyze_2 , levels=ord)
  }
  
  p <- ggplot(data = plot_df, aes(x=metadata_col_to_analyze_1, y=metadata_col_to_analyze_2, fill=mean_beta_diversity)) + 
    geom_raster() +
    theme_classic() +
    xlab("") + ylab("") +
    scale_fill_gradient(low = low_col, high = high_col) +
    theme(panel.grid = element_blank()) +
    labs(title=paste0("Mean ", dist_metric_name, " distance between pairs"))
  
  return(p)
}


#' Train machine learning classification models based on all features in the 
#'    feature table, to assess predictability of each binary metadata 
#'    feature. Model evaluation is carried out using cross-validated AUC. 
#'    
#' @param feat_table A feature table (taxa abundances, pathways, alpha-diversity 
#'    metrics, etc.) with rows as samples and columns as features. A 'sample id' 
#'    column is expected in order to merge the table with the metadata.
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `feat_table` and `metadata`
#' @param vars_to_analyze A vector of columns names from the metadata 
#'    table for which the analysis should be carried out. 
#' @param ml_method One of "glm" (logistic regression), or "rf" (random forest)
#' @param n_ml_repeats Auc threshold used to define significant findings 
#' @param n_ml_repeats Number of repeats for repeated cross-validation
#' @param n_ml_folds Number of folds for repeated cross-validation
#' @param n_ml_shuffle Number of repeats for the shuffle process.
#' @param rand_seed Random seed for reproducible results
#' @param stratify_by_col An optional categorical column defining sample groups,
#'    for which each group should be analyzed separately (i.e. stratification)
#' 
#' @return A list with two elements: 
#'    `results` contains summary results: AUC summarized over cv-folds and over 
#'    shuffle repeats.
#'    `plots` list of ROC curve plots (per each measure and per each 
#'    stratification group if exists)
#' @export
get_binary_predictability <- function(feat_table, 
                                      metadata, 
                                      sample_id_column, 
                                      vars_to_analyze, 
                                      ml_method = 'glm',
                                      auc_threshold = 0.7,
                                      n_ml_repeats = 10,
                                      n_ml_folds = 5,
                                      n_ml_shuffle = 10,
                                      rand_seed = 1111,
                                      quiet = FALSE,
                                      stratify_by_col = NULL) {
  # Required libraries loading
  require(rsample)
  require(dplyr)
  require(glmnet)
  require(ranger)
  require(pROC)
  require(ROCR)
  
  # Argument validations 
  if(! all(vars_to_analyze %in% names(metadata)))
    stop("Invalid 'vars_to_analyze' argument. Some of the variables provided are missing from the metadata table")
  if(is.null(vars_to_analyze))
    stop("Invalid 'vars_to_analyze' argument. Cannot be NULL")
  if(! ml_method %in% c('glm', 'rf'))
    stop("Invalid 'ml_method' argument. Only 'glm' and 'rf' supported.")

  if(!is.null(stratify_by_col)){
    if (!stratify_by_col %in% colnames(metadata))
      stop("stratify_by_col not found in metadata")
    if (n_distinct(metadata[,stratify_by_col]) > 10)
      message(paste0("Note that stratify_by_col: ", stratify_by_col, " have more then 10 distinct values, it might affect your FDR."))
  } else {
    # Add dummy stratification column.
    metadata[,"stratify_by_col"] = "no_strat"
    stratify_by_col = "stratify_by_col"
  }
  
  metadata[, "stratify_by_col"] <- metadata[, stratify_by_col]  
  metadata <- metadata %>% drop_na(stratify_by_col)
  
  # Fix names of features (causes errors in some models)
  names(feat_table) <- make.names(names(feat_table))
  vars_to_analyze <- make.names(vars_to_analyze)
  
  # Collect evaluation stats here
  cv_aucs <- data.frame()
  
  # Collect ROC plots
  roc_list <- list()
  
  # Iterate over stratification groups
  for (strat_val in unique(metadata$stratify_by_col)){
    if (strat_val != "no_strat" & !quiet) message('Analyzing group: ', strat_val)
    
    # Iterate over variables to predict
    for (measure in vars_to_analyze) {
      if (!quiet) message("Training model for variable: ", measure)
      
      # Prepare table for training (features + label only, and only samples in current stratification group)
      filtered_feat_table <- feat_table %>%
        left_join(metadata %>% select(all_of(c(sample_id_column, 'stratify_by_col', measure))) %>% rename(Label = 3), 
                  by = sample_id_column) %>%
        filter(stratify_by_col == strat_val) %>%
        filter(!is.na(Label)) %>%
        select(-any_of(c(sample_id_column, 'stratify_by_col'))) %>%
        mutate(Label = factor(Label))
      
      if (nrow(filtered_feat_table) < 30)
        warning('Stratifying to \'', strat_val, '\' group only, and taking only samples with non-missing \'', measure, '\', <30 samples remained. It is not recommended ot run the ml pipeline on such small sample sizes')
      
      # Place holders for predictions, for ROC plots     
      predictions_for_ROC <- data.frame()
      
      # Generate folds for cross-validation
      n_samples <- nrow(filtered_feat_table)
      set.seed(rand_seed)
      folds <- vfold_cv(
        data.frame(sample_num = 1:n_samples),
        v = n_ml_folds,
        repeats = n_ml_repeats
      )
      
      # Also train models on shuffled labels, as a null model (0 = no shuffling)
      for (shuffle_index in 0:n_ml_shuffle) {
        
        # Copy before shuffling label (if relevant)
        tmp <- filtered_feat_table 
        
        # Shuffle if needed
        if (shuffle_index > 0) tmp <- tmp %>% mutate(Label = sample(Label))
      
        # Cross validation loop
        for (i in 1:nrow(folds)) {
          
          fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
          
          # Get train/test samples of this fold
          train_samples <- folds$splits[[i]]$in_id
          train_data <- tmp[train_samples,]
          test_data <- tmp[-train_samples,]
          
          # When train/test data do not include examples of both categories, no model is trained
          if (n_distinct(train_data$Label) > 1 & n_distinct(test_data$Label) > 1) {
            
            # Train model, make predictions and record results
            preds <- train_predict(ml_method, train_data, test_data, n_samples)
            predictions_for_ROC <- bind_rows(
              predictions_for_ROC,
              data.frame(
                shuffle_index =  shuffle_index,
                fold_id = fold_id,
                preds = preds, 
                label = test_data$Label
              )
            )
            
            # Calculate and record AUC
            current_auc <- performance(prediction(preds, test_data$Label, label.ordering = levels(tmp$Label)), "auc")@y.values[[1]]
            
            cv_aucs <- bind_rows(
              cv_aucs,
              data.frame(
                shuffle_index =  shuffle_index,
                group = strat_val,
                measure = measure,
                fold_id = fold_id,
                auc = current_auc
              )
            )
          } 
        } # Completed CV
      } # Completed shuffling 
      
      # Prepare an ROC plot for the current analyzed variable
      # First, get ROC curves for true + shuffled (null) models
      if (nrow(predictions_for_ROC) > 0) {
        rocs <- lapply(0:n_ml_shuffle, function(i) {
          roc(
            response = predictions_for_ROC %>% filter(shuffle_index == i) %>% pull(label), 
            predictor	= predictions_for_ROC %>% filter(shuffle_index == i) %>% pull(preds), 
            levels = levels(filtered_feat_table$Label),
            direction = '<',
            ci = TRUE, plot = FALSE
          )})
        names(rocs) <- 0:n_ml_shuffle
        
        # Organize ROC curves with confidence intervals in a table
        df_roc_plot <- data.frame()
        for (i in 0:n_ml_shuffle) {
          ci_obj <- ci.se(rocs[[as.character(i)]], specificities = seq(0, 1, l=25), boot.n = 100)
          df_roc_plot <- bind_rows(
            df_roc_plot,
            data.frame(shuffle_index = i,
                       shuffled_flag = ifelse(i > 0, 'Null models', 'True model'),
                       x = as.numeric(rownames(ci_obj)),
                       lower = ci_obj[, 1],
                       median = ci_obj[, 2],
                       upper = ci_obj[, 3])
          )
        }
        
        # Define the plot title
        if (strat_val == "no_strat") {
          title <- measure
        } else {
          title <- paste0(measure, " (", strat_val, ")")
        }
        
        # Patch to make true model plotted last
        df_roc_plot <- df_roc_plot %>% mutate(shuffle_index = 100 - shuffle_index)
        
        # Build plot
        p <- ggplot(df_roc_plot, aes(x = x, group = shuffle_index)) +
          geom_path(aes(y = `median`, color = shuffled_flag)) +
          scale_color_manual(values = c('Null models' = 'darkgrey', 'True model' = 'black')) +
          theme_bw() + 
          geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + 
          coord_equal() + 
          scale_x_reverse() +
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = shuffled_flag, alpha = shuffled_flag)) +
          scale_alpha_manual(values = c('Null models' = 0.2, 'True model' = 0.4)) +
          scale_fill_manual(values = c('Null models' = 'lightgrey', 'True model' = 'steelblue')) +
          ggtitle(title) +
          labs(x = 'Specificity', y = 'Sensitivity') +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(legend.title = element_blank())
        
        roc_list[[title]] <- p
      } else {
          if(!quiet) message('To few samples - Cannot train and evaluate models')
      }
    } # Completed iterations over variables of interest
  } # Completed iterations over stratification groups
  
  # Summarise all AUC results over folds
  summary_aucs <- cv_aucs %>%
    mutate(is_shuffled = ifelse(shuffle_index > 0, 'shuffled_model', 'true_model')) %>%
    group_by(group, measure, is_shuffled) %>%
    summarise(n_folds = n(),
              mean_auc = mean(auc, na.rm = TRUE),
              sd_auc = sd(auc, na.rm = TRUE),
              .groups = 'drop') %>%
    tidyr::pivot_wider(id_cols = c(group, measure), 
                       names_from = is_shuffled,
                       values_from = c(mean_auc, sd_auc)) %>%
    mutate(Significant = mean_auc_true_model > auc_threshold)
  
  if (strat_val == "no_strat")
    summary_aucs <- summary_aucs %>% select(-group)
  
  # Report significant findings
  if (! any(summary_aucs$Significant)) {
    message('No significantly well-predicted metadata measures')
  } else {
    message('The following measures were well-predicted by the microbiome using a ', toupper(ml_method), ' model:')
    if (strat_val == "no_strat") {
      message(paste(summary_aucs$measure[summary_aucs$Significant], collapse = ', '))
    } else {
      message(paste0(summary_aucs$measure[summary_aucs$Significant], " (",
                     summary_aucs$group[summary_aucs$Significant], ")", collapse = ', '))
    }
  }
  return(list(results = summary_aucs, plots = roc_list))
}

train_predict <- function(ml_method, train_data, test_data, n_samples){
  # Train logistic regression / RF
  if (ml_method == 'glm') {
    lasso_model <- cv.glmnet(
      x = train_data %>% select(-Label) %>% as.matrix(), 
      y = train_data$Label, 
      grouped = (n_samples > 50), # grouped = FALSE relevant for very small sample sets only
      family = "binomial") 
  } else if (ml_method == 'rf') {
    rf_model <- ranger(Label ~ ., data = train_data, respect.unordered.factors=TRUE, probability = TRUE)
  }
  
  # Make prediction on held out samples
  if (ml_method == 'glm') {
    preds <- predict(lasso_model, 
                     test_data %>% select(-Label) %>% as.matrix(), 
                     s = "lambda.min") %>% unname()
    preds <- preds[, 1]
  } else if (ml_method == 'rf') {
    preds <- predict(rf_model, test_data %>% select(-Label))$predictions
    # Extract predictions of the "positive" (=1) class
    preds <- preds[, levels(train_data$Label)[2]]
  }
  return(preds)
}


#' Train machine learning classification models based on all features in the 
#'    feature table, to assess predictability of each binary metadata 
#'    feature. Model evaluation is carried out using cross-validated AUC. 
#'    
#' @param feat_table A feature table (taxa abundances, pathways, alpha-diversity 
#'    metrics, etc.) with rows as samples and columns as features. A 'sample id' 
#'    column is expected in order to merge the table with the metadata.
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `feat_table` and `metadata`
#' @param vars_to_analyze A vector of columns names from the metadata 
#'    table for which the analysis should be carried out. 
#' @param n_ml_repeats Auc threshold used to define significant findings 
#' @param n_ml_repeats Number of repeats for repeated cross-validation
#' @param n_ml_folds Number of folds for repeated cross-validation
#' @param n_ml_shuffle Number of repeats for the shuffle process.
#' @param rand_seed Random seed for reproducible results
#' @param stratify_by_col An optional categorical column defining sample groups,
#'    for which each group should be analyzed separately (i.e. stratification)
#' 
#' @return A list with two elements: 
#'    `results` contains summary results: AUC summarized over cv-folds and over 
#'    shuffle repeats.
#'    `plots` list of ROC curve plots (per each measure and per each 
#'    stratification group if exists)
#' @export
get_multiclass_predictability <- function(feat_table, 
                                      metadata, 
                                      sample_id_column, 
                                      vars_to_analyze, 
                                      accuracy_threshold = 0.7,
                                      n_ml_repeats = 10,
                                      n_ml_folds = 5,
                                      n_ml_shuffle = 10,
                                      rand_seed = 1111,
                                      quiet = FALSE,
                                      stratify_by_col = NULL) {
  # Required libraries loading
  require(rsample)
  require(dplyr)
  require(glmnet)
  require(ranger)
  require(pROC)
  require(ROCR)
  require(caret)
  
  # Argument validations 
  if(! all(vars_to_analyze %in% names(metadata)))
    stop("Invalid 'vars_to_analyze' argument. Some of the variables provided are missing from the metadata table")
  if(is.null(vars_to_analyze))
    stop("Invalid 'vars_to_analyze' argument. Cannot be NULL")

  if(!is.null(stratify_by_col)){
    if (!stratify_by_col %in% colnames(metadata))
      stop("stratify_by_col not found in metadata")
    if (n_distinct(metadata[,stratify_by_col]) > 10)
      message(paste0("Note that stratify_by_col: ", stratify_by_col, " have more then 10 distinct values, it might affect your FDR."))
  } else {
    # Add dummy stratification column.
    metadata[,"stratify_by_col"] = "no_strat"
    stratify_by_col = "stratify_by_col"
  }
  
  metadata[, "stratify_by_col"] <- metadata[, stratify_by_col]  
  metadata <- metadata %>% drop_na(stratify_by_col)
  
  # Fix names of features (causes errors in some models)
  names(feat_table) <- make.names(names(feat_table))
  vars_to_analyze <- make.names(vars_to_analyze)
  
  # Collect evaluation stats here
  cv_aucs <- data.frame()
  
  # Collect ROC plots
  roc_list <- list()
  
  # Iterate over stratification groups
  for (strat_val in unique(metadata$stratify_by_col)){
    
    # Iterate over variables to predict
    for (measure in vars_to_analyze) {
      if (!quiet) message("Training model for variable: ", measure)
      
      # Prepare table for training (features + label only, and only samples in current stratification group)
      filtered_feat_table <- feat_table %>%
        left_join(metadata %>% select(all_of(c(sample_id_column, 'stratify_by_col', measure))) %>% rename(Label = 3), 
                  by = sample_id_column) %>%
        filter(stratify_by_col == strat_val) %>%
        filter(!is.na(Label)) %>%
        select(-any_of(c(sample_id_column, 'stratify_by_col'))) %>%
        mutate(Label = factor(Label))
      
      # Place holders for predictions, for ROC plots     
      predictions_for_ROC <- data.frame()
      
      # Generate folds for cross-validation
      n_samples <- nrow(filtered_feat_table)
      set.seed(rand_seed)
      folds <- vfold_cv(
        data.frame(sample_num = 1:n_samples),
        v = n_ml_folds,
        repeats = n_ml_repeats
      )
      
      # Also train models on shuffled labels, as a null model (0 = no shuffling)
      for (shuffle_index in 0:n_ml_shuffle) {
        # Copy before shuffling label (if relevant)
        tmp <- filtered_feat_table 
        
        # Shuffle if needed
        if (shuffle_index > 0) tmp <- tmp %>% mutate(Label = sample(Label))
        
        # Cross validation loop
        for (i in 1:nrow(folds)) {
          fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
          
          # Get train/test samples of this fold
          train_samples <- folds$splits[[i]]$in_id
          train_data <- tmp[train_samples,]
          test_data <- tmp[-train_samples,]
          
          # Train true model, make predictions and record results
          preds <- train_predict_multiclass(train_data, test_data, n_samples)
          
          # Get confusion matrix and accuracy
          cm <- confusionMatrix(preds,test_data$Label)
          accuracy <- cm$overall[["Accuracy"]]
          
          cv_aucs <- bind_rows(
            cv_aucs,
            data.frame(
              shuffle_index =  shuffle_index,
              group = strat_val,
              measure = measure,
              fold_id = fold_id,
              accuracy = accuracy
            )
          )
        } # Completed CV
      } # Completed shuffling 
    } # Completed iterations over variables of interest
  } # Completed iterations over stratification groups
  
  # Summarise all AUC results over folds
  summary_aucs <- cv_aucs %>%
    mutate(is_shuffled = ifelse(shuffle_index > 0, 'shuffled_model', 'true_model')) %>%
    group_by(group, measure, is_shuffled) %>%
    summarise(n_folds = n(),
              mean_accuracy = mean(accuracy, na.rm = TRUE),
              sd_accuracy = sd(accuracy, na.rm = TRUE),
              .groups = 'drop') %>%
    tidyr::pivot_wider(id_cols = c(group, measure), 
                       names_from = is_shuffled,
                       values_from = c(mean_accuracy, sd_accuracy)) %>%
    mutate(Significant = mean_accuracy_true_model > accuracy_threshold)
  
  if (strat_val == "no_strat")
    summary_aucs <- summary_aucs %>% select(-group)
  
  # Report significant findings
  if (! any(summary_aucs$Significant)) {
    message('No significantly well-predicted metadata measures')
  } else {
    message('The following measures were well-predicted by the microbiome using a ', toupper(ml_method), ' model:')
    if (strat_val == "no_strat") {
      message(paste(summary_aucs$measure[summary_aucs$Significant], collapse = ', '))
    } else {
      message(paste0(summary_aucs$measure[summary_aucs$Significant], " (",
                     summary_aucs$group[summary_aucs$Significant], ")", collapse = ', '))
    }
  }
  return(list(results = summary_aucs))
}

train_predict_multiclass <- function(train_data, test_data, n_samples){
  # Train  RF
  rf_model <- ranger(Label ~ ., data = train_data, respect.unordered.factors=TRUE, probability = FALSE)
  
  # Make prediction on held out samples
  preds <- predict(rf_model, test_data %>% select(-Label))$predictions

  return(preds)
}


# Read tab-delimited files
read_tsv2 <- function(p, cmnt = "") { 
  return(read_delim(p, comment = cmnt, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)) 
}

# Convert Kraken table to convenient feature table
reorganize_kraken_table <- function(df) {
  require(dplyr)
  require(tibble)
  
  df %>%
    column_to_rownames('Taxon') %>%
    t() %>%
    data.frame() %>%
    rownames_to_column('sample_id')
}

# Convert metaphlan table to convenient feature table
# level should be one of: 'p','c','o','f','g','s'
reorganize_metaphlan_table <- function(df, level = 's') {
  # Regex for required taxonomic level
  rgx_in <- paste0('\\|',level,'__')
  
  # Regex to exclude (lower taxonomy levels)
  levels <- c('k'=1,'p'=2,'c'=3,'o'=4,'f'=5,'g'=6,'s'=7,'t'=8)
  level_out <- names(levels[1+levels[level]])
  rgx_out <- paste0('\\|',level_out,'__')
  
  df %>%
    filter(clade_name == 'UNCLASSIFIED' | grepl(rgx_in, clade_name)) %>%
    filter(! grepl(rgx_out, clade_name)) %>%
    mutate(clade_name = gsub(paste0('^.*',rgx_in), '', clade_name)) %>%
    column_to_rownames('clade_name') %>%
    t() %>%
    data.frame() %>%
    rownames_to_column('sample_id') %>%
    mutate(sample_id = gsub('_metaphlan4_bugs_list', '', sample_id))
}

# Convert to relative abundance
# Assumes first column is sample id's and all the rest are taxa abundances
convert_to_relab <- function(df) {
  require(vegan)
  relab_df <- decostand(df %>% select(-1), method = 'total', MARGIN = 1)
  relab_df <- bind_cols(df %>% select(1), relab_df)
  return(relab_df)
}

# Calculate rarefaction curves
# Note that vegan::rarecurve function returns its own plot (the option to hide it with `tidy` argument is not working)
get_rarefaction_curves <- function(kraken, rarefy_step = 1000000) {
  require(dplyr)
  require(tidyr)
  require(vegan)
  
  rarecurves <- rarecurve(kraken %>% select(-sample_id), step = rarefy_step)
  rarecurves <- bind_rows(rarecurves, .id = 'sample_id')
  rarecurves <- pivot_longer(
    data = rarecurves,
    cols = -sample_id,
    names_to = 'depth', 
    names_prefix = 'N', 
    names_transform = as.integer, 
    values_to = 'species_count', 
    values_drop_na = TRUE
    )
  return(rarecurves)
}

# Create raw read quality plots
plot_quality <- function(df, facet_by = 'direction') {
  # Reorganize table
  df <- df %>%
    group_by(across(all_of(c(facet_by, 'position')))) %>%
    summarise(N = n(),
              q9 = quantile(mean_quality, probs = c(0.09)),
              q25 = quantile(mean_quality, probs = c(0.25)),
              q50 = quantile(mean_quality, probs = c(0.5)),
              q75 = quantile(mean_quality, probs = c(0.75)),
              q91 = quantile(mean_quality, probs = c(0.91))) %>%
    rename(Position = position)
  
  ggplot(df, aes(x = Position, group = Position)) +
    geom_boxplot(aes(ymin = q9, 
                     lower = q25, 
                     middle = q50, 
                     upper = q75, 
                     ymax = q91),
                 stat = "identity", 
                 fill = 'cadetblue3', 
                 color = 'grey30') +
    geom_hline(yintercept = 30, linetype = 'dashed', color = 'darkred') +
    theme_classic() +
    facet_wrap(facet_by, nrow = 1) +
    ggtitle('Quality curves') +
    ylim(c(20,40)) +
    ylab('Quality score') +
    scale_x_continuous(expand = c(0,0)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(strip.background = element_rect(fill = "lightgrey"))
}

get_plot_taxonomy_data <- function(df, rare_otu_cutoff = 0.01){
  df <- df %>%
    tidyr::pivot_longer(cols = -sample_id, names_to = 'taxon', values_to = 'relab') %>%
    # Re-normalize into relative abundance
    group_by(sample_id) %>%
    mutate(relab = relab / sum(relab)) %>%
    ungroup()
  
  # Group rare otu's together for plot simplicity
  rare_otus <- df %>%
    group_by(taxon) %>%
    summarise(mean_relab = mean(relab), .groups = 'drop') %>%
    filter(mean_relab < rare_otu_cutoff) %>%
    pull(taxon)
  
  df <- df %>%
    mutate(otu_grouped = ifelse(taxon %in% rare_otus, 'Others', taxon)) %>%
    group_by(sample_id, otu_grouped) %>%
    summarise(relab = sum(relab), .groups = 'drop')
  
  # Slightly reorder to make the 'Others' category last
  taxa_list <- sort(unique(df$otu_grouped))
  taxa_list <- c(taxa_list[taxa_list != 'Others'], 'Others')
  df$otu_grouped <- factor(df$otu_grouped, levels = taxa_list)
  
  # Sort samples by the abundance of upmost taxon
  first_tax <- taxa_list[1]
  samp_order <- df %>%
    filter(otu_grouped == first_tax) %>%
    arrange(-relab) %>%
    pull(sample_id) 
  # Add samples with 0 abundance of this taxa
  samp_order <- c(samp_order, setdiff(unique(df$sample_id), samp_order))
  df$sample_id <- factor(df$sample_id, levels = samp_order)
  return(df)
}

plot_taxonomy <- function(df, rare_otu_cutoff = 0.01, plot_title = 'Taxonomy barplots', legend_title = 'Taxon', n_legend_cols = 1) {
  df <- get_plot_taxonomy_data(df, rare_otu_cutoff)
  tmp_palette <- colorRampPalette(RColorBrewer::brewer.pal(9,name = 'RdYlGn'))(n_distinct(df$otu_grouped)-1)
  tmp_palette <- c(tmp_palette, 'grey') # For the 'Others' category
  ggplot(df, aes(fill = otu_grouped, y = relab, x = sample_id)) + 
    geom_bar(position = "fill", stat = "identity", color = 'black') +
    scale_fill_manual(values = tmp_palette, name = legend_title) +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    ggtitle(plot_title) +
    xlab('Sample') +
    ylab('Relative abundance') +
    guides(fill = guide_legend(ncol = n_legend_cols)) +
    theme(axis.text.x = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom") +
    theme(legend.direction = "vertical") +
    theme(legend.text = element_text(size = 9)) +
    theme(legend.key.size = unit(9, "points"))
}

plot_pcoa_from_dist_mat <- function(beta_dist, metric_name, axis_x = 1, axis_y = 2) {
  beta_dist <- beta_dist %>%
    tibble::column_to_rownames(var = '...1')
  beta_dist <- as.dist(beta_dist)
  PCOA <- pcoa(beta_dist)
  pcoa_values <- PCOA$vectors[,1:5] %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('sample') %>% 
    mutate(subject = gsub('t[0-9]','',sample)) %>%
    mutate(time_point = gsub('.*t','',sample))
  
  ax_x_title <- paste0('Axis ', axis_x, ' (', round(PCOA$values$Relative_eig[axis_x]*100,1) , '% variance explained)')
  ax_y_title <- paste0('Axis ', axis_y, ' (', round(PCOA$values$Relative_eig[axis_y]*100,1) , '% variance explained)')
  
  # Plot
  tmp_palette <- colorRampPalette(RColorBrewer::brewer.pal(9,name = 'BrBG'))(n_distinct(pcoa_values$subject))
  ggplot(pcoa_values %>% arrange(time_point), 
         aes_string(x = paste0('Axis.', axis_x), 
                    y = paste0('Axis.', axis_y), 
                    fill = 'subject')) +
    geom_path(aes(group = subject), 
              color = 'darkgrey',
              arrow = arrow(length = unit(8, "points"), type = "closed")) +
    geom_point(color = 'black', alpha = 0.8, shape = 21, size = 3) +
    # geom_text(aes(label = time_point), size = 3) +
    scale_fill_manual(values = tmp_palette, name = 'Subject') +
    ggtitle(paste0('PCoA - ', metric_name, ' (Axes ', axis_x, ' vs. ', axis_y, ')')) +
    xlab(ax_x_title) +
    ylab(ax_y_title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}


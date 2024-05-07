# Extract subject ID and collection timepoint from each sample name
extract_metadata_from_sample_name <- function(data, sample_var = "sample"){
  data <- data %>% 
    mutate(subject = stringr::str_extract(!!sym(sample_var),"[0-9]+"),
           subject = paste0("subject_",subject),
           timepoint = stringr::str_extract(!!sym(sample_var), "PRE|POST"),
           timepoint = factor(tolower(timepoint),levels = c("pre", "post")))
  return(data)
}

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

# Convert to relative abundance
# Assumes first column is sample id's and all the rest are taxa abundances
convert_to_relab <- function(df) {
  require(vegan)
  relab_df <- decostand(df %>% select(-1), method = 'total', MARGIN = 1)
  relab_df <- bind_cols(df %>% select(1), relab_df)
  return(relab_df)
}





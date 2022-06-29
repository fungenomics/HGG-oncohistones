


#' Convert Ensembl gene IDs to gene symbols
ensembl2symbols_safe <- function(genes, ...) {
  
  sym <- lapply(genes, ensembl2symbols, ...)
  empty <- which(unlist(map(sym, ~ identical(unlist(.x), character(0)))))
  sym[empty] <- genes[empty]
  
  return(unlist(sym))
  
}

# bulk RNAseq pipeline helpers -------------------------------------------------

#' Based off of code from Nicolas De Jay. Extract counts for genes of interest from
#' the raptor bulk RNA-seq pipeline and wrangle them into tidy format
#'
#' @param path Character, path to counts file
#' @param goi Character, vector of one or more gene symbols for genes of interest
#' for which to retrieve counts
#'
#' @return Data frame with four columns: sample, gene_expression, gene_ensg, gene_symbol
extract_pipeline_counts <- function(path, goi, long = TRUE) {
  
  counts <- read.table(path, header = T, sep = "\t", check.names = FALSE)
  
  # Our pipeline labels genes using #{ENSG}:#{SYMBOL}
  ensg_to_symbol <- data.frame("gene_id" = rownames(counts) %>% as.character) %>%
    dplyr::mutate("gene_ensg"   = gene_id %>% as.character %>% strsplit(":") %>% sapply(`[[`, 1)) %>%
    dplyr::mutate("gene_symbol" = gene_id %>% as.character %>% strsplit(":") %>% sapply(`[[`, 2))
  
  # Convert all genes from #{SYMBOL} to #{ENSG}:#{SYMBOL}
  genes.symbol <- data.frame(genes = goi)
  colnames(genes.symbol) = "gene_symbol"
  
  genes.ensg <- genes.symbol %>% left_join(ensg_to_symbol) %>% dplyr::select(gene_id) %>% 
    filter(!is.na(gene_id))
  
  # Subset counts table
  counts.subset <- counts[genes.ensg %>% pull(gene_id) %>% as.character, , drop = F] %>%
    as.matrix()
  
  if (!long) return(counts.subset)
  
  counts.subset <- counts.subset %>%
    reshape2::melt() %>%
    setNames(c("gene_id", "sample", "gene_expression")) %>%
    left_join(ensg_to_symbol, by = "gene_id") %>%
    dplyr::select(-gene_id)
  
  return(counts.subset)
  
}









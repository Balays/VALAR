
### Normalize coverage
normalise_coverage <- function(mapped.cov,
                               method = 'sum', norm_to = 1e6,
                               contig_sizes = NULL, fasta.ref = NULL,
                               return = 'norm_cov') {
  
  # Calculate total counts per sample
  total_counts    <- mapped.cov[, .(total = sum(count)), by = sample]
  
  # Normalize coverage
  norm_mapped.cov <- mapped.cov[total_counts, on = 'sample']
  
  if (method == 'sum') {
    norm_mapped.cov[, norm_count := count / total * norm_to]
    
  } else if (method == 'genome') {
    
    # Check if contig_sizes is provided, otherwise use fasta.ref
    if (is.null(contig_sizes)) {
      if (is.null(fasta.ref)) {
        stop("Either contig_sizes or fasta.ref must be provided.")
      }
      
      # Load the reference FASTA file
      fasta  <- seqinr::read.fasta(fasta.ref)
      # Get contig sizes
      contig_sizes <- unlist(lapply(fasta, length))
    }
    contig_sizes <- data.table(seqnames=names(contig_sizes), contig_size=contig_sizes)
    
    norm_mapped.cov <- merge(norm_mapped.cov, contig_sizes, by='seqnames')
    
    norm_mapped.cov[, norm_count := count / (total / contig_size)]
    
  } else {
    stop("Invalid method. Currently, only 'sum' is supported.")
    
  }
  
  if (return == 'norm_cov') {
    norm_mapped.cov[, contig_size := NULL]
    
    # Drop the total column and return the result
    norm_mapped.cov[, total := NULL]
    norm_mapped.cov[, count := norm_count]
    norm_mapped.cov[, norm_count := NULL]
    
    return(norm_mapped.cov)
  } else if (return == 'summary') {
    
    colstouni <- setdiff(colnames(norm_mapped.cov),
                         c('strand', 'count', 'norm_count', 'pos'))
    norm_cov_summary <- unique(norm_mapped.cov[,
                                               ..colstouni])
    norm_cov_summary[, average_coverage := total / contig_size]
  }
}


### Associate each position with a window

window.cov <- function(mapped_cov, contig_sizes = NULL, fasta.ref = NULL, window_size, window_step) {
  require(data.table)
  require(Biostrings)
  
  # Check if contig_sizes is provided, otherwise use fasta.ref
  if (is.null(contig_sizes)) {
    if (is.null(fasta.ref)) {
      stop("Either contig_sizes or fasta.ref must be provided.")
    }
    
    # Load the reference FASTA file
    fasta  <- seqinr::read.fasta(fasta.ref)
    # Get contig sizes
    contig_sizes <- unlist(lapply(fasta, length))
  }
  
  mapped_cov <- mapped_cov[,.(sample, seqnames, strand, pos, count)] %>%
    spread(sample, count,  fill=0)  ;  setDT(mapped_cov)
  
  # Create a data.table with all positions
  all_windows <- data.table(seqnames= names(contig_sizes), pos=1:contig_sizes)
  all_windows[, window_start := floor((pos - 1) / window_step) * window_step + 1]
  
  # Calculate window end
  all_windows[, window_end := pmin(window_start + window_size - 1, contig_sizes[seqnames])]
  
  all_windows <- rbind(
    data.table(all_windows, strand = '+'), 
    data.table(all_windows, strand = '-') 
  ) 
  #setDT(all_windows)
  
  # window start
  mapped_cov[, window_start := floor((pos - 1) / window_step) * window_step + 1]
  mapped_cov[, window_end := pmin(window_start + window_size - 1, contig_sizes[seqnames])]
  
  # Merge with mapped_cov
  merged_cov <- merge(all_windows, mapped_cov,
                      by = c("seqnames", "strand", "pos", "window_start", "window_end"),
                      all = TRUE, allow.cartesian=F)
  
  return(merged_cov)
}


### Summarise the coverages in each window (mean, sd, sum, etc) in each sample

win_coverage_summary <- function(merged_cov,  stat_by=metacols) {
  
  stat_by <- c('sample', 'seqnames', 'strand', 'window_start', "window_end")
  
  win.cov.sum <-
    merged_cov[,
               .(sum_coverage    = sum (count,   na.rm = TRUE),
                 mean_coverage   = mean(count,   na.rm = TRUE),
                 sd_coverage     = sd  (count,   na.rm = TRUE),
                 median_coverage = median(count, na.rm = TRUE),
                 min_coverage    = min (count,   na.rm = TRUE),
                 max_coverage    = max (count,   na.rm = TRUE)),
               by = .(sample, seqnames, strand, window_start, window_end)]
  
  win.cov.sum[,varcoeff_coverage := sd_coverage / mean_coverage]
  
  return(win.cov.sum)
  
}

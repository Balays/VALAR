
library(data.table)
library(fitdistrplus)
library(sn)  # For skew-normal fitting
library(ggplot2)
library(pracma)  # For peak detection
library(BiocParallel)  # For parallel processing

# Convert TSSr processed matrix to data.table
#TSS.dt <- as.data.table(myTSSr@TSSprocessedMatrix)

## use the raw counts instead of TPM normalization - we'll carry out the fitting on each sample group so no normalization is neccessary
TSS.dt <- fread('TSS.raw.dt.txt', na.strings='')


## Consensus Clusters from TSSr
filtered_clusters <- fread(paste0(outdir, '/TSSr.filtered_clusters.tsv'), na.strings = '')


# Define parallel processing parameters
bp_param <- MulticoreParam(workers = parallel::detectCores() - 1)  # Use all but one core

# Define parallel processing parameters
bp_param <- MulticoreParam(workers = parallel::detectCores() - 1)  # Use all but one core

# Option to enable or disable KDE-based peak detection
use_KDE <- FALSE  # Set to FALSE to disable KDE peak separation

# Function to process a single cluster
process_cluster <- function(sample, cl) {
  message(paste("Processing cluster", cl, "in sample", sample))
  
  cluster_data <- filtered_clusters[group == sample & cluster == cl]
  
  region_data <- TSS.dt[chr == cluster_data$seqnames & 
                          strand == cluster_data$strand & 
                          pos >= cluster_data$cluster_start & 
                          pos <= cluster_data$cluster_end, .(chr, pos, strand, count = get(sample))]
  
  if (nrow(region_data) <= 1) {
    message(paste("Skipping empty cluster", cl, "in sample", sample))
    return(NULL)
  }
  
  # Expand positions based on read counts
  pos_vec <- tryCatch({
    unlist(lapply(1:nrow(region_data), function(i) rep(region_data$pos[i], times = as.integer(round(region_data$count[i])))))
  }, error = function(e) {
    message(paste("Error in expanding positions for cluster", cl, "in sample", sample, "-", e$message))
    return(NULL)
  })
  
  if (is.null(pos_vec) || length(pos_vec) < 5) {
    message(paste("Skipping small cluster", cl, "in sample", sample))
    return(NULL)
  }
  
  # Detect peaks using Kernel Density Estimation (KDE) if enabled
  peaks <- NULL
  if (use_KDE) {
    peaks <- tryCatch({
      density_data <- density(pos_vec, bw = "nrd0")
      findpeaks(density_data$y, nups = 1, ndowns = 1, npeaks = Inf, sortstr = TRUE)
    }, error = function(e) {
      message(paste("Error in peak detection for cluster", cl, "in sample", sample, "-", e$message))
      return(NULL)
    })
  }
  
  # If KDE is disabled or no peaks are found, use the dominant position as a single peak
  if (!use_KDE || is.null(peaks) || nrow(peaks) == 0) {
    peaks <- matrix(nrow = 1, ncol = 2)
    peaks[1, 2] <- which.max(tabulate(match(pos_vec, unique(pos_vec))))  # Most frequent position
  }
  
  # Process each peak separately
  peak_results <- lapply(1:nrow(peaks), function(peak_idx) {
    peak_pos <- tryCatch({
      if (use_KDE) density_data$x[peaks[peak_idx, 2]] else unique(pos_vec)[peaks[peak_idx, 2]]
    }, error = function(e) {
      message(paste("Error extracting peak position for cluster", cl, "in sample", sample, "-", e$message))
      return(NA)
    })
    
    if (is.na(peak_pos)) return(NULL)
    
    # Define local window around the peak (e.g., ±30 nt)
    window_start <- max(cluster_data$start, round(peak_pos - 30))
    window_end <- min(cluster_data$end, round(peak_pos + 30))
    
    # Subset data within this window
    sub_region_data <- region_data[pos >= window_start & pos <= window_end]
    sub_pos_vec <- tryCatch({
      unlist(lapply(1:nrow(sub_region_data), function(i) rep(sub_region_data$pos[i], times = as.integer(round(sub_region_data$count[i])))))
    }, error = function(e) {
      message(paste("Error in expanding sub-region data for cluster", cl, "in sample", sample, "-", e$message))
      return(NULL)
    })
    
    if (is.null(sub_pos_vec) || length(sub_pos_vec) < 5) return(NULL)  # Skip small peak regions
    
    dist_vec <- sub_pos_vec - peak_pos  # Normalize distances around peak
    
    # Initialize fit variables
    gamma_fit <- NULL
    gamma_TSS <- NA
    sn_fit <- NULL
    sn_TSS <- NA
    
    # Shift distribution for Gamma fitting
    shift_val <- abs(min(dist_vec)) + 1
    dist_vec_shifted <- dist_vec + shift_val
    
    # Fit gamma distribution
    if (length(dist_vec_shifted) > 5) {
      tryCatch({
        gamma_fit <- fitdist(dist_vec_shifted, "gamma")
        gamma_mode <- (gamma_fit$estimate["shape"] - 1) / gamma_fit$estimate["rate"]
        gamma_TSS <- round(peak_pos + (gamma_mode - shift_val))
      }, error = function(e) {
        message(paste("Error in Gamma fitting for cluster", cl, "in sample", sample, "-", e$message))
        gamma_fit <- NULL
        gamma_TSS <- NA
      })
    }
    
    # Fit skew-normal distribution
    if (length(dist_vec) > 5) {
      tryCatch({
        sn_fit <- selm(dist_vec ~ 1, family = "SN")
        sn_params <- coef(sn_fit)
        sn_mode <- sn_params[1]
        sn_TSS <- round(peak_pos + sn_mode)
      }, error = function(e) {
        message(paste("Error in Skew-Normal fitting for cluster", cl, "in sample", sample, "-", e$message))
        sn_fit <- NULL
        sn_TSS <- NA
      })
    }
    
    # Return results as a data.table
    return(data.table(
      sample = sample,
      cluster_id = cl,
      peak_id = peak_idx,
      chr = cluster_data$seqnames,
      strand = cluster_data$strand,
      cluster_start = window_start,
      cluster_end = window_end,
      peak_pos = peak_pos,
      gamma_TSS = gamma_TSS,
      gamma_shape = if (!is.null(gamma_fit)) gamma_fit$estimate["shape"] else NA,
      gamma_rate = if (!is.null(gamma_fit)) gamma_fit$estimate["rate"] else NA,
      sn_TSS = sn_TSS,
      sn_location = if (!is.null(sn_fit)) sn_params[1] else NA,
      sn_scale = if (!is.null(sn_fit)) sn_params[2] else NA,
      sn_skew = if (!is.null(sn_fit)) sn_params[3] else NA
    ))
  })
  
  return(rbindlist(peak_results, fill = TRUE))
}

# Prepare cluster information for parallel processing
cluster_info_list <- filtered_clusters[, .(group, cluster)]

# Run parallel processing
message("Starting parallel processing...")
result_list <- bplapply(seq_len(nrow(cluster_info_list)), function(i) {
  process_cluster(cluster_info_list$group[i], cluster_info_list$cluster[i])
}, BPPARAM = bp_param)
message("Parallel processing completed.")

# Convert list of data.tables to a single data.table safely
TSS_results <- rbindlist(result_list, fill = TRUE)

# Save refined TSS annotations
fwrite(TSS_results, "TSS_refined_annotations.tsv", sep = "\t")

# Print results
message("TSS annotation refinement completed.")
print(TSS_results)






stop()

TSS_results <- fread( "TSS_refined_annotations.tsv", sep = "\t")
TSS_results[,cluster_width := cluster_end - cluster_start + 1]

TSS_results[, .N, by = .(peak_id)]



TSS_results_comb <- merge(TSS_results, tagClusters, 
                          by.y=c('group', 'chr', 'strand', 'start', 'end', 'cluster'),
                          by.x=c('sample',  'chr', 'strand', 'cluster_start', 'cluster_end', 'cluster_id'),
                          all=T)
                          
ggplot(TSS_results_comb, aes(x = peak_pos  , y = gamma_TSS )) +
  geom_point( aes(fill = sample), shape=21) +
  labs(title = "Histogram of Skew-Normal Estimated TSS Positions", x = "Skew-Normal TSS Position", y = "Count")
 

# Print results
print(TSS_results)

# Visualization: Histogram of Gamma and Skew-Normal Estimated TSS
ggplot(TSS_results, aes(x = gamma_TSS)) +
  geom_histogram(binwidth = 5, fill = "blue", alpha = 0.5) +
  labs(title = "Histogram of Gamma Estimated TSS Positions", x = "Gamma TSS Position", y = "Count")

ggplot(TSS_results, aes(x = sn_TSS)) +
  geom_histogram(binwidth = 5, fill = "red", alpha = 0.5) +
  labs(title = "Histogram of Skew-Normal Estimated TSS Positions", x = "Skew-Normal TSS Position", y = "Count")

ggplot(TSS_results, aes(x = gamma_shape, y = gamma_rate)) +
  geom_point( aes(fill = sample), shape=21) +
  labs(title = "Histogram of Skew-Normal Estimated TSS Positions", x = "Skew-Normal TSS Position", y = "Count")


# Interpretation Table
interpretation_table <- data.table(
  Parameter = c("gamma_TSS", "gamma_shape", "gamma_rate", "sn_TSS", "sn_location", "sn_scale", "sn_skew"),
  Meaning = c(
    "Refined TSS position from Gamma distribution",
    "Shape parameter of Gamma distribution (higher = sharper peak)",
    "Rate parameter of Gamma distribution (higher = narrower peak)",
    "Refined TSS position from Skew-Normal distribution",
    "Location parameter (approximate mode) of Skew-Normal distribution",
    "Scale parameter (spread) of Skew-Normal distribution",
    "Skewness of Skew-Normal distribution (negative = left-skew, positive = right-skew)"
  )
)


TSS_results[,.N,gamma_shape]

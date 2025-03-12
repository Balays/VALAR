library(data.table)
library(fitdistrplus)
library(sn)  # For skew-normal fitting

# Convert TSSr processed matrix to data.table
TSS.dt <- as.data.table(myTSSr@TSSprocessedMatrix)
TSS.dt <- TSS.raw.dt
tagClusters <- as.data.table(tagClusters)

# List to store refined TSS annotations
result_list <- list()

# Iterate over each sample separately
for (sample in colnames(TSS.dt)[-(1:3)]) {  # Skip first 3 columns (chr, pos, strand)
  
  # Subset data for this sample
  sample_data <- TSS.dt[, .(chr, pos, strand, count = get(sample))]
  sample_data <- sample_data[count > 0]  # Remove zero-count positions
  
  # Iterate over each cluster in this sample
  for (cl in tagClusters[group == sample]$cluster) {
    
    # Extract cluster region
    cluster_data <- tagClusters[group == sample & cluster == cl]
    
    # Extract positions within cluster boundaries
    region_data <- sample_data[chr == cluster_data$chr & 
                                 strand == cluster_data$strand & 
                                 pos >= cluster_data$start & 
                                 pos <= cluster_data$end]
    
    if (nrow(region_data) == 0) next  # Skip empty clusters
    
    # Expand positions based on read counts
    pos_vec <- unlist(lapply(1:nrow(region_data), function(i) rep(region_data$pos[i], times = as.integer(round(region_data$count[i])))))
    
    if (length(pos_vec) < 5) next  # Skip clusters with insufficient data
    
    # Reference point: use dominant_tss as 0, measure distances
    peak_pos <- cluster_data$dominant_tss
    dist_vec <- pos_vec - peak_pos  # Negative = upstream, Positive = downstream
    
    # Initialize fit variables
    gamma_fit <- NULL
    gamma_TSS <- NA
    sn_fit <- NULL
    sn_TSS <- NA
    
    # Fit gamma distribution (for right-skewed TSS peaks)
    if (length(dist_vec[dist_vec >= 0]) > 5) {  # Ensure enough data for fitting
      tryCatch({
        gamma_fit <- fitdist(dist_vec[dist_vec >= 0], "gamma")
        gamma_mode <- (gamma_fit$estimate["shape"] - 1) / gamma_fit$estimate["rate"]
        gamma_TSS <- round(peak_pos + gamma_mode)  # Convert back to genomic position
      }, error = function(e) {
        gamma_fit <- NULL
        gamma_TSS <- NA
      })
    }
    
    # Fit skew-normal (for flexible peak fitting)
    if (length(dist_vec) > 5) {
      tryCatch({
        sn_fit <- selm(dist_vec ~ 1, family = "SN")  # Skew-normal fit
        sn_params <- coef(sn_fit)  # Extract parameters
        sn_mode <- sn_params[1]  # Location parameter (approximate mode)
        sn_TSS <- round(peak_pos + sn_mode)
      }, error = function(e) {
        sn_fit <- NULL
        sn_TSS <- NA
      })
    }
    
    # Store results as a named list
    result_list[[paste0(sample, "_", cl)]] <- data.table(
      sample = sample,
      cluster_id = cl,
      chr = cluster_data$chr,
      strand = cluster_data$strand,
      cluster_start = cluster_data$start,
      cluster_end = cluster_data$end,
      peak_pos = peak_pos,
      gamma_TSS = gamma_TSS,
      gamma_shape = if (!is.null(gamma_fit)) gamma_fit$estimate["shape"] else NA,
      gamma_rate = if (!is.null(gamma_fit)) gamma_fit$estimate["rate"] else NA,
      sn_TSS = sn_TSS,
      sn_location = if (!is.null(sn_fit)) sn_params[1] else NA,
      sn_scale = if (!is.null(sn_fit)) sn_params[2] else NA,
      sn_skew = if (!is.null(sn_fit)) sn_params[3] else NA
    )
  }
}

# Convert list of data.tables to a single data.table safely
TSS_results <- rbindlist(result_list, fill = TRUE)


TSS_results[,cluster_width := cluster_end - cluster_start + 1]

# Save refined TSS annotations
fwrite(TSS_results, "TSS_refined_annotations.tsv", sep = "\t")

TSS_results <- fread( "TSS_refined_annotations.tsv", sep = "\t")


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

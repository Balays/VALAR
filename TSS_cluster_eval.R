library(data.table)
library(fitdistrplus)  # For distribution fitting
library(sn)            # For skew-normal fitting
library(KernSmooth)    # For KDE smoothing


cluster_dt <- all_clusters[cell_line == 'PK-15' & hpi == '4h']
cluster_dt[, cluster := .GRP, by=.(strand, cluster)]
#cluster_dt[,cluster_id := paste0('cluster_', cluster)]

# Make sure data is sorted
setorder(cluster_dt, seqnames, strand, position)

# Store results
# Create a list to store results for each cluster
result_list <- list()

# Iterate over each unique cluster
for (cl in unique(cluster_dt$cluster)) {
  # Subset data for this cluster
  region_data <- cluster_dt[cluster == cl]
  
  # Expand positions based on read counts
  pos_vec <- rep(region_data$position, times = region_data$count)
  
  # Reference point: use cluster peak as 0, measure distances
  peak_pos <- unique(region_data$cluster_peak)
  dist_vec <- pos_vec - peak_pos  # Negative = upstream, Positive = downstream
  
  # Fit gamma distribution (for right-skewed TSS peaks)
  if (length(dist_vec[dist_vec >= 0]) > 5) {  # Ensure enough data for fitting
    gamma_fit <- fitdist(dist_vec[dist_vec >= 0], "gamma")
    gamma_mode <- (gamma_fit$estimate["shape"] - 1) / gamma_fit$estimate["rate"]
    gamma_TSS <- round(peak_pos + gamma_mode)  # Convert back to genomic position
  } else {
    gamma_fit <- NULL
    gamma_TSS <- NA
  }
  
  # Fit skew-normal (if we expect both upstream & downstream skew)
  if (length(dist_vec) > 5) {
    sn_fit <- selm(dist_vec ~ 1, family = "SN")  # Skew-normal fit
    sn_params <- coef(sn_fit)  # Extract parameters
    sn_mode <- sn_params[1]  # Location parameter (approximate mode)
    sn_TSS <- round(peak_pos + sn_mode)
  } else {
    sn_fit <- NULL
    sn_TSS <- NA
  }
  
  # Store results in a list
  result_list[[as.character(cl)]] <- list(
    cluster_id = cl,
    peak_pos = peak_pos,
    gamma_TSS = gamma_TSS,
    gamma_fit = gamma_fit,
    sn_TSS = sn_TSS,
    sn_fit = sn_fit
  )
}

# Convert list to data.table
TSS_results <- rbindlist(lapply(result_list, as.data.table), fill = TRUE)

# Display final refined TSS annotations
print(TSS_results)


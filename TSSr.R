library(TSSr)

TssData <- dcast.data.table(prime5.counts[correct_tss == T], 
                            seqnames + pos + strand ~ sample, 
                            value.var = 'count')
colnames(TssData)[1] <- 'chr'

TssData <- TSSr::createTSSobj(TssData)
# TssData <- clusterTSS(TssData, ...)


library(fitdistrplus)
result_list <- list()

for (cl in unique(clusters$cluster_id)) {
  region_data <- dt[clusters$cluster_id == cl]
  # Expand positions by coverage
  # This step might be big if coverage is large, but let's show the idea:
  pos_vec <- rep(region_data$position, times=region_data$count)
  
  # Optionally shift to a reference, e.g. subtract min(pos_vec) to make it 0-based if you want gamma
  # Or center around the cluster peak or median:
  peak_pos <- median(pos_vec)   # or region_data$peak
  dist_vec <- pos_vec - peak_pos
  
  # If you want to fit a skew-normal (sn package) directly:
  # library(sn)
  # sn_fit <- selm(dist_vec ~ 1, family="SN")
  # param_est <- coef(sn_fit)
  
  # If you want gamma with only non-negative data, ensure dist_vec >= 0
  dist_vec_gamma <- dist_vec[dist_vec >= 0]
  if (length(dist_vec_gamma) > 5) {
    gamma_fit <- fitdist(dist_vec_gamma, "gamma")
  } else {
    gamma_fit <- NULL
  }
  
  # Evaluate goodness-of-fit or store results
  result_list[[cl]] <- list(
    cluster_id = cl,
    peak_pos = peak_pos,
    gamma_fit = gamma_fit
  )
}

# Then parse the 'gamma_fit' object to get shape, rate, etc.


library(data.table)
library(fitdistrplus)
library(sn)  # For skew-normal fitting
library(ggplot2)
library(pracma)  # For peak detection
library(BiocParallel)  # For parallel processing
library(parallel)  # Load parallel for detectCores

# Convert TSSr processed matrix to data.table
TSS.dt <- as.data.table(myTSSr@TSSprocessedMatrix)
filtered_clusters <- as.data.table(filtered_clusters)  # Use filtered_clusters instead of tagClusters

# Define parallel processing parameters
bp_param <- MulticoreParam(workers = parallel::detectCores() - 1)  # Use all but one core


# Function to fit distributions and visualize fit for a single cluster
fit_and_plot_cluster <- function(cl,
                                 sample  = NULL, 
                                 use_KDE = F # Option to enable or disable KDE-based peak detection. Set to FALSE to disable
                                 ) {
  message(paste("Fitting and plotting cluster", cl, "in sample", sample))
  
  cluster_data <- filtered_clusters[cluster == cl]
  if(!is.null(sample)) { 
    cluster_data <- cluster_data[group == sample] 
  }
  
  region_data <- TSS.dt[chr      == cluster_data$seqnames & 
                          strand == cluster_data$strand & 
                          pos >= cluster_data$cluster_start & 
                          pos <= cluster_data$cluster_end, .(chr, pos, strand, count = get(sample))]
  
  if (nrow(region_data) <= 1) {
    message(paste("Skipping empty or small cluster", cl, "in sample", sample))
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
  
  # Select the first peak for fitting
  peak_pos <- if (use_KDE) density_data$x[peaks[1, 2]] else unique(pos_vec)[peaks[1, 2]]
  
  # Normalize distances around peak
  dist_vec <- pos_vec - peak_pos
  
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
      gamma_fit  <- fitdist(dist_vec_shifted, "gamma")
      gamma_mode <- (gamma_fit$estimate["shape"] - 1) / gamma_fit$estimate["rate"]
      gamma_TSS  <- round(peak_pos + (gamma_mode - shift_val))
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
  
  # Create a plot
  p <- ggplot(region_data, aes(x = pos, y = count)) +
    geom_bar(stat = "identity", fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = peak_pos, color = "black", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = gamma_TSS, color = "red", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = sn_TSS, color = "green", linetype = "dotted", linewidth = 1) +
    labs(title = paste("TSS Fit for Cluster", cl, "in Sample", sample),
         x = "Genomic Position", y = "Read Count",
         subtitle = paste("Gamma TSS:", gamma_TSS, "| Skew-Normal TSS:", sn_TSS)) +
    theme_minimal()
  
  print(p)
  
}

# Example usage
fit_and_plot_cluster(cl=2229, sample="PK_15_6h" , F) + scale_x_continuous(n.breaks = 10, limits=(c(126900, 127000)))

###############################################################################
# Load libraries
###############################################################################
library(data.table)
library(DECIPHER)
library(Biostrings)
library(BiocParallel)
library(future.apply)
library(pwalign)



# Adapter sequence
adapter_seq <- "GCTGATATTGCTGGG"
min_adapter_length <- 10  # User-defined threshold for the length of the aligned adapter

if(!exists('check_in_match') ) { check_in_match <- 10 }

###############################################################################
# Define Adapter Alignment Functions
###############################################################################

# 1) DECIPHER-based alignment (AlignSeqs)
decipher_adapter_alignment <- function(seq_to_align, strand, adapter_seq, min_adapter_length = 10) { ### alignment parameters !!!
  
  # Reverse complement the adapter if strand == '-'
  adapter_to_align <- if (strand == "-") {
    as.character(reverseComplement(DNAString(adapter_seq)))
  } else {
    adapter_seq
  }
  
  # Perform alignment with DECIPHER::AlignSeqs
  
  ### put the alignment parameters into the main functino's arguments !!!
  aligned <- AlignSeqs(verbose = F,
    c(DNAStringSet(seq_to_align), DNAStringSet(adapter_to_align)),
   # gapPower = -1
    gapOpening = -6, # -3, # 
    gapExtension = -6, # -3, # 
    misMatch = -2, # -3, # 
    perfectMatch = 3, # 2 
  )
  
  aln_score <- ScoreAlignment(aligned,
                              gapOpening = -6, # -3, # 
                              gapExtension = -6, # -3, # 
                              misMatch = -2, # -3, # 
                              perfectMatch = 3, # 2 
                 )
  
  aligned_adapter    <- as.character(aligned[[2]])
  non_gap_characters <- sum(strsplit(aligned_adapter, "")[[1]] %in% c("A", "T", "G", "C"))
  
  # Default results
  adapter_found <- "F"
  distance <- NA_integer_
  
  # Check if enough nucleotides aligned
  if (non_gap_characters >= min_adapter_length) {
    adapter_found <- "T"
    # Distance logic
    if (strand == "+") {
      last_nucleotide_pos <- max(which(strsplit(aligned_adapter, "")[[1]] %in% c("A", "T", "G", "C")))
      distance <- nchar(aligned_adapter) - last_nucleotide_pos
    } else {
      first_nucleotide_pos <- min(which(strsplit(aligned_adapter, "")[[1]] %in% c("A", "T", "G", "C")))
      distance <- first_nucleotide_pos - 1
    }
  }
  
  list(found = adapter_found, dist = distance, aln_score = aln_score)
  
}

# 2) Biostrings-based alignment (pairwiseAlignment)
biostrings_adapter_alignment <- function(seq_to_align, strand, adapter_seq, min_adapter_length = 10) {
  
  # Reverse complement the adapter if strand == '-'
  adapter_to_align <- if (strand == "-") {
    as.character(reverseComplement(DNAString(adapter_seq)))
  } else {
    adapter_seq
  }
  
  # Build a substitution matrix (example scoring: match=2, mismatch=-3)
  subs_mat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -3)
  
  # Local alignment with Biostrings::pairwiseAlignment
  palign <- pairwiseAlignment(
    pattern = DNAString(adapter_to_align),
    subject = DNAString(seq_to_align),
    type = "local",
    substitutionMatrix = subs_mat,
    gapOpening = -6,
    gapExtension = -6
  )
  
  aln_score <- palign@score
  
  #aligned_pattern <- as.character(alignedPattern(palign))
  
  n_nucl  <- palign@subject@range@width
  #n_nucl <- sum(strsplit(aligned_pattern, "")[[1]] %in% c("A", "T", "G", "C"))
  
  # Default results
  adapter_found <- "F"
  distance <- NA_integer_
  
  #if(palign@subject@range@width != palign@pattern@range@width) {
    
    # Check if enough nucleotides aligned
    if (n_nucl >= min_adapter_length) {
      adapter_found <- "T"
      
      # Distance logic
      if (strand == "+") {
        # position of last nucleotide in aligned pattern
        last_nucleotide_pos <- (n_nucl + palign@subject@range@start - 1)
        distance            <- nchar(seq_to_align) - last_nucleotide_pos
      } else {
        # position of first nucleotide
        first_nucleotide_pos <- palign@subject@range@start
        distance             <- first_nucleotide_pos - 1
      }
      
    }
    
  #}
  
  list(found = adapter_found, dist = distance, aln_score = aln_score)
}

###############################################################################
# Parallel Wrapper Functions (BiocParallel vs. future.apply)
###############################################################################

# -- DECIPHER + BiocParallel
run_DECIPHER_BiocParallel <- function(dt, adapter_seq = "GCTGATATTGCTGGG", min_adapter_length = 10, workers = 4, check_in_match=10) {
  register(MulticoreParam(workers = workers)) 
  res_list <- bplapply(seq_len(nrow(dt)), function(i) {
    row_data <- dt[i]
    decipher_adapter_alignment(row_data$sc_seq, row_data$strand, adapter_seq, min_adapter_length)
  })
  # Combine results
  dt[, adapter_found    := sapply(res_list, `[[`, "found")]
  dt[, adapter_distance := sapply(res_list, `[[`, "dist")]
  
  dt[, adapter_found    := factor(adapter_found, levels = c("F", "T"))]
  dt[, adapter_distance := adapter_distance - check_in_match]
  
  dt[, adapt_aln_score  := sapply(res_list, `[[`, "aln_score")]
  dt
}

# -- DECIPHER + future
run_DECIPHER_future <- function(dt, adapter_seq = "GCTGATATTGCTGGG", min_adapter_length = 10, workers = 4, check_in_match=10) {
  plan(multicore, workers = workers)  # or multisession on Windows
  res_list <- future_lapply(seq_len(nrow(dt)), function(i) {
    row_data <- dt[i]
    decipher_adapter_alignment(row_data$sc_seq, row_data$strand, adapter_seq, min_adapter_length)
  })
  # Combine results
  dt[, adapter_found    := sapply(res_list, `[[`, "found")]
  dt[, adapter_distance := sapply(res_list, `[[`, "dist")]
  
  dt[, adapter_found := factor(adapter_found, levels = c("F", "T"))]
  dt[, adapter_distance := adapter_distance - check_in_match]
  
  dt[, adapt_aln_score  := sapply(res_list, `[[`, "aln_score")]
  dt
}

# -- Biostrings + BiocParallel
run_Biostrings_BiocParallel <- function(dt, adapter_seq = "GCTGATATTGCTGGG", min_adapter_length = 10, workers = 4, check_in_match=10) {
  register(MulticoreParam(workers = workers)) 
  res_list <- bplapply(seq_len(nrow(dt)), function(i) {
    row_data <- dt[i]
    biostrings_adapter_alignment(row_data$sc_seq, row_data$strand, adapter_seq, min_adapter_length)
  })
  # Combine results
  dt[, adapter_found    := sapply(res_list, `[[`, "found")]
  dt[, adapter_distance := sapply(res_list, `[[`, "dist")]
  
  dt[, adapter_found := factor(adapter_found, levels = c("F", "T"))]
  dt[, adapter_distance := adapter_distance - check_in_match]
  
  dt[, adapt_aln_score  := sapply(res_list, `[[`, "aln_score")]
  dt
}

# -- Biostrings + future
run_Biostrings_future <- function(dt, adapter_seq = "GCTGATATTGCTGGG", min_adapter_length = 10, workers = 4, check_in_match=10) {
  plan(multicore, workers = workers)  # or multisession on Windows
  res_list <- future_lapply(seq_len(nrow(dt)), function(i) {
    row_data <- dt[i]
    biostrings_adapter_alignment(row_data$sc_seq, row_data$strand, adapter_seq, min_adapter_length)
  })
  # Combine results
  dt[, adapter_found    := sapply(res_list, `[[`, "found")]
  dt[, adapter_distance := sapply(res_list, `[[`, "dist")]
  
  dt[, adapter_found := factor(adapter_found, levels = c("F", "T"))]
  dt[, adapter_distance := adapter_distance - check_in_match]
  
  dt[, adapt_aln_score  := sapply(res_list, `[[`, "aln_score")]
  dt
}

###############################################################################
# Sample data (assume 'bam.correct' is your full data.table)
###############################################################################
# We'll take a random 5000-row subset for testing
bam.sub <- bam.sc.correct[sample(.N, 500)]

###############################################################################
# Run Timings
###############################################################################
# We'll make copies of the 5000-row subset so each method operates on a "fresh" dt.
# (Otherwise, columns get added/overwritten in-place.)

bam.decipher.bioc  <- copy(bam.sub)
bam.decipher.fut   <- copy(bam.sub)
bam.bio.bp         <- copy(bam.sub)
bam.bio.fut        <- copy(bam.sub)

adapter_seq <- "GCTGATATTGCTGGG"  # example adapter
workers     <- 14                  # adjust as desired
min_adapter_length <- 5


# 1) DECIPHER + BiocParallel
t1 <- system.time({
  bam.decipher.bioc <- run_DECIPHER_BiocParallel(bam.decipher.bioc, adapter_seq, min_adapter_length = min_adapter_length, workers = workers, check_in_match=check_in_match)
})

# 2) DECIPHER + future
t2 <- system.time({
  bam.decipher.fut <- run_DECIPHER_future(bam.decipher.fut, adapter_seq, min_adapter_length = min_adapter_length, workers = workers, check_in_match=check_in_match)
})

# 3) Biostrings + BiocParallel
t3 <- system.time({
  bam.bio.bp <- run_Biostrings_BiocParallel(bam.bio.bp, adapter_seq, min_adapter_length = min_adapter_length, workers = workers, check_in_match=check_in_match)
})

# 4) Biostrings + future
t4 <- system.time({
  bam.bio.fut <- run_Biostrings_future(bam.bio.fut, adapter_seq, min_adapter_length = min_adapter_length, workers = workers, check_in_match=check_in_match)
})

###############################################################################
# Print Timings
###############################################################################
cat("DECIPHER + BiocParallel:\n"); print(t1)
cat("DECIPHER + future:\n");       print(t2)
cat("Biostrings + BiocParallel:\n"); print(t3)
cat("Biostrings + future:\n");       print(t4)


###############################################################################
# The best result was DECIPHER with future. Run this on all the data.
###############################################################################

adapter_seq <- "GCTGATATTGCTGGG"  # example adapter
workers     <- 14                  # adjust as desired
min_adapter_length <- 5

# Increase limit to, e.g., 10 GB
options(future.globals.maxSize = 10 * 1024^3)  # 10 GiB


bam.sc.correct.adapt <- copy(bam.sc.correct[,])
bam.sc.correct.adapt <- run_Biostrings_future(bam.sc.correct.adapt, adapter_seq, min_adapter_length = min_adapter_length, workers = workers, check_in_match=check_in_match)



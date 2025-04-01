library(data.table)

# Provided foverlaps2 function
foverlaps2 <- function(DTx, DTy,
                        by = c('seqnames', 'strand', 'start', 'end'),
                        by.x = by, by.y = by,
                        type = 'any',
                        minoverlap = 20,
                        maxgap = 0) {
  setkeyv(DTx, by.x)
  setkeyv(DTy, by.y)
  
  DTxy <- foverlaps(DTx, DTy, minoverlap = 1, maxgap = maxgap, type = type)
  
  DTxy[, width_x := end - start + 1]
  DTxy[, width_y := i.end - i.start + 1]
  DTxy[, overlap_size := pmax(0, pmin(end, i.end) - pmax(start, i.start) + 1)]
  
  DTxy <- DTxy[overlap_size >= minoverlap]
  return(DTxy)
}

### 1. Prepare CDS UTR intervals with corrected boundaries

# Assume an extension of 1000 bp (adjust as needed)
extension <- extension
# Make a copy of CDS.dt for the UTR definitions
CDS_UTR <- copy(CDS.dt)
CDS_UTR[, `:=`(
  start = ifelse(strand == "+", CDS_end + 1, CDS_start - 1 - extension),
  end   = ifelse(strand == "+", CDS_end + 1 + extension, CDS_start - 1)
)]

### 2. Prepare the TES clusters (ConsClusters)

ConsClusters <- rbindlist(myTSSr@consensusClusters, idcol = 'group', use.names = T)
ConsClusters[, cluster_width := end - start +1]
# 
ConsClusters  <- merge(ConsClusters, unique(prime3.counts[,.(group, hpi, cell_line, Time)]), by='group', all.x=T)

ConsClusters <- as.data.table(ConsClusters)
# Rename chromosome column if needed
if("chr" %in% names(ConsClusters)) {
  setnames(ConsClusters, "chr", "seqnames")
}
# Save original cluster coordinates before merging

# Add a unique identifier
#ConsClusters[, cluster_id := .I]

### 3. First foverlaps2: assign TES clusters (in the UTR) to CDSs

mergedResult <- foverlaps2(
  DTx = ConsClusters,
  DTy = CDS_UTR,
  by = c("seqnames", "strand", "start", "end"),
  by.x = c("seqnames", "strand", "start", "end"),
  by.y = c("seqnames", "strand", "start", "end"),
  type = "any",
  minoverlap = 1,
  maxgap = 0
)

setnames(mergedResult, c('start' , 'end'), c('3-UTR_start', '3-UTR-end'))
setnames(mergedResult, c('i.start' , 'i.end'), c('start', 'end'))
setnames(mergedResult, c('width_x' , 'width_y', 'overlap_size'), c('3-UTR-width', 'width', 'TES_UTR_ov_width'))
mergedResult[, width := NULL]

# Merge back those TES clusters that did not overlap with a 3-UTR of a gene

mergedResult <- merge(mergedResult, 
                      ConsClusters, #[, .(seqnames, strand, group, , 
                      by =  colnames(ConsClusters), #c('seqnames', 'strand', 'gene', 'cluster', 'group' ),
                      all = T
                      )

# (Optional) For each gene and group, keep only the cluster with the highest 'tags' -->> canonic TES
merged_best <- mergedResult[, .SD[which.max(tags)], by = .(cluster, gene, group)]


# 


### 4. Second foverlaps2: use original cluster coordinates to assign the inCoding CDS

# Build a CDS intervals table using the coding region boundaries from CDS.dt
CDS_intervals <- CDS.dt[, .(seqnames, strand, gene, start = CDS_start, end = CDS_end)]

# Use the original TES cluster positions (orig_start/orig_end) as the query interval
coding_overlaps <- foverlaps2(
  DTx = merged_best, #[, .(cluster_id, seqnames, strand, start = orig_start, end = orig_end)],
  DTy = CDS_intervals,
  by = c("seqnames", "strand", "start", "end"),
  type = "any",
  minoverlap = 1,   # Minimal overlap requirement; adjust if needed
  maxgap = 0
)

setnames(coding_overlaps, c('start' , 'end'), c('CDS_Overlap_start', 'CDS_Overlap_end'))
setnames(coding_overlaps, c('i.start' , 'i.end'), c('start', 'end'))
setnames(coding_overlaps, c('gene', 'i.gene', 'overlap_size'), c('CDS_Overlap', 'gene', 'CDS_UTR_ov_width'))

coding_overlaps <- unique(coding_overlaps[,.(seqnames, strand, gene, cluster, start, end, CDS_Overlap, group, group, hpi, Time, cell_line, CDS_UTR_ov_width )])

# If multiple CDS overlap a cluster, choose the one with the maximum overlap
coding_overlaps <- coding_overlaps[, .SD[which.max(CDS_UTR_ov_width)], by = .(seqnames, strand, gene, cluster, start, end, group, hpi, Time, cell_line)]

coding_overlaps[,inCoding := CDS_Overlap]
coding_overlaps[,CDS_Overlap := NULL]

# Merge the inCoding annotation (the CDS 'gene' from CDS_intervals) back into merged_best
merged_best <- merge(merged_best, 
                     coding_overlaps,
                     by =  c('seqnames', 'strand',  'gene', 'cluster', 'start', 'end', 'group',  'hpi', 'Time', 'cell_line' ),
                     all = T)



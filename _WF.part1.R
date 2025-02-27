
##
#### WF Part 1. ####

require(BiocParallel)
safeBPParam <- function(nworkers) {
  if (.Platform$OS.type=="windows") {
    BiocParallel::SerialParam(bpparam())
  } else {
    BiocParallel::MulticoreParam(workers=nworkers, tasks = 20, progressbar = TRUE)
  }
}
if (.Platform$OS.type!="windows") {safeBPParam(config$nproc)}

#### Filter bamfiles to virus only ?
if (config$filter_bams) {
  source('filt_bamfiles_to_contig.R')
  
  bam_filt_outdir <- config$bam_filt_outdir
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  bamfiles <- list.files(config$bamdir, pattern = "\\.bam$", full.names = TRUE)
  contig_name <- config$genome
  
  filtered_bamfiles <- bplapply(
    bamfiles,
    filt_bamfiles_to_contig,
    contig_name = contig_name,
    output_dir = bam_filt_outdir
  )
  
  bamfiles <- filtered_bamfiles
}

## Import (filtered?) bamfiles
if (.Platform$OS.type=="windows") {
  bam.all.list <- purrr::map(
    bamfiles[],
    dt.from.bam,
    pattern=config$pattern,
    flag=config$flag,
    what = c("rname", "qname", "qwidth", "flag", "pos", "mapq", "cigar", "strand", "seq"),
    is.lortia = config$is.lortia,
    crop.na.cigar = TRUE,
    rm.gaps.in.aln = config$rm.gaps.in.aln,
    add.primes=TRUE
  )
} else {
  bam.all.list <- bplapply(
    bamfiles[],
    dt.from.bam,
    pattern=config$pattern,
    flag=config$flag,
    what = c("rname", "qname", "qwidth", "flag", "pos", "mapq", "cigar", "strand", "seq"),
    is.lortia = config$is.lortia,
    crop.na.cigar = TRUE,
    rm.gaps.in.aln = config$rm.gaps.in.aln,
    add.primes=TRUE
  )
}

bam.all <- rbindlist(bam.all.list)

## Overwrite alignment ID
bam.all[,aln_ID:=paste0(qname, '_', aln_nr)]



### Modify seqnames if multiple host genomes
if (config$rename_host_contigs) {
  bam.all[!grepl('PK', sample) & !is.na(seqnames), seqnames:=paste0(seqnames, '_Rnor')]
  bam.all[ grepl('PK', sample) & !is.na(seqnames), seqnames:=paste0(seqnames, '_Sscrofa')]
  bam.all[,.N,by=.(sample, seqnames, strand)]
}


### Fix viral contig
if (config$fix.viral.contigs) {
  bam.all[grepl(config$genome, seqnames),seqnames:=config$genome]
}


### calculate mapped and non-mapped read counts
reads <- rbind( unique(bam.all[!is.na(seqnames),.(sample, qname)])[,is.na:=F],
                unique(bam.all[ is.na(seqnames),.(sample, qname)])[,is.na:=T])
readcounts <- reads[, .(count = .N), by=.(sample, is.na)]
readcounts <- dcast(readcounts, sample ~ is.na, value.var = 'count')

readcounts[order(sample),]

print(readcounts)

#### ####
##


#### Write outputs ####
fwrite(readcounts, paste0(outdir, '/readcounts.tsv'), sep = '\t')

if (config$write.all) {
  fwrite(bam.all,    paste0(config$outdir, '/bam.all.tsv'), sep = '\t')
}

#### ####
##

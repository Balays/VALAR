
## -->> this chunk is the content of the '_WF.part0.R' file - now applied for monkeypox virus

##
##### Import libraries and functions
library(rtracklayer)
library(Rsamtools)
library(ggsci)
library(seqinr)
library(data.table)
library(stringi)
library(ggplot2)
library(tidygenomics)
library(tidyr)
library(dplyr)
library(fuzzyjoin)
library(future.apply)
library(purrr)
library(BiocParallel)

### Own functions
misc.dir <- config$misc_dir
minitax.dir <- config$minitax_dir
functions_dir <- 'functions'

if (.Platform$OS.type!="windows") {
  misc.dir <- paste0('/mnt/', gsub(':', '/', tolower(gsub('/.*', '', misc.dir))),
                     stri_replace_first_regex(misc.dir, '.*:/', ''))
  minitax.dir <- paste0('/mnt/', gsub(':', '/', tolower(gsub('/.*', '', minitax.dir))),
                        stri_replace_first_regex(minitax.dir, '.*:/', ''))
  
}


for(f in list.files(misc.dir, '*.R', full.names = T)) { try({source(f)}) }
for(f in list.files(minitax.dir, '*.R', full.names = T)) { try({source(f)}) }
for(f in list.files(functions_dir, '*.R', full.names = T)) { try({source(f)}) }

bam.flags <- fread(paste0(misc.dir, '/bam.flags.tsv'))
gff_compare_classes <- fread(paste0(misc.dir, '/gff_compare.txt'))
nproc <- config$nproc

#### ####

##### Metadata ####
metadata_file <- config$metadata_file
metadata_from_bamfiles <- config$metadata_from_bamfiles
bamfiles <- config$bamfiles
pattern  <- config$pattern

if(metadata_from_bamfiles) {
  bamfiles_df <- data.table(bamfile=bamfiles, sample = gsub(pattern, '', gsub('.*/', '', bamfiles)))
  bamfiles_df[, libtype := ifelse(grepl('dRNA', sample), 'dRNA', 'dcDNA')]
  bamfiles_df[, hpi     := stri_extract_first_regex(sample, 'C6_.*h|PC-12_.*h|PK-15_.*h|NRK_.*h')]
  bamfiles_df[, hpi     := gsub('C6_|PC-12_|PK-15_|NRK_', '', hpi)]
  bamfiles_df[, Time    := as.numeric(gsub('h', '', hpi))]
  
  bamfiles_df[, rep     := as.integer(stri_extract_last_regex(sample, '[1-3]')), by=.(sample)]
  
  bamfiles_df[, cell_line := fifelse(grepl('PC-12', sample), 'PC-12', 
                                     fifelse(grepl('PK-15', sample), 'PK-15', 
                                             fifelse(grepl('C6', sample), 'C6', 
                                                     fifelse(grepl('NRK', sample), 'NRK', 
                                                             'other'))))]
  
  bamfiles_df[, group   := paste(na.omit(c(cell_line, libtype, hpi)), collapse = '_'), by=.(sample)]
  
  bamfiles_df[, Infected  := 'inf']
  
  bamfiles_df$basecaller <- 'dorado'
  bamfiles_df$basecall_pass <- 'pass'
  
  bamfiles_df[, sample_name := paste(na.omit(c(group, rep)), collapse = '_'), by=.(sample)]
  
  metadata <- as.data.frame(bamfiles_df)
} else {
  metadata <- fread(metadata_file)
  setDF(metadata)
}

## check if metadata contains the bamfiles
stopifnot(length(c(
  setdiff(bamfiles_df$sample, gsub(pattern, '', gsub('.*/', '', bamfiles))),
  setdiff(gsub(pattern, '', gsub('.*/', '', bamfiles)), bamfiles_df$sample)
)) == 0)

##
metacols <- config$metacols
nmeta <- length(metacols)

metafilt <- metadata[,metacols]

#### ####

##### Genome and annotation ####
by <- c('seqnames', 'start', 'end')

genome <- config$genome
fasta.ref <- config$fasta_ref
fasta <- seqinr::read.fasta(fasta.ref)
l_genome <- length(fasta[[1]])
virus <- config$virus
viruses <- virus

### Annotation -->> This part should be adapted to each annotation !
 
if (config$create.ann.from.gff) {
 
  ## virus gff
  gff          <- as.data.frame(rtracklayer::import.gff(config$gff_file))
  CDS.df       <- gff[gff$type == 'CDS', c("seqnames", "start", "end", "strand", "type", "Name", "gene", "part")]
  CDS.df$ID    <- CDS.df$gene
  CDS.df[is.element(CDS.df$gene, dup(CDS.df$gene)), 'ID'] <- paste0(
    CDS.df[is.element(CDS.df$gene, dup(CDS.df$gene)), 'Name'], '_',
    CDS.df[is.element(CDS.df$gene, dup(CDS.df$gene)), 'part']
  )
  CDS.df$ID <- gsub('_NA', '', CDS.df$ID)
  luniq(CDS.df$ID) == nrow(CDS.df)
  
  feature.df <- CDS.df[,c("seqnames", "start", "end", "strand", "type", "gene", 'part', 'ID')]
  
  ## differentiate multicopy genes
  feature.dt <- data.table(feature.df)
  feature.dt[,copy_number := 1:n_distinct(strand), by=.(gene)]
  feature.dt[copy_number > 1, gene := paste0(gene, '_', copy_number)]
  feature.df <- as.data.frame(feature.dt)[,c("seqnames", "start", "end", "strand", "type", "gene", 'part', 'ID')]
  
  ## 
  feature.colname <- 'gene'
  by <- c('seqnames', 'start', 'end')
  
  ## Import gene clusters
  
  gene.clusters.all <- fread('refgenome/PRV.genes.TSS.TES.txt')
  gene.clusters.all <- gene.clusters.all[order(CDS.start),]
  ORFs  <- unique(gene.clusters.all$gene)
  genes <- unique(gene.clusters.all$gene)
  
  ## add non-coding genes
  feature.nc <- as.data.frame(gene.clusters.all[type=='NC-gene',.(seqnames, start=CDS.start, end=CDS.end, strand, type, gene, part=1, ID=gene)])
  feature.nc$gene_name <- feature.nc$gene
  feature.nc$gene_name[feature.nc$gene == 'NOIR']  <- 'NOIR1'
  feature.nc$gene_name[feature.nc$gene == 'NOIR2'] <- 'NOIR1'
  feature.nc$gene_name[feature.nc$gene == 'NOIR2-2'] <- 'NOIR2'
  feature.nc$gene_name[feature.nc$gene == 'NOIR-2'] <- 'NOIR2'
  
  feature.df$gene_name <- feature.df$gene
  feature.df$gene_name <- gsub('_2', '', feature.df$gene_name)
  feature.df <- rbind(feature.df, feature.nc)
 
  
  write.table(feature.df, config$feature.df.file, sep = '\t', row.names = F, quote = F)
  
} else {
  feature.df <- read.delim(config$feature.df.file)
  feature.colname <- 'ID'
}

#### ####

## -->> end of chunk to be included in the '_WF.part0.R' file 

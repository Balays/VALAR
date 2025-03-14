

## TSSr Results
#CAGE_counts <- rbind(
#  data.table(readxl::read_xlsx('ALL.samples.TSS.orient.xlsx', 1)),
#  data.table(readxl::read_xlsx('ALL.samples.TSS.orient.xlsx', 2))
#)


##
bamdir   <- "D:/data/PRV_3cell/CAGE/Gmail"
if (.Platform$OS.type!="windows") {
  bamdir <- paste0('/mnt/', gsub(':', '/', tolower(gsub('/.*', '', bamdir))),
                   stri_replace_first_regex(bamdir, '.*:/', ''))
}

pattern  <- '_mappedAligned.out.bam$'
bamfiles <- list.files(bamdir, pattern, recursive = T, full.names = T)


#pattern  <- ''
#bamfiles <- colnames(CAGE_counts)[c(4,6,7,9,10,11)]


### Config file creation
CAGE_config <- list(

  ###
  ## input bamfiles
  bamdir = bamdir,
  pattern = pattern,
  bamfiles = bamfiles,
  
  ## filter bamfiles for virus ?
  bam_filt_outdir = paste0(bamdir, "_filtered"),
  filter_bams = FALSE,
  
  ## results output - project folder
  outdir = outdir,
  
  ## Neccessary packages and functions
  misc_dir = 'C:/GitHub/Rlyeh/R',
  minitax_dir = 'C:/GitHub/minitax/R',
  
  ## Metadata
  metadata_from_bamfiles = F,
  metadata_file = 'metadata.cage.tsv',
  metacols = c('sample', 'hpi', 'Time', 'libtype', 'group', 'cell_line', 'sample_name'),
  
  ## Genome and annotation
  genome = 'LT934125.1',
  fasta_ref = 'refgenome/LT934125.1.fasta',
  virus = 'PRV-MDBIO',
  gff_file = 'refgenome/LT934125.1.gff3',
  create.ann.from.gff = T,
  feature.df.file = "refgenome/LT934125.1_feature.df.tsv",
  ## reference transcripts
  viral.ref.gff = "refgenome/LT934125.1_Torma_et_al_corrected.gff3",
  
  ## Other settings
  nproc = 48,
  
  write.all = TRUE,
  save.images = TRUE,
  
  rename_host_contigs = FALSE,
  fix.viral.contigs = TRUE,
  multiCellLines <- F,
  
  make.plots = TRUE,
  
  ## CAGE ?
  include.cage = T,
  
  ## Coverage analysis settings
  window_size = 50,
  window_step = 50,
  
  ## Alignment analysis settings
  is.lortia = F,
  rm.gaps.in.aln = TRUE,
  
  ## Coverage and alignment Analysis settings
  flag  = scanBamFlag(isSupplementaryAlignment=FALSE),
  param = ScanBamParam(what=scanBamWhat(), flag=scanBamFlag(isSupplementaryAlignment=FALSE)),
  
  ## GFF-compare settings
  TR.reads.gfffile = paste0(outdir, '/TR.reads.gff2'),
  
  thresh.eq.prime5 = 10,
  thresh.eq.prime3 = 10,
  thresh.eq.junc   = 2
  
)


saveRDS(CAGE_config, file = 'CAGE_config.rds')


if(CAGE_config$metadata_from_bamfiles) {
  
  meta_cage <- data.table(bamfile=bamfiles, sample = gsub(pattern, '', gsub('.*/', '', bamfiles)))
  meta_cage[,libtype := 'CAGE']
  meta_cage[,group := paste(cell_line, libtype, hpi, sep='_')]
  meta_cage[,sample_name := group]
  #meta_cage[, hpi  := stri_extract_first_regex(sample, '[0-9]*h')]
  #meta_cage[, Time := as.integer(gsub('h', '', hpi))]
  #meta_cage[, cell_line  := gsub('_.*h', '', sample)]
  #meta_cage[,libtype := 'CAGE']
  
  #metadata <- meta_cage
  fwrite(meta_cage, 'metadata.cage.tsv', sep='\t')
    
} else {
  
  meta_cage <- fread(CAGE_config$metadata_file)
  setDF(meta_cage)
  metadata <- meta_cage
  metafilt <- meta_cage
  
}

bamfiles <- grep(paste0(meta_cage$sample, collapse = '|'), bamfiles, value = T) 



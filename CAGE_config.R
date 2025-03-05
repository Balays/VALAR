

## TSSr Results
CAGE_counts <- rbind(
  data.table(readxl::read_xlsx('ALL.samples.TSS.orient.xlsx', 1)),
  data.table(readxl::read_xlsx('ALL.samples.TSS.orient.xlsx', 2))
)


##
bamdir   <- "D:/data/PRV_3cell/CAGE/Gmail"
pattern  <- 'Aligned.out.bam'
bamfiles <- grep('.bai',
                 list.files(bamdir, pattern, recursive = T, full.names = T),
                 invert = T, value = T)

pattern  <- ''
bamfiles <- colnames(CAGE_counts)[c(4,6,7,9,10,11)]


### Config file creation
CAGE_config <- list(
  bamdir = bamdir,
  pattern = pattern,
  outdir = outdir,
  misc_dir = 'C:/GitHub/Rlyeh/R',
  minitax_dir = 'C:/GitHub/minitax/R',
  genome = 'MPXV.MDBIO.1.0.0',
  fasta_ref = '../refgenome/MPOX_MD-NOVA.fasta',
  virus = 'MPXV',
  gff_file = '../refgenome/MPOX_MDBIO_USA.gff3',
  nproc = 48,
  metadata_from_bamfiles = F,
  metacols = c('sample', 'rep', 'hpi', 'Time', 'libtype', 'group', 'cell_line', 'sample_name'),
  write.all = TRUE,
  rename_host_contigs = FALSE,
  fix.viral.contigs = TRUE,
  make.plots = TRUE,
  save.images = TRUE,
  include.cage = FALSE,
  filter_bams = FALSE,
  is.lortia = TRUE,
  rm.gaps.in.aln = TRUE,
  window_size = 50,
  window_step = 50
)


saveRDS(CAGE_config, file = 'CAGE_config.rds')


metadata_from_bamfiles <- CAGE_config$metadata_from_bamfiles
if(metadata_from_bamfiles) {
  
  meta_cage <- data.table(sample=bamfiles)
  meta_cage[, hpi  := stri_extract_first_regex(sample, '[0-9]*h')]
  meta_cage[, Time := as.integer(gsub('h', '', hpi))]
  meta_cage[, cell_line  := gsub('_.*h', '', sample)]
  meta_cage[,libtype := 'CAGE']
  
  metadata <- meta_cage
  fwrite(meta_cage, 'CAGE_metadata.tsv', sep='\t')
    
} else {
  
  meta_cage <- fread('CAGE_metadata.tsv')
  metadata <- meta_cage
  
  
}




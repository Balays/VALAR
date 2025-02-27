project_config <- data.table(outdir = outdir)

bamdir   <- "D:/data/PRV_3cell/CAGE/Gmail"
pattern  <- 'Aligned.out.bam'
bamfiles <- grep('.bai',
                 list.files(bamdir, pattern, recursive = T, full.names = T),
                 invert = T, value = T)

CAGE_config <- CAGE_

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
  metadata_from_bamfiles = TRUE,
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

config

saveRDS(CAGE_config, file = 'CAGE_config.rds')


metadata_from_bamfiles <- config$metadata_from_bamfiles
if(metadata_from_bamfiles) {
  
  meta.cage <- data.frame(sample = gsub('.*\/', '', gsub(pattern, '', bamfiles)))
  
  meta.cage$group <- gsub('_202.*', '', meta.cage$sample)
  meta.cage$group <- gsub('_cage', '', meta.cage$group)
  
  meta.cage$cell_line <- 'RK-13'
  meta.cage$rep       <- c(1,1,1,2,2,2,3,3,3)

  metadata <- fwrite('CAGE_metadata.tsv', sep='\t')
    
} else {
  metadata <- fread('CAGE_metadata.tsv')
}




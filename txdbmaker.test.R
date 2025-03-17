
library(txdbmaker)

## TESTING GFF3
gffFile <- system.file("extdata", "GFF3_files", "a.gff3", package="txdbmaker")
gff     <- setDT(as.data.frame(import.gff3(gffFile)))
gff     <- gff[,.(seqnames, strand, start, end, source, phase, score, Name, ID)]
export.gff3(as.data.frame(gff), 'refgenome/test.gff3')

txdb <- makeTxDbFromGFF('refgenome/test.gff3',
                        dataSource="partial gtf file for Tomatoes for testing",
                        organism="Solanum lycopersicum")
ref  <- setDT(as.data.frame(genes(txdb)))


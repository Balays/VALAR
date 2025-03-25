

## OR !! use ref_mRNAs to accept 3-primes?
TR.ref[,transcript_prime3 := fifelse(strand == '+', transcript_end, transcript_start)]
TR.ref[,transcript_prime5 := fifelse(strand == '-', transcript_end, transcript_start)]

valid.prime3 <- unique(TR.ref[,.(seqnames, strand, transcript_prime3)])


fwrite(valid.prime3, 'valid.prime3.ori.tsv', sep='\t')



#### Import 
viral.ref <- fread(paste0(config$viral.ref.gff,'.tsv'), skip = 3 )
viral.ref <- data.table(separate(viral.ref, 9, into = c('transcript_id', 'Name', 'ID', 'gene_id', 'ORF_id'), sep = ';'))
viral.ref[, Name := gsub('Name=', '', Name)]
viral.ref[, transcript_id := gsub('transcript_id=', '', transcript_id)]
viral.ref[, ID := gsub('ID=', '', ID)]
viral.ref[, gene_id := gsub('gene_id=', '', gene_id)]
viral.ref[, ORF_id := gsub('ORF_id=', '', gene_id)]
viral.ref[type != 'exon', type := 'transcript']

## One canonic per gene?
length(unique(viral.ref[is.canonic == 'canonic', gene])) == length(unique(viral.ref[, gene]))
## OK

## start is smaller than end?
viral.ref[start > end]
## OK

## modified TES or TSS
viral.ref[, TSS := fifelse(strand == '+', start, end)]
viral.ref[, TES := fifelse(strand == '+', end, start)]

viral.corrected <- viral.ref[TSS != `TSS eredeti` | TES != `TES eredeti`]

tr.miss <- unique(viral.corrected[,transcript_id])

## gene annotation
viral.ref[gene != gene_id]
## OK


## export into .gff3 format
viral.ref <- viral.ref[,.(seqnames, source, type, start, end, strand, phase, score, transcript_id, ID, gene_id, ORF_id, is.canonic)]

export.gff3(viral.ref, config$viral.ref.gff )

## import back

#### Import 
viral.ref <- data.table(as.data.frame(rtracklayer::import.gff(config$viral.ref.gff, version='3')))
viral.ref[type=='exon', exon_number := 1:.N, by = .(transcript_id)]
viral.mrna <- viral.ref[,.(seqnames, source, type, phase, score, strand, start, end, transcript_id)]




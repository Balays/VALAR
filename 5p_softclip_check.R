library(Biostrings)
library(data.table)
library(GenomicRanges)
library(stringi)
library(Rsamtools)
library(BiocManager)
library(DECIPHER)


############## Part 1 ##############
### Get softclip seq, calculate adapter distance ###

# read data
#bam.all <- fread('dorado_MD-NOVA_long5p/bam.all.tsv')


# create test data with one sample
# bam.24h.A <- bam.all[
#  sample == "MPXV_dcDNA_24h_A_MPXV.MDBIO.1.0.0_sorted", 
#  .(seqnames, qname, strand, start, end, cigar, seq, tag.l5, tag.r5)
#]


# Use only the first exon (in some cases there were introns -> check !!)
# (The CIGAR and the Seq is the same in every exon of an alignment)
bam.all.filtered <- bam.all[ ex_nr == 1, 
                             .(seqnames, qname, strand, start, end, cigar, seq, tag.l5, tag.r5, sample) ] # 

# Create subset of correct 5' reads  ??????
# CONSIDER STRANDS !!!
bam.all.correct <- bam.all.filtered[
 # (strand == '+' & grepl("correct", tag.l5) ) | 
 # (strand == '-' & grepl("correct", tag.r5) )
]

# Number of correct 5-prime adapters per strand, according to LoRTIA
bam.all.correct[,tag.l5.type := gsub('.*,', '', tag.l5)]
bam.all.correct[,tag.r5.type := gsub('.*,', '', tag.r5)]
bam.all.correct[,.N,by=.(strand, tag.r5.type, tag.l5.type)]

# get SoftClip length
bam.all.correct[, softclip_length := fifelse(
  strand == "+",
  as.integer(stri_extract_first_regex(cigar, '[0-9]+(?=S)')),
  as.integer(stri_extract_last_regex(cigar, '[0-9]+(?=S)'))
)]

# Filter for those that have a 5-prime SoftClip
bam.sc.correct <- bam.all.correct[!is.na(softclip_length), ]
                                   

# get sc seq (+check_in_match from mapped [10 or 15]) for these
check_in_match <- 10
bam.sc.correct <- bam.sc.correct[,
                                   sc_seq := fifelse(
                                     strand == "+" & softclip_length > 0,
                                     paste0(
                                       stri_split_boundaries(seq, type = 'character', simplify = TRUE)[1:min(softclip_length + check_in_match, nchar(seq))], 
                                       collapse = ""
                                     ),
                                     fifelse(
                                       strand == "-" & softclip_length > 0,
                                       paste0(
                                         stri_split_boundaries(seq, type = 'character', simplify = TRUE)[(nchar(seq) - min(softclip_length + check_in_match, nchar(seq)) + 1):nchar(seq)], 
                                         collapse = ""
                                       ),
                                       NA_character_
                                     )
                                   ), 
                                   by = .I]


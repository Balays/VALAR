#!/bin/bash

main_dir=/mnt/d/data/PRV_3cell
reference=/mnt/d/data/genomes/Sus_scrofa.Sscrofa11.1.and.LT934125.1.fasta
#contig=LT934125.1


# PK-15
fastq_dir=${main_dir}/rebasecall/fastq_pass_PK-15
bamoutdir=${main_dir}/rebasecall/mapped_v6_virus/bam/PK-15
#filt_out=${main_dir}/rebasecall/mapped_v6_virus/filtered_bam/PK-15
out_sam=${main_dir}/rebasecall/mapped_v6_whost/sam/PK-15
outdir=${main_dir}/rebasecall/mapped_v6_whost/LoRTIA/PK-15

bash /mnt/k/programs/minimap2_LoRTIA.sh $fastq_dir $reference $bamoutdir 16 20 .fastq splice

#bash /mnt/k/programs/filter_bams.sh $bamoutdir $contig $filt_out 24

bash /mnt/d/programs/samfrombam.sh $bamoutdir $out_sam 24

bash /mnt/d/programs/run_LoRTIA_v0.sh $out_sam $outdir $reference

cd ${main_dir}/rebasecall


# C6
fastq_dir=${main_dir}/rebasecall/fastq_pass_C6
bamoutdir=${main_dir}/rebasecall/mapped_v6_virus/bam/C6
filt_out=${main_dir}/rebasecall/mapped_v6_virus/filtered_bam/C6
filt_out_sam=${main_dir}/rebasecall/mapped_v6_virus/filtered_sam/C6
outdir=${main_dir}/rebasecall/mapped_v6_virus/LoRTIA_virus/C6

bash /mnt/k/programs/minimap2_LoRTIA.sh $fastq_dir $reference $bamoutdir 56 20 .fastq splice

bash /mnt/k/programs/filter_bams.sh $bamoutdir $contig $filt_out

bash /mnt/k/programs/samfrombam.sh $filt_out $filt_out_sam

bash /mnt/k/programs/run_LoRTIA_v0.sh $filt_out_sam $outdir $reference

cd ${main_dir}/rebasecall


# PC-12
fastq_dir=${main_dir}/rebasecall/fastq_pass_PC-12
bamoutdir=${main_dir}/rebasecall/mapped_v6_virus/bam/PC-12
filt_out=${main_dir}/rebasecall/mapped_v6_virus/filtered_bam/PC-12
filt_out_sam=${main_dir}/rebasecall/mapped_v6_virus/filtered_sam/PC-12
outdir=${main_dir}/rebasecall/mapped_v6_virus/LoRTIA_virus/PC-12

bash /mnt/k/programs/minimap2_LoRTIA.sh $fastq_dir $reference $bamoutdir 56 20 .fastq splice

bash /mnt/k/programs/filter_bams.sh $bamoutdir $contig $filt_out 24

bash /mnt/k/programs/samfrombam.sh $filt_out $filt_out_sam 24

bash /mnt/k/programs/run_LoRTIA_v0.sh $filt_out_sam $outdir $reference

cd ${main_dir}/rebasecall


# NRK
fastq_dir=${main_dir}/NRK/merged.pass.fastq.gz
bamoutdir=${main_dir}/rebasecall/mapped_v6_virus/bam/NRK
filt_out=${main_dir}/rebasecall/mapped_v6_virus/filtered_bam/NRK
filt_out_sam=${main_dir}/rebasecall/mapped_v6_virus/filtered_sam/NRK
outdir=${main_dir}/rebasecall/mapped_v6_virus/LoRTIA_virus/NRK

bash /mnt/d/programs/minimap2_LoRTIA.sh $fastq_dir $reference $bamoutdir 12 20 .fastq.gz splice

bash /mnt/d/programs/filter_bams.sh $bamoutdir $contig $filt_out 12

bash /mnt/d/programs/samfrombam.sh $filt_out $filt_out_sam 12

bash /mnt/d/programs/run_LoRTIA_v0.sh $filt_out_sam $outdir $reference

cd ${main_dir}/rebasecall

#!/bin/bash

main_dir=/mnt/c/data/PRV_4cell
reference=/mnt/c/data/PRV_4cell/refgenome/LT934125.1.fasta
contig=LT934125.1

# C6
fastq_dir=${main_dir}/fastq_pass_C6
bamoutdir=${main_dir}/mapped_v6_virus/bam/C6
filt_out=${main_dir}/mapped_v6_virus/filtered_bam/C6
filt_out_sam=${main_dir}/mapped_v6_virus/filtered_sam/C6
outdir=${main_dir}/mapped_v6_virus/LoRTIA_virus/C6

bash /mnt/c/GitHub/HeadBash/minimap2_dRNA_LoRTIA.sh $fastq_dir $reference $bamoutdir 16 20 .fastq.gz splice

bash /mnt/c/GitHub/HeadBash/filter_bams.sh $bamoutdir $contig $filt_out 16

bash /mnt/c/GitHub/HeadBash/samfrombam.sh $filt_out $filt_out_sam 16

#bash /mnt/c/GitHub/HeadBash/run_LoRTIA_v0.sh $filt_out_sam $outdir $reference

cd ${main_dir}


# PK-15
fastq_dir=${main_dir}/fastq_pass_PK-15
bamoutdir=${main_dir}/mapped_v6_virus/bam/PK-15
filt_out=${main_dir}/mapped_v6_virus/filtered_bam/PK-15
filt_out_sam=${main_dir}/mapped_v6_virus/filtered_sam/PK-15
outdir=${main_dir}/mapped_v6_virus/LoRTIA_virus/PK-15

bash /mnt/c/GitHub/HeadBash/minimap2_LoRTIA.sh $fastq_dir $reference $bamoutdir 16 20 .fastq.gz splice

bash /mnt/c/GitHub/HeadBash/filter_bams.sh $bamoutdir $contig $filt_out 16

bash /mnt/c/GitHub/HeadBash/samfrombam.sh $filt_out $filt_out_sam 16

#bash /mnt/c/GitHub/HeadBash/run_LoRTIA_v0.sh $filt_out_sam $outdir $reference

cd ${main_dir}


# PC-12
fastq_dir=${main_dir}/fastq_pass_PC-12
bamoutdir=${main_dir}/mapped_v6_virus/bam/PC-12
filt_out=${main_dir}/mapped_v6_virus/filtered_bam/PC-12
filt_out_sam=${main_dir}/mapped_v6_virus/filtered_sam/PC-12
outdir=${main_dir}/mapped_v6_virus/LoRTIA_virus/PC-12

bash /mnt/c/GitHub/HeadBash/minimap2_LoRTIA.sh $fastq_dir $reference $bamoutdir 16 20 .fastq.gz splice

bash /mnt/c/GitHub/HeadBash/filter_bams.sh $bamoutdir $contig $filt_out 16

bash /mnt/c/GitHub/HeadBash/samfrombam.sh $filt_out $filt_out_sam 16

filt_out_sam=${main_dir}/mapped_v6_virus/filtered_sam
outdir=${main_dir}/mapped_v6_virus/LoRTIA_virus

bash /mnt/c/GitHub/HeadBash/run_LoRTIA_v0.sh $filt_out_sam $outdir $reference

cd ${main_dir}


# NRK
fastq_dir=${main_dir}/NRK/merged.pass.fastq.gz
bamoutdir=${main_dir}/mapped_v6_virus/bam/NRK
filt_out=${main_dir}/mapped_v6_virus/filtered_bam/NRK
filt_out_sam=${main_dir}/mapped_v6_virus/filtered_sam/NRK
outdir=${main_dir}/mapped_v6_virus/LoRTIA_virus/NRK

bash /mnt/d/programs/minimap2_LoRTIA.sh $fastq_dir $reference $bamoutdir 12 20 .fastq.gz splice

bash /mnt/d/programs/filter_bams.sh $bamoutdir $contig $filt_out 12

bash /mnt/d/programs/samfrombam.sh $filt_out $filt_out_sam 12

bash /mnt/d/programs/run_LoRTIA_v0.sh $filt_out_sam $outdir $reference

cd ${main_dir}

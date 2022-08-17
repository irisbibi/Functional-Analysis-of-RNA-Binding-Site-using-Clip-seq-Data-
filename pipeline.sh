#!/bin/bash
#SBATCH -N 2
#SBATCH -p RM-shared
#SBATCH -t 00-30:00:00
#SBATCH --ntasks-per-node=128
#SBATCH -o flow.sh_%j.out
#SBATCH -e flow.sh_%j.err


#$1, $2: genome fasta file and gff file
#$3, $4: experiment read1, read2
#$5, $6: control read1, read2
#$7: barcodecollapses.py
#$8: hg38_acc_sizes.txt

#Do the first quality control using FastQC.
module load FastQC
mkdir step1_first_quality_control

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "First quality control starts at $tic\n\n";

fastqc $3 $4 $5 $6 -o step1_first_quality_control -f fastq

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "First quality control finishes at $toc\n\n";


#Trim the adapter.
tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Adapter trimming starts at $tic\n\n";

mkdir step2_adapter_trimming

module load cutadapt
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o step2_adapter_trimming/exptrimmed_round1_1.fastq -p step2_adapter_trimming/exptrimmed_round1_2.fastq $3 $4

cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o step2_adapter_trimming/exptrimmed_round2_1.fastq -p step2_adapter_trimming/exptrimmed_round2_2.fastq step2_adapter_trimming/exptrimmed_round1_1.fastq step2_adapter_trimming/exptrimmed_round1_2.fastq

cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o step2_adapter_trimming/contrimmed_round1_1.fastq -p step2_adapter_trimming/contrimmed_round1_2.fastq $5 $6

cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o step2_adapter_trimming/contrimmed_round2_1.fastq -p step2_adapter_trimming/contrimmed_round2_2.fastq step2_adapter_trimming/contrimmed_round1_1.fastq step2_adapter_trimming/contrimmed_round1_2.fastq

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Adapter trimming finishes at $toc\n\n";


#Generate genome indice.
tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Genome indice generation starts at $tic\n\n";

mkdir step3_genome_indice

module load STAR/2.7.6a
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir step3_genome_indice --genomeFastaFiles $1 --sjdbGTFtagExonParentTranscript $2

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Genome indice generation finishes at $toc\n\n";


#Use STAR to map experiment reads to the genome.
tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "STAR alignment for replicate1 starts at $tic\n\n";

mkdir -p step4_alignment/experiment

STAR  --outSAMtype BAM SortedByCoordinate --runThreadN 10 --genomeDir step3_genome_indice --readFilesIn step2_adapter_trimming/exptrimmed_round2_1.fastq step2_adapter_trimming/exptrimmed_round2_2.fastq --outFileNamePrefix step4_alignment/experiment/exp --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "STAR alignment for replicate1 finishes at $toc\n\n";


#Use STAR to map control reads to the genome.
tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "STAR alignment for control starts at $tic\n\n";

mkdir -p step4_alignment/control

STAR  --outSAMtype BAM SortedByCoordinate --runThreadN 10 --genomeDir step3_genome_indice --readFilesIn step2_adapter_trimming/contrimmed_round2_1.fastq step2_adapter_trimming/contrimmed_round2_2.fastq --outFileNamePrefix step4_alignment/control/con --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "STAR alignment for control finishes at $toc\n\n";


#Sort and index the bam file after alignment.
tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Sorting and indexing bam file starts at $tic\n\n";

module load samtools

samtools view -b -F 8 step4_alignment/experiment/expAligned.sortedByCoord.out.bam > step4_alignment/experiment/exp_mapped.bam

samtools sort step4_alignment/experiment/exp_mapped.bam -o step4_alignment/experiment/expsorted.bam

samtools sort -n -o step4_alignment/experiment/expnamesorted.bam step4_alignment/experiment/exp_mapped.bam

samtools index step4_alignment/experiment/expsorted.bam

samtools view -b -F 8 step4_alignment/control/conAligned.sortedByCoord.out.bam > step4_alignment/control/con_mapped.bam

samtools sort step4_alignment/control/con_mapped.bam -o step4_alignment/control/consorted.bam

samtools sort -n -o step4_alignment/control/connamesorted.bam step4_alignment/control/con_mapped.bam

samtools index step4_alignment/control/consorted.bam

mv step4_alignment/experiment/expsorted.bam.bai step4_alignment/experiment/expnamesorted.bam.bai

mv step4_alignment/control/consorted.bam.bai step4_alignment/control/connamesorted.bam.bai

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Sorting and indexing bam file finishes at $toc\n\n";


#Eliminate the duplications.
tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Deduplication starts at $tic\n\n";

mkdir step5_deduplication

module load python
pip install pysam
python $7 -o  step5_deduplication/expnamesorteddup.bam -m  step5_deduplication/expnamesorteddup.metrics -b step4_alignment/experiment/expnamesorted.bam

python $7 -o  step5_deduplication/connamesorteddup.bam -m  step5_deduplication/connamesorteddup.metrics -b step4_alignment/control/connamesorted.bam

module load samtools

samtools sort step5_deduplication/expnamesorteddup.bam -o step5_deduplication/expdupsorted.bam
samtools index step5_deduplication/expdupsorted.bam

samtools sort step5_deduplication/connamesorteddup.bam -o step5_deduplication/condupsorted.bam
samtools index step5_deduplication/condupsorted.bam

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Deduplication finishes at $toc\n\n";


#Do the second quality control using FastQC.
module load FastQC
mkdir step6_second_quality_control

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Second quality control starts at $tic\n\n";

fastqc step5_deduplication/expdupsorted.bam step5_deduplication/condupsorted.bam -o step6_second_quality_control -f bam

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Second quality control finishes at $toc\n\n";


#Read coverage analysis.
tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Read coverage analysis starts at $tic\n\n";

mkdir step7_coverage_analysis

module load anaconda3
conda create --name correlation
conda activate correlation

conda install -c bioconda deeptools

plotFingerprint --bamfiles step5_deduplication/expdupsorted.bam step5_deduplication/condupsorted.bam --extendReads 110  --binSize=1000 --plotFile step7_coverage_analysis/fingerprint.png --labels HepG2input HepGcontrol

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Read coverage analysis finishes at $toc\n\n";


#Correlation analysis.
tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Correlation analysis starts at $tic\n\n";

mkdir -p step8_correlation_analysis/multiBamSummary

multiBamSummary bins --bamfiles step5_deduplication/expdupsorted.bam step5_deduplication/condupsorted.bam -bs 1000 --labels HepG2input HepG2control -o step8_correlation_analysis/multiBamSummary/multiBamSummaryresults.npz

plotCorrelation -in step8_correlation_analysis/multiBamSummary/multiBamSummaryresults.npz \
 -c spearman -p heatmap -o step8_correlation_analysis/correlation.png --skipZeros -l HepG2input HepG2control

conda deactivate correlation
toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Correlation analysis finishes at $toc\n\n";


#Peak calling.
mkdir step9_peak_calling

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Peak calling starts at $tic\n\n";

#Install peakachu package for peak calling
module load anaconda3/2020.07

#add necessary channels for PEAKachu installation
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

cd PEAKachu
make readme_rst
make package

conda create -n peakachu peakachu python=3

conda activate peakachu 

cd ..

peakachu adaptive\
 -M 200\
 -m 0.0\
 -Q 0.05\
 -c step5_deduplication/condupsorted.bam\
 -t step5_deduplication/expdupsorted.bam\
 -o step9_peak_calling

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "Peak calling finishes at $toc\n\n";

module unload anaconda3/2020.07


#MEME step
tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "meme starts at $tic\n\n";
#Column 2 is starting position, 3 is ending position, column 8 is fold_change, column 4 is peak_strand.
mkdir step10_meme_chip

awk 'NR>1{if ($2 < $3) {print $1,$2,$3,"clip_peak_",NR-1,$8,$4;}else {print $1,$3,$2,"clip_peak_",NR-1,$8,$4;}}' step9_peak_calling/initial_peaks.csv > step10_meme_chip/ini_peak2.bed

#Make it tab-delimited.
cat step10_meme_chip/ini_peak2.bed | tr ' ' '\t' > step10_meme_chip/ini_peak2_tab.bed

#load bedtools
module load bedtools

#hg38 should be accession id labeled instead of chromosome number labeled
bedtools slop -i step10_meme_chip/ini_peak2_tab.bed -g $8 -b 20 >  step10_meme_chip/peak_add20.bed

#align peak_add 20 to reference genome (hg38)
bedtools getfasta -fi $1  -bed step10_meme_chip/peak_add20.bed > step10_meme_chip/primary_seq.fna

# de-duplicate IDs
awk '/^>/{f=!d[$1];d[$1]=1}f' step10_meme_chip/primary_seq.fna > step10_meme_chip/primary_seq_dedup.fa

#run MEME
module load MEME-suite
meme-chip -minw 3 -maxw 20 step10_meme_chip/primary_seq_dedup.fa

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "meme finishes at $toc\n\n";

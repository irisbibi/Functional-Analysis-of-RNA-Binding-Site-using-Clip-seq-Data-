#!/bin/bash
#SBATCH -N 2
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node=128
#SBATCH -o /ocean/projects/mcb180074p/li6/cutadapt/cutadapt.sh_%j.out
#SBATCH -e /ocean/projects/mcb180074p/li6/cutadapt/cutadapt.sh_%j.err

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "start at $tic\n\n";

module load cutadapt
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o /ocean/projects/mcb180074p/li6/cutadapt/trimmed_round1_1.fastq -p /ocean/projects/mcb180074p/li6/cutadapt/trimmed_round1_2.fastq /ocean/projects/mcb180074p/li6/encodedata/TAF15/replicate1/R11.fastq /ocean/projects/mcb180074p/li6/encodedata/TAF15/replicate1/R12.fastq

cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o /ocean/projects/mcb180074p/li6/cutadapt/trimmed_round2_1.fastq -p /ocean/projects/mcb180074p/li6/cutadapt/trimmed_round2_2.fastq /ocean/projects/mcb180074p/li6/cutadapt/trimmed_round1_1.fastq /ocean/projects/mcb180074p/li6/cutadapt/trimmed_round1_2.fastq


toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "finish at $toc\n\n";




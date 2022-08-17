#!/bin/bash
#SBATCH -N 2
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node=128
#SBATCH -o /ocean/projects/mcb180074p/li6/star/experiment/sort.sh_%j.out
#SBATCH -e /ocean/projects/mcb180074p/li6/star/experiment/sort.sh_%j.err

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "start at $tic\n\n";

module load samtools
samtools sort -n -o /ocean/projects/mcb180074p/li6/star/experiment/expAligned.sortedByCoord.outSo.bam /ocean/projects/mcb180074p/li6/star/experiment/expAligned.sortedByCoord.out.bam


toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "finish at $toc\n\n";




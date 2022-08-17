#!/bin/bash
#SBATCH -N 2
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node=128
#SBATCH -o /ocean/projects/mcb180074p/biyuz/correlation/correlation.sh_%j.out
#SBATCH -e /ocean/projects/mcb180074p/biyuz/correlation/correlation.sh_%j.err

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "start at $tic\n\n";
#module load anaconda3
#cd $PROJECT
#conda create --name correlation
#conda activate correlation
#conda install -c bioconda deeptools

multiBamSummary bins --bamfiles /ocean/projects/mcb180074p/biyuz/resultAligned.sortedByCoord.out.bam /ocean/projects/mcb180074p/biyuz/control.out.bam -bs 1000 --labels HepG2input HepG2control -o /ocean/projects/mcb180074p/biyuz/correlation/multiBamSummary/multiBamSummaryresults.npz


plotCorrelation -in /ocean/projects/mcb180074p/biyuz/correlation/multiBamSummary/multiBamSummaryresults.npz \
 -c spearman -p heatmap -o /ocean/projects/mcb180074p/biyuz/correlation/correlation.png --skipZeros -l HepG2input HepG2control

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "finish at $toc\n\n";
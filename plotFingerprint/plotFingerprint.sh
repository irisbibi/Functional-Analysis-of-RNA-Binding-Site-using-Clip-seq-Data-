#!/bin/bash
#SBATCH -N 2
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node=128
#SBATCH -o /ocean/projects/mcb180074p/biyuz/plotFingerprint/plotFingerprint.sh_%j.out
#SBATCH -e /ocean/projects/mcb180074p/biyuz/plotFingerprint/plotFingerprint.sh_%j.err

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "start at $tic\n\n";

plotFingerprint --bamfiles /ocean/projects/mcb180074p/biyuz/resultAligned.sortedByCoord.out.bam /ocean/projects/mcb180074p/biyuz/control.out.bam --extendReads 110  --binSize=1000 --plotFile /ocean/projects/mcb180074p/biyuz/plotFingerprint/fingerprint.png --labels HepG2input HepGcontrol

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "finish at $toc\n\n";
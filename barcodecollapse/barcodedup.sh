#!/bin/bash
#SBATCH -N 2
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node=128
#SBATCH -o /ocean/projects/mcb180074p/li6/barcodedup/barcodedup.sh_%j.out
#SBATCH -e /ocean/projects/mcb180074p/li6/barcodedup/barcodedup.sh_%j.err

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "start at $tic\n\n";

module load python
pip install pysam
python /jet/home/li6/03708/barcodedup/barcodecollapsepe.py -o /ocean/projects/mcb180074p/li6/barcodedup/expAligned.sortedByCoord.outSo.rmdup.bam -m /ocean/projects/mcb180074p/li6/barcodedup/expAligned.sortedByCoord.outSo.rmdup.metrics -b /ocean/projects/mcb180074p/li6/star/experiment/expAligned.sortedByCoord.outSo.bam

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "finish at $toc\n\n";




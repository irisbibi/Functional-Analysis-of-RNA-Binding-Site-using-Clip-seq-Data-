mkdir step11_downstream_analysis

git clone https://github.com/vkkodali/cthreepo.git
cd cthreepo
python setup.py install
cd ..
cthreepo -i $1 -if rs -it uc -f gff3 -m h38 -o step11_downstream_analysis/ucsc_genome.gff3
cthreepo -i step10_meme_chip/peak_add20.bed -if rs -it uc -f bed -m h38 -o step11_downstream_analysis/peak_ucsc.bed

R CMD BATCH $2



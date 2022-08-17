# Functional-Analysis-of-RNA-Binding-Site-using-Clip-seq-Data-

##
### Sequence Motif Discovery Final Report

May 6th, 2022

### Introduction

A critical step in unraveling genome processing is to study the _in vivo_ regulatory networks of RNA-binding proteins (RBPs). In recent years, there has been more awareness of the function of RBPs as the number of RBPs that are associated with human diseases are being discovered. The purpose of this project is to build a pipeline that utilizes enhanced Cross-linking immunoprecipitation coupled with high-throughput (eCLIP) data to discover possible motifs in the human genome, specifically, on the TAF15 target (TATA-Box Binding Protein Associated Factor 15).

CLIP-seq data is often purposed to resolve binding sites in a given target gene sequence, allowing the discovery of sequence motifs recognized by specific RBPs and lending insights into binding functionality. The eCLIP data, a set of RNA elements recognized by RBPs, are collected from the Encyclopedia of DNA Elements (ENCODE) project. For this project, we used the eCLIP data in HepG2 dataset on the TAF15 target.

TAF15 is a protein coding gene that is associated with numerous diseases such as Chondrosarcoma, Extraskeletal Myxoid and Myxoid Chondrosarcoma. Several pathways are related to TAF15, such as RNA polymerase II transcription initiation and promoter clearance and apoptotic pathways in synovial fibroblasts3.

The workflow of the motif discovery is inspired by Blue et al. paper and the Galaxy Project that showed a comprehensive data processing procedure for eCLIP data. Our workflow is composed of mainly four steps: quality control using FastQC, reads alignment using RNA STAR, peak calling with PEAKachu and motif discovery using MEME-ChIP. The whole pipeline was constructed based on _bash_ scripts and is tailored for running on the supercomputer Bridges-2 from PSC.

In the end, we have discovered a decent amount of motifs with the top motif very rich in thymine. An additional downstream analysis is also provided along with a modified _RNA Centric Annotation System_ R script for specific usage. This additional part generates functional annotation and GO term enrichment analysis for further biology research.

### Material and Methods

## Material

Experimental pair-end eCLIP data for TAF15: [experiment data](https://www.encodeproject.org/experiments/ENCSR841EQA/?fbclid=IwAR1ZP4e__lwZyVRS1XRbGIhc7DyxwrCVIf9V8RHTQMSYJnSMs6X0Tfdxn3A)

Mock input pair-end eCLIP data for HepG2: [mock input data](https://www.encodeproject.org/experiments/ENCSR716AKC/)

Reference genome for GRCh38: [NCBI genome](https://www.ncbi.nlm.nih.gov/genome/guide/human/?fbclid=IwAR1_m3hC04evFzwCali2u5FBKmAFKdNQNt0SZqt3r_dL7-hDybbcf2N7aWo)

## Methods

In this project, a pipeline for finding RNA-binding proteins' binding motifs was built. It was tested on the Pittsburgh SuperComputer (PSC). Experimental and control pair-end data files were firstly downloaded from the ENCODE website: two reads fastq files for TAF15 assay and two mock input reads files for HepG2 cell line. In addition, reference genome files were retrieved from the NCBI human genome resource website.

In the beginning, data processing steps were performed. As the first step of the pipeline, _ **FastQC** _ was used to examine the duplication level from PCR of the four raw read files. PCR duplication occurs commonly for eCLIP experiments so a quality control process is a crucial step before executing any further steps. The results from this step will be compared with those from the second quality control step, which is conducted after a deduplication step to check if the duplication level has been correctly reduced. Then, specific adapters which were added to the sequences during PCR amplification needed to be removed before being aligned to the reference genome. The specific adapter sequences were obtained from Yeo's github page and _ **Cutadapt** _ was used to remove them. Two rounds of adapter trimming were performed as double ligation at the 5' end occurs sometimes in the first read file.

After the adapters were trimmed off from the experimental fastq files and the control fastq files, _ **STAR** _ was used to produce reference genome index and align the reads to the reference hg38. Two bam files were generated here, where one was the alignment file for the experiment group and the other was for the mock input (control) group. Once the alignment was done, a deduplication step was conducted to eliminate duplications generated during PCR using _ **barcodecollapsepe.py** _ from Yeo's github. This step required a preliminary step. _ **Samtools** _ was used to pre-process the two output bam files to make them valid input to the deduplication tool. Firstly, because in some reads, one of the pair-end reads did not match the reference genome, a flag _-f 8_ in _ **Samtools** _ was used to filter out these unmapped reads. Then, the bam files were sorted by name using _ **Samtools** _ so that the two reads with the same name were placed next to each other. In addition, the original _ **barcodecollapsepe.py** _ was modified: "itertools.izip" method in the code was replaced by "zip" because "itertools.zip" was no longer used in python3.

A second quality control was applied at this point. _ **FastQC** _ was again used for examining the duplication level after deduplication. _ **plotFingerprint** _ and _ **plotCorrelation** _ were used to check if the enrichment of this eCLIP experiment was strong enough to separate the signal from background. After the second quality control, _ **Samtools** _ was used to sort and index the two bam files. Two index files with suffix bai were generated respectively for the experiment group and control group. These sorted and indexed bam files were now ready for peak calling, which was done using _ **PEAKachu** _.

Several peak files with different formats were generated after peak calling, among which a peak csv file was further treated for _ **meme-chip** _ motif discovery tool. It was a csv file with several columns, where the second and third columns were the starting and ending position of the peak. For example, the second column value of 157 and the third column value of 345 indicated a peak starting from 157 and ending at 345. However, sometimes a program might return a result with a starting point at 345 and an ending point at 157. Therefore, a check point was set using _ **awk** _ command lines in linux to ensure that the starting position was always smaller than the ending position. One adjusted file was generated in bed format after the _ **meme-chip** _ process. Moreover, the searching window of this bed file for _ **meme-chip** _ was manually enlarged using _ **bedtools** _. After all of these were done, _ **meme-chip** _ took the processed file as input and returned found motifs.

After the motif discovery step, an additional downstream analysis was added to the pipeline; sequence structure annotation and GO term enrichment analysis were applied on the peak bed file generated in the last step. As a preparatory step before the functional analysis, a python tool _ **cthreepo** _ was used to convert the NCBI-formatted chromosome name (NC\_) to UCSC-formatted name (chr\_) in both the peak bed file and our reference genome GFF file. This step was necessary because the _**RNA Centric Annotation System (RCAS)**_ R package required a UCSC-format chromosome name system in the input file. _ **RCAS** _ produced two figures, one of which showed the annotated structure of the RNA peak sequences (eg. protein coding region) and the other presented significant GO terms that indicated the potential biological functions of the sequences where our motifs came from. The downstream analysis could not be run on PSC because the installation of _ **cthreepo** _ required system authority.

For all the steps in the bash script, time consumption tracking commands were used to record the time usage for each part of the pipeline, which was important for measuring the pipeline efficiency.

## Tools

### FastQC

_FastQC_ package is a pre-installed tool in Bridges-2. This package is used for evaluating the quality of the raw sequencing data file by measuring the duplication level of reads with different lengths. PCR duplicates need to be erased before the following steps in the pipeline as they can give false positive results to the analysis.

For usage on Bridges-2, refer to: [FastQC | PSC](https://www.psc.edu/resources/software/fastqc)

For documentations, refer to: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### Cutadapt

_Cutadapt_ package is a pre-installed tool in Bridges-2, and is used for finding and removing adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

For usage on Bridges-2, refer to: [Cutadapt | PSC](https://www.psc.edu/resources/software/cutadapt/)

For documentations, refer to: [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)

### STAR Aligner

_STAR Aligner_ package is a pre-installed tool in Bridges-2, which is used for mapping reads from the spliced transcripts alignment to a reference with support for splice-junction and fusion read detection.

For usage on Bridges-2, refer to: [SAMtools | PSC](https://www.psc.edu/resources/software/samtools/)

For documentations, refer to: [SAMtools](http://samtools.sourceforge.net/)

### SAMtools

_SAMtools_ package is a pre-installed tool in Bridges-2, which provides various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.

For usage on Bridges-2, refer to: [STAR Aligner | PSC](https://www.psc.edu/resources/software/star-aligner/)

For documentations, refer to: [STAR Aligner](https://github.com/alexdobin/STAR)

### barcodecollapsepe.py

barcodecollapsepe.py is a python script from yeo's lab Github homepage used for deduplication of pair-end eclip data. For single-end eCLIP data, use UMI-tools for deduplication.

barcodecollapsepe.py can be retrieved from [yeo's lab Github](https://github.com/yeolab/eclip).

The script is under the bin directory.

### Anaconda

Anaconda environment is a pre-installed tool in Bridges-2. It is a data science platform which includes Python and R. Multiple versions of Anaconda are available on Bridges-2. For this pipeline, Anaconda3 is installed for deepTools and PEAKachu packages. To add the environment on Bridges-2, refer to the Adding to an environment section: [Anaconda | PSC](https://www.psc.edu/resources/software/anaconda). Then install the following packages.

### deepTools

deepTools is a suite of python tools particularly developed for the efficient analysis of high-throughput sequencing data, such as ChIP-seq, RNA-seq or MNase-seq. Installation is through the [Bioconda](https://bioconda.github.io/) channel in the Anaconda environment. We use _plotFingerprint_, _multiBamSummary_, and _plotCorrelation_ tools to perform the second quality control for the mapped reads.

For installation details, refer to: [Installation — deepTools](https://deeptools.readthedocs.io/en/develop/content/installation.html)

For more information, refer to: [deepTools](https://deeptools.readthedocs.io/en/develop/)

### PEAKachu

PEAKachu is a peak calling tool for CLIP-seq data. Installation is through the [Bioconda](https://bioconda.github.io/) channel in the Anaconda environment. It takes input in BAM format and identifies regions of statistically significant read enrichment. PEAKachu uses signal and control libraries to detect binding sites. It implements two peak calling approaches.

For installation details, refer to: [Peakachu](https://anaconda.org/bioconda/peakachu)

For usages, refer to: [PEAKachu](https://github.com/tbischler/PEAKachu/blob/master/bin/peakachu)

### MEME Suite

MEME Suite is a collection of pre-installed tools in Bridges-2, which is used for the discovery and analysis of sequence motifs. We use the MEME-ChIP tool to discover novel motifs. The Position-Specific Weight Matrix (PSWM) is the core of this algorithm, which is a probability matrix of bases at different positions of a sequence. A motif can be considered as a pattern among multiple sequences, so that it is possible to distinguish between motif sequences and background sequences based on their different position-specific weight. Normally, the parameters of this probabilistic model are estimated via estimation-maximization algorithm, where in the estimation step we randomly assign parameters and in the maximization step we update parameters to maximize the expected likelihood.

For usage on Bridges-2, refer to: [MEME-Suite | PSC](https://www.psc.edu/resources/software/meme-suite/)

For usages of MEME-ChIP, refer to: [MEME-ChIP](https://meme-suite.org/meme/doc/meme-chip.html?man_type=web)

### Cthreepo

Cthreepo is a python script that helps convert reads whose chromosome name in NCBI format to the chromosome name used by the UCSC system. It supports various input file formats, such as GFF, GTF, bed and so on.
 To install Cthreepo, refer to: [Cthreepo-Github](https://github.com/vkkodali/cthreepo)

### RNA Centric Annotation System

RNA Centric Annotation System (RCAS) is a R package and it conducts functional analysis on transcriptome-wide regions. It provides various analysis, including annotation summaries, GO term enrichment analysis, gene set enrichment analysis and motif discovery.
 To learn more details, refer to the original paper: [RCAS](https://academic.oup.com/nar/article/45/10/e91/3038237)
 For installation and usage instructions, refer to: [RCAS user guideline](https://www.bioconductor.org/packages/release/bioc/vignettes/RCAS/inst/doc/RCAS.vignette.html)

## Using the Pipeline

All the commands the pipeline uses are summarized into one command line as explained in the manual.

### Required Input Files

#### Part 1: Motif Discovery (bridges2 supported)

The reference data file can be downloaded by the users from the public database: [Human Genome Resources at NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/?fbclid=IwAR1_m3hC04evFzwCali2u5FBKmAFKdNQNt0SZqt3r_dL7-hDybbcf2N7aWo).

The experimental data can be generated either by the users or downloaded from the public ENCODE database. Link to the website that contains the files used in our pipeline: [Experiment summary for ENCSR841EQA](https://www.encodeproject.org/experiments/ENCSR841EQA/?fbclid=IwAR1ZP4e__lwZyVRS1XRbGIhc7DyxwrCVIf9V8RHTQMSYJnSMs6X0Tfdxn3A).

1. Reference genome fasta file: serves as reference genome file for STAR Aligner.
2. Reference genome annotation GFF or GTF file: genome annotation file, which can be in either GFF or GTF format.
3. Two pair-end eCLIP experiment fasta files: experimental data to build our own transcripts. These include experiment\_read1.fastq and experiment\_read2.fastq.

NOTICE: Due to different versions of the python you have installed on your computer, you may need to manually replace "itertools.izip" statement in barcodecollapsepe.py you download from yeo's lab github with "zip". In newer python3. The itertools.izip is no longer used and replaced by the built-in "zip."

#### Part 2: Additional downstream analysis (locally supported)

Downstream\_analysis.r R script for sequence annotation and GO analysis.

NOTICE: You should have R installed on your computer. Test it with "R —version" in your terminal.

### Usage

#### Part 1: MOTIF DISCOVERY

As a bash script, it can be simpy run with one command in the terminal:

_sbatch pipeline.sh_ _ **$1 $2 $3 $4 $5 $6 $7 $8** _

(on bridges)

$1: the absolute path to the reference genome fasta file

$2: the absolute path to the reference genome annotation GFF/GTF file

$3: the absolute path to the experiment pair-end eclip read\_1 fastq file

$4: the absolute path to the experiment pair-end eclip read\_2 fastq file

$5: the absolute path to the control pair-end eclip read\_1 fastq file

$6: the absolute path to the control pair-end eclip read\_2 fastq file

$7: the absolute path to the barcodecollapsepe.py python file

$8: the absolute path to the hg38\_acc\_sizes.txt

#### Part 2: DOWNSTREAM ANALYSIS (Additional)

The downstream analysis cannot be run on bridges, but it can be done locally as an additional analysis based on the discovered peaks. There are two scripts used for this part, including one bash script and one R script. To run this, give the following command in your terminal:

_bash downstream\_analysis.sh_ _ **$1 $2** _

$1: the absolute path to the reference genome annotation GFF/GTF file

$2: the absolute path to the downstream\_analysis.r R script (in case you change its directory)

### Discussion

RNA binding proteins (RBP) play important roles in the control of gene expression, the post-transcriptional regulation and the initiation of translation. Given that a single RBP can have multiple roles in a cell and could have an effect on tens of thousands of target mRNA transcripts, each RBP possesses significant research potentials. However, identification of RBPs and transcripts are difficult to predict as the interactions depend on the sequence as well as RNA structure. Though finding a comprehensive list of RBP binding sites is challenging, it is possible to find candidate motifs using tools relevant to studying and analyzing RBPs. Our focus was to build a user-friendly pipeline that identifies RBP binding motifs. It is designed so that it can be run using any experimental eCLIP data from the ENCODE website, allowing the user to conduct motif discovery on their RBP of interest found in ENCODE website.

The pipeline was run on Bridges-2 PSC using the HepG2 dataset on the TAF15 protein as the experimental data. Conducting a quality control step using _fastQC_ was an essential preprocessing step. CLIP-Seq data is prone to many PCR duplicates because CLIP-seq experiment starts with sparse materials that go through several PCR cycles. The results of the first quality control show that the experiment data has a pretty low duplication rate even before deduplication. However, even a few duplicated reads could lead to errors further downstream, so deduplication was conducted regardless. As anticipated, the duplication level of experiment data did not change much. In the control files, we observed some duplication of the reads and saw a significant drop in the duplication level after.

Before the official peak calling step, a series of quality control steps were performed to check the quality of our mapped reads given by the alignment and to see if the samples are correlated or not. Large data such as CLIP-seq is required to first check if some samples encompass major quality problems. _plotFingerprint_ exhibited that the CLIP signals were fairly good compared to the control signal. The experimental line had a prominent and steep rise at the end of the graph which is a sign of a very specific and strong CLIP enrichment. The negative control showed similar trends with the same sharp slope at the end as the CLIP experiment and this could be a sign that the CLIP experiment and the control are closely related. Hence, we checked the correlation between the experiment and the control using _plotCorrelation_. The results suggested that the plot had two disparate clusters, one for the biological replicates of the CLIP-Seq experiment and one for the replicates of control. Despite the result of _plotFingerprint_, the correlation between experiment and control is not nearly as significant as we thought. These results implied good signs for us to move on to the peak calling step.

PEAKachu was used for peak calling analysis because it allowed us to incorporate control data in finding binding regions that are significantly enriched in comparison to the experimental data. The identified peaks are shown in the plot as the read dots. The low-expressed peaks were shown as the points in a straight, diagonal line as these reads can change the fold-change quite drastically. This shows that the peak analysis is quite successful with both high-expressed and low-expressed depicted in our experiment data. We then used the peak table generated which contains 26,319 peaks to perform motif discovery using MEME Suite. MEME-ChIp discovered in total 45 significant peaks, with E-values all fewer than 0.05. Specifically there is a unique motif with many Ts with a very small E-value. Low E value means that this sequence is a significant motif, telling us that finding the motif in our peak file just by chance is very low. Other sequences in the graph have some similar repeats which have higher E-values.

After the motif discovery, we found out what the possible functions of the RBP induced genes could be. We performed a functional analysis using the RNA Centric Annotation System in the R programming language locally. Gene Ontology enrichment analysis was also performed to see the biological processes of motifs. Many pathways involved in cellular metabolic processes with all the p-values as zeros. Having p-values of zeros may be due to the fact that we did not specify the entire motifs, i.e. the statistical background of this analysis, but only the identified significantly enriched motifs as an enrichment analysis must always be done by specifying explicitly the set of genes measured. Also, the p-values of zeros could be due to failing to correct for multiple comparisons. Since ontologies have hierarchical relationships, it is important for us to apply appropriate correction factors to minimize errors. As our focus was on building a pipeline for RBP data available on ENCODE website, we conducted no further analysis on the TAF15 RBP.

We would like to mention some caveats in using our pipeline. The current pipeline is built upon a lot of other external packages which might cause problems in version compatibility. If one of the packages updated or made changes, our pipeline risks not running smoothly as expected. The pipeline requires a large amount of computational power, and so it is ideal to be run on PSC. However, because PSC has a distinct system environment than many other computation systems, the current pipeline might not be compatible with all computing environments. In the future, we may want to explore using R programming packages to run on PSC for downstream analysis and use better parameters and statistical tests to perform gene ontology enrichment analysis.

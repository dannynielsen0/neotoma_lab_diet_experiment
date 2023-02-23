#This qiime2 workflow will produce an OTU table, taxa table, and phylogenetic tree that can be then used in R 

#working directory for getting to this project

	/Volumes/GoogleDrive/My Drive/chapter 3/analyses/microbiome/C_FG_16S


# Can use conda to activate and deactivate the qiime environment, if you do not already have qiime2 installed

conda activate qiime2-2022.2 #activate environment within working dir. that has the files

#to deactivate
conda deactivate


#step one of qiime
#import data to qiime

#we will use the paired end type and input format

 qiime tools import \
 --type "SampleData[PairedEndSequencesWithQuality]" \
 --input-format PairedEndFastqManifestPhred33V2 \
 --input-path C_FG_manifest.tsv \
 --output-path qiime_output/demux_seqs.qza
  
#visualize the above  
qiime demux summarize \
  --i-data qiime_output/demux_seqs.qza \
  --o-visualization qiime_output/demux_seqs.qzv
  
#denoise with dada2

# use cutadapt to remove primer seqs first, then denoise with dada2 - seems more accurate?
#GTGCCAGCCGCCGCGGTA
##Forward primer: GTGYCAGCMGCCGCGGTAA
##Reverse primer: GGACTACNVGGGTWTCTAAT

qiime cutadapt trim-paired \
 --i-demultiplexed-sequences qiime_output/demux_seqs.qza \
 --p-front-f GTGYCAGCMGCCGCGGTAA \
 --p-front-r GGACTACNVGGGTWTCTAAT \
 --p-error-rate 0 \
 --o-trimmed-sequences qiime_output/trimmed-seqs.qza
	
#visualize the above  
qiime demux summarize \
  --i-data qiime_output/trimmed-seqs.qza \
  --o-visualization qiime_output/cutadapt_demux_seqs.qzv
  
#Then, use this if taking trimmed seqs from cutadapt, since primers were F - 19 and R - 20 BPs respectively, we can set the trunc F and R to 231 and 230...

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime_output/trimmed-seqs.qza \
  --p-trunc-len-f 231 \
  --p-trunc-len-r 230 \
  --o-table qiime_output/table.qza \
  --o-representative-sequences qiime_output/rep-seqs.qza \
  --p-n-threads 0 \
  --o-denoising-stats qiime_output/denoising-stats.qza
  
  
#generate summary visualization files
  
  qiime feature-table summarize \
  --i-table qiime_output/table.qza \
  --o-visualization qiime_output/table.qzv \
  --m-sample-metadata-file C_FG_metadata.txt

qiime feature-table tabulate-seqs \
  --i-data qiime_output/rep-seqs.qza \
  --o-visualization qiime_output/rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file qiime_output/denoising-stats.qza \
  --o-visualization qiime_output/denoising-stats.qzv


  
#phylo tree step
  
#alignment 

qiime alignment mafft \
  --i-sequences qiime_output/rep-seqs.qza \
  --o-alignment qiime_output/aligned-rep-seqs.qza
  
#mask
qiime alignment mask \
  --i-alignment qiime_output/aligned-rep-seqs.qza \
  --o-masked-alignment qiime_output/masked-aligned-rep-seqs.qza
  
#construct phylogeny #https://docs.qiime2.org/2021.8/tutorials/phylogeny/ 

qiime phylogeny fasttree \
--i-alignment qiime_output/masked-aligned-rep-seqs.qza \
--p-n-threads auto \
--o-tree qiime_output/fasttree-tree.qza


#root the phylogeny

qiime phylogeny midpoint-root \
--i-tree qiime_output/fasttree-tree.qza \
--o-rooted-tree qiime_output/fasttree-tree-rooted.qza

  
  
#download pre-trained Naive Bayes classifier

wget \
-O "silva-138-99-515-806-nb-classifier.qza" \
 "https://data.qiime2.org/2021.8/common/silva-138-99-515-806-nb-classifier.qza"

#generate taxonomy and visualize 
qiime feature-classifier classify-sklearn \
--i-classifier silva-138-99-515-806-nb-classifier.qza \
--i-reads qiime_output/rep-seqs.qza \
--o-classification qiime_output/taxonomy.qza


# remove mitochondria, chloroplasts, unassigned, 
#and those that aren't ID below kingdom level

qiime taxa filter-table \
  --i-table qiime_output/table.qza \
  --i-taxonomy qiime_output/taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Eukaryota \
  --p-include p_ \
  --o-filtered-table qiime_output/mito-chloro-filt-table.qza
  
#visualize taxonomy list  
qiime metadata tabulate \
--m-input-file qiime_output/taxonomy.qza \
--o-visualization qiime_output/taxonomy.qzv
  
  
#taxa barplot 
qiime taxa barplot \
--i-table qiime_output/mito-chloro-filt-table.qza \
--i-taxonomy qiime_output/taxonomy.qza \
--m-metadata-file C_FG_metadata.txt \
--o-visualization qiime_output/taxa-bar-plots.qzv




################
#Picrust 
################

#need to create compatible qiime enviornment

wget https://data.qiime2.org/distro/core/qiime2-2021.11-py38-osx-conda.yml
conda env create -n qiime2-2021.11 --file qiime2-2021.11-py38-osx-conda.yml

# OPTIONAL CLEANUP
rm qiime2-2021.11-py38-osx-conda.yml

#activate this environment
	conda activate qiime2-2021.11

#activate the picrust2 plugin
	conda install q2-picrust2=2021.11 -c conda-forge -c bioconda -c gavinmdouglas 



# run PICRUSt2 full pipeline
qiime picrust2 full-pipeline \
   --i-table qiime_output/mito-chloro-filt-table.qza \
   --i-seq qiime_output/rep-seqs.qza \
   --output-dir picrust2_output \
   --p-threads 8 \
   --p-hsp-method mp \
   --p-max-nsti 2 \
   --verbose

# Summarize the PICRUSt2 pathways both pathway abundance and ko_metagenome

qiime feature-table summarize \
   --i-table picrust2_output/pathway_abundance.qza \
   --o-visualization picrust2_output/pathway_abundance.qzv

qiime feature-table summarize \
   --i-table picrust2_output/ko_metagenome.qza \
   --o-visualization picrust2_output/ko_metagenome.qzv









  
  


###Start with demultiplexed fastq files

###run these qiime2 steps for foregut and caecum 16S fastq processing


conda activate qiime2-2021.8 #activate environment within working dir. that has the files

#conda deactivate


#step one of qiime
#import data to qiime

#we will use the paired end type and input format

 qiime tools import \
 --type "SampleData[PairedEndSequencesWithQuality]" \
 --input-format PairedEndFastqManifestPhred33V2 \
 --input-path qiime_input/manifest.tsv \
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
  --m-sample-metadata-file qiime_input/metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data qiime_output/rep-seqs.qza \
  --o-visualization qiime_output/rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file qiime_output/denoising-stats.qza \
  --o-visualization qiime_output/denoising-stats.qzv
  

#download ref database; this used in greengenes 13 8 SEPPw

wget \
  -O qiime_input/"sepp-refs-gg-13-8.qza" \
  "https://data.qiime2.org/2021.8/common/sepp-refs-gg-13-8.qza"
  
  
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
-O qiime_input/"gg-13-8-99-515-806-nb-classifier.qza" \
 "https://data.qiime2.org/2021.8/common/gg-13-8-99-515-806-nb-classifier.qza"

#generate taxonomy and visualize 
qiime feature-classifier classify-sklearn \
--i-classifier qiime_input/gg-13-8-99-515-806-nb-classifier.qza \
--i-reads qiime_output/rep-seqs.qza \
--o-classification qiime_output/taxonomy.qza


qiime metadata tabulate \
--m-input-file qiime_output/taxonomy.qza \
--o-visualization qiime_output/taxonomy.qzv
  
  
#taxa barplot 
qiime taxa barplot \
--i-table qiime_output/table.qza \
--i-taxonomy qiime_output/taxonomy.qza \
--m-metadata-file qiime_input/metadata.tsv \
--o-visualization qiime_output/taxa-bar-plots.qzv
  
# Woodrat experimental diet analysis

## This includes example code used for processing 16S rRNA sequences derived from woodrat foregut and caecum contents 

### The QIIME2 portion of the workflow includes commands used in shell environment

### QIIME2 artifacts are then imported into R for statistical analysis


### To generate a manifest file for qiime, run the make_manifest bash script in the directory where the fastq files are, or provide the directory as an argument when calling the script

```
bash make_manifest.sh . #here I am providing the directory with '.' as this is the directory with the fastq files
```

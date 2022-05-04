
rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/manuscript_working_directory/FG_C_qiime2/r_analysis")

library("ape")          # read tree file
library("Biostrings")   # read fasta file
library("phyloseq")     # filtering and utilities for such objects
library("biomformat") # perhaps unnecessary 
library(ggplot2)
library(picante) 
library(decontam)
library(WGCNA)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

# import QIIME2 data to phyloseq

install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
library(qiime2R)

physeq <- qza_to_phyloseq(features="../qiime_output/table.qza",tree="../qiime_output/fasttree-tree-rooted.qza",
                          taxonomy = "../qiime_output/taxonomy.qza", metadata = "../qiime_input/metadata.tsv")


#save physeq object for FG
physeq_FG <- subset_samples(physeq, physeq@sam_data$GI_region == "Foregut")
save(physeq_FG, file = "physeq_FG.rdata")

#save physeq object for C
physeq_C <- subset_samples(physeq, physeq@sam_data$GI_region == "Caecum")
save(physeq_C, file = "physeq_C.rdata")



#remove BLANKS
physeq <- subset_samples(physeq, physeq@sam_data$Sample.Type != "BLANK")
physeq_blank <- subset_samples(physeq, physeq@sam_data$Sample.Type == "BLANK")

#remove singletons 
physeq_no_singles <- filter_taxa(physeq, function (x) {sum(x>0) > 1}, prune = TRUE)



#see what is in the blank and remove potential contaminants
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample.Type == "BLANK"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant)

contamdf.prev05 <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

physeq_decon <- prune_taxa(contamdf.prev05$contaminant!="TRUE", physeq_noBlank)
physeq<- physeq_decon

#convert to RRA

physeq <- transform_sample_counts(physeq, function(x) x/sum(x)) #converts to relative abundance

##Ordinate

PCoA_cnvrg_unifrac <- ordinate(physeq, method = "PCoA", distance = "unifrac") #ordination using bray-curtis distances
PCoA_cnvrg_wunifrac <- ordinate(physeq, method = "PCoA", distance = "wunifrac") #ordination using bray-curtis distances
PCoA_cnvrg_bray <- ordinate(physeq, method = "PCoA", distance = "bray") #ordination using bray-curtis distances


PCoA_plot1 <- plot_ordination(physeq, PCoA_cnvrg_bray, color = "Species", shape = "Diet_treatment", axes = 1:2)

plot <- PCoA_plot1 + geom_point(size = 5) + facet_wrap(~GI_region) +
  scale_color_manual(values= c("forestgreen", "red")) + 
  scale_shape_manual(values=c(17,15)) + theme_bw() 


#PERMANOVA

#First, we'll make a subset of data each for FG and C
physeq_FG <- subset_samples(physeq, GI_region=="Foregut")
physeq_C <- subset_samples(physeq, GI_region=="Caecum")

#PERMANOVA on FG
metadata.FG <- as(sample_data(physeq_FG), "data.frame")
distmat.FG_bray <- phyloseq::distance(physeq_FG, method="bray")
distmat.FG_unifrac <- phyloseq::distance(physeq_FG, method="unifrac")
distmat.FG_wunifrac <- phyloseq::distance(physeq_FG, method="wunifrac")

#PERMANOVA
permanova_diet_trial_FG_bray <- adonis(distmat.FG_bray ~ Species + Diet_treatment + Species*Diet_treatment, data=metadata.FG, permutations = 999)

permanova_diet_trial_FG_unifrac <- adonis(distmat.FG_unifrac ~ Species + Diet_treatment + Species*Diet_treatment, data=metadata.FG, permutations = 999)

permanova_diet_trial_FG_wunifrac <- adonis(distmat.FG_wunifrac ~ Species + Diet_treatment + Species*Diet_treatment, data=metadata.FG, permutations = 999)


#PERMANOVA on Caecum
metadata.C <- as(sample_data(physeq_C), "data.frame")
distmat.C_bray <- phyloseq::distance(physeq_C, method="bray")
distmat.C_unifrac <- phyloseq::distance(physeq_C, method="unifrac")
distmat.C_wunifrac <- phyloseq::distance(physeq_C, method="wunifrac")

#PERMANOVA
permanova_diet_trial_C_bray <- adonis(distmat.C_bray ~ Species + Diet_treatment + Species*Diet_treatment, data=metadata.C, permutations = 999)

permanova_diet_trial_C_unifrac <- adonis(distmat.C_unifrac ~ Species + Diet_treatment + Species*Diet_treatment, data=metadata.C, permutations = 999)

permanova_diet_trial_C_wunifrac <- adonis(distmat.C_wunifrac ~ Species + Diet_treatment + Species*Diet_treatment, data=metadata.C, permutations = 999)




#############
#Filter, normalize, and visualize some bar charts
#############
physeq@sam_data$Species_gut_diet <- paste(physeq@sam_data$Species, physeq@sam_data$GI_region, 
                                          physeq@sam_data$Diet_treatment, sep="_")

#first, we will do some filtering of low abundance taxa
physeq_filt <- filter_taxa(physeq, function(x) sum(x > 5) > (0.2*length(x)), TRUE) #filter any otu with less than 5 reads in 20% of samples

#Next, take the 50 most abundant taxa and order them into new dataset
physeq_filt_50 <- prune_taxa(names(sort(taxa_sums(physeq_filt), decreasing = T)[1:50]), physeq)

family_genus <- paste(tax_table(physeq)[ ,"Family"],tax_table(physeq)[ ,"Genus"], sep = "_")
tax_table(physeq) <- cbind(tax_table(physeq), family_genus)


y1 <- tax_glom(physeq, taxrank = 'family_genus') # agglomerate taxa
y2 = merge_samples(y1, "Species_gut_diet") # merge samples on sample variable of interest
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #convert to RRA (compositional)
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$family_genus <- as.character(y4$family_genus) #convert to character
y4$family_genus[y4$Abundance < 0.025] <- "Less than 2.5% abund." #rename genera with < 1% abundance

y4 <- cbind(y4, read.table(text=y4$Sample, sep="_", header=FALSE, col.names = paste0("col", 1:3), stringsAsFactors=FALSE))
y4$geno_diet <- paste(y4$col1, y4$col3, sep="_")

y4$col1 <- factor(y4$col1, levels=c("N. bryanti", "N. lepida"))
names <- list(unique(y4$family_genus))
y4$family_genus<- factor(y4$family_genus, levels=rev(unique(y4$family_genus)))
colnames(y4)[21] <- "Family_Genus"
y4$col2 <- factor(y4$col2, levels=c("Foregut", "Caecum"))


#Then plot
library(RColorBrewer)

mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(length(unique(y4$family_genus)))

plot1 <- ggplot(y4, aes(x=geno_diet, y=Abundance, fill=family_genus, order=family_genus)) + coord_flip() +
  geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=mycolors) +
  facet_wrap(~col2) + guides(fill = guide_legend(reverse = TRUE)) +
  ylab("Percentage of Sequences")

plot1

ggsave(filename = "feeding_trials_gut_16S.jpg", plot1 , width = 12, height = 6, dpi = 300)

p <- ggplot(data=y4, aes(x=col1, y=Abundance, fill=Microbial_Family, order = Microbial_Family)) +
  geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=mycolors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="bottom", legend.title=element_text(size=16), legend.text=element_text(size=14)) + guides(fill=guide_legend(nrow=5)) + coord_flip() + facet_wrap(~col2) + 
  scale_x_discrete(limits=rev(levels(as.factor(y4$col1)))) + theme(text = element_text(size=18), strip.text = element_text(size = 20)) +
  guides(fill = guide_legend(reverse = TRUE)) + 
  ylab("16S relative read abundance") + xlab("")
p

ggsave(filename = "experimental_gut_16S.jpg", plot1 , width = 6.5, height = 3.6, dpi = 300)


###DESeq

library(DESeq2)


#combine at genus level
physeq_C <- tax_glom(physeq_C, taxrank = "Family")

physeq_C@sam_data$Species_diet <- paste(physeq_C@sam_data$Species, physeq_C@sam_data$Diet_treatment, sep="_")

Caecum_diet = phyloseq_to_deseq2(physeq_C, ~Species_diet)
Caecum_diet = DESeq(Caecum_diet, test="Wald", fitType ="parametric")


res = results(Caecum_diet, cooksCutoff = FALSE)
alpha = 0.05
alpha = 1
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_C)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab


library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Family), levels=names(x))

Caecum_plot <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip()

sigtab<-arrange(sigtab, desc(log2FoldChange))

write.csv(sigtab, "Lep_bry_deseq.csv")

masilia <- subset_taxa(physeq_pure_diet, Rank6=="Prevotella")
plot_bar(masilia, facet_grid=~Habitat)



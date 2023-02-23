
rm(list=ls())

setwd("")

library("ape")          # read tree file
library("Biostrings")   # read fasta file
library("phyloseq")     # filtering and utilities for such objects
library("biomformat") # perhaps unnecessary 
library(ggplot2)
library(picante) 
library(decontam)
library(WGCNA)


# import QIIME2 data to phyloseq

install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
library(qiime2R)

physeq <- qza_to_phyloseq(features="../qiime_output/table.qza",tree="../qiime_output/fasttree-tree-rooted.qza",
                          taxonomy = "../qiime_output/taxonomy.qza", metadata = "../qiime_input/metadata.tsv")


#save physeq object for Foregut
physeq_FG <- subset_samples(physeq, physeq@sam_data$GI_region == "Foregut")
save(physeq_FG, file = ".rdata")

#save physeq object for Caecum
physeq_C <- subset_samples(physeq, physeq@sam_data$GI_region == "Caecum")
save(physeq_C, file = ".rdata")



#remove BLANKS
physeq <- subset_samples(physeq, physeq@sam_data$Sample.Type != "BLANK")
physeq_blank <- subset_samples(physeq, physeq@sam_data$Sample.Type == "BLANK")

#remove singletons 
physeq_no_singles <- filter_taxa(physeq, function (x) {sum(x>0) > 1}, prune = TRUE)



#convert to RRA

physeq_RRA <- transform_sample_counts(physeq, function(x) x/sum(x)) #converts to relative abundance

##Ordinate

PCoA_cnvrg_unifrac <- ordinate(physeq_RRA, method = "PCoA", distance = "unifrac") #ordination using bray-curtis distances
PCoA_cnvrg_wunifrac <- ordinate(physeq_RRA, method = "PCoA", distance = "wunifrac") #ordination using bray-curtis distances
PCoA_cnvrg_bray <- ordinate(physeq_RRA, method = "PCoA", distance = "bray") #ordination using bray-curtis distances


PCoA_plot1 <- plot_ordination(physeq_RRA, PCoA_cnvrg_unifrac, color = "Species", shape = "Diet_treatment", axes = 1:2)

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

write.csv(sigtab, ".csv")

#sanity check plotting
masilia <- subset_taxa(physeq_pure_diet, Rank6=="")
plot_bar(masilia, facet_grid=~Habitat)



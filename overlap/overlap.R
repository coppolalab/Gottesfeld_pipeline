library(Biobase)
library(annotate)
library(WGCNA)
library(biomaRt)
library(broom)
library(RRHO)

library(stringr)
library(readr)
library(openxlsx)

library(rlist)
library(reshape2)
library(dplyr)
library(purrr)
library(magrittr)
library(tidyr)
library(vadr)

library(ggplot2)
library(Cairo)

GetRRHO <- function(colname1, colname2, dataset1, dataset2, symbolname, output.suffix, stepsize = 100) {
    print(colname1)
    print(colname2)
    subset1 <- select_(dataset1, symbolname, str_c("Log.Pvalue.", colname1))
    subset2 <- select_(dataset2, symbolname, str_c("Log.Pvalue.", colname2))
    output.dir <- str_c("./rrho", output.suffix, sep = "_")
    dir.create(output.dir, showWarnings = FALSE)
    rrho.out <- RRHO(subset1, subset2, alternative = "enrichment", stepsize = stepsize, labels = c(colname1, colname2), plots = TRUE, outputdir = output.dir, BY = TRUE)
    return(rrho.out)
}

source("../../FRDA project/common_functions.R")
batch1.toptable <- ReadRDSgz("../batch1_DE/save/batch1.toptable.rda")
batch2.toptable <- ReadRDSgz("../batch2_DE/save/batch2.toptable.rda")
batch3.toptable <- ReadRDSgz("../batch3_DE/save/batch3.toptable.rda")

batch1.comparisons <- ReadRDSgz("../batch1_DE/save/batch1.comparisons.rda")
batch2.comparisons <- ReadRDSgz("../batch2_DE/save/batch2.comparisons.rda")
batch3.comparisons <- ReadRDSgz("../batch3_DE/save/batch3.comparisons.rda")

batch1.comp.format <- c("FRDA vs. Control", "FRDA vs. Carrier", "Carrier vs. Control")
batch2.comp.format <- c("FRDA vs. Control (DMSO)", "109 vs. DMSO (Control)", "109 vs. DMSO (FRDA)")
batch3.comp.format <- c("FRDA vs. Control (sensory neuron)", "FRDA vs. Control (IPSC)", "IPSC vs sensory neuron (FRDA)", "IPSC vs sensory neuron (Control)")

#Overlap
batch1.voom.reduce <- select(batch1.toptable, Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
batch1.voom.reduce$Log.Pvalue.pco <- -log10(batch1.voom.reduce$P.Value.pco) * sign(batch1.voom.reduce$logFC.pco)
batch1.voom.reduce$Log.Pvalue.pca <- -log10(batch1.voom.reduce$P.Value.pca) * sign(batch1.voom.reduce$logFC.pca)
batch1.voom.reduce$Log.Pvalue.cc <- -log10(batch1.voom.reduce$P.Value.cc) * sign(batch1.voom.reduce$logFC.cc)

batch2.voom.reduce <- select(batch2.toptable, Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
batch2.voom.reduce$Log.Pvalue.pco <- -log10(batch2.voom.reduce$P.Value.pco) * sign(batch2.voom.reduce$logFC.pco)
batch2.voom.reduce$Log.Pvalue.109c <- -log10(batch2.voom.reduce$P.Value.109c) * sign(batch2.voom.reduce$logFC.109c)
batch2.voom.reduce$Log.Pvalue.109p <- -log10(batch2.voom.reduce$P.Value.109p) * sign(batch2.voom.reduce$logFC.109p)

batch3.voom.reduce <- select(batch3.toptable, Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
batch3.voom.reduce$Log.Pvalue.pcos <- -log10(batch3.voom.reduce$P.Value.pcos) * sign(batch3.voom.reduce$logFC.pcos)
batch3.voom.reduce$Log.Pvalue.pcoi <- -log10(batch3.voom.reduce$P.Value.pcoi) * sign(batch3.voom.reduce$logFC.pcoi)
batch3.voom.reduce$Log.Pvalue.pis <- -log10(batch3.voom.reduce$P.Value.pis) * sign(batch3.voom.reduce$logFC.pis)
batch3.voom.reduce$Log.Pvalue.cis <- -log10(batch3.voom.reduce$P.Value.cis) * sign(batch3.voom.reduce$logFC.cis)

#IPSC
ipsc.groups <- "Healthy.vs.FRDA"
ipsc <- ReadRDSgz("../../ipsc/save/top.object.frdah")

ipsc.reduce <- select(ipsc, logFC, P.Value)
ipsc.reduce$Log.Pvalue <- -(log10(ipsc.reduce$P.Value)) * sign(ipsc.reduce$logFC)
colnames(ipsc.reduce) <- str_c(colnames(ipsc.reduce), "ipsc", sep = ".")
ipsc.reduce$Symbol <- rownames(ipsc.reduce)

ipsc.overlap <- intersect(batch1.voom.reduce$Symbol, ipsc.reduce$Symbol)
batch1.ipsc <- filter(batch1.voom.reduce, Symbol %in% ipsc.overlap)
ipsc.human <- filter(ipsc.reduce, Symbol %in% ipsc.overlap)

ipsch.rrho <- map(batch1.comparisons, GetRRHO, "ipsc", batch1.ipsc, ipsc.human, "Symbol", "batch1_vs_ipsc", 100)
ipsch.rrho.logpval <- map(ipsch.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
ipsch.rrho.pval <- exp(-ipsch.rrho.logpval)

ipsc2.overlap <- intersect(batch2.voom.reduce$Symbol, ipsc.reduce$Symbol)
batch2.ipsc <- filter(batch2.voom.reduce, Symbol %in% ipsc2.overlap)
ipsc2.human <- filter(ipsc.reduce, Symbol %in% ipsc2.overlap)

ipsc2.rrho <- map(batch2.comparisons, GetRRHO, "ipsc", batch2.ipsc, ipsc2.human, "Symbol", "batch2_vs_ipsc", 100)
ipsc2.rrho.logpval <- map(ipsc2.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
ipsc2.rrho.pval <- exp(-ipsc2.rrho.logpval)

ipsc3.overlap <- intersect(batch3.voom.reduce$Symbol, ipsc.reduce$Symbol)
batch3.ipsc <- filter(batch3.voom.reduce, Symbol %in% ipsc3.overlap)
ipsc3.human <- filter(ipsc.reduce, Symbol %in% ipsc3.overlap)

ipsc3.rrho <- map(batch3.comparisons, GetRRHO, "ipsc", batch3.ipsc, ipsc3.human, "Symbol", "batch3_vs_ipsc", 100)
ipsc3.rrho.logpval <- map(ipsc3.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
ipsc3.rrho.pval <- exp(-ipsc3.rrho.logpval)

#Overlap between 3 and current data
batch2.batch1.overlap <- intersect(batch2.voom.reduce$Symbol, batch1.voom.reduce$Symbol)
batch2.batch1 <- filter(batch2.voom.reduce, Symbol %in% batch2.batch1.overlap) %>% filter(!duplicated(Symbol))
batch1.batch2 <- filter(batch1.voom.reduce, Symbol %in% batch2.batch1.overlap) %>% filter(!duplicated(Symbol))

batch2.batch1.rrho <- map(batch1.comparisons, mkchain( map(batch2.comparisons, GetRRHO, .,  batch2.batch1, batch1.batch2, "Symbol", "batch2_vs_batch1", 100))) %>% flatten
batch2.batch1.rrho.logpval <- map(batch2.batch1.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
batch2.batch1.rrho.pval <- exp(-batch2.batch1.rrho.logpval) %>% signif(3)
dim(batch2.batch1.rrho.pval) <- c(3,3)
rownames(batch2.batch1.rrho.pval) <- batch2.comp.format
colnames(batch2.batch1.rrho.pval) <- batch1.comp.format
write.csv(batch2.batch1.rrho.pval[c(1,2),], "overlap_batch2_batch1.csv")

batch3.batch1.overlap <- intersect(batch3.voom.reduce$Symbol, batch1.voom.reduce$Symbol)
batch3.batch1 <- filter(batch3.voom.reduce, Symbol %in% batch3.batch1.overlap) %>% filter(!duplicated(Symbol))
batch1.batch3 <- filter(batch1.voom.reduce, Symbol %in% batch3.batch1.overlap) %>% filter(!duplicated(Symbol))

batch3.batch1.rrho <- map(batch1.comparisons, mkchain( map(batch3.comparisons, GetRRHO, .,  batch3.batch1, batch1.batch3, "Symbol", "batch3_vs_batch1", 100))) %>% flatten
batch3.batch1.rrho.logpval <- map(batch3.batch1.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
batch3.batch1.rrho.pval <- exp(-batch3.batch1.rrho.logpval) %>% signif(3)
dim(batch3.batch1.rrho.pval) <- c(4,3)
rownames(batch3.batch1.rrho.pval) <- batch3.comp.format
colnames(batch3.batch1.rrho.pval) <- batch1.comp.format
write.csv(batch3.batch1.rrho.pval[c(1,2),], "overlap_batch3_batch1.csv")

batch3.batch2.overlap <- intersect(batch3.voom.reduce$Symbol, batch2.voom.reduce$Symbol)
batch3.batch2 <- filter(batch3.voom.reduce, Symbol %in% batch3.batch2.overlap) %>% filter(!duplicated(Symbol))
batch2.batch3 <- filter(batch2.voom.reduce, Symbol %in% batch3.batch2.overlap) %>% filter(!duplicated(Symbol))

batch3.batch2.rrho <- map(batch2.comparisons, mkchain( map(batch3.comparisons, GetRRHO, .,  batch3.batch2, batch2.batch3, "Symbol", "batch3_vs_batch2", 100))) %>% flatten
batch3.batch2.rrho.logpval <- map(batch3.batch2.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
batch3.batch2.rrho.pval <- exp(-batch3.batch2.rrho.logpval) %>% signif(3)
dim(batch3.batch2.rrho.pval) <- c(4,3)
rownames(batch3.batch2.rrho.pval) <- batch3.comp.format
colnames(batch3.batch2.rrho.pval) <- batch2.comp.format
write.csv(batch3.batch2.rrho.pval[c(1,2),], "overlap_batch3_batch2.csv")

#Shared genes
batch1.filter.up.strict <- filter(batch1.toptable, nchar(Symbol) > 0 & adj.P.Val.pco < 0.01 & logFC.pco > 1.5)
batch1.filter.down.strict <- filter(batch1.toptable, nchar(Symbol) > 0 & adj.P.Val.pco < 0.01 & logFC.pco < -1.5)

batch2.filter.up.strict <- filter(batch2.toptable, nchar(Symbol) > 0 & adj.P.Val.pco < 0.01)
batch2.filter.down.strict <- filter(batch2.toptable, nchar(Symbol) > 0 & adj.P.Val.pco < 0.01 & logFC.pco < -1.5)

batch3.filter.up.strict <- filter(batch3.toptable, nchar(Symbol) > 0 & adj.P.Val.pcos < 0.01 & logFC.pcos > 1.5)
batch3.filter.down.strict <- filter(batch3.toptable, nchar(Symbol) > 0 & adj.P.Val.pcos < 0.01 & logFC.pcos < -1.5)

batch1.batch2.up.strict <- intersect(batch1.filter.up.strict$Symbol, batch2.filter.up.strict$Symbol) %>% sort
batch1.batch2.down.strict <- intersect(batch1.filter.down.strict$Symbol, batch2.filter.down.strict$Symbol) %>% sort

batch1.batch3.up.strict <- intersect(batch1.filter.up.strict$Symbol, batch3.filter.up.strict$Symbol) %>% sort
batch1.batch3.down.strict <- intersect(batch1.filter.down.strict$Symbol, batch3.filter.down.strict$Symbol) %>% sort

batch2.batch3.up.strict <- intersect(batch2.filter.up.strict$Symbol, batch3.filter.up.strict$Symbol) %>% sort
batch2.batch3.down.strict <- intersect(batch2.filter.down.strict$Symbol, batch3.filter.down.strict$Symbol) %>% sort

all.up.strict <-  intersect(batch1.batch2.up.strict, batch1.batch3.up.strict) %>% sort
all.down.strict <-  intersect(batch1.batch2.down.strict, batch1.batch3.down.strict) %>% sort
all.strict <- c(all.up.strict, all.down.strict)
shared.ensembl <- filter(batch1.toptable, Symbol %in% all.strict)$Ensembl.ID

#Load expression
pheno.batch1 <- ReadRDSgz("../batch1_DE/save/pheno.batch1.rda")
pheno.batch2 <- ReadRDSgz("../batch2_DE/save/pheno.batch2.rda")
pheno.batch3 <- ReadRDSgz("../batch3_DE/save/pheno.batch3.rda")

batch1.voom <- ReadRDSgz("../batch1_DE/save/batch1.voom.rda")
batch2.voom <- ReadRDSgz("../batch2_DE/save/batch2.voom.rda")
batch3.voom <- ReadRDSgz("../batch3_DE/save/batch3.voom.rda")

batch1.shared.expr <- batch1.voom[shared.ensembl,]
batch2.shared.expr <- batch2.voom[shared.ensembl,]
batch3.shared.expr <- batch3.voom[shared.ensembl,]

rownames(batch1.shared.expr) <- all.strict
rownames(batch2.shared.expr) <- all.strict
rownames(batch3.shared.expr) <- all.strict

GenePlot <- function(gene.symbol, voom.object, pheno.object, pheno.col, levels.vector, prefix) {
    gene.expr <- voom.object[gene.symbol,]
    gene.df <- data.frame(Sample = pheno.object[[pheno.col]], Expression = gene.expr)
    gene.df$Sample %<>% factor(levels = levels.vector)

    p <- ggplot(gene.df, aes(x = Sample, y = Expression, color = Sample)) + geom_boxplot() + geom_jitter() + theme_bw()
    p <- p + theme(legend.position = "none")
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
    p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    #p <- p + ggtitle(maintitle)
    CairoPDF(str_c(prefix, gene.symbol, sep = "_"), height = 4, width = 4, bg = "transparent")
    print(p)
    dev.off()
}

batch1.levels <- c("normal", "carrier", "FRDA")
batch2.levels <- c("normal.DMSO", "FRDA.DMSO", "normal.109", "FRDA.109")
batch3.levels <- c("FRDA_sensory_neurons", "Control_sensory_neurons", "FRDA_IPSCs", "Control_IPSCs")

map(all.strict, GenePlot, batch1.shared.expr, pheno.batch1, "Disease", batch1.levels, "batch1")
map(all.strict, GenePlot, batch2.shared.expr, pheno.batch2, "Combined", batch2.levels, "batch2")
map(all.strict, GenePlot, batch3.shared.expr, pheno.batch3, "Combined", batch3.levels, "batch3")

#Less strict
batch1.filter.up <- filter(batch1.toptable, nchar(Symbol) > 0 & adj.P.Val.pco < 0.01 & logFC.pco > 0)
batch1.filter.down <- filter(batch1.toptable, nchar(Symbol) > 0 & adj.P.Val.pco < 0.01 & logFC.pco < 0)

batch2.filter.up <- filter(batch2.toptable, nchar(Symbol) > 0 & adj.P.Val.pco < 0.01 & logFC.pco > 0)
batch2.filter.down <- filter(batch2.toptable, nchar(Symbol) > 0 & adj.P.Val.pco < 0.01 & logFC.pco < 0)

batch3.filter.up <- filter(batch3.toptable, nchar(Symbol) > 0 & adj.P.Val.pcos < 0.01 & logFC.pcos > 0)
batch3.filter.down <- filter(batch3.toptable, nchar(Symbol) > 0 & adj.P.Val.pcos < 0.01 & logFC.pcos < 0)

batch1.batch2.up <- intersect(batch1.filter.up$Symbol, batch2.filter.up$Symbol) %>% sort
batch1.batch2.down <- intersect(batch1.filter.down$Symbol, batch2.filter.down$Symbol) %>% sort

batch1.batch3.up <- intersect(batch1.filter.up$Symbol, batch3.filter.up$Symbol) %>% sort
batch1.batch3.down <- intersect(batch1.filter.down$Symbol, batch3.filter.down$Symbol) %>% sort

batch2.batch3.up <- intersect(batch2.filter.up$Symbol, batch3.filter.up$Symbol) %>% sort
batch2.batch3.down <- intersect(batch2.filter.down$Symbol, batch3.filter.down$Symbol) %>% sort

all.up <-  intersect(batch1.batch2.up, batch1.batch3.up) %>% sort
all.down <-  intersect(batch1.batch2.down, batch1.batch3.down) %>% sort
write_lines(c(all.up, all.down), "shared_genes.txt")

#PBMC 2011
pbmc.groups <- str_c(c("pca", "pco", "cc"), "pbmc", sep = ".")
pbmc.comp.formats <- c("Patient vs. Carrier (PBMC)", "Patient vs. Control (PBMC)", "Carrier vs. Control (PBMC)")
    
pco.pbmc.genes <- ReadRDSgz("../../pbmc/save/top.object.pco.rda")
colnames(pco.pbmc.genes) <- str_c(colnames(pco.pbmc.genes), "pco.pbmc", sep = ".")
pco.pbmc.genes$Log.Pvalue.pco.pbmc <- (-log10(pco.pbmc.genes$P.Value.pco.pbmc) * sign(pco.pbmc.genes$logFC.pco.pbmc))

pca.pbmc.genes <- ReadRDSgz("../../pbmc/save/top.object.pca.rda")
colnames(pca.pbmc.genes) <- str_c(colnames(pca.pbmc.genes), "pca.pbmc", sep = ".")
pca.pbmc.genes$Log.Pvalue.pca.pbmc <- (-log10(pca.pbmc.genes$P.Value.pca.pbmc) * sign(pca.pbmc.genes$logFC.pca.pbmc))

cc.pbmc.genes <- ReadRDSgz("../../pbmc/save/top.object.cc.rda")
colnames(cc.pbmc.genes) <- str_c(colnames(cc.pbmc.genes), "cc.pbmc", sep = ".")
cc.pbmc.genes$Log.Pvalue.cc.pbmc <- (-log10(cc.pbmc.genes$P.Value.cc.pbmc) * sign(cc.pbmc.genes$logFC.cc.pbmc))

all.pbmc.genes <- cbind(pca.pbmc.genes, pco.pbmc.genes, cc.pbmc.genes)
all.pbmc.genes$Symbol <- rownames(all.pbmc.genes)

pbmc.reduce <- select(all.pbmc.genes, Symbol, dplyr::contains("Log.Pvalue"))

pbmc.overlap <- intersect(batch1.voom.reduce$Symbol, pbmc.reduce$Symbol)
batch1.pbmc <- filter(batch1.voom.reduce, Symbol %in% pbmc.overlap) %>% filter(!duplicated(Symbol))
pbmc.batch1 <- filter(pbmc.reduce, Symbol %in% pbmc.overlap)

pbmc1.rrho <- map(pbmc.groups, mkchain( map(batch1.comparisons, GetRRHO, .,  batch1.pbmc, pbmc.batch1, "Symbol", "batch1_vs_pbmc", 100))) %>% flatten
pbmc1.rrho.logpval <- map(pbmc1.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
pbmc1.rrho.pval <- exp(-pbmc1.rrho.logpval)
dim(pbmc1.rrho.pval) <- c(3,3)
rownames(pbmc1.rrho.pval) <- batch1.comp.format
colnames(pbmc1.rrho.pval) <- pbmc.comp.formats

pbmc.overlap2 <- intersect(batch2.voom.reduce$Symbol, pbmc.reduce$Symbol)
batch2.pbmc <- filter(batch2.voom.reduce, Symbol %in% pbmc.overlap2) %>% filter(!duplicated(Symbol))
pbmc.batch2 <- filter(pbmc.reduce, Symbol %in% pbmc.overlap2)

pbmc2.rrho <- map(pbmc.groups, mkchain( map(batch2.comparisons, GetRRHO, .,  batch2.pbmc, pbmc.batch2, "Symbol", "batch2_vs_pbmc", 100))) %>% flatten
pbmc2.rrho.logpval <- map(pbmc2.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
pbmc2.rrho.pval <- exp(-pbmc2.rrho.logpval)
dim(pbmc2.rrho.pval) <- c(3,3)
rownames(pbmc2.rrho.pval) <- batch2.comp.format
colnames(pbmc2.rrho.pval) <- pbmc.comp.formats

pbmc.overlap3 <- intersect(batch3.voom.reduce$Symbol, pbmc.reduce$Symbol)
batch3.pbmc <- filter(batch3.voom.reduce, Symbol %in% pbmc.overlap3) #%>% filter(!duplicated(Symbol))
pbmc.batch3 <- filter(pbmc.reduce, Symbol %in% pbmc.overlap3) #%>% filter(!duplicated(Symbol))

pbmc3.rrho <- map(pbmc.groups, mkchain( map(batch3.comparisons, GetRRHO, .,  batch3.pbmc, pbmc.batch3, "Symbol", "batch3_vs_pbmc", 100))) %>% flatten
pbmc3.rrho.logpval <- map(pbmc3.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
pbmc3.rrho.pval <- exp(-pbmc3.rrho.logpval)
dim(pbmc3.rrho.pval) <- c(4,3)
rownames(pbmc3.rrho.pval) <- batch3.comp.format
colnames(pbmc3.rrho.pval) <- pbmc.comp.formats

#PBMC 2016
human.groups <- str_c(c("pca", "pco", "cc"), "human", sep = ".")
pco.human.genes <- ReadRDSgz("../../FRDA\\ project/baseline_lumi/save/toptable.pco.rda")
#colnames(pco.human.genes) <- str_c(colnames(pco.human.genes), "human", sep = ".")
pco.human.genes$Log.Pvalue.pco <- (-log10(pco.human.genes$P.Value.pco) * sign(pco.human.genes$logFC.pco))

pca.human.genes <- ReadRDSgz("../../FRDA\\ project/baseline_lumi/save/toptable.pca.rda")
#colnames(pca.human.genes) <- str_c(colnames(pca.human.genes), "human", sep = ".")
pca.human.genes$Log.Pvalue.pca <- (-log10(pca.human.genes$P.Value.pca) * sign(pca.human.genes$logFC.pca))

cc.human.genes <- ReadRDSgz("../../FRDA\\ project/baseline_lumi/save/toptable.cc.rda")
#colnames(cc.human.genes) <- str_c(colnames(cc.human.genes), "human", sep = ".")
cc.human.genes$Log.Pvalue.cc <- (-log10(cc.human.genes$P.Value.cc) * sign(cc.human.genes$logFC.cc))

all.human.genes <- left_join(pca.human.genes, pco.human.genes) %>% left_join(cc.human.genes)

human.reduce <- select(all.human.genes, Symbol, dplyr::contains("Log.Pvalue"))
colnames(human.reduce)[2:ncol(human.reduce)] %<>% str_c("human", sep = ".")

human.overlap <- intersect(batch1.voom.reduce$Symbol, human.reduce$Symbol)
batch1.human <- filter(batch1.voom.reduce, Symbol %in% human.overlap) %>% filter(!duplicated(Symbol))
human.batch1 <- filter(human.reduce, Symbol %in% human.overlap)

human1.rrho <- map(human.groups, mkchain( map(batch1.comparisons, GetRRHO, .,  batch1.human, human.batch1, "Symbol", "batch1_vs_human", 100))) %>% flatten
human1.rrho.logpval <- map(human1.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
human1.rrho.pval <- exp(-human1.rrho.logpval)
dim(human1.rrho.pval) <- c(3,3)
rownames(human1.rrho.pval) <- batch1.comparisons
colnames(human1.rrho.pval) <- human.groups

human.overlap <- intersect(batch2.voom.reduce$Symbol, human.reduce$Symbol)
batch2.human <- filter(batch2.voom.reduce, Symbol %in% human.overlap) %>% filter(!duplicated(Symbol))
human.batch2 <- filter(human.reduce, Symbol %in% human.overlap)

human2.rrho <- map(human.groups, mkchain( map(batch2.comparisons, GetRRHO, .,  batch2.human, human.batch2, "Symbol", "batch2_vs_human", 100))) %>% flatten
human2.rrho.logpval <- map(human2.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
human2.rrho.pval <- exp(-human2.rrho.logpval)
dim(human2.rrho.pval) <- c(3,3)
rownames(human2.rrho.pval) <- batch2.comparisons
colnames(human2.rrho.pval) <- human.groups

human3.overlap <- intersect(batch3.voom.reduce$Symbol, human.reduce$Symbol)
batch3.human <- filter(batch3.voom.reduce, Symbol %in% human3.overlap) %>% filter(!duplicated(Symbol))
human.batch3 <- filter(human.reduce, Symbol %in% human3.overlap)

human3.rrho <- map(human.groups, mkchain( map(batch3.comparisons, GetRRHO, .,  batch3.human, human.batch3, "Symbol", "batch3_vs_human", 100))) %>% flatten
human3.rrho.logpval <- map(human3.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
human3.rrho.pval <- exp(-human3.rrho.logpval)
dim(human3.rrho.pval) <- c(4,3)
rownames(human3.rrho.pval) <- batch3.comparisons
colnames(human3.rrho.pval) <- human.groups

#Vijay mouse
vijay.drg.doxnd <- ReadRDSgz("../../Vijay_mouse/baseline_drg/save/top.object.doxnd.rda") 
vijay.drg.doxnd$Symbol <- rownames(vijay.drg.doxnd)
vijay.drg.doxnd %<>% select(Symbol, logFC, P.Value)
vijay.drg.doxnd$Log.Pvalue <- -log10(vijay.drg.doxnd$P.Value) * sign(vijay.drg.doxnd$logFC)
colnames(vijay.drg.doxnd)[2:4] %<>% str_c(".doxnd")

vijay.drg.tgwt <- ReadRDSgz("../../Vijay_mouse/baseline_drg/save/top.object.tgwt.rda")
vijay.drg.tgwt$Symbol <- rownames(vijay.drg.tgwt)
vijay.drg.tgwt %<>% select(Symbol, logFC, P.Value)
vijay.drg.tgwt$Log.Pvalue <- -log10(vijay.drg.tgwt$P.Value) * sign(vijay.drg.tgwt$logFC)
colnames(vijay.drg.tgwt)[2:4] %<>% str_c(".tgwt")

vijay.drg.rescue <- ReadRDSgz("../../Vijay_mouse/baseline_drg/save/top.object.rescue.rda")
vijay.drg.rescue$Symbol <- rownames(vijay.drg.rescue)
vijay.drg.rescue %<>% select(Symbol, logFC, P.Value)
vijay.drg.rescue$Log.Pvalue <- -log10(vijay.drg.rescue$P.Value) * sign(vijay.drg.rescue$logFC)
colnames(vijay.drg.rescue)[2:4] %<>% str_c(".rescue")

vijay.drg.all <- left_join(vijay.drg.doxnd, vijay.drg.tgwt) %>% left_join(vijay.drg.rescue)

vijay.cerebellum.doxnd <- ReadRDSgz("../../Vijay_mouse/baseline_cerebellum/save/top.object.doxnd.rda")
vijay.cerebellum.doxnd$Symbol <- rownames(vijay.cerebellum.doxnd)
vijay.cerebellum.doxnd %<>% select(Symbol, logFC, P.Value)
vijay.cerebellum.doxnd$Log.Pvalue <- -log10(vijay.cerebellum.doxnd$P.Value) * sign(vijay.cerebellum.doxnd$logFC)
colnames(vijay.cerebellum.doxnd)[2:4] %<>% str_c(".doxnd")

vijay.cerebellum.tgwt <- ReadRDSgz("../../Vijay_mouse/baseline_cerebellum/save/top.object.tgwt.rda")
vijay.cerebellum.tgwt$Symbol <- rownames(vijay.cerebellum.tgwt)
vijay.cerebellum.tgwt %<>% select(Symbol, logFC, P.Value)
vijay.cerebellum.tgwt$Log.Pvalue <- -log10(vijay.cerebellum.tgwt$P.Value) * sign(vijay.cerebellum.tgwt$logFC)
colnames(vijay.cerebellum.tgwt)[2:4] %<>% str_c(".tgwt")

vijay.cerebellum.rescue <- ReadRDSgz("../../Vijay_mouse/baseline_cerebellum/save/top.object.rescue.rda")
vijay.cerebellum.rescue$Symbol <- rownames(vijay.cerebellum.rescue)
vijay.cerebellum.rescue %<>% select(Symbol, logFC, P.Value)
vijay.cerebellum.rescue$Log.Pvalue <- -log10(vijay.cerebellum.rescue$P.Value) * sign(vijay.cerebellum.rescue$logFC)
colnames(vijay.cerebellum.rescue)[2:4] %<>% str_c(".rescue")

vijay.cerebellum.all <- left_join(vijay.cerebellum.doxnd, vijay.cerebellum.tgwt) %>% left_join(vijay.cerebellum.rescue)

vijay.comparisons <- c("doxnd", "tgwt", "rescue")
vijay.comp.format <- c("DOX vs. ND", "Tg vs. WT", "Rescue")
mouse.homology <- read_tsv("../../FRDA project/overlap/HOM_MouseHumanSequence.rpt.txt") %>% data.frame %>% select(HomoloGene.ID, NCBI.Taxon.ID, Symbol)
mouse.only <- filter(mouse.homology, NCBI.Taxon.ID == 10090)
human.only <- filter(mouse.homology, NCBI.Taxon.ID == 9606)

batch1.hom <- left_join(batch1.voom.reduce, human.only) %>% filter(!is.na(HomoloGene.ID))
vijay.drg.hom <- left_join(vijay.drg.all, mouse.only) %>% filter(!is.na(HomoloGene.ID))
vijay.cerebellum.hom <- left_join(vijay.cerebellum.all, mouse.only) %>% filter(!is.na(HomoloGene.ID))

vd.batch1.overlap <- intersect(batch1.hom$HomoloGene.ID, vijay.drg.hom$HomoloGene.ID)
batch1.vd <- filter(batch1.hom, HomoloGene.ID %in% vd.batch1.overlap) %>% filter(!duplicated(HomoloGene.ID))
vd.batch1 <- filter(vijay.drg.hom, HomoloGene.ID %in% vd.batch1.overlap) %>% filter(!duplicated(HomoloGene.ID))

batch1.vd.rrho <- map(batch1.comparisons, mkchain( map(vijay.comparisons, GetRRHO, .,  vd.batch1, batch1.vd, "HomoloGene.ID", "batch1_vs_vijay_drg", 50))) %>% flatten
batch1.vd.rrho.logpval <- map(batch1.vd.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
batch1.vd.rrho.pval <- exp(-batch1.vd.rrho.logpval) %>% signif(3)
dim(batch1.vd.rrho.pval) <- c(3,3)
rownames(batch1.vd.rrho.pval) <- vijay.comp.format
colnames(batch1.vd.rrho.pval) <- batch1.comp.format

vc.batch1.overlap <- intersect(batch1.hom$HomoloGene.ID, vijay.cerebellum.hom$HomoloGene.ID)
batch1.vc <- filter(batch1.hom, HomoloGene.ID %in% vc.batch1.overlap) %>% filter(!duplicated(HomoloGene.ID))
vc.batch1 <- filter(vijay.cerebellum.hom, HomoloGene.ID %in% vc.batch1.overlap) %>% filter(!duplicated(HomoloGene.ID))

batch1.vc.rrho <- map(batch1.comparisons, mkchain( map(vijay.comparisons, GetRRHO, .,  vc.batch1, batch1.vc, "HomoloGene.ID", "batch1_vs_vijay_cerebellum", 50))) %>% flatten
batch1.vc.rrho.logpval <- map(batch1.vc.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
batch1.vc.rrho.pval <- exp(-batch1.vc.rrho.logpval) %>% signif(3)
dim(batch1.vc.rrho.pval) <- c(3,3)
rownames(batch1.vc.rrho.pval) <- vijay.comp.format
colnames(batch1.vc.rrho.pval) <- batch1.comp.format

batch2.hom <- left_join(batch2.voom.reduce, human.only) %>% filter(!is.na(HomoloGene.ID))
vd.batch2.overlap <- intersect(batch2.hom$HomoloGene.ID, vijay.drg.hom$HomoloGene.ID)
batch2.vd <- filter(batch2.hom, HomoloGene.ID %in% vd.batch2.overlap) %>% filter(!duplicated(HomoloGene.ID))
vd.batch2 <- filter(vijay.drg.hom, HomoloGene.ID %in% vd.batch2.overlap) %>% filter(!duplicated(HomoloGene.ID))

batch2.vd.rrho <- map(batch2.comparisons, mkchain( map(vijay.comparisons, GetRRHO, .,  vd.batch2, batch2.vd, "HomoloGene.ID", "batch2_vs_vijay_drg", 50))) %>% flatten
batch2.vd.rrho.logpval <- map(batch2.vd.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
batch2.vd.rrho.pval <- exp(-batch2.vd.rrho.logpval) %>% signif(3)
dim(batch2.vd.rrho.pval) <- c(3,3)
rownames(batch2.vd.rrho.pval) <- vijay.comp.format
colnames(batch2.vd.rrho.pval) <- batch2.comp.format

vc.batch2.overlap <- intersect(batch2.hom$HomoloGene.ID, vijay.cerebellum.hom$HomoloGene.ID)
batch2.vc <- filter(batch2.hom, HomoloGene.ID %in% vc.batch2.overlap) %>% filter(!duplicated(HomoloGene.ID))
vc.batch2 <- filter(vijay.cerebellum.hom, HomoloGene.ID %in% vc.batch2.overlap) %>% filter(!duplicated(HomoloGene.ID))

batch2.vc.rrho <- map(batch2.comparisons, mkchain( map(vijay.comparisons, GetRRHO, .,  vc.batch2, batch2.vc, "HomoloGene.ID", "batch2_vs_vijay_cerebellum", 50))) %>% flatten
batch2.vc.rrho.logpval <- map(batch2.vc.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
batch2.vc.rrho.pval <- exp(-batch2.vc.rrho.logpval) %>% signif(3)
dim(batch2.vc.rrho.pval) <- c(3,3)
rownames(batch2.vc.rrho.pval) <- vijay.comp.format
colnames(batch2.vc.rrho.pval) <- batch2.comp.format

batch3.hom <- left_join(batch3.voom.reduce, human.only) %>% filter(!is.na(HomoloGene.ID))
vd.batch3.overlap <- intersect(batch3.hom$HomoloGene.ID, vijay.drg.hom$HomoloGene.ID)
batch3.vd <- filter(batch3.hom, HomoloGene.ID %in% vd.batch3.overlap) %>% filter(!duplicated(HomoloGene.ID))
vd.batch3 <- filter(vijay.drg.hom, HomoloGene.ID %in% vd.batch3.overlap) %>% filter(!duplicated(HomoloGene.ID))

batch3.vd.rrho <- map(batch3.comparisons, mkchain( map(vijay.comparisons, GetRRHO, .,  vd.batch3, batch3.vd, "HomoloGene.ID", "batch3_vs_vijay_drg", 50))) %>% flatten
batch3.vd.rrho.logpval <- map(batch3.vd.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
batch3.vd.rrho.pval <- exp(-batch3.vd.rrho.logpval) %>% signif(3)
dim(batch3.vd.rrho.pval) <- c(3,4)
rownames(batch3.vd.rrho.pval) <- vijay.comp.format
colnames(batch3.vd.rrho.pval) <- batch3.comp.format

vc.batch3.overlap <- intersect(batch3.hom$HomoloGene.ID, vijay.cerebellum.hom$HomoloGene.ID)
batch3.vc <- filter(batch3.hom, HomoloGene.ID %in% vc.batch3.overlap) %>% filter(!duplicated(HomoloGene.ID))
vc.batch3 <- filter(vijay.cerebellum.hom, HomoloGene.ID %in% vc.batch3.overlap) %>% filter(!duplicated(HomoloGene.ID))

batch3.vc.rrho <- map(batch3.comparisons, mkchain( map(vijay.comparisons, GetRRHO, .,  vc.batch3, batch3.vc, "HomoloGene.ID", "batch3_vs_vijay_cerebellum", 50))) %>% flatten
batch3.vc.rrho.logpval <- map(batch3.vc.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
batch3.vc.rrho.pval <- exp(-batch3.vc.rrho.logpval) %>% signif(3)
dim(batch3.vc.rrho.pval) <- c(3,4)
rownames(batch3.vc.rrho.pval) <- vijay.comp.format
colnames(batch3.vc.rrho.pval) <- batch3.comp.format

#Final tables
batch1.overlaps <- cbind("Patient vs. Control (IPSC)" = ipsch.rrho.pval, pbmc1.rrho.pval) 
batch1.overlaps %<>% apply(2, p.adjust, method = "fdr", n = nrow(batch1.overlaps) * ncol(batch1.overlaps)) %>% signif(3)
batch1.overlaps <- cbind(Group = rownames(batch1.overlaps), batch1.overlaps)
write.csv(batch1.overlaps, "./batch1_overlaps.csv", row.names = FALSE)

batch2.overlaps <- cbind("Patient vs. Control (IPSC)" = ipsc2.rrho.pval, pbmc2.rrho.pval)
batch2.overlaps %<>% apply(2, p.adjust, method = "fdr", n = nrow(batch2.overlaps) * ncol(batch2.overlaps)) %>% signif(3)
batch2.overlaps <- cbind(Group = rownames(batch2.overlaps), batch2.overlaps)
write.csv(batch2.overlaps, "./batch2_overlaps.csv", row.names = FALSE)

batch3.overlaps <- cbind("Patient vs. Control (IPSC)" = ipsc3.rrho.pval, pbmc3.rrho.pval)[1:2,]
batch3.overlaps %<>% apply(2, p.adjust, method = "fdr", n = nrow(batch3.overlaps) * ncol(batch3.overlaps)) %>% signif(3)
batch3.overlaps <- cbind(Group = rownames(batch3.overlaps), batch3.overlaps)
write.csv(batch3.overlaps, "./batch3_overlaps.csv", row.names = FALSE)

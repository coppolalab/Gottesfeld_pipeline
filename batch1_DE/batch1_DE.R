library(Biobase)
library(annotate)
library(WGCNA)
library(biomaRt)
library(broom)
library(RRHO)

library(edgeR)
library(limma)

library(stringr)
library(readr)
library(openxlsx)

library(rlist)
library(reshape2)
library(dplyr)
library(purrr)
library(magrittr)
library(tidyr)
library(matrixStats)

library(ggplot2)
library(Cairo)

DecidePlot <- function(file.name, decide.plot, y.lab, bar.padding = 100, pos.hjust = -0.4, neg.hjust = 1.3) {
    y.max <- max(decide.plot$Num.Genes) + nchar(max(decide.plot$Num.Genes)) * bar.padding
    y.min = min(decide.plot$Num.Genes) - nchar(abs(min(decide.plot$Num.Genes))) * bar.padding

    p <- ggplot() + geom_bar(data = filter(decide.plot, Direction == "positive"),  aes(x = Comparison, y = Num.Genes), stat = "identity", colour = "black", fill = "red", position = "dodge")
    p <- p + geom_text(data = filter(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Comparison, y = Num.Genes, hjust = pos.hjust, label = Num.Genes), position = position_dodge(width = 1))
    p <- p + geom_bar(data = filter(decide.plot, Direction == "negative"),  aes(x = Comparison, y = Num.Genes), stat = "identity", colour = "black", fill = "green", position = "dodge")
    p <- p + geom_text(data = filter(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Comparison, y = Num.Genes, hjust = neg.hjust, label = abs(Num.Genes)), position = position_dodge(width = 1))
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0)) + ylab(y.lab) + ylim(y.min, y.max)
    p <- p + theme(panel.border = element_rect(size = 1, color = "black")) + facet_grid(Threshold + Total.Genes ~ .) 
    CairoPDF(file.name, width = 7, height = 8, bg = "transparent")
    print(p)
    dev.off()
}

GetSizes <- function(p.val, pval.column, log.column, dataset) {
    dataset.sig <- filter_(dataset, str_c(pval.column, " < ", p.val))
    dataset.up <- filter_(dataset.sig, str_c(log.column, " > 0"))
    dataset.down <- filter_(dataset.sig, str_c(log.column, " < 0"))
    return(c(positive = dim(dataset.up)[1], negative = -(dim(dataset.down)[1]), Threshold = p.val))
}

GetThresholds <- function(pval.column, log.column, dataset, thresholds, pval.label) {
    nums.table <- map(thresholds, GetSizes, pval.column, log.column, dataset) %>% reduce(rbind) %>% data.frame
    nums.table$Threshold <- str_c(pval.label, " < ", nums.table$Threshold)
    return(nums.table)
}

GenWorkbook <- function(dataset, filename, pval.name, fdr.name, log.name) {
    pval.cols <- colnames(dataset) %>% str_detect(pval.name) %>% which
    adj.pval.cols <- colnames(dataset) %>% str_detect(fdr.name) %>% which
    logfc.cols <- colnames(dataset) %>% str_detect(log.name) %>% which
    description.cols <- colnames(dataset) %>% str_detect("Description") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.005", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = adj.pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = logfc.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 3:ncol(dataset), widths = "auto")
    setColWidths(wb, 1, cols = description.cols, widths = 45)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

PCAPlot <- function(filename, dataset, targetset, colorscheme = "none", variablename, plot.height = 6, plot.width = 6) {
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$SampleName, factor(targetset[[variablename]]))
    colnames(target.data) <- c("CellLine", variablename)
    colnames(dataset.plot) <- c("CellLine", "Component.1", "Component.2")

    dataset.plot <- left_join(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename, label = "CellLine")) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Component 1") + ylab("Component 2") + theme(panel.border = element_rect(size = 1, color = "black"))
    p <- p + theme(plot.background = element_blank(), legend.background = element_blank())
    CairoPDF(file = filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

VolcanoPlot <- function(top.table, filename, plot.name, cutoff = 0.05, cutoff.column = "P.Value", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Log.Pvalue") {
    top.table$Significant <- factor(top.table[[cutoff.column]] < cutoff)
    top.table$Log.Pvalue <- -log10(top.table[[cutoff.column]])
    p <- ggplot(top.table, aes_string(x = log.column, y = "Log.Pvalue")) + geom_point(aes(color = Significant))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(legend.position = "none", plot.background = element_blank())
    p <- p + xlab(xlabel) + ylab(ylabel) + theme(panel.border = element_rect(size = 1, color = "black"))
    p <- p + ggtitle(plot.name) + theme(plot.title = element_text(hjust = 0.5))
    CairoPDF(filename, width = 6, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

GetEnrichr <- function(comparison, submit.df, cutoff, logFC.cutoff, enrichr.terms) {
    comparison.pval <- str_c("adj.P.Val", comparison, sep = ".")
    comparison.logFC <- str_c("abs(logFC.", comparison, ")")
    filter.df <- filter_(submit.df, str_c(comparison.pval, "<", cutoff)) %>% filter_(str_c(comparison.logFC, " > ", logFC.cutoff)) 
    print(dim(filter.df))
    enrichr.data <- map(enrichr.terms, GetEnrichrData, filter.df, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), EnrichrWorkbook, enrichr.data, comparison)
}

EnrichrWorkbook <- function(subindex, full.df, comparison) {
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    sub.dir <- file.path("./enrichr", comparison)
    dir.create(sub.dir, showWarnings = FALSE, recursive = TRUE)
    filename = str_c(sub.dir, "/", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

EnrichrPlot <- function(enrichr.df, enrichr.expr, filename, plot.title, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Adj.P.value <- p.adjust(enrichr.df$P.value, method = "fdr")
    enrichr.df$Log.P.value <- log10(enrichr.df$Adj.P.value)
    log.column <- str_split_fixed(filename, ".", 2)[1]
    enrichr.updown <- map(enrichr.df$Genes, UpDown, enrichr.expr) %>% reduce(rbind)
    colnames(enrichr.updown) <- c("Down", "Up")
    enrichr.df <- cbind(enrichr.df, enrichr.updown)
    enrichr.df$Log.Up <- -enrichr.df$Log.P.value * enrichr.df$Up / enrichr.df$Gene.Count
    enrichr.df$Log.Down <- -enrichr.df$Log.P.value * enrichr.df$Down / enrichr.df$Gene.Count
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(-Log.P.value)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Up, Log.Down) %>% gather(Direction, Length, -Format.Name) 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Length, fill = Direction)) 
    p <- p + geom_bar(stat = "identity", size = 1) 
    #p <- p + geom_text(label = c(as.character(enrichr.df$Format.Name), rep("", nrow(enrichr.df))), hjust = "left", aes(y = 0.1)) 
    p <- p + scale_fill_discrete(name = "Direction", labels = c("Up", "Down")) 
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' P-Value')))
    p <- p + theme(plot.background = element_blank(), legend.background = element_blank(), axis.text.y = element_text(size = 12), panel.border = element_blank())
    p <- p + theme(axis.line.x = element_line(size = 1, color = "black"), panel.background = element_blank()) + geom_hline(color = "blue", yintercept = -log10(0.05)) 
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

UpDown <- function(filter.vector, enrichr.df) {
    split.vector <- str_split(filter.vector, ",")[[1]]
    enrichr.filter <- filter(enrichr.df, is.element(Symbol, split.vector))
    enrichr.vector <- c("Up" = length(which(sign(enrichr.filter$logFC) == 1)), "Down" = length(which(sign(enrichr.filter$logFC) == -1)))
    enrichr.vector
}

FilterEnrichr <- function(enrichr.df, size = 200) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    if (nrow(enrichr.df) > size) {
        enrichr.filter %<>% slice(1:size)
    }
    enrichr.filter
}

Top5Plot <- function(suffix, toptable.object, voom.object, pheno.object, pheno.col, levels.vector, maintitle, filename, direction.change = "up") {
    if(direction.change == "up") {
        col.name <- str_c("desc(logFC.", suffix, ")")
    } else {
        col.name <- str_c("logFC.", suffix)
    }

    top5.genes <- arrange_(toptable.object, col.name)$Ensembl.ID[1:5]
    top5.symbol <- arrange_(toptable.object, col.name)$Symbol[1:5]
    top5.expr <- t(voom.object$E[top5.genes,])
    colnames(top5.expr) <- top5.symbol
    top5.df <- data.frame(Sample = pheno.object[[pheno.col]], top5.expr) %>%
        gather_("Gene", "Expression", names(.)[-1])
    top5.df$Gene %<>% factor(levels = top5.symbol)
    top5.df$Sample %<>% factor(levels = levels.vector)

    p <- ggplot(top5.df, aes(x = Sample, y = Expression, color = Sample)) + geom_jitter() + theme_bw()
    p <- p + facet_wrap(~ Gene, ncol = 5, scales = "free") + theme(legend.position = "none")
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
    p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + ggtitle(maintitle)
    CairoPDF(filename, height = 4, width = 16, bg = "transparent")
    print(p)
    dev.off()
}

source("../../FRDA project/common_functions.R")
ensembl.counts <- read_csv("../raw/counts/batch1_metaReadCount.csv") %>% data.frame %>% filter(!(grepl("PAR_Y", gene_id)))
pheno.data <- read_csv("../targets.csv") 

pheno.batch1             <- filter(pheno.data, is.na(Condition2)) %>% filter(Condition != "KiPS") %>% select(-Condition2) %>% data.frame
SaveRDSgz(pheno.batch1, "./save/pheno.batch1.rda")

ensembl.batch1           <- select_(ensembl.counts, .dots = c("transcript_id", pheno.batch1$SampleName))
rownames(ensembl.batch1) <- str_replace(ensembl.batch1$transcript_id, "\\..*$", "") #Remove .## after main ENST id
batch1.counts            <- select(ensembl.batch1, -transcript_id)
batch1.filter            <- batch1.counts[rowSums(batch1.counts) > 1,]
#batch1.rowMads           <- rowMads(as.matrix(batch1.filter))
#batch1.filter            <- batch1.filter[batch1.rowMads != 0,]

#Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
vega    <- useMart("ENSEMBL_MART_VEGA", dataset = "hsapiens_gene_vega")

batch1.ensembl <- getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'ensembl_transcript_id', values = rownames(batch1.filter), mart = ensembl)
batch1.ensembl$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(batch1.ensembl) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")
batch1.ensembl %<>% filter(!duplicated(Ensembl.ID))

batch1.vega              <- getBM(attributes = c('enst_ident', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'enst_ident', values = rownames(batch1.filter), mart = vega)
colnames(batch1.vega)    <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

#Prepare collapsed version
batch1.ensembl.annot <- filter(batch1.ensembl, nchar(Symbol) > 0)
batch1.annot <- batch1.filter[batch1.ensembl.annot$Ensembl.ID,]
batch1.collapse <- collapseRows(batch1.annot, batch1.ensembl.annot$Symbol, rownames(batch1.annot))$datETcollapsed
batch1.ensembl.filter <- filter(batch1.ensembl.annot, !duplicated(Symbol))

#batch1.ensembl.keep <- filter(batch1.ensembl, nchar(Symbol) > 0)
#batch1.annot <- batch1.filter[rownames(batch1.filter) %in% batch1.ensembl.keep$Ensembl.ID, ]
batch1.dge.filter <- DGEList(batch1.collapse)
batch1.dgenorm    <- calcNormFactors(batch1.dge.filter)
batch1.design <- model.matrix(~ 0 + Disease, pheno.batch1)
colnames(batch1.design) %<>% str_replace_all("Disease", "")

#ShrinkBayes
#batch1.normfactors <- batch1.dgenorm$samples[,3]
#batch1.libsize <- colSums(batch1.annot)
#batch1.rellibsize <- batch1.libsize / exp(mean(log(batch1.libsize)))
#batch1.nf <- batch1.normfactors * batch1.rellibsize
#batch1.norm <- round(sweep(batch1.annot, 2, batch1.nf, "/"))
#SaveRDSgz(batch1.norm, "./save/batch1.norm.rda")

#disease <- pheno.batch1$Disease
#batch1.form <- ~ 1 + disease
#batch1.shrink <- ShrinkBayesWrap(batch1.norm, batch1.form)

#voom
batch1.voom <- voom(batch1.dgenorm, batch1.design, normalize = "quantile")
SaveRDSgz(batch1.voom$E, "./save/batch1.voom.rda")
batch1.variable <- batch1.voom[rowMads(batch1.voom$E) > 0, ]
batch1.variable$E <- batch1.variable$E + abs(min(batch1.variable$E))

#MDS
batch1.mds <- batch1.voom$E %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
PCAPlot("mds_batch1", batch1.mds, pheno.batch1, "none", "Disease") #label PCs by status

#limma
batch1.comparisons <- c("pco", "pca", "cc")
batch1.comp.format <- c("FRDA_vs._Control", "FRDA_vs._Carrier", "Carrier_vs._Control")
SaveRDSgz(batch1.comparisons, "./save/batch1.comparisons.rda")
SaveRDSgz(batch1.comp.format, "./save/batch1.comp.format.rda")

batch1.contrasts <- makeContrasts(FRDA - normal, FRDA - carrier, carrier - normal, levels = batch1.design)
batch1.fit    <- lmFit(batch1.variable, batch1.design) %>% contrasts.fit(batch1.contrasts)
batch1.ebayes <- eBayes(batch1.fit)

#Batch 1
toptable.pco <- topTable(batch1.ebayes, coef = 1, n = Inf) 
colnames(toptable.pco) %<>% str_c(".pco")
toptable.pco$Symbol <- rownames(toptable.pco)

toptable.pca <- topTable(batch1.ebayes, coef = 2, n = Inf) 
colnames(toptable.pca) %<>% str_c(".pca")
toptable.pca$Symbol <- rownames(toptable.pca)

toptable.cc <- topTable(batch1.ebayes, coef = 3, n = Inf) 
colnames(toptable.cc) %<>% str_c(".cc")
toptable.cc$Symbol <- rownames(toptable.cc)

VolcanoPlot(toptable.pco, filename = "volcano.pco", plot.name = "FRDA vs. Control", cutoff = 0.05, cutoff.column = "adj.P.Val.pco", log.column = "logFC.pco", ylabel = "-Log10 Adj. P-value")
VolcanoPlot(toptable.pca, filename = "volcano.pca", plot.name = "FRDA vs. Carrier", cutoff = 0.05, cutoff.column = "adj.P.Val.pca", log.column = "logFC.pca", ylabel = "-Log10 Adj. P-value")
VolcanoPlot(toptable.cc, filename = "volcano.cc", plot.name = "Carrier vs. Control", cutoff = 0.05, cutoff.column = "adj.P.Val.cc", log.column = "logFC.cc", ylabel = "-Log10 Adj. P-value")

batch1.toptable <- left_join(toptable.pco, toptable.pca) %>% 
    left_join(toptable.cc) %>% 
    left_join(batch1.ensembl.filter) %>% 
    select(Symbol, Description, Gene.Type, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val"), dplyr::contains("AveEXPR"), matches("^t\\."), matches("^B\\.")) 
SaveRDSgz(batch1.toptable, "./save/batch1.toptable.rda")
GenWorkbook(batch1.toptable, "./batch1_voom.xlsx", "P.Value", "adj.P.Val", "logFC")

thresholds     <- c(0.01, 0.005, 0.001)
thresholds.fdr <- c(0.05, 0.01)

#voom
voom.pval.batch1  <- str_c("P.Value", batch1.comparisons, sep = ".")
voom.fdr.batch1   <- str_c("adj.P.Val", batch1.comparisons, sep = ".")
voom.logfc.batch1 <- str_c("logFC", batch1.comparisons, sep = ".")

toptable.batch1.pval <- map2(voom.pval.batch1, voom.logfc.batch1, GetThresholds, batch1.toptable, thresholds, "P value")
toptable.batch1.fdr  <- map2(voom.fdr.batch1, voom.logfc.batch1, GetThresholds, batch1.toptable, thresholds.fdr, "FDR")

names(toptable.batch1.pval) <- batch1.comp.format
names(toptable.batch1.fdr)  <- batch1.comp.format

batch1.melt.pval.voom <- melt(toptable.batch1.pval)
batch1.melt.fdr.voom  <- melt(toptable.batch1.fdr)
batch1.melt.voom      <- rbind(batch1.melt.pval.voom, batch1.melt.fdr.voom)
colnames(batch1.melt.voom)[2:4] <- c("Direction", "Num.Genes", "Comparison")

batch1.genetotals.voom <- group_by(batch1.melt.voom, Threshold) %>% summarise(sum(abs(Num.Genes)))
colnames(batch1.genetotals.voom)[2] <- "Total.Genes"

batch1.toptable.plot <- left_join(batch1.melt.voom, batch1.genetotals.voom)
batch1.toptable.plot$Threshold %<>% factor
batch1.toptable.plot$Comparison %<>% str_replace_all("_", " ")
DecidePlot("batch1_toptable_thresholds", batch1.toptable.plot, "Differentially Expressed Genes", bar.padding = 200)

batch1.levels <- c("normal", "carrier", "FRDA")
batch1.top.filter <- filter(batch1.toptable, nchar(Symbol) > 0)
batch1.top.filter$Symbol %<>% str_replace("\\-", "")
Top5Plot("pco", batch1.top.filter, batch1.voom, pheno.batch1, "Disease", batch1.levels, "Upregulated", "top5.pco.up")
Top5Plot("pco", batch1.top.filter, batch1.voom, pheno.batch1, "Disease", batch1.levels, "Downregulated", "top5.pco.down", "down")
Top5Plot("pca", batch1.top.filter, batch1.voom, pheno.batch1, "Disease", batch1.levels, "Upregulated", "top5.pca.up")
Top5Plot("pca", batch1.top.filter, batch1.voom, pheno.batch1, "Disease", batch1.levels, "Downregulated", "top5.pca.down", "down")
Top5Plot("cc", batch1.top.filter, batch1.voom, pheno.batch1, "Disease", batch1.levels, "Upregulated", "top5.cc.up")
Top5Plot("cc", batch1.top.filter, batch1.voom, pheno.batch1, "Disease", batch1.levels, "Downregulated", "top5.cc.down", "down")

fxn.expr <- batch1.voom$E["ENST00000377270",]
fxn.df <- data.frame(Disease = pheno.batch1$Disease, Expression = fxn.expr)

p <- ggplot(fxn.df, aes(x = Disease, y = Expression, color = Disease)) + geom_point() + theme_bw() + ylab("voom-normalized expression")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) + ggtitle("Ensembl")
p <- p + theme(plot.background = element_blank(), legend.position = "none", panel.border = element_rect(size = 1, color = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + ggtitle("Batch1")
CairoPDF("fxn.gencode.batch1.pdf", bg = "transparent")
print(p)
dev.off()

#Enrichr
source("../../code/GO/enrichr.R")
source("../../FRDA project/common_functions.R")

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up", "GTEx_Tissue_Sample_Gene_Expression_Profiles_down") 
batch1.enrichr <- map(batch1.comparisons, GetEnrichr, batch1.top.filter, "0.01", "1.5", enrichr.terms)

#Enrichr plots
#Patient vs. Control
pco.toptable.enrichr <- filter(batch1.top.filter, abs(logFC.pco) > 1.5 & adj.P.Val.pco < 0.01) %>% select(Symbol, logFC.pco)
colnames(pco.toptable.enrichr)[2] <- "logFC"
pco.symbols <- pco.toptable.enrichr$Symbol

pco.gobiol.file <- "./enrichr/pco/GO_Biological_Process_2015.xlsx"
pco.gobiol <- read.xlsx(pco.gobiol.file) 
pco.gobiol.filter <- FilterEnrichr(pco.gobiol, 200)
GetKappaCluster(file_path_sans_ext(pco.gobiol.file), pco.gobiol.filter, pco.symbols)
pco.gobiol.final <- slice(pco.gobiol.filter, c(1, 10, 12, 32))
pco.gobiol.final$Database <- "GO Biological Process"

pco.gomole.file <- "./enrichr/pco/GO_Molecular_Function_2015.xlsx"
pco.gomole <- read.xlsx(pco.gomole.file) 
pco.gomole.filter <- FilterEnrichr(pco.gomole, 200)
GetKappaCluster(file_path_sans_ext(pco.gomole.file), pco.gomole.filter, pco.symbols)
pco.gomole.final <- slice(pco.gomole.filter, c(3, 5))
pco.gomole.final$Database <- "GO Molecular Function"

pco.reactome.file <- "./enrichr/pco/Reactome_2016.xlsx"
pco.reactome <- read.xlsx(pco.reactome.file) 
pco.reactome.filter <- FilterEnrichr(pco.reactome, 200)
GetKappaCluster(file_path_sans_ext(pco.reactome.file), pco.reactome.filter, pco.symbols)
pco.reactome.final <- slice(pco.reactome.filter, c(1, 8))
pco.reactome.final$Database <- "Reactome"

pco.kegg.file <- "./enrichr/pco/KEGG_2016.xlsx"
pco.kegg <- read.xlsx(pco.kegg.file) 
pco.kegg.filter <- FilterEnrichr(pco.kegg, 200)
GetKappaCluster(file_path_sans_ext(pco.kegg.file), pco.kegg.filter, pco.symbols)
pco.kegg.final <- slice(pco.kegg.filter, c(1))
pco.kegg.final$Database <- "KEGG"

pco.enrichr.final <- rbind(pco.gobiol.final, pco.gomole.final, pco.reactome.final, pco.kegg.final)
EnrichrPlot(pco.enrichr.final, pco.toptable.enrichr, "pco.enrichr")

#Patient vs. Carrier
pca.toptable.enrichr <- filter(batch1.top.filter, abs(logFC.pca) > 1.5 & adj.P.Val.pca < 0.01) %>% select(Symbol, logFC.pca)
colnames(pca.toptable.enrichr)[2] <- "logFC"
pca.symbols <- pco.toptable.enrichr$Symbol

pca.gobiol.file <- "./enrichr/pca/GO_Biological_Process_2015.xlsx"
pca.gobiol <- read.xlsx(pca.gobiol.file) 
pca.gobiol.filter <- FilterEnrichr(pca.gobiol)
GetKappaCluster(file_path_sans_ext(pca.gobiol.file), pca.gobiol.filter, pca.symbols)
pca.gobiol.final <- slice(pca.gobiol.filter, c(1, 3, 21))
pca.gobiol.final$Database <- "GO Biological Process"

pca.gomole.file <- "./enrichr/pca/GO_Molecular_Function_2015.xlsx"
pca.gomole <- read.xlsx(pca.gomole.file) 
pca.gomole.filter <- FilterEnrichr(pca.gomole)
GetKappaCluster(file_path_sans_ext(pca.gomole.file), pca.gomole.filter, pca.symbols)
pca.gomole.final <- slice(pca.gomole.filter, c(4, 2, 8))
pca.gomole.final$Database <- "GO Molecular Function"

pca.reactome.file <- "./enrichr/pca/Reactome_2016.xlsx"
pca.reactome <- read.xlsx(pca.reactome.file) 
pca.reactome.filter <- FilterEnrichr(pca.reactome)
GetKappaCluster(file_path_sans_ext(pca.reactome.file), pca.reactome.filter, pca.symbols)
pca.reactome.final <- slice(pca.reactome.filter, c(1, 19))
pca.reactome.final$Database <- "Reactome"

pca.kegg.file <- "./enrichr/pca/KEGG_2016.xlsx"
pca.kegg <- read.xlsx(pca.kegg.file) 
pca.kegg.filter <- FilterEnrichr(pca.kegg, 200)
GetKappaCluster(file_path_sans_ext(pca.kegg.file), pca.kegg.filter, pca.symbols)
pca.kegg.final <- slice(pca.kegg.filter, c(1))
pca.kegg.final$Database <- "KEGG"

pca.enrichr.final <- rbind(pca.gobiol.final, pca.gomole.final, pca.reactome.final, pca.kegg.final)
EnrichrPlot(pca.enrichr.final, pca.toptable.enrichr, "pca.enrichr", plot.width = 10)

#Carrier vs. Control
cc.toptable.enrichr <- filter(batch1.top.filter, abs(logFC.cc) > 1.5 & adj.P.Val.cc < 0.01) %>% select(Symbol, logFC.cc)
colnames(cc.toptable.enrichr)[2] <- "logFC"
cc.symbols <- pco.toptable.enrichr$Symbol

cc.gobiol.file <- "./enrichr/cc/GO_Biological_Process_2015.xlsx"
cc.gobiol <- read.xlsx(cc.gobiol.file) 
cc.gobiol.filter <- FilterEnrichr(cc.gobiol)
GetKappaCluster(file_path_sans_ext(cc.gobiol.file), cc.gobiol.filter, cc.symbols)
cc.gobiol.final <- slice(cc.gobiol.filter, c(1, 6, 22))
cc.gobiol.final$Database <- "GO Biological Process"

cc.gomole.file <- "./enrichr/cc/GO_Molecular_Function_2015.xlsx"
cc.gomole <- read.xlsx(cc.gomole.file) 
cc.gomole.filter <- FilterEnrichr(cc.gomole)
GetKappaCluster(file_path_sans_ext(cc.gomole.file), cc.gomole.filter, cc.symbols)
cc.gomole.final <- slice(cc.gomole.filter, c(1, 3, 12))
cc.gomole.final$Database <- "GO Molecular Function"

cc.reactome.file <- "./enrichr/cc/Reactome_2016.xlsx"
cc.reactome <- read.xlsx(cc.reactome.file) 
cc.reactome.filter <- FilterEnrichr(cc.reactome)
GetKappaCluster(file_path_sans_ext(cc.reactome.file), cc.reactome.filter, cc.symbols)
cc.reactome.final <- slice(cc.reactome.filter, c(1))
cc.reactome.final$Database <- "Reactome"

cc.kegg.file <- "./enrichr/cc/KEGG_2016.xlsx"
cc.kegg <- read.xlsx(cc.kegg.file) 
cc.kegg.filter <- FilterEnrichr(cc.kegg, 200)
GetKappaCluster(file_path_sans_ext(cc.kegg.file), cc.kegg.filter, cc.symbols)
cc.kegg.final <- slice(cc.kegg.filter, c(1, 2))
cc.kegg.final$Database <- "KEGG"

cc.enrichr.final <- rbind(cc.gobiol.final, cc.gomole.final, cc.reactome.final, cc.kegg.final)
EnrichrPlot(cc.enrichr.final, cc.toptable.enrichr, "cc.enrichr")

#FXN
#ensembl.fxn1 <- filter(ensembl.batch1, gene_name == "FXN") %>% select(-gene_name, -gene_id)
#counts.fxn1 <- data.frame(Sample = pheno.batch1$SampleName, Disease = pheno.batch1$Disease, FXN = as.vector(t(ensembl.fxn1)))

#batch1.refseq <- read_csv("./counts/rawCounts_refSeq.csv") %>% data.frame %>% select(Gene, F1_3816:F6_FA1, F10_E35:F12_E35)
#refseq.fxn1 <- filter(batch1.refseq, Gene == "FXN") %>% select(-Gene)
#refseq.fxn1.df <- data.frame(Sample = pheno.batch1$SampleName, Disease = pheno.batch1$Disease, FXN = as.vector(t(refseq.fxn1)))
#p <- ggplot(refseq.fxn1.df, aes(x = Disease, y = FXN, color = Disease)) + geom_boxplot(width = 0.5) + geom_point() + theme_bw() + ylab("Counts")
#p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) + ggtitle("RefSeq")
#p <- p + theme(plot.background = element_blank(), legend.position = "none", panel.border = element_rect(size = 1, color = "black"))
#CairoPDF("fxn.refseq.batch1.pdf", bg = "transparent")
#print(p)
#dev.off()


library(Biobase)
library(annotate)
library(WGCNA)
library(biomaRt)
library(broom)
library(matrixStats)

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

library(ggplot2)
library(Cairo)
library(RColorBrewer)

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
    target.data <- data.frame(targetset$CellLine, factor(targetset[[variablename]]))
    colnames(target.data) <- c("CellLine", variablename)
    colnames(dataset.plot) <- c("CellLine", "Component.1", "Component.2")

    dataset.plot <- merge(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename, label = "CellLine")) + geom_text()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if (colorscheme != "none") {
        p <- p + scale_fill_manual(values = colorscheme) 
    }
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
    enrichr.df$Format.Name <- str_c(enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
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

FilterEnrichr <- function(enrichr.df, size = 100) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    enrichr.filter %<>% slice(1:size)
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

pheno.batch3              <- read.xlsx("../IDs for Daniel.xlsx") 
colnames(pheno.batch3)[1] <- "CellLine"
pheno.batch3 %<>% filter(!grepl("ES4|ES5|ES6", CellLine))
pheno.batch3$cell.type    %<>% str_replace("human ", "")
pheno.batch3$Combined <- str_c(str_split_fixed(pheno.batch3$disease, " ", 4)[,2], pheno.batch3$cell.type, sep = "_") %>% str_replace("Unaffected", "Control") %>% str_replace_all(" ", "_")
SaveRDSgz(pheno.batch3, "./save/pheno.batch3.rda")

ensembl.batch3           <- read_csv("./batch3_raw/metaReadCount.csv") %>% data.frame %>% filter(!(grepl("PAR_Y", gene_id))) #PAR_Y are duplicates of genes found on X and Y chromosomes - these should be dropped, not combined
rownames(ensembl.batch3) <- str_replace(ensembl.batch3$transcript_id, "\\..*$", "") #Remove .## after main ENST id
batch3.counts            <- select(ensembl.batch3, ES1:ES3, ES7:ES18) #only counts columns
batch3.filter            <- batch3.counts[rowSums(batch3.counts) > 1,] #drop any transcript with 0 counts in all samples

#Retrieve BioMart annotations
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
vega    <- useMart("ENSEMBL_MART_VEGA", dataset = "hsapiens_gene_vega")

batch3.ensembl <- getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'ensembl_transcript_id', values = rownames(batch3.filter), mart = ensembl)
batch3.ensembl$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(batch3.ensembl) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")
batch3.ensembl %<>% filter(!duplicated(Ensembl.ID)) #Four miRNAs return multiple symbols and descriptions for the same ENST id

batch3.vega              <- getBM(attributes = c('enst_ident', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'enst_ident', values = rownames(batch3.filter), mart = vega)
colnames(batch3.vega)    <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

#Prepare collapsed version
batch3.ensembl.annot <- filter(batch3.ensembl, nchar(Symbol) > 0)
batch3.annot <- batch3.filter[batch3.ensembl.annot$Ensembl.ID,]
batch3.collapse <- collapseRows(batch3.annot, batch3.ensembl.annot$Symbol, rownames(batch3.annot))$datETcollapsed
batch3.ensembl.filter <- filter(batch3.ensembl.annot, !duplicated(Symbol))

batch3.comparisons <- c("pcos", "pcoi", "pis", "cis")
batch3.comp.format <- c("FRDA_vs._Control_(sensory_neuron)", "FRDA_vs._Control_(IPSC)", "IPSC_vs_sensory_neuron_(FRDA)", "IPSC_vs_sensory_neuron_(Control)")
SaveRDSgz(batch3.comparisons, "./save/batch3.comparisons.rda")
SaveRDSgz(batch3.comp.format, "./save/batch3.comp.format.rda")

thresholds     <- c(0.01, 0.005, 0.001)
thresholds.fdr <- c(0.05, 0.01)

#voom
voom.pval.batch3  <- str_c("P.Value", batch3.comparisons, sep = ".")
voom.fdr.batch3   <- str_c("adj.P.Val", batch3.comparisons, sep = ".")
voom.logfc.batch3 <- str_c("logFC", batch3.comparisons, sep = ".")

batch3.dge.filter <- DGEList(batch3.collapse)
batch3.dgenorm    <- calcNormFactors(batch3.dge.filter)

batch3.design <- model.matrix(~ 0 + Combined, pheno.batch3)
colnames(batch3.design) %<>% str_replace_all("Combined", "") 

batch3.voom <- voom(batch3.dgenorm, batch3.design, normalize = "quantile")
SaveRDSgz(batch3.voom$E, "./save/batch3.voom.rda")
batch3.variable <- batch3.voom[rowMads(batch3.voom$E) > 0, ]
batch3.variable$E <- batch3.variable$E + abs(min(batch3.variable$E))

batch3.mds <- batch3.voom$E %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
PCAPlot("mds_batch3", batch3.mds, pheno.batch3, "none", "Combined", plot.width = 8) #label PCs by status

batch3.contrasts <- makeContrasts(pcos = FRDA_sensory_neurons - Control_sensory_neurons, pcoi = FRDA_IPSCs - Control_IPSCs, pis = FRDA_sensory_neurons - FRDA_IPSCs, cis = Control_sensory_neurons - Control_IPSCs, levels = batch3.design)
batch3.fit    <- lmFit(batch3.variable, batch3.design) %>% contrasts.fit(batch3.contrasts)
batch3.ebayes <- eBayes(batch3.fit)

#Batch 3
toptable.pcos <- topTable(batch3.ebayes, coef = 1, n = Inf) 
colnames(toptable.pcos) %<>% str_c(".pcos")
toptable.pcos$Symbol <- rownames(toptable.pcos)

toptable.pcoi <- topTable(batch3.ebayes, coef = 2, n = Inf) 
colnames(toptable.pcoi) %<>% str_c(".pcoi")
toptable.pcoi$Symbol <- rownames(toptable.pcoi)

toptable.pis <- topTable(batch3.ebayes, coef = 3, n = Inf) 
colnames(toptable.pis) %<>% str_c(".pis")
toptable.pis$Symbol <- rownames(toptable.pis)

toptable.cis <- topTable(batch3.ebayes, coef = 4, n = Inf) 
colnames(toptable.cis) %<>% str_c(".cis")
toptable.cis$Symbol <- rownames(toptable.cis)

VolcanoPlot(toptable.pcos, filename = "volcano.pcos", plot.name = "FRDA vs. Control (sensory neurons)", cutoff = 0.05, cutoff.column = "adj.P.Val.pcos", log.column = "logFC.pcos", ylabel = "-Log10 Adj. P-value")
VolcanoPlot(toptable.pcoi, filename = "volcano.pcoi", plot.name = "FRDA vs. Control (IPSCs)", cutoff = 0.05, cutoff.column = "adj.P.Val.pcoi", log.column = "logFC.pcoi", ylabel = "-Log10 Adj. P-value")
VolcanoPlot(toptable.pis, filename = "volcano.pis", plot.name = "IPSCs vs. sensory neurons (FRDA)", cutoff = 0.05, cutoff.column = "adj.P.Val.pis", log.column = "logFC.pis", ylabel = "-Log10 Adj. P-value")
VolcanoPlot(toptable.cis, filename = "volcano.cis", plot.name = "IPSCs vs. sensory neurons (Control)", cutoff = 0.05, cutoff.column = "adj.P.Val.cis", log.column = "logFC.cis", ylabel = "-Log10 Adj. P-value")

batch3.toptable <- left_join(toptable.pcos, toptable.pcoi) %>% 
    left_join(toptable.pis) %>% 
    left_join(toptable.cis) %>% 
    left_join(batch3.ensembl.filter) %>% 
    select(Symbol, Description, Gene.Type, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val"), dplyr::contains("AveEXPR"), matches("^t\\."), matches("^B\\.")) 
SaveRDSgz(batch3.toptable, "./save/batch3.toptable.rda")

toptable.batch3.pval <- map2(voom.pval.batch3, voom.logfc.batch3, GetThresholds, batch3.toptable, thresholds, "P value")
toptable.batch3.fdr <- map2(voom.fdr.batch3, voom.logfc.batch3, GetThresholds, batch3.toptable, thresholds.fdr, "FDR")

names(toptable.batch3.pval) <- batch3.comp.format
names(toptable.batch3.fdr) <- batch3.comp.format

batch3.melt.pval.voom <- melt(toptable.batch3.pval)
batch3.melt.fdr.voom <- melt(toptable.batch3.fdr)
batch3.melt.voom <- rbind(batch3.melt.pval.voom, batch3.melt.fdr.voom) 
colnames(batch3.melt.voom)[2:4] <- c("Direction", "Num.Genes", "Comparison")

batch3.genetotals.voom <- group_by(batch3.melt.voom, Threshold) %>% summarise(sum(abs(Num.Genes)))
colnames(batch3.genetotals.voom)[2] <- "Total.Genes"

batch3.toptable.plot <- left_join(batch3.melt.voom, batch3.genetotals.voom)
batch3.toptable.plot$Threshold %<>% factor
batch3.toptable.plot$Comparison %<>% str_replace_all("_", " ")
DecidePlot("batch3_toptable_thresholds", batch3.toptable.plot, "Differentially Expressed Genes", bar.padding = 500)
GenWorkbook(batch3.toptable, "./batch3_voom.xlsx", "P.Value", "adj.P.Val", "logFC")

batch3.logFC <- filter(batch3.toptable, abs(logFC.pcos) > 1.5 & adj.P.Val.pcos < 0.01)
GenWorkbook(batch3.logFC, "./batch3_voom_pcos.xlsx", "P.Value", "adj.P.Val", "logFC")

batch3.top.filter <- filter(batch3.toptable, nchar(Symbol) > 0)
batch3.top.filter$Symbol %<>% str_replace_all("\\-", ".")
batch3.levels <- c("FRDA_sensory_neurons", "Control_sensory_neurons", "FRDA_IPSCs", "Control_IPSCs")
Top5Plot("pcos", batch3.top.filter, batch3.voom, pheno.batch3, "Combined", batch3.levels, "FRDA vs. Control (sensory neurons)", "top5.pcos.up")
Top5Plot("pcos", batch3.top.filter, batch3.voom, pheno.batch3, "Combined", batch3.levels, "FRDA vs. Control (sensory neurons)", "top5.pcos.down", "down")
Top5Plot("pcoi", batch3.top.filter, batch3.voom, pheno.batch3, "Combined", batch3.levels, "FRDA vs. Control (IPSC)")
Top5Plot("pis", batch3.top.filter, batch3.voom, pheno.batch3, "Combined", batch3.levels, "IPSC vs. sensory neurons (FRDA)")
Top5Plot("cis", batch3.top.filter, batch3.voom, pheno.batch3, "Combined", batch3.levels, "IPSC vs. sensory neurons (Control)", "top5.cis.up")
Top5Plot("cis", batch3.top.filter, batch3.voom, pheno.batch3, "Combined", batch3.levels, "IPSC vs. sensory neurons (Control)", "top5.cis.down", "down")

#FXN expression

fxn3.expr <- batch3.voom$E["ENST00000377270",]
fxn3.df <- data.frame(Combined = pheno.batch3$Combined, Counts = fxn3.expr)

p <- ggplot(fxn3.df, aes(x = Combined, y = Counts, color = Combined)) + geom_boxplot(width = 0.5) + geom_point() + theme_bw() + ylab("Voom-tranformed expression")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), legend.position = "none", panel.border = element_rect(size = 1, color = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5))
CairoPDF("fxn.gencode.batch3", bg = "transparent")
print(p)
dev.off()

#Enrichr
source("../../code/GO/enrichr.R")
source("../../FRDA project/common_functions.R")

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up", "GTEx_Tissue_Sample_Gene_Expression_Profiles_down") 
batch3.enrichr <- GetEnrichr("pcos", batch3.top.filter, "0.01", "1.5", enrichr.terms)

#Enrichr plots
#Batch 3
#FRDA vs. Control (sensory neuron)
pcos.symbols.final <- filter(batch3.top.filter, adj.P.Val.pcos < 0.01) %>% filter(abs(logFC.pcos) > 1.5) %>% select(Symbol) %>% unlist

gobiol.file.pcos <- "./enrichr/pcos/GO_Biological_Process_2015.xlsx"
gobiol.pcos <- read.xlsx(gobiol.file.pcos) 
gobiol.pcos$Num.Genes <- map(gobiol.pcos$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
gobiol.filter.pcos <- FilterEnrichr(gobiol.pcos, 200)
write.xlsx(gobiol.filter.pcos, "./enrichr/pcos/gobiol_all.xlsx")
GetKappaCluster(file_path_sans_ext(gobiol.file.pcos), gobiol.filter.pcos, pcos.symbols.final)
gobiol.final.pcos <- slice(gobiol.pcos, c(1, 3, 23, 6))

EnrichrPCA <- function(enrichr.output, gene.names) {
    num.genes <- length(gene.names)
    enrichr.list <- map(enrichr.output$Genes, str_split, ",") %>% map(getElement, 1) 
    enrichr.match <- map(enrichr.list, is.element, el = toupper(gene.names)) %>% reduce(rbind) 
    #rownames(enrichr.match) <- toupper(gene.names)
    #colnames(enrichr.match) <- enrichr.output$Term
    enrichr.match.char <- enrichr.match[rowSums(enrichr.match) > 0,] %>% apply(2, as.character)
    enrichr.match.df <- data.frame(enrichr.match.char) 
    enrichr.match.df
}

test <- EnrichrPCA(gobiol.filter.pcos, pcos.symbols.final) #%>% apply(2, as.character) %>% apply(2, as.factor)
test.mds <- MASS::mca(test)
p <- ggplot(data.frame(test.mds$rs), aes(X1, X2)) + geom_point()
plot(p)

gomole.file.pcos <- "./enrichr/pcos/GO_Molecular_Function_2015.xlsx"
gomole.pcos <- read.xlsx(gomole.file.pcos) 
gomole.pcos$Num.Genes <- map(gomole.pcos$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
gomole.filter.pcos <- FilterEnrichr(gomole.pcos, 200)
GetKappaCluster(file_path_sans_ext(gomole.file.pcos), gomole.filter.pcos, pcos.symbols.final)
gomole.final.pcos <- slice(gomole.pcos, c(1, 4, 20))

kegg.file.pcos <- "./enrichr/pcos/KEGG_2016.xlsx"
kegg.pcos <- read.xlsx(kegg.file.pcos) 
kegg.pcos$Num.Genes <- map(kegg.pcos$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
kegg.filter.pcos <- FilterEnrichr(kegg.pcos, 200)
GetKappaCluster(file_path_sans_ext(kegg.file.pcos), kegg.filter.pcos, pcos.symbols.final)
kegg.final.pcos <- slice(kegg.pcos, 5)

reactome.file.pcos <- "./enrichr/pcos/Reactome_2016.xlsx"
reactome.pcos <- read.xlsx(reactome.file.pcos) 
reactome.pcos$Num.Genes <- map(reactome.pcos$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
reactome.filter.pcos <- FilterEnrichr(reactome.pcos, 200)
GetKappaCluster(file_path_sans_ext(reactome.file.pcos), reactome.filter.pcos, pcos.symbols.final)
reactome.final.pcos <- slice(reactome.pcos, c(38, 29))
reactome.final.pcos$Term %<>% str_replace(" And.*$", "")
 
enrichr.final.pcos <- rbind(gobiol.final.pcos, reactome.final.pcos, gomole.final.pcos, kegg.final.pcos)
toptable.pcos.enrichr <- filter(batch3.toptable, !is.na(Symbol) & nchar(Symbol) > 0) %>% select(Symbol, logFC.pcos)
colnames(toptable.pcos.enrichr)[2] <- "logFC"
EnrichrPlot(enrichr.final.pcos, toptable.pcos.enrichr, "pcos.enrichr")

#Read in tables for Young and Flegel paper
young.genes <- read.delim('./Young 2014 genes.txt', sep = '\n', header = TRUE) 
young.join.pis <- left_join(young.genes, batch3.top.filter) %>% filter(adj.P.Val.pis < 0.01)
young.join.cis <- left_join(young.genes, batch3.top.filter) %>% filter(adj.P.Val.cis < 0.01 & abs(logFC.cis) > 1.5)
size.diff <- filter(batch3.top.filter, adj.P.Val.cis < 0.01 & abs(logFC.cis) > 1.5) %>% nrow
young.hyper <- phyper(nrow(young.join.cis), size.diff, nrow(batch3.top.filter) - size.diff, nrow(young.genes))

flegel.gpcrs <- read.delim('./flegel_gpcrs.txt', sep = '\n', header = TRUE, fileEncoding = 'UTF-16')
flegel.gpcrs.join.pis <- left_join(flegel.gpcrs, batch3.top.filter) %>% filter(adj.P.Val.pis < 0.01)
flegel.gpcrs.join.cis <- left_join(flegel.gpcrs, batch3.top.filter) %>% filter(adj.P.Val.cis < 0.01)

flegel.ion_channels <- read.delim('./flegel_ion_channels.txt', sep = '\n', header = TRUE, fileEncoding = 'UTF-16')
flegel.ion_channels.join.pis <- left_join(flegel.ion_channels, batch3.top.filter) %>% filter(adj.P.Val.pis < 0.01)
flegel.ion_channels.join.cis <- left_join(flegel.ion_channels, batch3.top.filter) %>% filter(adj.P.Val.cis < 0.01)

combined.sensory.genes <- c(young.join.cis$Ensembl.ID, flegel.gpcrs.join.cis$Ensembl.ID, flegel.ion_channels.join.cis$Ensembl.ID) 
combined.sensory.symbol <- c(young.join.cis$Symbol, flegel.gpcrs.join.cis$Symbol, flegel.ion_channels.join.cis$Symbol) 
pheno.cis <- filter(pheno.batch3, grepl('Control_sensory_neurons|Control_IPSCs', Combined))
sensory.expr <- batch3.voom$E[combined.sensory.genes,pheno.cis$CellLine]
rownames(sensory.expr) <- combined.sensory.symbol
sensory.expr.df <- cbind(pheno.cis, t(sensory.expr)) %>% select(CellLine, Combined, RET:KCND1) %>% gather(Gene, Expression, -CellLine, -Combined)

sensory.plot <- ggplot(sensory.expr.df, aes(x = Combined, y = Expression, color = Combined)) + geom_boxplot() + geom_jitter() + theme_bw()
sensory.plot <- sensory.plot + facet_wrap(~ Gene, ncol = 3, scales = "free") + theme(legend.position = "none")
sensory.plot <- sensory.plot + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
sensory.plot <- sensory.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
sensory.plot <- sensory.plot + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
sensory.plot <- sensory.plot + theme(plot.title = element_text(hjust = 0.5))
#sensory.plot <- sensory.plot + ggtitle(maintitle)
CairoPDF('sensory.genes', height = 12, width = 12, bg = "transparent")
print(sensory.plot)
dev.off()

#Not significant
CutoffOverlap <- function(cutoff, top.table, col.name, genes) {
    table.filter <- arrange_(top.table, col.name) %>% slice(1:cutoff)
    filter.symbols <- table.filter$Symbol
    overlap.size <- intersect(filter.symbols, genes$Symbol) %>% length
    overlap.size
}

cutoffs.vector <- map(seq(50, 2500, 50), CutoffOverlap, batch3.top.filter, 'P.Value.cis', young.genes)

#Heatmap of putative markers
liz.genes <- read.xlsx("./heat map genes.xlsx")
liz.genes$Symbol %<>% str_replace(" ", "")
liz.ensembl <- left_join(liz.genes, batch3.ensembl)

match.liz <- intersect(liz.ensembl$Ensembl.ID, rownames(batch3.voom$E)) 
labels.df <- filter(liz.ensembl, Ensembl.ID %in% match.liz)
liz.filter <- filter(liz.genes, Symbol %in% labels.df$Symbol)
labels.sort <- labels.df[match(liz.filter$Symbol, labels.df$Symbol),] 

expr.liz <- batch3.voom$E[labels.sort$Ensembl.ID,]
rownames(expr.liz) <- labels.sort$Symbol

#colnames(expr.liz) <- str_c(pheno.batch3$CellLine, " (", pheno.batch3$Combined, ")")
annotate.df <- select(pheno.batch3, Combined)
rownames(annotate.df) <- pheno.batch3$CellLine
annotate.df$Combined %<>% str_replace_all("_", " ")

CairoPDF("genes_heatmap", width = 10, height = 10, bg = "transparent")
    pheatmap(expr.liz, cluster_rows = FALSE, annotation_col = annotate.df)
    #heatmap.2(expr.liz, Rowv = NA, dendrogram = "column", scale = "none", trace = "none", col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
dev.off()

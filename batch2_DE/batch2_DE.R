library(Biobase)
library(annotate)
library(WGCNA)
library(biomaRt)
library(broom)
library(RRHO)
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
    target.data <- data.frame(targetset$SampleName.format, factor(targetset[[variablename]]))
    colnames(target.data) <- c("CellLine", variablename)
    colnames(dataset.plot) <- c("CellLine", "Component.1", "Component.2")

    dataset.plot <- merge(dataset.plot, target.data)
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

    p <- ggplot(top5.df, aes(x = Sample, y = Expression, color = Sample)) + geom_boxplot() + geom_jitter() + theme_bw()
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
ensembl.counts <- read_csv("../raw/counts/batch2_metaReadCount.csv") %>% data.frame %>% filter(!(grepl("PAR_Y", gene_id)))
pheno.data <- read_csv("../targets.csv") 

#Load raw data and extract and filter counts
pheno.batch2          <- filter(pheno.data, !is.na(Condition2)) %>% data.frame
pheno.batch2$Combined <- str_c(pheno.batch2$Disease, pheno.batch2$Condition2, sep = ".") %>% factor
pheno.batch2$Condition   %<>% factor %>% droplevels
pheno.batch2$SampleName.format <- str_replace_all(pheno.batch2$SampleName, "-", "_") 
pheno.batch2$SampleName.format <- str_c("X", pheno.batch2$SampleName.format)
SaveRDSgz(pheno.batch2, "./save/pheno.batch2.rda")

#ensembl.batch2           <- data.frame(ensembl.counts[,c("gene_id", "gene_name", as.character(pheno.batch2$CellLine))])
ensembl.batch2           <- select_(ensembl.counts, .dots = c("transcript_id", str_c(pheno.batch2$SampleName.format)))
rownames(ensembl.batch2) <- str_replace(ensembl.batch2$transcript_id, "\\..*$", "") #Remove .## after main ENST id
batch2.counts            <- select(ensembl.batch2, -transcript_id)
batch2.filter            <- batch2.counts[rowSums(batch2.counts) > 1,]

#Retrieve BioMart annotations
ensembl        <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
vega           <- useMart("ENSEMBL_MART_VEGA", dataset = "hsapiens_gene_vega")

batch2.ensembl           <- getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'ensembl_transcript_id', values = rownames(batch2.filter), mart = ensembl)
batch2.ensembl$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(batch2.ensembl) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")
batch2.ensembl %<>% filter(!duplicated(Ensembl.ID))

batch2.vega              <- getBM(attributes = c('enst_ident', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'enst_ident', values = rownames(batch2.filter), mart = vega)
colnames(batch2.vega)    <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

batch2.ensembl.annot <- filter(batch2.ensembl, nchar(Symbol) > 0)
batch2.annot <- batch2.filter[batch2.ensembl.annot$Ensembl.ID,]
batch2.collapse <- collapseRows(batch2.annot, batch2.ensembl.annot$Symbol, rownames(batch2.annot))$datETcollapsed
batch2.ensembl.filter <- filter(batch2.ensembl.annot, !duplicated(Symbol))

batch2.comparisons <- c("pco", "109c", "109p")
batch2.comp.format <- c("FRDA_vs._Control_(DMSO)", "109_vs._DMSO_(Control)", "109_vs._DMSO_(FRDA)")
SaveRDSgz(batch2.comparisons, "./save/batch2.comparisons.rda")
SaveRDSgz(batch2.comp.format, "./save/batch2.comp.format.rda")

thresholds     <- c(0.01, 0.005, 0.001)
thresholds.fdr <- c(0.05, 0.01)

#voom
voom.pval.batch2  <- str_c("P.Value", batch2.comparisons, sep = ".")
voom.fdr.batch2   <- str_c("adj.P.Val", batch2.comparisons, sep = ".")
voom.logfc.batch2 <- str_c("logFC", batch2.comparisons, sep = ".")

batch2.dge.filter <- DGEList(batch2.collapse)
batch2.dgenorm    <- calcNormFactors(batch2.dge.filter)
SaveRDSgz(batch2.dgenorm, "./save/batch2.dgenorm.rda")

batch2.design <- model.matrix(~ 0 + Combined, pheno.batch2)
colnames(batch2.design) %<>% str_replace_all("Combined", "")

batch2.voom <- voom(batch2.dgenorm, batch2.design, normalize = "quantile")
batch2.variable <- batch2.voom[rowMads(batch2.voom$E) > 0, ]
batch2.variable$E <- batch2.variable$E + abs(min(batch2.variable$E))
SaveRDSgz(batch2.voom$E, "./save/batch2.voom.rda")

batch2.mds <- batch2.voom$E %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
PCAPlot("mds_batch2", batch2.mds, pheno.batch2, "none", "Combined", plot.width = 7) #label PCs by status

pheno.batch2$Line.Condition <- factor(str_c(pheno.batch2$Condition, pheno.batch2$Condition2, sep = "."))
batch2.dupe <- duplicateCorrelation(batch2.variable, block = pheno.batch2$Line.Condition)
batch2.contrasts <- makeContrasts(FRDA.DMSO - normal.DMSO, normal.109 - normal.DMSO, FRDA.109 - FRDA.DMSO, levels = batch2.design)
batch2.fit  <- lmFit(batch2.variable, batch2.design, block = pheno.batch2$Line.Condition, correlation = batch2.dupe$consensus.correlation) %>% 
    contrasts.fit(batch2.contrasts)
batch2.fit  <- lmFit(batch2.variable, batch2.design) %>% contrasts.fit(batch2.contrasts)
batch2.ebayes <- eBayes(batch2.fit, robust = TRUE)

#Batch 2
toptable.pco2 <- topTable(batch2.ebayes, coef = 1, n = Inf) 
colnames(toptable.pco2) %<>% str_c(".pco")
toptable.pco2$Symbol <- rownames(toptable.pco2)

toptable.109c <- topTable(batch2.ebayes, coef = 2, n = Inf) 
colnames(toptable.109c) %<>% str_c(".109c")
toptable.109c$Symbol <- rownames(toptable.109c)

toptable.109p <- topTable(batch2.ebayes, coef = 3, n = Inf) 
colnames(toptable.109p) %<>% str_c(".109p")
toptable.109p$Symbol <- rownames(toptable.109p)

VolcanoPlot(toptable.pco2, filename = "volcano.pco2", plot.name = "FRDA vs. Control (DMSO)", cutoff = 0.05, cutoff.column = "adj.P.Val.pco", log.column = "logFC.pco", ylabel = "-Log10 Adj. P-value")
VolcanoPlot(toptable.109c, filename = "volcano.109c", plot.name = "109 vs. DMSO (Control)", cutoff = 0.05, cutoff.column = "adj.P.Val.109c", log.column = "logFC.109c", ylabel = "-Log10 Adj. P-value")
VolcanoPlot(toptable.109p, filename = "volcano.109p", plot.name = "109 vs. DMSO (FRDA)", cutoff = 0.05, cutoff.column = "adj.P.Val.109p", log.column = "logFC.109p", ylabel = "-Log10 Adj. P-value")

batch2.toptable <- left_join(toptable.pco2, toptable.109c) %>% 
    left_join(toptable.109p) %>% 
    left_join(batch2.ensembl.filter) %>% 
    select(Symbol, Description, Gene.Type, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val"), dplyr::contains("AveExpr"), matches("^t\\."), matches("^B\\.")) 
SaveRDSgz(batch2.toptable, "./save/batch2.toptable.rda")

toptable.batch2.pval <- map2(voom.pval.batch2, voom.logfc.batch2, GetThresholds, batch2.toptable, thresholds, "P value")
toptable.batch2.fdr <- map2(voom.fdr.batch2, voom.logfc.batch2, GetThresholds, batch2.toptable, thresholds.fdr, "FDR")

names(toptable.batch2.pval) <- batch2.comp.format
names(toptable.batch2.fdr) <- batch2.comp.format

batch2.melt.pval.voom <- melt(toptable.batch2.pval)
batch2.melt.fdr.voom <- melt(toptable.batch2.fdr)
batch2.melt.voom <- rbind(batch2.melt.pval.voom, batch2.melt.fdr.voom) 
colnames(batch2.melt.voom)[2:4] <- c("Direction", "Num.Genes", "Comparison")

batch2.genetotals.voom <- group_by(batch2.melt.voom, Threshold) %>% summarise(sum(abs(Num.Genes)))
colnames(batch2.genetotals.voom)[2] <- "Total.Genes"

batch2.toptable.plot <- left_join(batch2.melt.voom, batch2.genetotals.voom)
batch2.toptable.plot$Threshold %<>% factor
batch2.toptable.plot$Comparison %<>% str_replace_all("_", " ")

DecidePlot("batch2_toptable_thresholds", batch2.toptable.plot, "Differentially Expressed Genes", bar.padding = 200)
GenWorkbook(batch2.toptable, "./batch2_voom.xlsx", "P.Value", "adj.P.Val", "logFC")

batch2.toptable.filter <- filter(batch2.toptable, !is.na(Symbol) & nchar(Symbol) > 0) 
batch2.toptable.filter$Symbol %<>% str_replace_all("\\-", ".")
colnames(batch2.toptable.filter) %<>% str_replace_all("pco", "pco2")
batch2.levels <- c("normal.DMSO", "FRDA.DMSO", "normal.109", "FRDA.109")

Top5Plot("pco2", batch2.toptable.filter, batch2.voom, pheno.batch2, "Combined", batch2.levels, "Upregulated", "top5.pco2.up")
Top5Plot("pco2", batch2.toptable.filter, batch2.voom, pheno.batch2, "Combined", batch2.levels, "Downregulated", "top5.pco2.down", "down")
Top5Plot("109p", batch2.toptable.filter, batch2.voom, pheno.batch2, "Combined", batch2.levels, "109 vs. DMSO (FRDA)", "top5.109p.up")
Top5Plot("109p", batch2.toptable.filter, batch2.voom, pheno.batch2, "Combined", batch2.levels, "109 vs. DMSO (FRDA)", "top5.109p.down", "down")
Top5Plot("109c", batch2.toptable.filter, batch2.voom, pheno.batch2, "Combined", batch2.levels, "109 vs. DMSO (Control)", "top5.109c.up")
Top5Plot("109c", batch2.toptable.filter, batch2.voom, pheno.batch2, "Combined", batch2.levels, "109 vs. DMSO (Control)", "top5.109c.down", "down")

reversed.table.up.strict <- filter(batch2.toptable.filter, adj.P.Val.pco2 < 0.01 & logFC.pco2 > 1.5 & adj.P.Val.109p < 0.01 & logFC.109p < -1.5) %>% arrange(desc(logFC.pco2)) %>% filter(Gene.Type == "protein_coding")
reversed.table.down.strict <- filter(batch2.toptable.filter, adj.P.Val.pco2 < 0.01 & logFC.pco2 < -1.5 & adj.P.Val.109p < 0.01 & logFC.109p > 1.5) %>% arrange(logFC.pco2) %>% filter(Gene.Type == "protein_coding")

reversed.table.up <- filter(batch2.toptable.filter, adj.P.Val.pco2 < 0.01 & logFC.pco2 > 0 & adj.P.Val.109p < 0.01 & logFC.109p < 0) %>% arrange(desc(logFC.pco2)) %>% filter(Gene.Type == "protein_coding")
reversed.table.down <- filter(batch2.toptable.filter, adj.P.Val.pco2 < 0.01 & logFC.pco2 < 0 & adj.P.Val.109p < 0.01 & logFC.109p > 0) %>% arrange(logFC.pco2) %>% filter(Gene.Type == "protein_coding")

reversed.all.up <- filter(batch2.toptable.filter, adj.P.Val.pco2 < 0.01 & logFC.pco2 > 0 & logFC.109p < 0) %>% arrange(desc(logFC.pco2)) %>% filter(Gene.Type == "protein_coding")
reversed.all.down <- filter(batch2.toptable.filter, adj.P.Val.pco2 < 0.01 & logFC.pco2 < 0 & logFC.109p > 0) %>% arrange(logFC.pco2) %>% filter(Gene.Type == "protein_coding")

pco2.all.up <- filter(batch2.toptable.filter, adj.P.Val.pco2 < 0.01 & logFC.pco2 > 0) %>% arrange(desc(logFC.pco2)) %>% filter(Gene.Type == "protein_coding")
pco2.all.down <- filter(batch2.toptable.filter, adj.P.Val.pco2 < 0.01 & logFC.pco2 < 0) %>% arrange(logFC.pco2) %>% filter(Gene.Type == "protein_coding")

counts.col <- c(nrow(reversed.table.up), nrow(reversed.all.up) - nrow(reversed.table.up), nrow(pco2.all.up) - nrow(reversed.all.up), nrow(reversed.table.down), nrow(reversed.all.down) - nrow(reversed.table.down), nrow(pco2.all.down) - nrow(reversed.all.down))
direction.col <- c(rep("Upregulated", 3), rep("Downregulated", 3))
description.col <- rep(c("Reversed (Significant)", "Reversed (Not Significant)", "Not Reversed"), 2)

reverse.df <- data.frame(Num.Genes = counts.col, Direction = direction.col, Description = description.col)

p <- ggplot(reverse.df, aes(x = Direction, y = Num.Genes, fill = Description, label = Num.Genes)) + geom_bar(stat = "identity") + geom_text(position = position_stack(vjust = 0.5))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + ylab("Number of Genes") + theme(axis.title.x = element_blank(), panel.border = element_rect(size = 1, color = "black"))
p <- p + theme(axis.text.x = element_text(hjust = 1, angle = 45), plot.background = element_blank(), legend.background = element_blank())
CairoPDF(file = "reversed_genes", height = 4, width = 4, bg = "transparent")
print(p)
dev.off()

reversed.top5.up.genes <- reversed.table.up.strict$Ensembl.ID[1:5]
reversed.top5.up.symbol <- reversed.table.up.strict$Gene.Name[1:5]
reversed.top5.up.expr <- t(batch2.voom$E[reversed.top5.up.genes,])
colnames(reversed.top5.up.expr) <- make.names(reversed.top5.up.symbol)
reversed.top5.up.df <- data.frame(Combined = pheno.batch2$Combined, reversed.top5.up.expr) %>% gather(Gene, Expression, -Combined)
reversed.top5.up.df$Gene %<>% factor(levels = make.names(reversed.top5.up.symbol))
reversed.top5.up.df$Combined %<>% factor(levels = c("normal.DMSO", "FRDA.DMSO", "normal.109", "FRDA.109"))

p <- ggplot(reversed.top5.up.df, aes(x = Combined, y = Expression, color = Combined)) + geom_boxplot() + geom_point() + theme_bw()
p <- p + facet_wrap(~ Gene, ncol = 5, scales = "free") + theme(legend.position = "none")
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
CairoPDF("reversed.top5.up", height = 4, width = 16, bg = "transparent")
print(p)
dev.off()

reversed.top5.down.genes <- reversed.table.down$Ensembl.ID[1:5]
reversed.top5.down.symbol <- reversed.table.down$Gene.Name[1:5]
reversed.top5.down.expr <- t(batch2.voom$E[reversed.top5.down.genes,])
colnames(reversed.top5.down.expr) <- make.names(reversed.top5.down.symbol)
reversed.top5.down.df <- data.frame(Combined = pheno.batch2$Combined, reversed.top5.down.expr) %>% gather(Gene, Expression, -Combined)
reversed.top5.down.df$Gene %<>% factor(levels = make.names(reversed.top5.down.symbol))
reversed.top5.down.df$Combined %<>% factor(levels = c("normal.DMSO", "FRDA.DMSO", "normal.109", "FRDA.109"))

p <- ggplot(reversed.top5.down.df, aes(x = Combined, y = Expression, color = Combined)) + geom_boxplot() + geom_point() + theme_bw()
p <- p + facet_wrap(~ Gene, ncol = 5, scales = "free") + theme(legend.position = "none")
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
CairoPDF("reversed.top5.down", height = 4, width = 16, bg = "transparent")
print(p)
dev.off()

#LogFC plot
logfc.reverse <- data.frame(LogFC.pco = arrange(toptable.pco2, Ensembl.ID)$logFC.pco, LogFC.109p = arrange(toptable.109p, Ensembl.ID)$logFC.109p, LogFC.109c = arrange(toptable.109c, Ensembl.ID)$logFC.109c)
reverse.cor <- bicor(logfc.reverse$LogFC.pco, logfc.reverse$LogFC.109p)
reverse.cor <- bicor(logfc.reverse$LogFC.pco, logfc.reverse$LogFC.109p)
reverse.cor <- bicor(logfc.reverse$LogFC.pco, logfc.reverse$LogFC.109p)
reverse.pval <- corPvalueStudent(reverse.cor, nrow(logfc.reverse))

p <- ggplot(logfc.reverse, aes(LogFC.pco, LogFC.109p)) + geom_point()
p <- p + geom_smooth(method = "lm", se = FALSE)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
CairoPDF("logfc.cor", height = 6, width = 6)
print(p)
dev.off()


fxn2.expr <- batch2.voom$E["ENST00000377270",]
fxn2.df <- data.frame(Combined = pheno.batch2$Combined, Expression = fxn2.expr)

p <- ggplot(fxn2.df, aes(x = Combined, y = Expression, color = Combined)) + geom_boxplot() + geom_point() + theme_bw() + ylab("voom-normalized expression")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) + ggtitle("Batch 2")
p <- p + theme(plot.background = element_blank(), legend.position = "none", panel.border = element_rect(size = 1, color = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5))
CairoPDF("fxn.gencode.batch2.pdf", bg = "transparent")
print(p)
dev.off()

#Enrichr
source("../../code/GO/enrichr.R")
source("../../FRDA project/common_functions.R")

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up", "GTEx_Tissue_Sample_Gene_Expression_Profiles_down") 
batch2.enrichr <- map(c("pco2", "109p"), GetEnrichr, batch2.toptable.filter, "0.01", "1.5", enrichr.terms)

#Enrichr plots
#Patient vs. Control Batch 2
pco2.toptable.enrichr <- filter(batch2.top.filter, abs(logFC.pco) > 1.5 & adj.P.Val.pco < 0.01) %>% select(Symbol, logFC.pco)
colnames(pco2.toptable.enrichr)[2] <- "logFC"
pco2.symbols <- pco2.toptable.enrichr$Symbol

pco2.gobiol.file <- "./enrichr/pco2/GO_Biological_Process_2015.xlsx"
pco2.gobiol <- read.xlsx(pco2.gobiol.file) 
pco2.gobiol.filter <- FilterEnrichr(pco2.gobiol)
GetKappaCluster(file_path_sans_ext(pco2.gobiol.file), pco2.gobiol.filter, pco2.symbols)
pco2.gobiol.final <- slice(pco2.gobiol.filter, c(1, 3, 11, 12))
pco2.gobiol.final$Database <- "GO Biological Process"

pco2.gomole.file <- "./enrichr/pco2/GO_Molecular_Function_2015.xlsx"
pco2.gomole <- read.xlsx(pco2.gomole.file) 
pco2.gomole.filter <- FilterEnrichr(pco2.gomole)
GetKappaCluster(file_path_sans_ext(pco2.gomole.file), pco2.gomole.filter, pco2.symbols)
pco2.gomole.final <- slice(pco2.gomole.filter, c(1, 19))
pco2.gomole.final$Database <- "GO Molecular Function"

pco2.reactome.file <- "./enrichr/pco2/Reactome_2016.xlsx"
pco2.reactome <- read.xlsx(pco2.reactome.file) 
pco2.reactome.filter <- FilterEnrichr(pco2.reactome)
GetKappaCluster(file_path_sans_ext(pco2.reactome.file), pco2.reactome.filter, pco2.symbols)
pco2.reactome.final <- slice(pco2.reactome.filter, 1)
pco2.reactome.final$Database <- "Reactome"

pco2.kegg.file <- "./enrichr/pco2/KEGG_2016.xlsx"
pco2.kegg <- read.xlsx(pco2.kegg.file) 
pco2.kegg.filter <- FilterEnrichr(pco2.kegg)
GetKappaCluster(file_path_sans_ext(pco2.kegg.file), pco2.kegg.filter, pco2.symbols)
pco2.kegg.final <- slice(pco2.kegg.filter, 1)
pco2.kegg.final$Database <- "KEGG"

pco2.enrichr.final <- rbind(pco2.gobiol.final, pco2.gomole.final, pco2.reactome.final, pco2.kegg.final)
EnrichrPlot(pco2.enrichr.final, pco2.toptable.enrichr, "pco2.enrichr")

#109p
toptable.enrichr.109p <- filter(batch2.top.filter, abs(logFC.109p) > 1.5 & adj.P.Val.109p < 0.01) %>% select(Symbol, logFC.109p)
colnames(toptable.enrichr)[2] <- "logFC"
symbols.109p <- toptable.enrichr.109p$Symbol

gobiol.file.109p <- "./enrichr/109p/GO_Biological_Process_2015.xlsx"
gobiol.109p <- read.xlsx(gobiol.file.109p) 
gobiol.filter.109p <- FilterEnrichr(gobiol.109p)
GetKappaCluster(file_path_sans_ext(gobiol.file.109p), gobiol.filter.109p, symbols.109p)
gobiol.final.109p <- slice(gobiol.filter.109p, c(1, 32))
gobiol.final.109p$Database <- "GO Biological Process"

gomole.file.109p <- "./enrichr/109p/GO_Molecular_Function_2015.xlsx"
gomole.109p <- read.xlsx(gomole.file.109p) 
gomole.filter.109p <- FilterEnrichr(gomole.109p)
GetKappaCluster(file_path_sans_ext(gomole.file.109p), gomole.filter.109p, symbols.109p)
gomole.final.109p <- slice(gomole.filter.109p, c(1, 2))
gomole.final.109p$Database <- "GO Molecular Function"

reactome.file.109p <- "./enrichr/109p/Reactome_2016.xlsx"
reactome.109p <- read.xlsx(reactome.file.109p) 
reactome.filter.109p <- FilterEnrichr(reactome.109p)
GetKappaCluster(file_path_sans_ext(reactome.file.109p), reactome.filter.109p, symbols.109p)
reactome.final.109p <- slice(reactome.filter.109p, 2)
reactome.final.109p$Database <- "Reactome"

kegg.file.109p <- "./enrichr/109p/KEGG_2016.xlsx"
kegg.109p <- read.xlsx(kegg.file.109p) 
kegg.filter.109p <- FilterEnrichr(kegg.109p)
GetKappaCluster(file_path_sans_ext(kegg.file.109p), kegg.filter.109p, symbols.109p)
kegg.final.109p <- slice(kegg.filter.109p, 1)
kegg.final.109p$Database <- "KEGG"

enrichr.final.109p <- rbind(gobiol.final.109p, gomole.final.109p, reactome.final.109p, kegg.final.109p)
EnrichrPlot(enrichr.final.109p, toptable.enrichr.109p, "109p.enrichr")

#FXN
#ensembl.fxn <- filter(ensembl.batch2, gene_name == "FXN") %>% select(-gene_name, -gene_id)
#counts.fxn <- data.frame(Sample = pheno.batch2$SampleName, Combined = str_replace(pheno.batch2$Combined, "\\.", " "), FXN = as.vector(t(ensembl.fxn)))

#batch2.refseq <- read_csv("./counts/rawCounts_refSeq.csv") %>% data.frame %>% select(Gene, X2_4078_DMSO:X9_3816_DMSO, X11_3816_DMSO:X32_8333_109)
#refseq.fxn <- filter(batch2.refseq, Gene == "FXN") %>% select(-Gene)
#refseq.fxn.df <- data.frame(Sample = pheno.batch2$SampleName, Combined = pheno.batch2$Combined, FXN = as.vector(t(refseq.fxn)))
#p <- ggplot(refseq.fxn.df, aes(x = Combined, y = FXN, color = Combined)) + geom_boxplot(width = 0.5) + geom_point() + theme_bw() + ylab("Counts")
#p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) + ggtitle("RefSeq")
#p <- p + theme(plot.background = element_blank(), legend.position = "none", panel.border = element_rect(size = 1, color = "black"))
#CairoPDF("fxn.refseq.batch2.pdf", bg = "transparent")
#print(p)
#dev.off()

##Old RPKM stuff
#rpkm.batch2 <- read.xlsx("./RPKM_NS8081.xlsx") 
#colnames(rpkm.batch2) <- c("Symbol", make.names(colnames(batch2.counts)))
#rownames(rpkm.batch2) <- rpkm.batch2$Symbol
#rpkm.batch2 %<>% select(-Symbol)
#rpkm.df <- data.frame(Combined = pheno.batch2$Combined, FXN = t(rpkm.batch2["FXN",-1]))

#p <- ggplot(rpkm.df, aes(x = Combined, y = FXN)) + geom_point()
#CairoPDF("fxn.refseq.old.pdf")
#print(p)
#dev.off()

##Total gene counts
#old.refseq.sums <- rowSums(rpkm.batch2) %>% sort(decreasing = TRUE)
#old.refseq.sorted <- rpkm.batch2[names(old.refseq.sums),]

#batch2.refseq.nodup <- batch2.refseq[!duplicated(batch2.refseq$Gene),]
#new.refseq.counts <- select(batch2.refseq.nodup, -Gene)
#rownames(new.refseq.counts) <- batch2.refseq.nodup$Gene
#new.refseq.sums <- rowSums(new.refseq.counts) %>% sort(decreasing = TRUE)
#new.refseq.sorted <- new.refseq.counts[names(new.refseq.sums),]

#new.ensembl.counts <- select(ensembl.batch2, -gene_name, -gene_id)
#rownames(new.ensembl.counts) <- ensembl.batch2$gene_name
#new.ensembl.sums <- rowSums(new.ensembl.counts) %>% sort(decreasing = TRUE)
#new.ensembl.sorted <- new.ensembl.counts[names(new.ensembl.sums),]


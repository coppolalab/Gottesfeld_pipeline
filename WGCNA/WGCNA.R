#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()
library(biomaRt)

#For baseline processing
library(limma)
library(sva)
library(R.utils)
library(tools)
library(Biobase)
library(matrixStats)
library(broom)
library(PMCMR)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(vadr)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(igraph)

#Reading and writing tables
library(readr)
library(openxlsx)
library(parallel)

PCAPlot <- function(filename, dataset, facet.bool, size.height, size.width) {
    colnames(dataset)[2] <- "Module"
    dataset$Module %<>% str_replace("ME", "") 
    p <- ggplot(dataset, aes(x = as.numeric(x), y = value, fill = Module, color = Module)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("First Principal Component")
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(dataset$x)))
    p <- p + scale_color_manual(values = sort(unique(dataset$Module)))
    if (facet.bool == TRUE)
    {
        p <- p + facet_wrap(~ Module)
        p <- p + theme(legend.position = "none")
    } 
    CairoPDF(filename, height = size.height, width = size.width)
    print(p)
    dev.off()
}

Heatmap <- function(dataset, ME.genes) {
    color <- as.character(unique(dataset$module.colors))
    dataset %<>% select(-module.colors) %>% scale
    max.dataset <- max(abs(dataset))
    print(dim(dataset))
    CairoPDF(paste("./modules/", color, sep = ""), width = 21, height = 12)
    par(mar = c(3.5,3,2,3))
    par(oma = c(4,0,2,0))
    plotMat(dataset, zlim = c(-max.dataset, max.dataset), main = paste(color, " (", nrow(dataset), ")", sep = ""))

    ME.genes.plot <- select(ME.genes, Sample.ID, matches(color))
    p <- ggplot(ME.genes.plot, aes_string(x = "Sample.ID", y = color))
    p <- p + geom_bar(stat = "identity") + xlab("Eigengene Expression")#+ ylim(c(-6, 16)) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 2))  
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    print(p)
    dev.off()
}

EnrichrSubmit <- function(index, full.df, enrichr.terms, use.weights = FALSE) {
    dataset <- filter(full.df, Module == index)
    dir.create(file.path("./enrichr", index), recursive = TRUE, showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, GetEnrichrData, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), EnrichrWorkbook, enrichr.data, index)
}

EnrichrWorkbook <- function(subindex, full.df, index) {
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#plot.eigencor <- function(module.traits.pval, status.col, status.vector) {
    #sig.col <- paste(status.col, ".p.value", sep = "")
    #cor.status.labeled <- data.frame(Color = rownames(module.traits.pval), select_(data.frame(module.traits.pval), sig.col))
    #filter.cond <- paste(sig.col, "< 0.05")
    #status.sig <- filter_(cor.status.labeled, filter.cond)
    #me.genes.status <- select(ME.genes, one_of(as.character(status.sig$Color)))
    #me.genes.status$Status <- status.vector
    #split.cols <- str_split(status.col, "\\.")[[1]]
    #me.status.melt <- melt(me.genes.status, id.vars = "Status") %>% filter(Status == split.cols[1] | Status == split.cols[2])
    #colnames(me.status.melt)[2] <- "Module"

    #sum.fun <- function(data.vector){ data.frame(ymin = min(data.vector), ymax = max(data.vector), y = mean(data.vector)) }
    #me.status.melt$Module %<>% as.character
    #p <- ggplot(me.status.melt, aes(x = factor(Status), y = value, col = Module)) + geom_point(position = "jitter")
    #p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    #p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
    #p <- p + theme(axis.title.x = element_blank()) + stat_summary(aes(group = 1), fun.y = mean, geom = "line", col = "black", position = position_dodge(width = 0.9))
    #p <- p + scale_color_manual(values = sort(unique(me.status.melt$Module)))
    #p <- p + facet_wrap(~ Module, scales = "free_y")
    #p <- p + theme(legend.position = "none")

    #filename <- paste(split.cols[1], "_", split.cols[2], "_eigengenes_05", sep = "")
    #CairoPDF(filename, height = 13, width = 20)
    #print(p)
    #dev.off()
#}

EigengeneModel <- function(ME.vector, trait.vector, contrasts.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.contrasts <- posthoc.kruskal.dunn.test(ME ~ Trait, trait.df, p.adjust.method = "none") 
    trait.pvals <- trait.contrasts$p.value
    pval.pco <- trait.pvals["normal.DMSO","FRDA.DMSO"]
    pval.109c <- trait.pvals["normal.DMSO","normal.109"]
    pval.109p <- trait.pvals["FRDA.DMSO","FRDA.109"]
    pvals.subset <- p.adjust(c(pval.pco, pval.109c, pval.109p), method = "fdr")

    trait.medians <- group_by(trait.df, Trait) %>% summarise(median(ME)) %>% data.frame
    colnames(trait.medians)[2] <- "Median"
    rownames(trait.medians) <- trait.medians$Trait
    diff.pco <- trait.medians["FRDA.DMSO","Median"] - trait.medians["normal.DMSO","Median"]
    diff.109c <- trait.medians["normal.109","Median"] - trait.medians["normal.DMSO","Median"]
    diff.109p <- trait.medians["FRDA.109","Median"] - trait.medians["FRDA.DMSO","Median"]
    diffs.subset <- c(diff.pco, diff.109c, diff.109p)
    anova.return <- data.frame(Diff = diffs.subset, P.value = pvals.subset)
    rownames(anova.return) <- contrasts.vector

    return(anova.return)
}

EigengeneAnova <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- kruskal.test(ME ~ Trait, trait.df) %>% tidy
    return(trait.anova$p.value)
}

ModuleWorkbook <- function(dataset, filename) {
    pval.detect <- colnames(dataset) %>% str_detect("pvalue") 
    pval.cols <- which(pval.detect)
    cor.detect <- colnames(dataset) %>% str_detect("MM.*") 
    cor.cols <- which(cor.detect & !(pval.detect))
    description.cols <- colnames(dataset) %>% str_detect("Description") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.005", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = cor.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 3:ncol(dataset), widths = "auto")
    setColWidths(wb, 1, cols = description.cols, widths = 45)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

EnrichrPlot <- function(enrichr.df, filename, plot.title, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Adj.P.value <- p.adjust(enrichr.df$P.value, method = "fdr")
    enrichr.df$Log.P.value <- -log10(enrichr.df$Adj.P.value)
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.P.value)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)

    p <- ggplot(enrichr.df, aes(Format.Name, Log.P.value)) + geom_bar(stat = "identity", fill = "mediumblue", size = 1, color = "black") 
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' P-Value')))
    p <- p + theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), axis.text.y = element_text(size = 12), panel.border = element_blank(), axis.line.x = element_line(size = 1, color = "black"))
    p <- p + geom_hline(color = "red", yintercept = -log10(0.05), size = 1)
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

voom.raw <- ReadRDSgz("../batch2_DE/save/batch2.voom.rda")
voom.top <- ReadRDSgz("../batch2_DE/save/batch2.toptable.rda") %>% filter(nchar(Symbol) > 0)
voom.annot <- voom.raw[voom.top$Ensembl.ID,]
voom.mads <- rowMads(voom.annot)
voom.expr <- voom.annot[voom.mads != 0,]
voom.top.expr <- filter(voom.top, Ensembl.ID %in% rownames(voom.expr))
voom.collapse <- collapseRows(voom.expr, voom.top.expr$Symbol, rownames(voom.expr))$datETcollapsed

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

sft <- pickSoftThreshold(t(voom.collapse), powerVector = powers, verbose = 5, networkType = "signed")
sft.bicor <- pickSoftThreshold(t(voom.collapse), powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.spearman <- pickSoftThreshold(t(voom.collapse), powerVector = powers, verbose = 5, corFnc = cor, corOptions = list(method = "spearman"), networkType = "signed")
sft.df <- sft.bicor$fitIndices
SaveRDSgz(sft.bicor, file = "./save/sft.bicor.rda")

#Plot scale indendence and mean connectivity as functions of power
sft.df$multiplied <- sft.df$SFT.R.sq * -sign(sft.df$slope)
p <- ggplot(sft.df, aes(x = Power,  y = multiplied, label = Power)) + geom_point() + geom_text(vjust = -0.6, size = 4, col = "red")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(aes(yintercept = 0.9))
p <- p + xlab("Soft Threshold") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale Independence")
CairoPDF(file = "./scaleindependence", width = 6, height = 6)
print(p)
dev.off()

p <- ggplot(sft.df, aes(x = Power,  y = mean.k.)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Soft Threshold") + ylab("Mean Connectivity") + ggtitle("Mean Connectivity")
CairoPDF(file = "./meanconnectivity", width = 6, height = 6)
print(p)
dev.off()

softPower <- 7
adjacency.expr <- adjacency(t(voom.collapse), power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
#SaveRDSgz(adjacency.expr, file = "./save/adjacency.rda")
gc()

TOM <- TOMsimilarity(adjacency.expr, verbose = 5)
dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()
#SaveRDSgz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
SaveRDSgz(geneTree, file = "./save/gene.tree.rda")
gc()

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 20

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
SaveRDSgz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(t(voom.collapse), colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
SaveRDSgz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50_tree", height = 10, width = 15, bg = "transparent")
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(t(voom.collapse), dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering", height = 10, width = 15, bg = "transparent")
plotDendroAndColors(geneTree, cbind("Original Modules" = dynamic.colors, "Merged Modules" = merged.colors), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
SaveRDSgz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
SaveRDSgz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
SaveRDSgz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 8, bg = "transparent")
par(bg = "transparent")
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

#modules.out <- data.frame(Symbol = colnames(t(voom.collapse)), Module = module.colors)
#write.xlsx(modules.out, "modules_out.xlsx")

pheno.import <- read_csv("../targets.csv") %>% filter(!is.na(Condition2))
pheno.import$Combined <- str_c(pheno.import$Disease, pheno.import$Condition2, sep = ".") %>% factor
kw.model <- map_dbl(ME.genes, EigengeneAnova, pheno.import$Combined) %>% p.adjust("fdr") %>% signif(3)

combined.design <- model.matrix(~ 0 + Combined, data = pheno.import)
colnames(combined.design) %<>% str_replace("Combined", "")
anova.contrasts <- makeContrasts(FRDA.DMSO - normal.DMSO, normal.109 - normal.DMSO, FRDA.109 - normal.109, levels = combined.design)
contrasts(pheno.import$Combined) <- anova.contrasts

color.values <- unique(module.colors)
cor.status <- map(ME.genes, EigengeneModel, pheno.import$Combined, colnames(anova.contrasts))
status.diff <- map(cor.status, select, Diff) %>% map(t) %>% reduce(rbind)
rownames(status.diff) <- names(cor.status)
colnames(status.diff) %<>% str_replace(" - ", "_vs._")
status.pval <- map(cor.status, select, P.value) %>% map(t) %>% reduce(rbind) 
rownames(status.pval) <- names(cor.status)
colnames(status.pval) %<>% str_replace(" - ", "_vs._")

pval.adjust <- map(data.frame(status.pval), p.adjust, method = "fdr", n = length(color.values)) %>% reduce(cbind) %>% data.frame
rownames(pval.adjust) <- rownames(status.pval)
colnames(pval.adjust) <- str_c(colnames(status.pval), ".pval")

text.matrix.traits <- paste(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(pval.adjust), 1), ')', sep = '')
dim(text.matrix.traits) = dim(status.diff)

#TextHeatmap(as.matrix(status.diff), text.matrix.traits, colnames(status.diff), colnames(ME.genes), "", "module-trait relationships", heatmap.range)

heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
width.dynamic <- 3 + (1 * ncol(text.matrix.traits))
CairoPDF("module_trait_relationships", width = width.dynamic, height = 10, bg = "transparent")
par(mar = c(10, 8, 3, 3))
labeledHeatmap(Matrix = as.matrix(status.diff), xLabels = colnames(status.diff), yLabels = colnames(ME.genes), ySymbols = colnames(ME.genes), yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix.traits, setStdMargins = F, zlim = heatmap.range, main = "")
dev.off()

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description', 'gene_biotype'), filters = 'hgnc_symbol', values = as.character(rownames(voom.collapse)), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "") 
colnames(bm.table) <- c("Symbol", "Description", "Gene.Type")
bm.table %<>% filter(!duplicated(Symbol))


all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Symbol = rownames(all.degrees), Module = module.colors, all.degrees, stringsAsFactors = FALSE) %>% arrange(Module)
gene.info$kscaled <- by(gene.info, gene.info$Module, select, kWithin) %>% map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)
gene.info[3:ncol(gene.info)] <- signif(gene.info[3:ncol(gene.info)], 3)
gene.info.annot <- arrange(gene.info, Module, desc(kscaled)) %>% left_join(bm.table) %>% select(Symbol, Description, Module:kscaled)

gene.module.membership <- data.frame(bicor(t(voom.collapse), ME.genes, maxPOutliers = 0.05)) %>% signif(3)
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(t(voom.collapse)))) %>% signif(3)
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

module.membership <- left_join(gene.info.annot, gene.module.membership) %>% left_join(module.membership.pvalue, by = "Symbol")

#Annotate gene table
ModuleWorkbook(module.membership, "./module_membership.xlsx")

#gene.info <- data.frame(Symbol = rownames(all.degrees), module.color = module.colors, all.degrees)
#gene.info$kscaled <- by(gene.info, gene.info$module.color, select, kWithin) %>% map(function(x) { x / max (x) }) %>% reduce(c)

#gene.module.membership <- as.data.frame(bicor(t(voom.collapse), ME.genes, maxPOutliers = 0.05))
#module.membership.pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene.module.membership), ncol(voom.collapse)))
#names(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
#colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")

#module.membership <- cbind(gene.info, gene.module.membership, module.membership.pvalue)
#colnames(module.membership)[1] <- "Symbol"

#final.genetable <- join(module.membership, voom.bm.table) %>% 
    #dplyr::select(Symbol, Description, Gene.Type, module.color, kTotal:kscaled, matches("MM.*")) %>%
    #arrange(module.color, desc(kscaled))

#SaveRDSgz(gene.info, file = "./save/gene.info.rda")
#write_csv(data.frame(table(module.colors)), path = "./final_eigengenes.csv") 
#write_csv(module.membership, "module_membership.csv")

#colnames(gene.module.membership) %<>% str_replace("MM.", "")
#colnames(module.membership.pvalue) %<>% str_replace("MM.pvalue.", "")
#gene.module.membership$Symbol <- rownames(gene.module.membership)
#module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

#gene.module.membership.long <- gather_(data = gene.module.membership, "module.comparison", "correlation", colnames(ME.genes))
#module.membership.pvalue.long <- gather_(data = module.membership.pvalue, "module.comparison", "p.value", colnames(ME.genes))
#membership.join <- join(gene.module.membership.long, module.membership.pvalue.long)
#eigengene.connectivity.long <- join(membership.join, gene.info) %>% select(Symbol, module.color:kscaled, module.comparison:p.value)
#write_csv(eigengene.connectivity.long, "eigengene_connectivity_long.csv")


all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y")
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$x <- as.factor(1:nrow(smooth.df))
smooth.plot <- melt(smooth.df, id.vars = "x")

PCAPlot("all_principal_components", smooth.plot, FALSE, 10, 15)
PCAPlot("facet_principal_components", smooth.plot, TRUE, 13, 25)
rownames(ME.genes) <- rownames(expr.collapse)

sample.ids <- factor(rownames(expr.collapse), levels = rownames(expr.collapse))
colnames(ME.genes) %<>% str_replace("ME", "")
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.collapse), module.colors)
cluster <- makeForkCluster(8)
split(expr.data.plot, expr.data.plot$module.colors) %>% map(Heatmap, ME.genes.plot)
#by(expr.data.plot, expr.data.plot$module.colors, Heatmap, ME.genes.plot)

#modules.out <- select(gene.info, Symbol, module.color)
#targets.final.known$Sample.Name %<>% str_replace(" ", "")

source("../../code/GO/enrichr.R")
source("../../FRDA project/common_functions.R")

modules.filter <- select(module.membership, Symbol, Module)
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up", "GTEx_Tissue_Sample_Gene_Expression_Profiles_down") 
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names, EnrichrSubmit, modules.filter, enrichr.terms, FALSE)

ME.genes.plot <- mutate(ME.genes, Combined = pheno.import$Combined) %>%
    gather(Module.Color, Eigengene, -Combined) 
ME.genes.plot$Module.Color %<>% str_replace("ME", "")

p <- ggplot(ME.genes.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + scale_fill_manual(values = sort(unique(ME.genes.plot$Module.Color)))
p <- p + facet_wrap(~ Module.Color, ncol = 4, scales = "free") + theme(plot.background = element_blank())
CairoPDF("eigengene_plots", height = 9, width = 16, bg = "transparent")
print(p)
dev.off()

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#Enrichr plots
FilterEnrichr <- function(enrichr.df, size = 100) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    if (nrow(enrichr.df) > size) {
        enrichr.filter %<>% slice(1:size)
    }
    enrichr.filter
}

symbols.cleaned <- na.omit(final.genetable$Symbol)
symbols.length <- map(symbols.cleaned, nchar) %>% reduce(c)
symbols.final <- symbols.cleaned[symbols.length > 0]

#Turquoise
turquoise.symbols <- filter(module.membership, Module == "turquoise")$Symbol
turquoise.gobiol.file <- "./enrichr/turquoise/turquoise_GO_Biological_Process_2015.xlsx"
turquoise.gobiol <- read.xlsx(turquoise.gobiol.file) 
turquoise.gobiol.filter <- FilterEnrichr(turquoise.gobiol)
GetKappaCluster(file_path_sans_ext(turquoise.gobiol.file), turquoise.gobiol.filter, turquoise.symbols)
turquoise.gobiol.final <- slice(turquoise.gobiol.filter, c(1, 2, 3))
turquoise.gobiol.final$Database <- "GO BP"

turquoise.gomole.file <- "./enrichr/turquoise/turquoise_GO_Molecular_Function_2015.xlsx"
turquoise.gomole <- read.xlsx(turquoise.gomole.file) 
turquoise.gomole.filter <- FilterEnrichr(turquoise.gomole)
GetKappaCluster(file_path_sans_ext(turquoise.gomole.file), turquoise.gomole.filter, turquoise.symbols)
turquoise.gomole.final <- slice(turquoise.gomole.filter, c(2))
turquoise.gomole.final$Database <- "GO MF"

turquoise.reactome.file <- "./enrichr/turquoise/turquoise_Reactome_2016.xlsx"
turquoise.reactome <- read.xlsx(turquoise.reactome.file) 
turquoise.reactome.filter <- FilterEnrichr(turquoise.reactome)
GetKappaCluster(file_path_sans_ext(turquoise.reactome.file), turquoise.reactome.filter, turquoise.symbols)
turquoise.reactome.final <- slice(turquoise.gomole.filter, c(3))
turquoise.reactome.final$Database <- "Reactome"

turquoise.enrichr.final <- rbind(turquoise.gobiol.final, turquoise.gomole.final, turquoise.reactome.final)
EnrichrPlot(turquoise.enrichr.final, "turquoise.enrichr", "", plot.width = 9)

#Brown module
brown.symbols <- filter(module.membership, Module == "brown")$Symbol
brown.gobiol.file <- "./enrichr/brown/brown_GO_Biological_Process_2015.xlsx"
brown.gobiol <- read.xlsx(brown.gobiol.file) 
brown.gobiol.filter <- FilterEnrichr(brown.gobiol, size = 200)
GetKappaCluster(file_path_sans_ext(brown.gobiol.file), brown.gobiol.filter, brown.symbols)
brown.gobiol.final <- slice(brown.gobiol.filter, c(1, 3))
brown.gobiol.final$Database <- "GO BP"

brown.gomole.file <- "./enrichr/brown/brown_GO_Molecular_Function_2015.xlsx"
brown.gomole <- read.xlsx(brown.gomole.file) 
brown.gomole.filter <- FilterEnrichr(brown.gomole)
GetKappaCluster(file_path_sans_ext(brown.gomole.file), brown.gomole.filter, brown.symbols)
brown.gomole.final <- slice(brown.gomole.filter, c(1, 3))
brown.gomole.final$Database <- "GO MF"

brown.reactome.file <- "./enrichr/brown/brown_Reactome_2016.xlsx"
brown.reactome <- read.xlsx(brown.reactome.file) 
brown.reactome.filter <- FilterEnrichr(brown.reactome)
GetKappaCluster(file_path_sans_ext(brown.reactome.file), brown.reactome.filter, brown.symbols)
brown.reactome.final <- slice(brown.reactome.filter, 1)
brown.reactome.final$Database <- "Reactome"

brown.kegg.file <- "./enrichr/brown/brown_KEGG_2016.xlsx"
brown.kegg <- read.xlsx(brown.kegg.file) 
brown.kegg.filter <- FilterEnrichr(brown.kegg)
GetKappaCluster(file_path_sans_ext(brown.kegg.file), brown.kegg.filter, brown.symbols)
brown.kegg.final <- slice(brown.kegg.filter, c(1, 3))
brown.kegg.final$Database <- "KEGG"

brown.enrichr.final <- rbind(brown.gobiol.final, brown.gomole.final, brown.reactome.final, brown.kegg.final)
EnrichrPlot(brown.enrichr.final, "brown.enrichr", plot.width = 7)

#Light cyan
pink.symbols <- filter(module.membership, Module == "pink")$Symbol
pink.gobiol.file <- "./enrichr/pink/pink_GO_Biological_Process_2015.xlsx"
pink.gobiol <- read.xlsx(pink.gobiol.file) 
pink.gobiol.filter <- FilterEnrichr(pink.gobiol)
GetKappaCluster(file_path_sans_ext(pink.gobiol.file), pink.gobiol.filter, pink.symbols)
pink.gobiol.final <- slice(pink.gobiol.filter, c(5))
pink.gobiol.final$Database <- "GO BP"

pink.gomole.file <- "./enrichr/pink/pink_GO_Molecular_Function_2015.xlsx"
pink.gomole <- read.xlsx(pink.gomole.file) 
pink.gomole.filter <- FilterEnrichr(pink.gomole)
GetKappaCluster(file_path_sans_ext(pink.gomole.file), pink.gomole.filter, pink.symbols)
pink.gomole.final <- slice(pink.gomole.filter, c(1))
pink.gomole.final$Database <- "GO MF"

pink.reactome.file <- "./enrichr/pink/pink_Reactome_2016.xlsx"
pink.reactome <- read.xlsx(pink.reactome.file) 
pink.reactome.filter <- FilterEnrichr(pink.reactome)
GetKappaCluster(file_path_sans_ext(pink.reactome.file), pink.reactome.filter, pink.symbols)
pink.reactome.final <- slice(pink.reactome.filter, c(3))
pink.reactome.final$Database <- "Reactome"

pink.enrichr.final <- rbind(pink.gobiol.final, pink.gomole.final, pink.reactome.final)
EnrichrPlot(pink.enrichr.final, "pink.enrichr", plot.height = 3, plot.width = 6)

#Magenta
magenta.symbols <- filter(module.membership, Module == "magenta")$Symbol
magenta.gobiol.file <- "./enrichr/magenta/magenta_GO_Biological_Process_2015.xlsx"
magenta.gobiol <- read.xlsx(magenta.gobiol.file) 
magenta.gobiol.filter <- FilterEnrichr(magenta.gobiol, size = 200)
GetKappaCluster(file_path_sans_ext(magenta.gobiol.file), magenta.gobiol.filter, magenta.symbols)
magenta.gobiol.final <- slice(magenta.gobiol.filter, c(1, 7, 10))
magenta.gobiol.final$Database <- "GO BP"

magenta.gomole.file <- "./enrichr/magenta/magenta_GO_Molecular_Function_2015.xlsx"
magenta.gomole <- read.xlsx(magenta.gomole.file) 
magenta.gomole.filter <- FilterEnrichr(magenta.gomole)
GetKappaCluster(file_path_sans_ext(magenta.gomole.file), magenta.gomole.filter, magenta.symbols)
magenta.gomole.final <- slice(magenta.gomole.filter, c(3, 9))
magenta.gomole.final$Database <- "GO MF"

magenta.reactome.file <- "./enrichr/magenta/magenta_Reactome_2016.xlsx"
magenta.reactome <- read.xlsx(magenta.reactome.file) 
magenta.reactome.filter <- FilterEnrichr(magenta.reactome)
GetKappaCluster(file_path_sans_ext(magenta.reactome.file), magenta.reactome.filter, magenta.symbols)
magenta.reactome.final <- slice(magenta.reactome.filter, c(3, 5))
magenta.reactome.final$Database <- "Reactome"

magenta.enrichr.final <- rbind(magenta.gobiol.final, magenta.gomole.final, magenta.reactome.final)
EnrichrPlot(magenta.enrichr.final, "magenta.enrichr", "")

#PPI
genetable.symbol <- filter(module.membership, !is.na(Symbol) & nchar(Symbol) > 0)
turquoise.only <- filter(genetable.symbol, Module == "turquoise")$Symbol
brown.only <- filter(genetable.symbol, Module == "brown")$Symbol
pink.only <- filter(genetable.symbol, Module == "pink")$Symbol
magenta.only <- filter(genetable.symbol, Module == "magenta")$Symbol

#PPI
turquoise.ppi <-  GetPPI(turquoise.only)
brown.ppi <-  GetPPI(brown.only)
pink.ppi <-  GetPPI(pink.only)
magenta.ppi <-  GetPPI(magenta.only)

turquoise.all.graph <- graph_from_edgelist(as.matrix(turquoise.ppi))
turquoise.incident <- map(V(turquoise.all.graph), incident, graph = turquoise.all.graph) %>% map_int(length) %>% sort
turquoise.ppi.df <- data.frame(Symbol = names(turquoise.incident), PPI = turquoise.incident)

brown.all.graph <- graph_from_edgelist(as.matrix(brown.ppi))
brown.incident <- map(V(brown.all.graph), incident, graph = brown.all.graph) %>% map_int(length) %>% sort
brown.ppi.df <- data.frame(Symbol = names(brown.incident), PPI = brown.incident)

pink.all.graph <- graph_from_edgelist(as.matrix(pink.ppi))
pink.incident <- map(V(pink.all.graph), incident, graph = pink.all.graph) %>% map_int(length) %>% sort
pink.ppi.df <- data.frame(Symbol = names(pink.incident), PPI = pink.incident)

magenta.all.graph <- graph_from_edgelist(as.matrix(magenta.ppi))
magenta.incident <- map(V(magenta.all.graph), incident, graph = magenta.all.graph) %>% map_int(length) %>% sort
magenta.ppi.df <- data.frame(Symbol = names(magenta.incident), PPI = magenta.incident)

#Overlap with DE data
voom.pco.up <- filter(voom.top, adj.P.Val.pco < 0.01 & logFC.pco > 1.5) %>% select(Symbol, adj.P.Val.pco, logFC.pco)
voom.pco.down <- filter(voom.top, adj.P.Val.pco < 0.01 & logFC.pco < -1.5) %>% select(Symbol, adj.P.Val.pco, logFC.pco)

turquoise.pco.up <- left_join(voom.pco.up, filter(module.membership, Module == "turquoise")) %>% filter(!is.na(Module)) %>% arrange(desc(kscaled))
brown.pco.down <- left_join(voom.pco.down, filter(module.membership, Module == "brown")) %>% filter(!is.na(Module)) %>% arrange(desc(kscaled))
pink.pco.up <- left_join(voom.pco.up, filter(module.membership, Module == "pink")) %>% filter(!is.na(Module)) %>% arrange(desc(kscaled))
magenta.pco.down <- left_join(voom.pco.down, filter(module.membership, Module == "magenta")) %>% filter(!is.na(Module)) %>% arrange(desc(kscaled))

turquoise.ppi.up <- left_join(voom.pco.up, turquoise.ppi.df) %>% filter(!is.na(PPI)) %>% arrange(desc(PPI))
brown.ppi.down <- left_join(voom.pco.down, brown.ppi.df) %>% filter(!is.na(PPI)) %>% arrange(desc(PPI))
pink.ppi.up <- left_join(voom.pco.up, pink.ppi.df) %>% filter(!is.na(PPI)) %>% arrange(desc(PPI))
magenta.ppi.down <- left_join(voom.pco.down, magenta.ppi.df) %>% filter(!is.na(PPI)) %>% arrange(desc(PPI))

Top5Plot <- function(rank.column, toptable.object, voom.object, pheno.object, pheno.col, levels.vector, maintitle, filename) {
    col.name <- str_c("desc(", rank.column, ")")

    top5.symbol <- arrange_(toptable.object, col.name)$Symbol[1:5]
    top5.expr <- t(voom.object[top5.symbol,])
    colnames(top5.expr) <- top5.symbol
    top5.df <- data.frame(Sample = pheno.object[[pheno.col]], top5.expr) %>%
        gather_("Gene", "Expression", names(.)[-1])
    top5.df$Gene %<>% factor(levels = make.names(top5.symbol))
    top5.df$Sample %<>% factor(levels = make.names(levels.vector))

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

batch2.levels <- c("normal.DMSO", "FRDA.DMSO", "normal.109", "FRDA.109")
Top5Plot("kscaled", turquoise.pco.up, voom.collapse, pheno.import, "Combined", batch2.levels, "FRDA vs. Control (upregulated, turquoise module)", "top5.turquoise.up")
Top5Plot("kscaled", brown.pco.down, voom.collapse, pheno.import, "Combined", batch2.levels, "FRDA vs. Control (downregulated, brown module)", "top5.brown.down")
Top5Plot("kscaled", pink.pco.up, voom.collapse, pheno.import, "Combined", batch2.levels, "FRDA vs. Control (upregulated, pink module)", "top5.pink.up")
Top5Plot("kscaled", magenta.pco.down, voom.collapse, pheno.import, "Combined", batch2.levels, "FRDA vs. Control (downregulated, magenta module)", "top5.magenta.down")

Top5Plot("PPI", turquoise.ppi.up, voom.collapse, pheno.import, "Combined", batch2.levels, "FRDA vs. Control (upregulated, turquoise module)", "top5.turquoise.ppi.up")
Top5Plot("PPI", brown.ppi.down, voom.collapse, pheno.import, "Combined", batch2.levels, "FRDA vs. Control (downregulated, brown module)", "top5.brown.ppi.down")
Top5Plot("PPI", magenta.ppi.down, voom.collapse, pheno.import, "Combined", batch2.levels, "FRDA vs. Control (downregulated, magenta module)", "top5.magenta.ppi.down")

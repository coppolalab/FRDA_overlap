#String operations
library(stringr)

library(openxlsx)

library(GEOquery)
library(lumi)
library(lumiMouseIDMapping)
library(lumiMouseAll.db)
library(annotate)
library(ggplot2)
library(Cairo)
library(heatmap.plus)
library(WGCNA)
library(limma)
library(readr)
library(gplots)
library(siggenes)

#Functional programming
library(magrittr)
library(purrr)

#Data arrangement
library(dplyr)
library(tidyr)

saveRDS.gz <- function(object, file, threads=parallel::detectCores()) 
{
    con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
    saveRDS(object, file = con)
    close(con)
}

readRDS.gz <- function(file, threads=parallel::detectCores()) 
{
    con <- pipe(paste0("pigz -d -c -p",threads," ",file))
    object <- readRDS(file = con)
    close(con)
    return(object)
}

#Boxplot function
gen.boxplot <- function(filename, lumi.object, maintext, ylabtext)
{
    #dataset %<>% t %>% data.frame
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object))
    dataset.m <- melt(dataset.addvars, id = "Sample.Status")
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value)) + geom_boxplot() + theme_bw()
    p <- p + scale_fill_manual(values = colorscheme)
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Oxygen", width = 20 , height = 8)
}

#Heatmap bar function
gen.heatmapbars <- function(genotype.colorscheme, tissue.colorscheme, targetset)
{
    genotype.heatmap <- data.frame("Genotype" = levels(factor(targetset$Genotype)), "Genotype.Color" = genotype.colorscheme)
    tissue.heatmap <- data.frame("Tissue" = levels(factor(targetset$Tissue)), "Tissue.Color" = tissue.colorscheme)
    colorscheme <- data.frame("Genotype" = targetset$Genotype, "Tissue" = targetset$Tissue) %>% join(genotype.heatmap) %>% join(tissue.heatmap)
    colorscheme <- as.matrix(select(colorscheme, Genotype.Color, Tissue.Color))
    return(colorscheme)
}

#Heatmap function
gen.heatmap <- function(filename, lumi.object, maintitle)
{
    intensities1.cor <- cor(exprs(lumi.object))
    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(intensities1.cor, col = heat.colors(40), ColSideColors = cbind(lumi.object$Genotype.Color, lumi.object$Tissue.Color), scale = "none", main = maintitle)
    dev.off()
}

#IAC detection of outliers vix 
gen.IACcluster <- function(filename, dataset, maintitle)
{
    IAC = corFast(dataset, use = "p")
    cluster1 = flashClust::hclust(as.dist(1 - IAC))
    CairoPDF(filename, width = 13, height = 10)
    plot(cluster1, main = paste(maintitle, " (no = ", dim(IAC)[2], ")"))
    dev.off()
    return(IAC)
}

#Create plot of standard deviations of all interarray correlations.  
gen.sdplot <- function(filename, dataset, maintitle)
{
    meanIAC <- apply(dataset, 2, mean)
    sdCorr <- sd(meanIAC)
    numbersd <- (meanIAC - mean(meanIAC)) / sdCorr
    numbersd.plot <- data.frame(Sample.Num = 1:ncol(dataset), Z.score = numbersd, Sample.Status = colnames(dataset))

    p <- ggplot(numbersd.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Status) )
    p <- p + geom_text(size = 4, colour = "red")
    p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
    CairoPDF(filename, width = 10, height = 10)
    print(p)
    dev.off()
    return(numbersd)
}

#Median absolute deviation standardization function
standardize <- function(dataset)
{
    rowmed <- apply(dataset, 1, median)
    rowmad <- apply(dataset, 1, mad)
    rv <- sweep(dataset, 1, rowmed)
    rv <- sweep(rv, 1, rowmad, "/")
    return(rv)
}

#Create heatmap of top genes
gen.topgenes <- function(filename, lumi.object, maintitle, rowmads, num.genes)
{
    top.genes.names <- rowmads[1:num.genes]
    dataset <- exprs(lumi.object)
    top.genes.intensities <- dataset[top.genes.names,]
    top.genes.dist <- dist(t(standardize(top.genes.intensities)))
    top.genes.clust <- flashClust::hclust(top.genes.dist)
    top.genes.matrix <- as.matrix(top.genes.dist)

    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(top.genes.matrix, col = rev(heat.colors(75)), distfun = function (x) as.dist(x), main = maintitle, scale = "none", ColSideColors = cbind(lumi.object$Genotype.Color, lumi.object$Tissue.Color))
    dev.off()
    return(top.genes.dist)
}

#MDS function - may need work
gen.pca <- function(filename, dataset, targetset, colorscheme, variablename)
{
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$GSE.ID, factor(targetset[[variablename]]))
    colnames(target.data) <- c("Sample.Status", variablename)
    colnames(dataset.plot) <- c("Sample.Status", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + scale_fill_manual(values = colorscheme) 
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename)
    CairoPDF(file = filename, height = 7, width = 7)
    print(p)
    dev.off()
}

#Run statistical cutoff tests
gen.decide <- function(test, fit.object, write.results = FALSE, suffix = "")
{
    results <- decideTests(fit.object, adjust.method = test[1], p = as.numeric(test[2])) #Run test at specified cutoff
    if(write.results == TRUE)
    {
        write.fit(file = paste("./fit_", test[1], suffix, ".tsv", sep = ""), fit.object, adjust = test[1], results = results)
    }
    num.genes <- length(which(apply(results, 1, function (x) any(x, na.rm = T))))  #Make this better
    mysum <- summary(results)[-2,] %>% data.frame#Eliminate the row for no change in expression
    mysum[2,] <- -(mysum[2,])
    colnames(mysum) <- "KIKO_vs_WT"
    mysum <- data.frame("Test" = paste(test[1], " p<", test[2], sep = ""), "Num" = paste(num.genes, "Genes", sep = " "), "Direction" = c("positive", "negative"), mysum)
    return(mysum)
}

#Plot statistical cutoff tests
gen.decideplot <- function(filename, decide.plot, width.plot = 6, height.plot = 7)
{
    decide.plot$variable <- str_replace_all(decide.plot$variable, "_", " ")
    p <- ggplot()
    p <- p + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = variable, y = value), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = variable, y = value, ymax = max(value) + 110, label = value), hjust = -0.3, position = position_dodge(width = 1))
    p <- p + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = variable, y = value), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = variable, y = value, ymax = min(value) - 110, label = abs(value)), hjust = 1.3, position = position_dodge(width = 1))
    if (length(unique(decide.plot$Test)) > 1)
    {
        p <- p + facet_grid(Test ~ .) 
        #p <- p + ggtitle("Threshold Selection")
    }
    #else
    #{
        #p <- p + ggtitle(paste(decide.plot$Test, "\n", decide.plot$Num))
    #}
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0))# + ylab("Differentially Expressed Genes")
    CairoPDF(filename, width = width.plot, height = height.plot)
    print(p)
    dev.off()
}

gen.pval.hist <- function(filename, fit.pvals)
{
    colnames(fit.pvals) <- c("Tg.DOX_vs_Tg.ND", "Tg.DOX_vs_WT.DOX", "Tg.DOXR_vs_Tg.DOX")
    fit.pvals.plot <- melt(fit.pvals)
    fit.pvals.plot$Contrasts <- str_replace_all(fit.pvals.plot$Contrasts, "_", " ")
    p <- ggplot(fit.pvals.plot, aes(x = value)) + geom_histogram(binwidth = 1/80) + facet_grid(. ~ Contrasts)
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("P-value distribution across contrasts") + theme(axis.title.y = element_blank()) + xlab("p-value")
    CairoPDF(filename, height = 7, width = 21)
    print(p)
    dev.off()
}

gen.venndiagram <- function(filename, results)
{
    sumV <- colSums(summary(results)[-2,])
    v <- paste(c("Tg.DOX vs Tg.ND", "Tg.DOX vs WT.DOX", "Tg.DOXR vs Tg.DOX"), " (", sumV, ")", sep = "")
    CairoPDF(filename, width = 6, height = 6)
    vennDiagram(results, names = v, main = "", include = c("up", "down"), counts.col=c(2,3), cex = 0.8)
    dev.off()
}

#Calculate ratios.  
gen.ratios <- function(lumi.final)
{
    all.kiko <- exprs(lumi.final[,lumi.final$Genotype == "KIKO"])
    all.wt <- exprs(lumi.final[,lumi.final$Genotype == "WT"])
    all.wt.means <- rowMeans(all.wt)

    coef.kiko.wt <- all.kiko - all.wt.means
    colnames(coef.kiko.wt) <- paste(colnames(coef.kiko.wt), ".kikowt", sep = "")

    coefs <- data.frame("Symbol" = rownames(coef.kiko.wt), coef.kiko.wt)
    all.samples <- data.frame("Symbol" = rownames(lumi.final), exprs(lumi.final))
    colnames(all.samples)[2:ncol(all.samples)] <- paste(colnames(all.samples[2:ncol(all.samples)]), "expr", sep = ".") 
    ratio.exp <- merge(coefs, all.samples)
    return(ratio.exp)
}

gen.fit <- function(dataset, model.design)
{
    fit <- lmFit(dataset, design = model.design)
    contrasts.anova <- makeContrasts(KIKO_vs_WT = KIKO - WT, levels = model.design)
    fit2.anova <- contrasts.fit(fit, contrasts.anova)
    fitb <- eBayes(fit2.anova)
}

gen.workbook <- function(dataset, filename)
{
    pval.cols <- colnames(dataset) %>% str_detect("p.value") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("Coef") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:7, widths = "auto")
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Create genelists
gen.tables <- function(dataset, lumi.object, ratio.exp, suffix)
{
    treat.de <- data.frame("Symbol" = rownames(dataset), dataset)
    
    fitsel.ratio.all <- merge(treat.de, ratio.exp)

    fitsel.return.all <- select(fitsel.ratio.all, Symbol, contains("Coef"), contains("p.value"), contains("_vs_"), matches("^t$"), A, contains("GSM")) %>% arrange(desc(t))

    anovalist <- select(fitsel.return.all, contains("KIKO_vs_WT")) %>% apply(1, any, na.rm = T) %>% which
    fitsel.return <- fitsel.return.all[anovalist,]
    coef.cols <- colnames(fitsel.return) %>% str_detect("Coef") %>% which
    gen.workbook(fitsel.return, paste("./significant_geneList_", suffix, ".xlsx", sep = ""))

    write.csv(fitsel.return.all, paste("./complete_genelist_time1_", suffix, ".csv", sep = ""), row.names = FALSE)
    return(fitsel.return)
}

#Make anova objects
gen.anova <- function(dataset, suffix)
{
    plot.kikowt <- filter(dataset, KIKO_vs_WT != 0) %>% select(matches("kikowt"))
    gen.anova.heatmap(paste("./anova_heatmap_kiko_vs_wt", suffix, sep = "_"), plot.kikowt, "KIKO_vs_WT") 
    return(plot.kikowt)
}

#Generate anova heatmaps
gen.anova.heatmap <- function(filename, dataset, maintitle)
{ 
    CairoPDF(filename, width = 8, height = 8)
    heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 0.3, labRow = "", labCol = "", keysize = 0.9)
    dev.off()
}

#3 tissues, wt and KIKO
GSE15848.raw <- getGEO("GSE15848", destdir = ".")[[1]]
GSE15848.pdata <- pData(GSE15848.raw) %>% select(title, characteristics_ch1, characteristics_ch1.1)

kiko.pdata <- select(pData(GSE15848.raw), geo_accession, characteristics_ch1, characteristics_ch1.1, description)
colnames(kiko.pdata) <- c("GSE.ID", "Genotype", "Tissue", "Slide.ID")
kiko.pdata$Genotype %<>% str_replace("genotype: ", "")
kiko.pdata$Tissue %<>% str_replace("tissue: ", "")
kiko.supplemental <- getGEOSuppFiles('GSE15848', baseDir = '.')

kiko.raw <- lumiR("./GSE15848/GSE15848_raw_data.txt", lib.mapping = 'lumiMouseIDMapping')
pData(kiko.raw) <- kiko.pdata

kiko.norm <- lumiT(kiko.raw, method = "log2")

kiko.heart <- kiko.norm[,kiko.norm$Tissue == "Heart"]
kiko.heart.norm <- lumiN(kiko.heart, method = "rsn") #Normalize with robust spline regression
kiko.heart.cutoff <- detectionCall(kiko.heart.norm) #Get the count of probes which passed the detection threshold per sample
kiko.heart.expr <- kiko.heart.norm[which(kiko.heart.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.kiko.heart <- getSYMBOL(rownames(kiko.heart.expr), 'lumiMouseAll.db') %>% is.na #Determine which remaining probes are unannotated
kiko.heart.annot <- kiko.heart.expr[!symbols.kiko.heart,] #Drop any probe which is not annotated
saveRDS.gz(kiko.heart.annot, file = "kiko.heart.annot.rda")

gen.boxplot("heart_intensity_norm.jpg", kiko.heart.annot, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
genotype.colors <- c("red", "blue")
tissue.colors <- "black"
heatmap.bars <- gen.heatmapbars(genotype.colors, tissue.colors, pData(kiko.heart.annot)) %>% data.frame
kiko.heart.annot$Genotype.Color <- as.character(heatmap.bars$Genotype.Color)
kiko.heart.annot$Tissue.Color <- as.character(heatmap.bars$Tissue.Color)
gen.heatmap("heart_heatmap", kiko.heart.annot, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.heart <- gen.IACcluster("heart_IAC", exprs(kiko.heart.annot), "Outlier removed")
sd.raw.norm.heart <- gen.sdplot("heart_sd", IAC.norm.heart, "Outlier removed")
saveRDS.gz(IAC.norm.heart, file = "./save/IAC.norm.heart.rda")
saveRDS.gz(sd.raw.norm.heart, file = "./save/sd.raw.norm.heart.rda")

kiko.heart.mads <- apply(exprs(kiko.heart.annot), 1, mad)
kiko.heart.ordered <- order(kiko.heart.mads, decreasing = TRUE)
top500.dist.norm.heart <- gen.topgenes("heart_heatmap_500", kiko.heart.annot, "Clustering Based on the Top 500 Most Variable Genes", kiko.heart.ordered, 500)
top1000.dist.norm.heart <- gen.topgenes("heart_heatmap_1000", kiko.heart.annot, "Clustering Based on the Top 1000 Most Variable Genes", kiko.heart.ordered, 1000)

cm1.heart <- cmdscale(top1000.dist.norm.heart, eig = TRUE)
gen.pca("heart_mds_genotype", cm1.heart, phenoData(kiko.heart.annot), genotype.colors, "Genotype")

#Batch correction skipped
#No covariates to remove

fdata.heart <- fData(kiko.heart.annot) #Get featureData from lumi object
fdata.heart$SYMBOL <- getSYMBOL(rownames(fdata), "lumiMouseAll.db")
heart.expr.collapse <- collapseRows(exprs(kiko.heart.annot), factor(fdata$SYMBOL), rownames(kiko.heart.annot))$datETcollapsed #collapseRows by symbol
colnames(heart.expr.collapse) <- sampleNames(kiko.heart.annot)
kiko.heart.collapse <- ExpressionSet(assayData = heart.expr.collapse, phenoData = phenoData(kiko.heart.annot))
saveRDS.gz(kiko.heart.collapse, file = "./save/kiko.heart.collapse.rda")

heart.test <- GSE15848.raw[,GSE15848.raw$characteristics_ch1.1 == "tissue: Heart"]
heart.test$Genotype <- str_replace_all(heart.test$characteristics_ch1, "genotype: ", "")

model.heart <- model.matrix( ~ 0 + kiko.heart.collapse$Genotype)
colnames(model.heart) <- c("KIKO", "WT")

fit.object.heart <- gen.fit(exprs(kiko.heart.collapse), model.heart)
saveRDS.gz(fit.object.heart, file = "./save/fit.object.heart.rda")
top.object.heart <- topTable(fit.object.heart, coef = 1, n = Inf)
saveRDS.gz(top.object.heart, file = "./save/top.object.heart.rda")

ebam.heart <- limma2ebam(fit.object.heart, coef = 1)
ebam.heart.df <- data.frame(Z.score = ebam.heart@z, Posterior = ebam.heart@posterior)
ebam.heart.df$Significant <- ebam.heart.df$Posterior > 0.9
ebam.heart.df$Symbol <- rownames(ebam.heart.df)
saveRDS.gz(ebam.heart.df, "ebam.heart.df.rda")

ratio.exp.heart <- gen.ratios(kiko.heart.collapse)
saveRDS.gz(ratio.exp.heart, file = "./save/ratio.exp.heart.rda")

decide <- list(c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot.heart <- map_df(decide, gen.decide, fit.object.heart) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection_heart", decide.plot.heart, 10, 10) #Plot different significance cutoffs

decide.final <- gen.decide(c("none", 0.005), fit.object.heart, TRUE, "_heart") %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute signifance cutoff p < 0.001, no FDR adjustment
gen.decideplot("./selected_threshold", decide.final) #Plot cutoff

de.object.heart <- read_tsv("./fit_none_heart.tsv") #Read in unadjusted fit object
rownames(de.object.heart) <- featureNames(kiko.heart.collapse)

fit.selection.heart <- gen.tables(de.object.heart, kiko.heart.collapse, ratio.exp.heart, "pLess005_heart") #create differential expression table for unadjusted fit
saveRDS.gz(fit.selection.heart, file = "./save/fit.selection.heart.rda")

clust.none <- gen.anova(fit.selection.heart, "none_heart")

kiko.muscle <- kiko.norm[,kiko.norm$Tissue == "Skeletal Muscle"]
kiko.muscle.norm <- lumiN(kiko.muscle, method = "rsn") #Normalize with robust spline regression
kiko.muscle.cutoff <- detectionCall(kiko.muscle.norm) #Get the count of probes which passed the detection threshold per sample
kiko.muscle.expr <- kiko.muscle.norm[which(kiko.muscle.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.kiko.muscle <- getSYMBOL(rownames(kiko.muscle.expr), 'lumiMouseAll.db') %>% is.na #Determine which remaining probes are unannotated
kiko.muscle.annot <- kiko.muscle.expr[!symbols.kiko.muscle,] #Drop any probe which is not annotated
saveRDS.gz(kiko.muscle.annot, file = "kiko.muscle.annot.rda")

gen.boxplot("muscle_intensity_norm.jpg", kiko.muscle.annot, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
genotype.colors <- c("red", "blue")
tissue.colors <- "black"
heatmap.bars <- gen.heatmapbars(genotype.colors, tissue.colors, pData(kiko.muscle.annot)) %>% data.frame
kiko.muscle.annot$Genotype.Color <- as.character(heatmap.bars$Genotype.Color)
kiko.muscle.annot$Tissue.Color <- as.character(heatmap.bars$Tissue.Color)
gen.heatmap("muscle_heatmap", kiko.muscle.annot, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.muscle <- gen.IACcluster("muscle_IAC", exprs(kiko.muscle.annot), "Outlier removed")
sd.raw.norm.muscle <- gen.sdplot("muscle_sd", IAC.norm.muscle, "Outlier removed")
saveRDS.gz(IAC.norm.muscle, file = "./save/IAC.norm.muscle.rda")
saveRDS.gz(sd.raw.norm.muscle, file = "./save/sd.raw.norm.muscle.rda")

kiko.muscle.mads <- apply(exprs(kiko.muscle.annot), 1, mad)
kiko.muscle.ordered <- order(kiko.muscle.mads, decreasing = TRUE)
top500.dist.norm.muscle <- gen.topgenes("muscle_heatmap_500", kiko.muscle.annot, "Clustering Based on the Top 500 Most Variable Genes", kiko.muscle.ordered, 500)
top1000.dist.norm.muscle <- gen.topgenes("muscle_heatmap_1000", kiko.muscle.annot, "Clustering Based on the Top 1000 Most Variable Genes", kiko.muscle.ordered, 1000)

cm1.muscle <- cmdscale(top1000.dist.norm.muscle, eig = TRUE)
gen.pca("muscle_mds_genotype", cm1.muscle, phenoData(kiko.muscle.annot), genotype.colors, "Genotype")

#Batch correction skipped
#No covariates to remove

fdata.muscle <- fData(kiko.muscle.annot) #Get featureData from lumi object
fdata.muscle$SYMBOL <- getSYMBOL(rownames(fdata.muscle), "lumiMouseAll.db")
muscle.expr.collapse <- collapseRows(exprs(kiko.muscle.annot), factor(fdata.muscle$SYMBOL), rownames(kiko.muscle.annot))$datETcollapsed #collapseRows by symbol
colnames(muscle.expr.collapse) <- sampleNames(kiko.muscle.annot)
kiko.muscle.collapse <- ExpressionSet(assayData = muscle.expr.collapse, phenoData = phenoData(kiko.muscle.annot))
saveRDS.gz(kiko.muscle.collapse, file = "./save/kiko.muscle.collapse.rda")

model.muscle <- model.matrix( ~ 0 + kiko.muscle.collapse$Genotype)
colnames(model.muscle) <- c("KIKO", "WT")

fit.object.muscle <- gen.fit(exprs(kiko.muscle.collapse), model.muscle)
saveRDS.gz(fit.object.muscle, file = "./save/fit.object.muscle.rda")
top.object.muscle <- topTable(fit.object.muscle, coef = 1, n = Inf)
saveRDS.gz(top.object.muscle, file = "./save/top.object.muscle.rda")

ebam.muscle <- limma2ebam(fit.object.muscle, coef = 1)
ebam.muscle.df <- data.frame(Z.score = ebam.muscle@z, Posterior = ebam.muscle@posterior)
ebam.muscle.df$Significant <- ebam.muscle.df$Posterior > 0.9
ebam.muscle.df$Symbol <- rownames(ebam.muscle.df)
saveRDS.gz(ebam.muscle.df, "ebam.muscle.df.rda")

ratio.exp.muscle <- gen.ratios(kiko.muscle.collapse)
saveRDS.gz(ratio.exp.muscle, file = "./save/ratio.exp.muscle.rda")

decide <- list(c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot.muscle <- map_df(decide, gen.decide, fit.object.muscle) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection_muscle", decide.plot.muscle, 10, 10) #Plot different significance cutoffs

decide.final <- gen.decide(c("none", 0.005), fit.object.muscle, TRUE, "_muscle") %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute signifance cutoff p < 0.001, no FDR adjustment
gen.decideplot("./selected_threshold", decide.final) #Plot cutoff

de.object.muscle <- read_tsv("./fit_none_muscle.tsv") #Read in unadjusted fit object
rownames(de.object.muscle) <- featureNames(kiko.muscle.collapse)

fit.selection.muscle <- gen.tables(de.object.muscle, kiko.muscle.collapse, ratio.exp.muscle, "pLess005_muscle") #create differential expression table for unadjusted fit
saveRDS.gz(fit.selection.muscle, file = "./save/fit.selection.muscle.rda")

clust.none <- gen.anova(fit.selection.muscle, "none_muscle")

kiko.liver <- kiko.norm[,kiko.norm$Tissue == "liver"]
kiko.liver.norm <- lumiN(kiko.liver, method = "rsn") #Normalize with robust spline regression
kiko.liver.cutoff <- detectionCall(kiko.liver.norm) #Get the count of probes which passed the detection threshold per sample
kiko.liver.expr <- kiko.liver.norm[which(kiko.liver.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.kiko.liver <- getSYMBOL(rownames(kiko.liver.expr), 'lumiMouseAll.db') %>% is.na #Determine which remaining probes are unannotated
kiko.liver.annot <- kiko.liver.expr[!symbols.kiko.liver,] #Drop any probe which is not annotated
saveRDS.gz(kiko.liver.annot, file = "kiko.liver.annot.rda")

gen.boxplot("liver_intensity_norm.jpg", kiko.liver.annot, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
genotype.colors <- c("red", "blue")
tissue.colors <- "black"
heatmap.bars <- gen.heatmapbars(genotype.colors, tissue.colors, pData(kiko.liver.annot)) %>% data.frame
kiko.liver.annot$Genotype.Color <- as.character(heatmap.bars$Genotype.Color)
kiko.liver.annot$Tissue.Color <- as.character(heatmap.bars$Tissue.Color)
gen.heatmap("liver_heatmap", kiko.liver.annot, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.liver <- gen.IACcluster("liver_IAC", exprs(kiko.liver.annot), "Outlier removed")
sd.raw.norm.liver <- gen.sdplot("liver_sd", IAC.norm.liver, "Outlier removed")
saveRDS.gz(IAC.norm.liver, file = "./save/IAC.norm.liver.rda")
saveRDS.gz(sd.raw.norm.liver, file = "./save/sd.raw.norm.liver.rda")

kiko.liver.mads <- apply(exprs(kiko.liver.annot), 1, mad)
kiko.liver.ordered <- order(kiko.liver.mads, decreasing = TRUE)
top500.dist.norm.liver <- gen.topgenes("liver_heatmap_500", kiko.liver.annot, "Clustering Based on the Top 500 Most Variable Genes", kiko.liver.ordered, 500)
top1000.dist.norm.liver <- gen.topgenes("liver_heatmap_1000", kiko.liver.annot, "Clustering Based on the Top 1000 Most Variable Genes", kiko.liver.ordered, 1000)

cm1.liver <- cmdscale(top1000.dist.norm.liver, eig = TRUE)
gen.pca("liver_mds_genotype", cm1.liver, phenoData(kiko.liver.annot), genotype.colors, "Genotype")

#Batch correction skipped
#No covariates to remove

fdata.liver <- fData(kiko.liver.annot) #Get featureData from lumi object
fdata.liver$SYMBOL <- getSYMBOL(rownames(fdata.liver), "lumiMouseAll.db")
liver.expr.collapse <- collapseRows(exprs(kiko.liver.annot), factor(fdata.liver$SYMBOL), rownames(kiko.liver.annot))$datETcollapsed #collapseRows by symbol
colnames(liver.expr.collapse) <- sampleNames(kiko.liver.annot)
kiko.liver.collapse <- ExpressionSet(assayData = liver.expr.collapse, phenoData = phenoData(kiko.liver.annot))
saveRDS.gz(kiko.liver.collapse, file = "./save/kiko.liver.collapse.rda")

model.liver <- model.matrix( ~ 0 + kiko.liver.collapse$Genotype)
colnames(model.liver) <- c("KIKO", "WT")

fit.object.liver <- gen.fit(exprs(kiko.liver.collapse), model.liver)
saveRDS.gz(fit.object.liver, file = "./save/fit.object.liver.rda")
top.object.liver <- topTable(fit.object.liver, coef = 1, n = Inf)
saveRDS.gz(top.object.liver, file = "./save/top.object.liver.rda")

ebam.liver <- limma2ebam(fit.object.liver, coef = 1)
ebam.liver.df <- data.frame(Z.score = ebam.liver@z, Posterior = ebam.liver@posterior)
ebam.liver.df$Significant <- ebam.liver.df$Posterior > 0.9
ebam.liver.df$Symbol <- rownames(ebam.liver.df)
saveRDS.gz(ebam.liver.df, "ebam.liver.df.rda")

ratio.exp.liver <- gen.ratios(kiko.liver.collapse)
saveRDS.gz(ratio.exp.liver, file = "./save/ratio.exp.liver.rda")

decide <- list(c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot.liver <- map_df(decide, gen.decide, fit.object.liver) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection_liver", decide.plot.liver, 10, 10) #Plot different significance cutoffs

decide.final <- gen.decide(c("none", 0.005), fit.object.liver, TRUE, "_liver") %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute signifance cutoff p < 0.001, no FDR adjustment
gen.decideplot("./selected_threshold", decide.final) #Plot cutoff

de.object.liver <- read_tsv("./fit_none_liver.tsv") #Read in unadjusted fit object
rownames(de.object.liver) <- featureNames(kiko.liver.collapse)

fit.selection.liver <- gen.tables(de.object.liver, kiko.liver.collapse, ratio.exp.liver, "pLess005_liver") #create differential expression table for unadjusted fit
saveRDS.gz(fit.selection.liver, file = "./save/fit.selection.liver.rda")

clust.none <- gen.anova(fit.selection.liver, "none_liver")

decide.plot.heart$variable <- "Heart"
decide.plot.muscle$variable <- "Muscle"
decide.plot.liver$variable <- "Liver"
plot.all <- rbind(decide.plot.heart, decide.plot.muscle, decide.plot.liver)
gen.decideplo("./threshold_selection", plot.all)

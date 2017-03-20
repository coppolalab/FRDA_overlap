#String operations
library(stringr)

library(openxlsx)
library(readr)

library(GEOquery)
library(lumi)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(annotate)
library(ggplot2)
library(Cairo)
library(heatmap.plus)
library(WGCNA)
library(limma)
library(gplots)
library(siggenes)

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
library(doBy)

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

#Boxplot function - needs no change
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

#Heatmap bar function - must be updated 
gen.heatmapbars <- function(status.colorscheme, targetset)
{
    status.heatmap <- data.frame("Status" = levels(factor(targetset$Status)), "Status.Color" = status.colorscheme)
    colorscheme <- data.frame("Status" = targetset$Status) %>% join(status.heatmap)
    colorscheme <- as.matrix(select(colorscheme, Status.Color))
    return(colorscheme)
}

#Heatmap function - needs no change
gen.heatmap <- function(filename, lumi.object, maintitle)
{
    intensities1.cor <- cor(exprs(lumi.object))
    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(intensities1.cor, col = heat.colors(40), ColSideColors = as.matrix(data.frame(lumi.object$Status.Color, "white")), scale = "none", main = maintitle)
    dev.off()
}

#IAC detection of outliers - needs no change 
gen.IACcluster <- function(filename, dataset, maintitle)
{
    IAC = corFast(dataset, use = "p")
    cluster1 = flashClust::hclust(as.dist(1 - IAC))
    CairoPDF(filename, width = 13, height = 10)
    plot(cluster1, main = paste(maintitle, " (no = ", dim(IAC)[2], ")"))
    dev.off()
    return(IAC)
}

#Create plot of standard deviations of all interarray correlations - needs no change
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

#Median absolute deviation standardization function - needs no change
standardize <- function(dataset)
{
    rowmed <- apply(dataset, 1, median)
    rowmad <- apply(dataset, 1, mad)
    rv <- sweep(dataset, 1, rowmed)
    rv <- sweep(rv, 1, rowmad, "/")
    return(rv)
}

#Create heatmap of top genes - needs no change
gen.topgenes <- function(filename, lumi.object, maintitle, rowmads, num.genes)
{
    top.genes.names <- rowmads[1:num.genes]
    dataset <- exprs(lumi.object)
    top.genes.intensities <- dataset[top.genes.names,]
    top.genes.dist <- dist(t(standardize(top.genes.intensities)))
    top.genes.clust <- flashClust::hclust(top.genes.dist)
    top.genes.matrix <- as.matrix(top.genes.dist)

    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(top.genes.matrix, col = rev(heat.colors(75)), distfun = function (x) as.dist(x), main = maintitle, scale = "none", ColSideColors = as.matrix(data.frame(lumi.object$Status.Color, "white")))
    dev.off()
    return(top.genes.dist)
}

#MDS function - may need work - needs no change
gen.pca <- function(filename, dataset, targetset, colorscheme, variablename)
{
    dataset.plot <- data.frame(Sample.Name = rownames(dataset$points), dataset$points)
    print(dataset.plot)
    target.data <- data.frame(targetset$Sample.Name, factor(targetset[[variablename]]))
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

#Run statistical cutoff tests - needs no change
gen.decide <- function(test, fit.object, write.results = FALSE)
{
    results <- decideTests(fit.object, adjust.method = test[1], p = as.numeric(test[2])) #Run test at specified cutoff
    if(write.results == TRUE)
    {
        write.fit(file = paste("./fit_", test[1], ".tsv", sep = ""), fit.object, adjust = test[1], results = results)
    }
    num.genes <- length(which(apply(results, 1, function (x) any(x, na.rm = T))))  #Make this better
    mysum <- summary(results)[-2,] %>% data.frame#Eliminate the row for no change in expression
    mysum[2,] <- -(mysum[2,])
    colnames(mysum) <- c("FRDA_vs_Carrier", "FRDA_vs_Normal", "Carrier_vs_Normal")
    mysum <- data.frame("Test" = paste(test[1], " p<", test[2], sep = ""), "Num" = paste(num.genes, "Genes", sep = " "), "Direction" = c("positive", "negative"), mysum)
    return(mysum)
}

#Plot statistical cutoff tests - needs no change
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
        p <- p + facet_grid(Num + Test ~ .) 
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

#Calculate ratios - needs many changes 
gen.ratios <- function(lumi.final)
{
    all.patient <- exprs(lumi.final[,lumi.final$Status == "FRDA"])
    all.carrier <- exprs(lumi.final[,lumi.final$Status == "Carrier"])
    all.control <- exprs(lumi.final[,lumi.final$Status == "Normal"])
    carrier.means <- rowMeans(all.carrier)
    control.means <- rowMeans(all.control)

    coef.pca <- all.patient - carrier.means
    coef.pco <- all.patient - control.means
    coef.cc <- all.carrier - control.means
    colnames(coef.pca) <- paste(colnames(coef.pca), ".pca", sep = "")
    colnames(coef.pco) <- paste(colnames(coef.pco), ".pco", sep = "")
    colnames(coef.cc) <- paste(colnames(coef.cc), ".cc", sep = "")

    coefs <- data.frame("Symbol" = rownames(coef.pca), coef.pca, coef.pco, coef.cc)
    all.samples <- data.frame("Symbol" = rownames(lumi.final), exprs(lumi.final))
    colnames(all.samples)[2:ncol(all.samples)] <- paste(colnames(all.samples[2:ncol(all.samples)]), "expr", sep = ".") 
    ratio.exp <- merge(coefs, all.samples)
    return(ratio.exp)
}

#Needs new levels
gen.fit <- function(dataset, model.design) 
{
    fit <- lmFit(dataset, design = model.design)
    contrasts.anova <- makeContrasts(FRDA.vs.Healthy = FRDA - Healthy, levels = model.design)
    fit2.anova <- contrasts.fit(fit, contrasts.anova)
    fitb <- eBayes(fit2.anova)

    top.object.frdah <- topTable(fitb, coef = 1, n = Inf)
    saveRDS.gz(top.object.frdah, "./save/top.object.frdah.rda")
    return(fitb)
}

#doesn't need update
gen.workbook <- function(dataset, filename)
{
    pval.cols <- colnames(dataset) %>% str_detect("p.value.") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")
    #dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

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

#needs updated colnames
gen.tables <- function(dataset, lumi.object, ratio.exp, suffix)
{
    treat.de <- data.frame("Symbol" = rownames(dataset), dataset)
    
    fitsel.ratio.all <- merge(treat.de, ratio.exp)
    fitsel.return.all <- select(fitsel.ratio.all, Symbol, contains("Coef."), contains("p.value."), contains("Res."), contains("t."), F, F.p.value, A, matches("pca|pco|cc|expr")) %>% arrange(desc(F))

    anovalist <- select(fitsel.return.all, contains("Res.")) %>% apply(1, any, na.rm = T) %>% which
    fitsel.return <- fitsel.return.all[anovalist,]
    coef.cols <- colnames(fitsel.return) %>% str_detect("Coef.") %>% which
    colnames(fitsel.return)[coef.cols] <- c("Coef.FRDA vs. Carrier", "Coef.FRDA vs. Normal", "Coef.Carrier vs. Normal")
    gen.workbook(fitsel.return, paste("./significant_geneList_", suffix, ".xlsx", sep = ""))

    write.csv(fitsel.return.all, paste("./complete_genelist_", suffix, ".csv", sep = ""), row.names = FALSE)
    return(fitsel.return)
}

#no need to update
gen.small.workbook <- function(dataset, filename)
{
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:4, widths = "auto")
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}
    
#Make anova objects - must be updated
gen.anova <- function(dataset, suffix)
{
    plot.pca <- filter(dataset, Res.FRDA.vs.Carrier != 0) %>% select(matches("pca"))
    plot.pco <- filter(dataset, Res.FRDA.vs.Normal != 0) %>% select(matches("pco"))
    plot.cc <- filter(dataset, Res.Carrier.vs.Normal != 0) %>% select(matches("cc"))
    gen.anova.heatmap(paste("./anova_heatmap_frda_vs_carrier", suffix, sep = "_"), plot.pca, "FRDA vs. Carrier") 
    gen.anova.heatmap(paste("./anova_heatmap_frda_vs_normal", suffix, sep = "_"), plot.pco, "FRDA vs. Normal")
    gen.anova.heatmap(paste("./anova_heatmap_carrier_vs_normal", suffix, sep = "_"), plot.cc, "Carrier vs. Normal")
    return(list(plot.pca, plot.pco, plot.cc))
}

#Generate anova heatmap - no need to update 
gen.anova.heatmap <- function(filename, dataset, maintitle)
{ 
    CairoPDF(filename, width = 8, height = 8)
    heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 0.3, labRow = "", labCol = "", keysize = 0.9)
    dev.off()
}

#3 statuses, 5 treatments
GSE65399.raw <- getGEO("GSE65399", destdir = ".")[[1]]
GSE65399.pdata <- pData(GSE65399.raw) %>% select(title)

ipsc.supplemental <- getGEOSuppFiles('GSE65399', baseDir = '.')
ipsc.raw <- lumiR('./GSE65399/GSE65399_HT12v4_RAW.txt', lib.mapping = 'lumiHumanIDMapping')

ipsc.diff <- ipsc.raw[,str_detect(sampleNames(ipsc.raw), "neuronal.diff")]
ipsc.annotation <- read.xlsx("./GSE65399_fixed.xlsx") %>% filter(str_detect(title, "neuronal.diff")) %>% select(title, X2)
ipsc.annotation$title %<>% str_replace_all(" ", "_")

ipsc.pdata <- merge(pData(ipsc.diff), ipsc.annotation, by.x = "sampleID", by.y = "title")
ipsc.pdata$Status <- map(ipsc.pdata$X2, str_split, " ") %>% map(unlist) %>% map_chr(`[`, 1)
ipsc.pdata$sampleID %<>% str_replace("neuronal.diff_", "")
rownames(ipsc.pdata) <- ipsc.pdata$sampleID

pData(ipsc.diff) <- ipsc.pdata

ipsc.log2 <- lumiT(ipsc.diff, method = "log2")
ipsc.norm <- lumiN(ipsc.log2, method = "rsn") #Normalize with robust spline regression
ipsc.cutoff <- detectionCall(ipsc.norm) #Get the count of probes which passed the detection threshold per sample
ipsc.expr <- ipsc.norm[which(ipsc.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.ipsc <- getSYMBOL(rownames(ipsc.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
ipsc.expr.annot <- ipsc.expr[!symbols.ipsc,] #Drop any probe which is not annotated
saveRDS.gz(ipsc.expr.annot, file = "./save/ipsc.expr.annot.rda")

gen.boxplot("intensity_norm.jpg", ipsc.expr.annot, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
status.colors <- c("red", "blue")
heatmap.bars <- gen.heatmapbars(status.colors, pData(ipsc.expr.annot)) %>% data.frame
ipsc.expr.annot$Status.Color <- as.character(heatmap.bars$Status.Color)
gen.heatmap("heatmap", ipsc.expr.annot, "Clustering Based on Inter-Array Pearson Coefficient, RSN normalized")
IAC.norm <- gen.IACcluster("ipsc_IAC", exprs(ipsc.expr.annot), "")
sd.raw.norm <- gen.sdplot("ipsc_sd", IAC.norm, "")
saveRDS.gz(IAC.norm, file = "./save/IAC.norm.rda")
saveRDS.gz(sd.raw.norm, file = "./save/sd.raw.norm.rda")

ipsc.mads <- apply(exprs(ipsc.expr.annot), 1, mad)
ipsc.ordered <- order(ipsc.mads, decreasing = TRUE)
top500.dist.norm <- gen.topgenes("ipsc_heatmap_500", ipsc.expr.annot, "Clustering Based on the Top 500 Most Variable Genes", ipsc.ordered, 500)
top1000.dist.norm <- gen.topgenes("ipsc_heatmap_1000", ipsc.expr.annot, "Clustering Based on the Top 1000 Most Variable Genes", ipsc.ordered, 1000)

cm1 <- cmdscale(top1000.dist.norm, eig = TRUE)
gen.pca("mds_Status", cm1, pData(ipsc.expr.annot), genotype.colors, "Status")

fdata.heart <- fData(ipsc.expr.annot) #Get featureData from lumi object
fdata.heart$SYMBOL <- getSYMBOL(rownames(fdata.heart), "lumiHumanAll.db")
heart.expr.collapse <- collapseRows(exprs(ipsc.expr.annot), factor(fdata.heart$SYMBOL), rownames(ipsc.expr.annot))$datETcollapsed #collapseRows by symbol
colnames(heart.expr.collapse) <- sampleNames(ipsc.expr.annot)
ipsc.collapse <- ExpressionSet(assayData = heart.expr.collapse, phenoData = phenoData(ipsc.expr.annot))
saveRDS.gz(ipsc.collapse, file = "./save/ipsc.collapse.rda")

#heart.test <- GSE10745.raw
#heart.test <- heart.test[,str_detect(heart.test$characteristics_ch1, "Heart") & str_detect(heart.test$characteristics_ch1, "control")]
#heart.test$Genotype <- str_split_fixed(heart.test$characteristics_ch1, " ", 6)[,3]

model.ipsc <- model.matrix( ~ 0 + ipsc.collapse$Status) %>% data.frame
colnames(model.ipsc) <- c("Carrier", "FRDA", "Normal") 

fit.object <- gen.fit(exprs(ipsc.collapse), model.ipsc)
saveRDS.gz(fit.object, file = "./save/fit.object.rda")

#ebam.ipsc <- limma2ebam(fit.object, coef = 1)
ipsc.collapse$Status %<>% factor(levels = c("Healthy", "FA"))
a0.ipsc <- find.a0(exprs(ipsc.collapse), as.integer(ipsc.collapse$Status), B = 1000, rand = 12345)
ebam.ipsc <- ebam(a0.ipsc)
ebam.ipsc.df <- data.frame(Z.score = ebam.ipsc@z, Posterior = ebam.ipsc@posterior)
ebam.ipsc.df$Symbol <- rownames(ebam.ipsc.df)
ebam.ipsc.df %<>% arrange(desc(Posterior))
saveRDS.gz(ebam.ipsc.df, "ebam.ipsc.df.rda")

decide <- list(c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot <- map_df(decide, gen.decide, fit.object) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection", decide.plot, 6, 4) #Plot different significance cutoffs

ratio.exp <- gen.ratios(ipsc.collapse)
saveRDS.gz(ratio.exp, file = "./save/ratio.exp.rda")

decide.final <- gen.decide(c("none", 0.001), fit.object, TRUE) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute signifance cutoff p < 0.001, no FDR adjustment
gen.decideplot("./selected_threshold", decide.final) #Plot cutoff

de.object <- read_tsv("./fit_none.tsv") #Read in unadjusted fit object
rownames(de.object) <- featureNames(pbmc.collapse)

fit.selection <- gen.tables(de.object, pbmc.collapse, ratio.exp, "pLess001") #create differential expression table for unadjusted fit
saveRDS.gz(fit.selection, file = "./save/fit.selection.rda")

clust.none <- gen.anova(fit.selection, "none")

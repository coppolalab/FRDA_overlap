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
library(illuminaio)

match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|"))

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
    num.genes <- map(any, na.rm = TRUE)
    num.genes <- length(which(apply(results, 1, any, na.rm = TRUE)))  #Make this better
    mysum <- summary(results)[-2,] %>% data.frame#Eliminate the row for no change in expression
    mysum[2,] <- -(mysum[2,])
    colnames(mysum) <- "KIKI_vs_WT"
    mysum <- data.frame("Test" = paste(test[1], " p<", test[2], sep = ""), "Num" = paste(num.genes, "Genes", sep = " "), "Direction" = c("positive", "negative"), mysum)
    return(mysum)
}

#Plot statistical cutoff tests
gen.decideplot <- function(filename, decide.plot, width.plot = 6, height.plot = 7)
{
    print(decide.plot)
    decide.plot$variable <- str_replace_all(decide.plot$variable, "_", " ")
    p <- ggplot()
    p <- p + geom_bar(data = filter(decide.plot, Direction == "positive"),  aes(x = variable, y = value), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = filter(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = variable, y = value, ymax = max(value) + 110, label = value), hjust = -0.3, position = position_dodge(width = 1))
    p <- p + geom_bar(data = filter(decide.plot, Direction == "negative"),  aes(x = variable, y = value), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = filter(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = variable, y = value, ymax = min(value) - 110, label = abs(value)), hjust = 1.3, position = position_dodge(width = 1))
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
    all.kiki <- exprs(lumi.final[,lumi.final$Genotype == "KIKI"])
    all.wt <- exprs(lumi.final[,lumi.final$Genotype == "WT"])
    all.wt.means <- rowMeans(all.wt)

    coef.kiki.wt <- all.kiki - all.wt.means
    colnames(coef.kiki.wt) <- paste(colnames(coef.kiki.wt), ".kikiwt", sep = "")

    coefs <- data.frame("Symbol" = rownames(coef.kiki.wt), coef.kiki.wt)
    all.samples <- data.frame("Symbol" = rownames(lumi.final), exprs(lumi.final))
    colnames(all.samples)[2:ncol(all.samples)] <- paste(colnames(all.samples[2:ncol(all.samples)]), "expr", sep = ".") 
    ratio.exp <- merge(coefs, all.samples)
    return(ratio.exp)
}

gen.fit <- function(dataset, model.design)
{
    fit <- lmFit(dataset, design = model.design)
    contrasts.anova <- makeContrasts(KIKI_vs_WT = KIKI - WT, levels = model.design)
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

    fitsel.return.all <- select(fitsel.ratio.all, Symbol, contains("Coef"), contains("p.value"), contains("_vs_"), matches("^t$"), A, contains("kikiwt"), contains("expr")) %>% arrange(desc(t))

    anovalist <- select(fitsel.return.all, contains("KIKI_vs_WT")) %>% apply(1, any, na.rm = T) %>% which
    fitsel.return <- fitsel.return.all[anovalist,]
    coef.cols <- colnames(fitsel.return) %>% str_detect("Coef") %>% which
    gen.workbook(fitsel.return, paste("./significant_geneList_", suffix, ".xlsx", sep = ""))

    write.csv(fitsel.return.all, paste("./complete_genelist_time1_", suffix, ".csv", sep = ""), row.names = FALSE)
    return(fitsel.return)
}

#Make anova objects
gen.anova <- function(dataset, suffix)
{
    plot.kikiwt <- filter(dataset, KIKI_vs_WT != 0) %>% select(contains("kikiwt"))
    gen.anova.heatmap(paste("./anova_heatmap_kiki_vs_wt", suffix, sep = "_"), plot.kikiwt, "KIKI_vs_WT") 
    return(plot.kikiwt)
}

#Generate anova heatmaps
gen.anova.heatmap <- function(filename, dataset, maintitle)
{ 
    CairoPDF(filename, width = 8, height = 8)
    heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 0.3, labRow = "", labCol = "", keysize = 0.9)
    dev.off()
}

#3 tissues, wt and KIKI, treated and untreated
GSE10745.raw <- getGEO("GSE10745", destdir = ".")[[1]]
kiki.pdata.raw <- select(pData(GSE10745.raw), title, geo_accession, characteristics_ch1, description)
kiki.chardata <- map(kiki.pdata.raw$characteristics_ch1, str_split_fixed, ",", n = 3) %>% reduce(rbind) %>% data.frame
kiki.charsplit <- map(kiki.chardata$X1, str_split_fixed, " ", n = 4) %>% reduce(rbind) %>% data.frame
kiki.treatment <- map(kiki.chardata$X2, str_replace_all, "\\-treated|treated |with |^ ", "") %>% reduce(rbind) 
kiki.slideid <- map(kiki.pdata.raw$description, str_replace_all, ": .*$", "") %>% reduce(rbind)
kiki.pdata <- data.frame(Sample.Name = kiki.pdata.raw$title, GSE.ID = kiki.pdata.raw$geo_accession, Slide.ID = kiki.slideid, Tissue = kiki.charsplit$X1, Genotype = kiki.charsplit$X3, Treatment = kiki.treatment)
rownames(kiki.pdata) <- kiki.pdata$Sample.Name
                                                               
kiki.supplemental <- getGEOSuppFiles('GSE10745', baseDir = '.')

kiki.raw <- lumiR("./GSE10745/GSE10745_Microarray_raw_data.txt", convertNuID = FALSE, lib.mapping = 'lumiMouseIDMapping')
kiki.profile <- readBGX("./GSE10745/GPL6103_Illumina_MouseRef-8_V1_1_R1_11234312_A.bgx")$probes
kiki.match <- select(kiki.profile, Probe_Id, Obsolete_Probe_Id)
kiki.match$Obsolete_Probe_Id %<>% str_replace_all("S\\-.*$","S") %>% str_replace_all("I\\-.*$","I") %>% str_replace_all("A\\-.*$","A")
kiki.raw.names <- merge(data.frame(featureNames(kiki.raw)), kiki.match, all.x = TRUE, all.y = FALSE, by.x = "featureNames.kiki.raw.", by.y = "Obsolete_Probe_Id", sort = FALSE)
not.matched <- filter(kiki.raw.names, is.na(Probe_Id))$featureNames.kiki.raw. %>% as.character %>% match.exact
kiki.matched <- filter(kiki.raw.names, !is.na(Probe_Id) & !duplicated(featureNames.kiki.raw.)) 

kiki.raw.match <- kiki.raw[!grepl(not.matched, featureNames(kiki.raw)),]
featureNames(kiki.raw.match) <- kiki.matched$Probe_Id
kiki.nuIDs <- IlluminaID2nuID(IlluminaID = featureNames(kiki.raw.match), lib.mapping = 'lumiMouseIDMapping') %>% data.frame %>% select(nuID)
featureNames(kiki.raw.match) <- unlist(kiki.nuIDs) %>% as.character

kiki.included <- paste(kiki.pdata$Slide.ID, collapse = "|")
kiki.reduce <- kiki.raw.match[,str_detect(kiki.included, sampleNames(kiki.raw.match))]
pData(kiki.reduce) <- kiki.pdata
sampleNames(kiki.reduce) <- kiki.pdata$Sample.Name

kiki.norm <- lumiT(kiki.reduce)
kiki.heart <- kiki.norm[,(kiki.norm$Tissue == "Heart" & kiki.norm$Treatment == "control")]
kiki.cerebellum <- kiki.norm[,(kiki.norm$Tissue == "Cerebellum" & kiki.norm$Treatment == "control")]
kiki.brain <- kiki.norm[,(kiki.norm$Tissue == "Brain" & kiki.norm$Treatment == "control")]

kiki.heart.norm <- lumiN(kiki.heart, method = "rsn") #Normalize with robust spline regression
kiki.heart.cutoff <- detectionCall(kiki.heart.norm) #Get the count of probes which passed the detection threshold per sample
kiki.heart.expr <- kiki.heart.norm[which(kiki.heart.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.kiki.heart <- getSYMBOL(rownames(kiki.heart.expr), 'lumiMouseAll.db') %>% is.na #Determine which remaining probes are unannotated
kiki.heart.annot <- kiki.heart.expr[!symbols.kiki.heart,] #Drop any probe which is not annotated
saveRDS.gz(kiki.heart.annot, file = "kiki.heart.annot.rda")

gen.boxplot("heart_intensity_norm.jpg", kiki.heart.annot, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
genotype.colors <- c("red", "blue")
tissue.colors <- "black"
heatmap.bars <- gen.heatmapbars(genotype.colors, tissue.colors, pData(kiki.heart.annot)) %>% data.frame
kiki.heart.annot$Genotype.Color <- as.character(heatmap.bars$Genotype.Color)
kiki.heart.annot$Tissue.Color <- as.character(heatmap.bars$Tissue.Color)
gen.heatmap("heart_heatmap", kiki.heart.annot, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.heart <- gen.IACcluster("heart_IAC", exprs(kiki.heart.annot), "Outlier removed")
sd.raw.norm.heart <- gen.sdplot("heart_sd", IAC.norm.heart, "Outlier removed")
saveRDS.gz(IAC.norm.heart, file = "./save/IAC.norm.heart.rda")
saveRDS.gz(sd.raw.norm.heart, file = "./save/sd.raw.norm.heart.rda")

kiki.heart.mads <- apply(exprs(kiki.heart.annot), 1, mad)
kiki.heart.ordered <- order(kiki.heart.mads, decreasing = TRUE)
top500.dist.norm.heart <- gen.topgenes("heart_heatmap_500", kiki.heart.annot, "Clustering Based on the Top 500 Most Variable Genes", kiki.heart.ordered, 500)
top1000.dist.norm.heart <- gen.topgenes("heart_heatmap_1000", kiki.heart.annot, "Clustering Based on the Top 1000 Most Variable Genes", kiki.heart.ordered, 1000)

cm1.heart <- cmdscale(top1000.dist.norm.heart, eig = TRUE)
gen.pca("heart_mds_genotype", cm1.heart, phenoData(kiki.heart.annot), genotype.colors, "Genotype")

#Batch correction skipped
#No covariates to remove

fdata.heart <- fData(kiki.heart.annot) #Get featureData from lumi object
fdata.heart$SYMBOL <- getSYMBOL(rownames(fdata.heart), "lumiMouseAll.db")
heart.expr.collapse <- collapseRows(exprs(kiki.heart.annot), factor(fdata.heart$SYMBOL), rownames(kiki.heart.annot))$datETcollapsed #collapseRows by symbol
colnames(heart.expr.collapse) <- sampleNames(kiki.heart.annot)
kiki.heart.collapse <- ExpressionSet(assayData = heart.expr.collapse, phenoData = phenoData(kiki.heart.annot))
saveRDS.gz(kiki.heart.collapse, file = "./save/kiki.heart.collapse.rda")

heart.test <- GSE10745.raw
heart.test <- heart.test[,str_detect(heart.test$characteristics_ch1, "Heart") & str_detect(heart.test$characteristics_ch1, "control")]
heart.test$Genotype <- str_split_fixed(heart.test$characteristics_ch1, " ", 6)[,3]

model.heart <- model.matrix( ~ 0 + kiki.heart.collapse$Genotype)
colnames(model.heart) <- c("KIKI", "WT")

fit.object.heart <- gen.fit(exprs(kiki.heart.collapse), model.heart)
saveRDS.gz(fit.object.heart, file = "./save/fit.object.heart.rda")

decide <- list(c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot.heart <- map_df(decide, gen.decide, fit.object.heart) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection_heart", decide.plot.heart, 10, 10) #Plot different significance cutoffs

ratio.exp.heart <- gen.ratios(kiki.heart.collapse)
saveRDS.gz(ratio.exp.heart, file = "./save/ratio.exp.heart.rda")

decide.final <- gen.decide(c("none", 0.001), fit.object.heart, TRUE, "_heart") %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute signifance cutoff p < 0.001, no FDR adjustment
gen.decideplot("./selected_threshold_heart", decide.final) #Plot cutoff

de.object.heart <- read_tsv("./fit_none_heart.tsv") #Read in unadjusted fit object
rownames(de.object.heart) <- featureNames(kiki.heart.collapse)

fit.selection.heart <- gen.tables(de.object.heart, kiki.heart.collapse, ratio.exp.heart, "pLess001_heart") #create differential expression table for unadjusted fit
saveRDS.gz(fit.selection.heart, file = "./save/fit.selection.heart.rda")

clust.none <- gen.anova(fit.selection.heart, "none_heart")

kiki.cerebellum.norm <- lumiN(kiki.cerebellum, method = "rsn") #Normalize with robust spline regression
kiki.cerebellum.cutoff <- detectionCall(kiki.cerebellum.norm) #Get the count of probes which passed the detection threshold per sample
kiki.cerebellum.expr <- kiki.cerebellum.norm[which(kiki.cerebellum.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.kiki.cerebellum <- getSYMBOL(rownames(kiki.cerebellum.expr), 'lumiMouseAll.db') %>% is.na #Determine which remaining probes are unannotated
kiki.cerebellum.annot <- kiki.cerebellum.expr[!symbols.kiki.cerebellum,] #Drop any probe which is not annotated
saveRDS.gz(kiki.cerebellum.annot, file = "kiki.cerebellum.annot.rda")

gen.boxplot("cerebellum_intensity_norm.jpg", kiki.cerebellum.annot, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
genotype.colors <- c("red", "blue")
tissue.colors <- "black"
heatmap.bars <- gen.heatmapbars(genotype.colors, tissue.colors, pData(kiki.cerebellum.annot)) %>% data.frame
kiki.cerebellum.annot$Genotype.Color <- as.character(heatmap.bars$Genotype.Color)
kiki.cerebellum.annot$Tissue.Color <- as.character(heatmap.bars$Tissue.Color)
gen.heatmap("cerebellum_heatmap", kiki.cerebellum.annot, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.cerebellum <- gen.IACcluster("cerebellum_IAC", exprs(kiki.cerebellum.annot), "Outlier removed")
sd.raw.norm.cerebellum <- gen.sdplot("cerebellum_sd", IAC.norm.cerebellum, "Outlier removed")
saveRDS.gz(IAC.norm.cerebellum, file = "./save/IAC.norm.cerebellum.rda")
saveRDS.gz(sd.raw.norm.cerebellum, file = "./save/sd.raw.norm.cerebellum.rda")

kiki.cerebellum.mads <- apply(exprs(kiki.cerebellum.annot), 1, mad)
kiki.cerebellum.ordered <- order(kiki.cerebellum.mads, decreasing = TRUE)
top500.dist.norm.cerebellum <- gen.topgenes("cerebellum_heatmap_500", kiki.cerebellum.annot, "Clustering Based on the Top 500 Most Variable Genes", kiki.cerebellum.ordered, 500)
top1000.dist.norm.cerebellum <- gen.topgenes("cerebellum_heatmap_1000", kiki.cerebellum.annot, "Clustering Based on the Top 1000 Most Variable Genes", kiki.cerebellum.ordered, 1000)

cm1.cerebellum <- cmdscale(top1000.dist.norm.cerebellum, eig = TRUE)
gen.pca("cerebellum_mds_genotype", cm1.cerebellum, phenoData(kiki.cerebellum.annot), genotype.colors, "Genotype")

#Batch correction skipped
#No covariates to remove

fdata.cerebellum <- fData(kiki.cerebellum.annot) #Get featureData from lumi object
fdata.cerebellum$SYMBOL <- getSYMBOL(rownames(fdata.cerebellum), "lumiMouseAll.db")
cerebellum.expr.collapse <- collapseRows(exprs(kiki.cerebellum.annot), factor(fdata.cerebellum$SYMBOL), rownames(kiki.cerebellum.annot))$datETcollapsed #collapseRows by symbol
colnames(cerebellum.expr.collapse) <- sampleNames(kiki.cerebellum.annot)
kiki.cerebellum.collapse <- ExpressionSet(assayData = cerebellum.expr.collapse, phenoData = phenoData(kiki.cerebellum.annot))
saveRDS.gz(kiki.cerebellum.collapse, file = "./save/kiki.cerebellum.collapse.rda")

model.cerebellum <- model.matrix( ~ 0 + kiki.cerebellum.collapse$Genotype)
colnames(model.cerebellum) <- c("KIKI", "WT")

fit.object.cerebellum <- gen.fit(exprs(kiki.cerebellum.collapse), model.cerebellum)
saveRDS.gz(fit.object.cerebellum, file = "./save/fit.object.cerebellum.rda")

ratio.exp.cerebellum <- gen.ratios(kiki.cerebellum.collapse)
saveRDS.gz(ratio.exp.cerebellum, file = "./save/ratio.exp.cerebellum.rda")

decide <- list(c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot.cerebellum <- map_df(decide, gen.decide, fit.object.cerebellum) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection_cerebellum", decide.plot.cerebellum, 10, 10) #Plot different significance cutoffs

decide.final <- gen.decide(c("none", 0.001), fit.object.cerebellum, TRUE, "_cerebellum") %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute signifance cutoff p < 0.001, no FDR adjustment
gen.decideplot("./selected_threshold_cerebellum", decide.final) #Plot cutoff

de.object.cerebellum <- read_tsv("./fit_none_cerebellum.tsv") #Read in unadjusted fit object
rownames(de.object.cerebellum) <- featureNames(kiki.cerebellum.collapse)

fit.selection.cerebellum <- gen.tables(de.object.cerebellum, kiki.cerebellum.collapse, ratio.exp.cerebellum, "pLess001_cerebellum") #create differential expression table for unadjusted fit
saveRDS.gz(fit.selection.cerebellum, file = "./save/fit.selection.cerebellum.rda")

clust.none <- gen.anova(fit.selection.cerebellum, "none_cerebellum")

kiki.brain.norm <- lumiN(kiki.brain, method = "rsn") #Normalize with robust spline regression
kiki.brain.cutoff <- detectionCall(kiki.brain.norm) #Get the count of probes which passed the detection threshold per sample
kiki.brain.expr <- kiki.brain.norm[which(kiki.brain.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.kiki.brain <- getSYMBOL(rownames(kiki.brain.expr), 'lumiMouseAll.db') %>% is.na #Determine which remaining probes are unannotated
kiki.brain.annot <- kiki.brain.expr[!symbols.kiki.brain,] #Drop any probe which is not annotated
saveRDS.gz(kiki.brain.annot, file = "kiki.brain.annot.rda")

gen.boxplot("brain_intensity_norm.jpg", kiki.brain.annot, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
genotype.colors <- c("red", "blue")
tissue.colors <- "black"
heatmap.bars <- gen.heatmapbars(genotype.colors, tissue.colors, pData(kiki.brain.annot)) %>% data.frame
kiki.brain.annot$Genotype.Color <- as.character(heatmap.bars$Genotype.Color)
kiki.brain.annot$Tissue.Color <- as.character(heatmap.bars$Tissue.Color)
gen.heatmap("brain_heatmap", kiki.brain.annot, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.brain <- gen.IACcluster("brain_IAC", exprs(kiki.brain.annot), "Outlier removed")
sd.raw.norm.brain <- gen.sdplot("brain_sd", IAC.norm.brain, "Outlier removed")
saveRDS.gz(IAC.norm.brain, file = "./save/IAC.norm.brain.rda")
saveRDS.gz(sd.raw.norm.brain, file = "./save/sd.raw.norm.brain.rda")

kiki.brain.mads <- apply(exprs(kiki.brain.annot), 1, mad)
kiki.brain.ordered <- order(kiki.brain.mads, decreasing = TRUE)
top500.dist.norm.brain <- gen.topgenes("brain_heatmap_500", kiki.brain.annot, "Clustering Based on the Top 500 Most Variable Genes", kiki.brain.ordered, 500)
top1000.dist.norm.brain <- gen.topgenes("brain_heatmap_1000", kiki.brain.annot, "Clustering Based on the Top 1000 Most Variable Genes", kiki.brain.ordered, 1000)

cm1.brain <- cmdscale(top1000.dist.norm.brain, eig = TRUE)
gen.pca("brain_mds_genotype", cm1.brain, phenoData(kiki.brain.annot), genotype.colors, "Genotype")

#Batch correction skipped
#No covariates to remove

fdata.brain <- fData(kiki.brain.annot) #Get featureData from lumi object
fdata.brain$SYMBOL <- getSYMBOL(rownames(fdata.brain), "lumiMouseAll.db")
brain.expr.collapse <- collapseRows(exprs(kiki.brain.annot), factor(fdata.brain$SYMBOL), rownames(kiki.brain.annot))$datETcollapsed #collapseRows by symbol
colnames(brain.expr.collapse) <- sampleNames(kiki.brain.annot)
kiki.brain.collapse <- ExpressionSet(assayData = brain.expr.collapse, phenoData = phenoData(kiki.brain.annot))
saveRDS.gz(kiki.brain.collapse, file = "./save/kiki.brain.collapse.rda")

model.brain <- model.matrix( ~ 0 + kiki.brain.collapse$Genotype)
colnames(model.brain) <- c("KIKI", "WT")

fit.object.brain <- gen.fit(exprs(kiki.brain.collapse), model.brain)
saveRDS.gz(fit.object.brain, file = "./save/fit.object.brain.rda")

ratio.exp.brain <- gen.ratios(kiki.brain.collapse)
saveRDS.gz(ratio.exp.brain, file = "./save/ratio.exp.brain.rda")

decide <- list(c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot.brain <- map_df(decide, gen.decide, fit.object.brain) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection_brain", decide.plot.brain, 10, 10) #Plot different significance cutoffs

decide.final <- gen.decide(c("none", 0.001), fit.object.brain, TRUE, "_brain") %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute signifance cutoff p < 0.001, no FDR adjustment
gen.decideplot("./selected_threshold_brain", decide.final) #Plot cutoff

de.object.brain <- read_tsv("./fit_none_brain.tsv") #Read in unadjusted fit object
rownames(de.object.brain) <- featureNames(kiki.brain.collapse)

fit.selection.brain <- gen.tables(de.object.brain, kiki.brain.collapse, ratio.exp.brain, "pLess001_brain") #create differential expression table for unadjusted fit
saveRDS.gz(fit.selection.brain, file = "./save/fit.selection.brain.rda")

clust.none <- gen.anova(fit.selection.brain, "none_brain")

decide.plot.heart$variable <- "Heart"
decide.plot.cerebellum$variable <- "Cerebellum"
decide.plot.brain$variable <- "Brain"
plot.all <- rbind(decide.plot.heart, decide.plot.cerebellum, decide.plot.brain)
gen.decideplot("./threshold_selection", plot.all, width.plot = 10)


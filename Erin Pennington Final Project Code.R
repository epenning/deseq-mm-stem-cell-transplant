library(TCGAbiolinks)         # Bioconductor, searches TCGA
library(SummarizedExperiment) # Bioconductor, organizes TCGA search results
library(DESeq2)               # Bioconductor, performs DEG analysis
library(biomaRt)              # Bioconductor, annotates the DEG analysis
library(ggplot2)              # CRAN, for plotting
library(ggrepel)              # CRAN, allows for automatic placement of text labels in scatterplots
library(viridis)              # CRAN, color-blind friendly continuous color gradient
library(dplyr)

## Predefine functions for visualizations ######################################
# Credit to Brandon Burgman & Nolan Bentley for these functions

# Create a function to visualize values
countStats <- function(counts,
                       margin,
                       scaleBySum = F,
                       scaleMargin,
                       rmNa = F) {
    #Scale by the sum of the opposite margin and multiple by mean count
    if (scaleBySum) {
        counts <-
            apply(counts, scaleMargin, function(x) {
                x / sum(x, na.rm = rmNa)
            }) * mean(counts)
    }

    #Calculate vectors
    meanVec <- apply(counts, margin, mean, na.rm = rmNa)
    sdVec   <- apply(counts, margin, sd  , na.rm = rmNa)

    transMeanVec <-
        apply(log10(counts + 1), margin, mean, na.rm = rmNa)
    transSdVec   <-
        apply(log10(counts + 1), margin, sd  , na.rm = rmNa)

    out <- list(
        meanVec      = meanVec,
        sdVec        = sdVec,
        transMeanVec = transMeanVec,
        transSdVec   = transSdVec,
        margin       = margin
    )

    class(out) <- c(class(out), "countStats")

    return(out)
}

countPlot <- function(countStatOut, breaks, rmNa = F) {
    require("ggplot2")
    #Check if countStats class
    if (!"countStats" %in% class(countStatOut)) {
        stop("countStatOut needs to be countStats object")
    }

    #Create objects from countStats
    for (i in 1:length(countStatOut)) {
        assign(names(countStatOut)[i], countStatOut[[i]])
    }

    #Reset plot device
    while (!is.null(dev.list())) {
        dev.off()
    }

    #Color by order
    colorVec <- viridis::viridis(length(meanVec))
    #Randomize plotting order
    set.seed(1)
    orderVec <- sample(1:length(meanVec), size = length(meanVec))

    #Relationship between count and sd
    par(mfrow = c(3, 2))

    hist(meanVec[orderVec], breaks,
         main = "Mean count")
    plot(meanVec[orderVec],
         sdVec[orderVec],
         xlab = "Mean count",
         ylab = "SD of count",
         col = colorVec[orderVec])

    hist(log10(meanVec[orderVec] + 1),
         breaks = breaks,
         main   = "log10(Mean count + 1)")
    plot(
        log10(meanVec[orderVec] + 1),
        log10(sdVec[orderVec] + 1),
        xlab = "log10(Mean count + 1)",
        ylab = "log10(SD of count + 1)",
        col = colorVec[orderVec]
    )

    hist(transMeanVec[orderVec], breaks,
         main = "Mean log10(count+1)")
    plot(
        transMeanVec[orderVec],
        transSdVec[orderVec],
        xlab = "Mean log10(count+1)",
        ylab = "SD of log10(count+1)",
        col = colorVec[orderVec]
    )

    #2d histograms of these relationships
    nPlots  <- 2
    currDf  <- data.frame(mean = meanVec, sd = sdVec)[orderVec,]
    transDf <-
        data.frame(mean = transMeanVec, sd = transSdVec)[orderVec,]
    smMeth  <- "gam"
    smForm  <- y ~ s(x, bs = "cs")

    currBinDim <- c((
        max(currDf$mean, na.rm = rmNa) -
            min(currDf$mean, na.rm = rmNa)
    ) / (30 * 1),
    (max(currDf$sd  , na.rm = rmNa) -
         min(currDf$sd  , na.rm = rmNa)) / (30 * nPlots))
    p1 <- ggplot(currDf, aes(mean, sd)) +
        geom_bin_2d(binwidth = currBinDim) +
        scale_fill_viridis_c(trans = "log10") +
        geom_smooth(
            method = smMeth,
            formula = smForm,
            col = "grey50",
            se = F
        ) +
        labs(x = "mean(count)", y = "sd(count)") +
        theme_bw()

    currBinDim <- c((max(log10(currDf$mean + 1), na.rm = rmNa) -
                         min(log10(currDf$mean + 1), na.rm = rmNa)) / (30 *
                                                                           1),
                    (max(log10(currDf$sd   + 1), na.rm = rmNa) -
                         min(log10(currDf$sd   + 1), na.rm = rmNa)) / (30 *
                                                                           nPlots))
    p2 <- ggplot(currDf, aes(log10(mean + 1), log10(sd + 1))) +
        geom_bin_2d(binwidth = currBinDim) +
        scale_fill_viridis_c(trans = "log10") +
        geom_smooth(
            method = smMeth,
            formula = smForm,
            col = "grey50",
            se = F
        ) +
        labs(x = "log10( mean(count) + 1)", y = "log10( sd(count) + 1)") +
        theme_bw()

    currBinDim <- c((
        max(transDf$mean, na.rm = rmNa) -
            min(transDf$mean, na.rm = rmNa)
    ) / (30 * 1),
    (max(transDf$sd  , na.rm = rmNa) -
         min(transDf$sd  , na.rm = rmNa)) / (30 * nPlots))
    p3 <- ggplot(transDf, aes(mean, sd)) +
        geom_bin_2d(binwidth = currBinDim) +
        scale_fill_viridis_c(trans = "log10") +
        geom_smooth(
            method = smMeth,
            formula = smForm,
            col = "grey50",
            se = F
        ) +
        labs(x = "mean( log10(count + 1) )", y = "sd( log10(count + 1) )") +
        theme_bw()

    ##Combine them
    ##### margin is top, right, bottom, left
    combFormatting <-
        theme(plot.margin = margin(1, 1, 0, 0),
              legend.position = "none")
    pCow <- cowplot::plot_grid(
        p1 + combFormatting,
        p2 + combFormatting,
        p3 + combFormatting,
        align = "none",
        ncol = 3
    )
    return(pCow)
}

## Retrieve data ###############################################################

# Load multiple myeloma RNA-seq counts from TCGA
## Search for multiple myeloma gene expression profiles

# Project: Multiple Myeloma CoMMpass Study
query <- TCGAbiolinks::GDCquery(
    project       = "MMRF-COMMPASS",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts",
    legacy        = F
)

## Read count data using search query results

# Extract barcodes from query
barcodes <- query$results[[1]]$cases

#Search for the barcodes
query <- GDCquery(
    project       = c("MMRF-COMMPASS"),
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts",
    legacy        = F,
    barcode       = barcodes
)

#Download the current query files
GDCdownload(query, directory = "GDCdata/")

#Prepare the current query files
commpass.counts <- GDCprepare(query, directory = "GDCdata/")

#Clean up the environment to remove redundant large objects
rm(query)

## Collect stem cell transplant types performed on each sample #################

# Collect all of the treatments for every sample
# (sample is the same as "colnames" in compass.counts)
treatments <-
    bind_rows(colData(commpass.counts)$treatments, .id = "sample")

# Just look at the sample and treatment_type
# (treatment_type indicates stem cell transplants)
treatments <-
    treatments %>% select(sample, treatment = treatment_type)

# Move each treatment into a column with TRUE or FALSE for whether sample had it
processedTreatments <-
    treatments[!duplicated(treatments),] %>% mutate(val = TRUE) %>%
    tidyr::pivot_wider(names_from = treatment,
                       values_from = val,
                       values_fill = FALSE)

# Add new boolean sample annotations for whether they had each treatment
for (i in 3:ncol(processedTreatments)) {
    treatmentName <- colnames(processedTreatments)[i]
    colData(commpass.counts)[treatmentName] <-
        processedTreatments[, i]
}

colData(commpass.counts)$comp <-
    colData(commpass.counts)$`Stem Cell Transplantation, Autologous`

# Count genes and treatment levels before filtering
nrow(assay(commpass.counts))
table(colData(commpass.counts)$`Stem Cell Transplantation, Autologous`)

## Filtering ###################################################################

# Subset counts to exclude allogeneic stem cell transplants
commpass.subset <- commpass.counts %>%
    subset(select = colData(commpass.counts)$`Stem Cell Transplantation, Allogeneic` == FALSE)

# Visualize the counts to identify outliers
rowStats       <-
    countStats(assays(commpass.subset)[[1]], margin = 1)
colStats       <-
    countStats(assays(commpass.subset)[[1]], margin = 2)

pRow       <-
    countPlot(countStatOut = rowStats      , breaks = 1000)
pCol       <- countPlot(countStatOut = colStats      , breaks = 100)

plot(pRow)
plot(pCol)

# Removes samples with mean count sd of more than 100,000 (outliers)
sampleSds <- apply(assays(commpass.subset)[[1]], 2, sd)
hist(sampleSds, 1000)
abline(v = 100000, col = "red")

commpass.subset <- commpass.subset[, sampleSds < 100000]

# Ignore genes with low baseline expression (>=50% of samples have 0 reads)
hist(
    rowMeans(assay(commpass.subset) == 0),
    breaks = 100,
    main = "Proportion of samples with 0 reads",
    xlab = "Proportion no reads"
)
abline(v = 0.5, col = "red")

commpassAssay <- as.data.frame(assay(commpass.subset))
commpass.subset <-
    subset(commpass.subset, subset = rowMeans(commpassAssay == 0) < 0.5)

# Genes with count outliers? Leaving these in
geneMeans <- apply(assays(commpass.subset)[[1]], 1, mean)
mean(geneMeans)
head(sort(geneMeans, decreasing = T), n = 10)
hist(geneMeans, 1000, xlab = "Gene Count Mean", main = "Histogram of Gene Count Mean")
hist(log10(geneMeans), 1000, xlab = "log10(Gene Count Mean)", main = "Histogram of log10(Gene Count Mean)")

geneSds <- apply(assays(commpass.subset)[[1]], 1, sd)
mean(geneSds)
head(sort(geneSds, decreasing = T), n = 10)
hist(geneSds, 1000, xlab = "Gene Count SD", main = "Histogram of Gene Count SD")
hist(log10(geneSds), 1000, xlab = "log10(Gene Count SD)", main = "Histogram of log10(Gene Count SD)")

# Count genes and treatment levels after filtering
nrow(assay(commpass.subset))
table(colData(commpass.subset)$`Stem Cell Transplantation, Autologous`)

## Calculate Differential Expression ###########################################

# Use DESeq to analyze difference between treatments
ddsSE <- DESeq2::DESeqDataSet(commpass.subset, design = ~ comp)
dds <- DESeq(ddsSE)
res <- results(dds)


## Annotate genes with names for human interpretation ##########################

ensembl <- useMart(biomart = "ensembl",
                   dataset = "hsapiens_gene_ensembl",
                   verbose = T)

#Create a vector of ensembl ids to get information about
my.genes   <- rownames(res)

#Retrieve information from the BioMart database
ensemblout <- getBM(
    attributes = c(
        "ensembl_gene_id",
        "external_gene_name",
        "gene_biotype",
        "entrezgene_id"
    ),
    filters   = "ensembl_gene_id",
    values    = my.genes,
    mart      = ensembl
)

#Use match to extract the gene information in the order you asked for
ensemblout <-
    ensemblout[match(my.genes, ensemblout$ensembl_gene_id),]

#add gene names and biomart annotation to results output
res$symbol  <- rowData(dds)$external_gene_name
res$entrez  <- ensemblout$entrezgene_id
res$biotype <- ensemblout$gene_biotype

# Visualize results
head(res)

## Visualize results ###########################################################

## Define significance cut-offs

# Maximum false discovery rate adjusted p-value cutoff:
fdr.cut.off <- 0.0001
# Minimum Log fold DIFFERENCE (absolute value of change) cutoff to be considered
lfc.cut.off <- 1.50

## Visualize differentially expressed genes

ggdf <- as.data.frame(res)

# Determine significance by false discovery rate 0.0001
# Count number of genes significant / up / down regulated
ggdf$padjL <- ggdf$padj < fdr.cut.off
table(ggdf$padjL)

ggdf$upRegL <- ggdf$log2FoldChange > lfc.cut.off
ggdf$dwRegL <- ggdf$log2FoldChange < -lfc.cut.off
table(ggdf$dwRegL | ggdf$upRegL)
table(ggdf$dwRegL)
table(ggdf$upRegL)

ggdf$sig <- ifelse(ggdf$padjL, "Sig.", "Not sig.")
ggdf$sig <- ifelse(ggdf$upRegL, paste(ggdf$sig, "+"), ggdf$sig)
ggdf$sig <- ifelse(ggdf$dwRegL, paste(ggdf$sig, "-"), ggdf$sig)

# Extract the normalized count values
cnt.norm <- counts(dds, normalized = T)

#Isolate the groupings
##Check to make sure we can compare the objects
if (!all(commpass.subset$barcode == colnames(cnt.norm))) {
    stop("Oh, no! The barcode names don't match! Stopping...")
}

# Create logical vectors from the grouping column
groupF <- commpass.subset$comp == F
groupT <- commpass.subset$comp == T

# Create subsets of the normalized counts based on the grouping
cnt.norm.F <- cnt.norm[, groupF]
cnt.norm.T <- cnt.norm[, groupT]

# Calculate normalized counts for groups
all(rownames(ggdf) == rownames(cnt.norm))
normCountsMeanSd <- function(cnt.norm.group) {
    result <- c()
    result[1] <- mean(cnt.norm.group, na.rm = T)
    result[2] <- sd(cnt.norm.group, na.rm = T)
    return(result)
}

ggdf_less <- ggdf %>% arrange(desc(baseMean)) %>% head

# Normalized counts for untreated:
normCountsMeanSd(cnt.norm.F[ggdf$padjL])               # Significant
normCountsMeanSd(cnt.norm.F[ggdf$padjL & ggdf$upRegL]) # Upregulated
normCountsMeanSd(cnt.norm.F[ggdf$padjL &
                                ggdf$dwRegL])          # Downregulated

# Normalized counts for treated:
normCountsMeanSd(cnt.norm.T[ggdf$padjL])               # Significant
normCountsMeanSd(cnt.norm.T[ggdf$padjL & ggdf$upRegL]) # Upregulated
normCountsMeanSd(cnt.norm.T[ggdf$padjL &
                                ggdf$dwRegL])          # Downregulated

# Normalized counts for treated

# Create volcano plot
ggplot(ggdf, aes(
    x = log2FoldChange,
    y = -log10(padj),
    col = sig,
    shape = sig
)) +
    geom_point() +
    ggtitle("MMRF-COMMPASS, Differential Gene Expression with Autologous Stem Cell transplant") +
    xlab("log2(Fold Change)") +
    ylab("-log10(FDR-Adjusted p-values)") +
    geom_vline(xintercept = lfc.cut.off,
               linetype = 'dashed',
               col = "red") +
    geom_vline(xintercept = -lfc.cut.off,
               linetype = 'dashed',
               col = "blue") +
    geom_hline(
        yintercept = -log10(fdr.cut.off),
        linetype = 'dashed',
        col = "green"
    ) +
    geom_text_repel(
        data = head(ggdf[order(ggdf$pvalue),], 10),
        aes(label = symbol),
        show.legend = FALSE,
        size = 6
    ) +
    theme(text = element_text(size = 20))

# Count downregulated biotypes
table(ggdf[ggdf$dwRegL, c("biotype")])
# Count upregulated biotypes
table(ggdf[ggdf$upRegL, c("biotype")])

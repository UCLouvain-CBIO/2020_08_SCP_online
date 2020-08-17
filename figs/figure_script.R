library(export) ## from GitHub:tomwenseleers/export
library(ggVennDiagram) ## from fork GitHub:cvanderaa/ggVennDiagram
library(QFeatures)  ## from GitHub:rformassspectrometry/QFeatures
library(scp) ## from GitHub:UClouvain-CBIO/scp
library(scpdata) ## from GitHub:UClouvain-CBIO/scpdata
library(SingleCellExperiment)
setwd("~/PhD/2020_08_SCP_online/")

## QFeatures demo figure
## ---------------------

## Prepare the data 
data(hlpsms)
hl <- readQFeatures(hlpsms, ecol = 1:10, name = "psms")
hl <- aggregateFeatures(hl, "psms", "Sequence", name = "peptides", fun = colMeans)
hl <- aggregateFeatures(hl, "peptides", "ProteinGroupAccessions", name = "proteins", fun = colMeans)
hl$tag <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
            "130N", "130C", "131")
## Extract to proteins of interest STAT1 and STAT3
stat3 <- hl["P42227-2", , ]
stat3_df <- data.frame(longFormat(stat3))
stat <- hl[c("P42227-2", "P42225"), , ]
## Get the plot and export it to pdf
stat %>%
    longFormat %>%
    data.frame %>%
    mutate(stat3 = ifelse(rowname %in% stat3_df$rowname, "STAT3", "STAT1"),
           assay = factor(assay, levels = c("psms", "peptides", "proteins")),
           colname = sub("^X", "", colname)) %>%
    ggplot(aes(x = colname,
               y = value,
               group = rowname)) +
    geom_line() +
    geom_point() +
    facet_grid(stat3 ~ assay) +
    xlab("Channel") + ylab("Intensity") +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1)) -> 
    p
graph2pdf(p, file = "figs/QFeatures_data.pdf", width = 8, height = 4)

## Export colData table 
## --------------------

data("specht2019v2")

colData(specht2019v2) %>%
    data.frame %>%
    select(-SampleAnnotation) %>%
    as_tibble ->
    cd
table2tex(cd[1:6, ], file = "colData.tex", standAlone = FALSE)

## PSM filtering 
## -------------

load("../scpScripts/20200525-specht2019v2-replication/data/vignette_final.RData")

computeSCR(specht2019v2,
           i = 1:5,
           colDataCol = "SampleType",
           carrierPattern = "Carrier",
           samplePattern = ".") ->
    specht2019v2    

library(tidyverse)
specht2019v2 %>%
    rowDataToDF(i = 1:5, vars = "meanSCR") %>%
    data.frame %>%
    ggplot(aes(x = meanSCR)) +
    geom_histogram() +
    scale_x_log10() +
    geom_vline(xintercept = 0.1) ->
    p
graph2pdf(p, file = "figs/meanSCR.pdf", width = 5, height = 2.5)



## Benchmark the data replication 
## ------------------------------

## Data is obtained when running the SCoPE2 replication vignette
## TODO put file online
load("../scpScripts/20200525-specht2019v2-replication/data/vignette_final.RData")

## Overlap of selected peptides
ggVennDiagram(list(SCoPE2 = rownames(peptides), 
                   SCP = rownames(specht2019v2[["peptides_log"]]))) +
    theme(legend.position = "none") ->
    p
graph2pdf(p, file = "figs/Benchmark_pep_venn.pdf", width = 4, height = 2.5)

## The difference between SCoPE2 and scp peptide data
rows <- intersect(rownames(peptides), 
                  rownames(specht2019v2[["peptides_log"]]))
cols <- intersect(colnames(proteins), 
                  colnames(specht2019v2[["peptides_log"]]))
err <- assay(peptides)[rows, cols] - assay(specht2019v2[["peptides_log"]])[rows, cols]
data.frame(difference = as.vector(err[!is.na(err)])) %>%
    mutate(difference = abs(difference)) %>%
    ggplot() +
    geom_histogram(aes(x = difference)) +
    xlab("|SCoPE2 - scp|") ->
    p
graph2pdf(p, file = "figs/Benchmark_pep_err.pdf", width = 4, height = 2.5)

## Overlap of selected proteins
ggVennDiagram(list(SCoPE2 = rownames(proteins), 
                   SCP = rownames(specht2019v2[["proteins_batchC"]]))) +
    theme(legend.position = "none") ->
    p
graph2pdf(p, file = "figs/Benchmark_prot_venn.pdf", width = 4, height = 2.5)

## The difference between SCoPE2 and scp protein data
rows <- intersect(rownames(proteins), 
                  rownames(specht2019v2[["proteins_batchC"]]))
cols <- intersect(colnames(proteins), 
                  colnames(specht2019v2[["proteins_batchC"]]))
err <- assay(proteins)[rows, cols] - assay(specht2019v2[["proteins_batchC"]])[rows, cols]
data.frame(difference = as.vector(err[!is.na(err)])) %>%
    mutate(difference = abs(difference)) %>%
    ggplot() +
    geom_histogram(aes(x = difference)) +
    xlab("|SCoPE2 - scp|") ->
    p
graph2pdf(p, file = "figs/Benchmark_prot_err.pdf", width = 4, height = 2.5)


## Replicate the wPCA
## ------------------

plotwPCA <- function(sce) {
    ## Extract the protein expression matrix
    X <- assay(sce)
    ## Compute the weights
    w <- rowSums(cor(t(X))^2)
    ## Normalize the data
    X <- sweep(X, 2, colMedians(X, na.rm = TRUE), FUN = "-")
    X <- sweep(X, 1, rowMeans(X, na.rm = TRUE), FUN = "-")
    ## Compute the PCs (code taken from the SCoPE2 script)
    Xw <- diag(w) %*%  X
    Xcor <- cor(Xw)
    pcaRes <- eigen(Xcor)
    pcaPercentVar <- round(pcaRes$values[1:2] / sum(pcaRes$values) * 100)
    ## Start plotting
    data.frame(PC = pcaRes$vectors[, 1:2], 
               colData(sce)) %>%
        ggplot(aes(x = PC.1, y = PC.2, col = SampleType)) +
        geom_point(alpha = 0.5) +
        ## Annotate plot
        xlab(paste0("PC1 (", pcaPercentVar[1], "%)")) +
        ylab(paste0("PC2 (", pcaPercentVar[2], "%)")) +
        ## Adapt the visual style to match the preprint figure
        scale_color_manual(name = "", values = c("#048ABF","#FF5733"),
                           labels = c("Macrophages", "Monocytes")) +
        theme_minimal() +
        theme(legend.position = "top")
}


plotwPCA(proteins) ->
    p1
graph2pdf(p1, file = "figs/wPCA_SCoPE2.pdf", width = 4, height = 4)

specht2019v2 %>%
    transferColDataToAssay("proteins_batchC") %>%
    .[["proteins_batchC"]] %>%
    plotwPCA ->
    p2
graph2pdf(p2, file = "figs/wPCA_scp.pdf", width = 4, height = 4)

data("specht2019v2")
specht2019v2[,, 1:2]  %>%
    longFormat(colDataCols = "SampleType") %>%
    data.frame %>%
    filter(grepl("^M", SampleType)) %>%
    mutate(value = ifelse(is.na(value) | value == 0, 1, value)) %>%
    mutate(primary = paste("Cell", as.numeric(as.factor(primary)))) %>%
    ggplot(aes(x = value, fill = primary)) +
    geom_density(alpha = 0.5, adjust = 0.2) +
    scale_x_log10() +
    theme(text = element_text(size = 14)) +
    xlab("PSM intensity") ->
    p
export::graph2pdf(p, file = "figs/PSM_intensity.pdf", 
                  width = 7, height = 6)


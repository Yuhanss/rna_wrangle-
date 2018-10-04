# rna_wrangle-
To get an example file, type:
git clone https://github.com/tgrn510/RNASeqExample.git

---
title: "RNA_Wrangle"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# Assignment 6: RNA Wrangle
## Load the libraries: Sample data & RNA expression result
```{r}
library(dplyr)
 setwd("/Users/lundai/Desktop/RNASeqExample/")
 samples <- read.csv('sample_info.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
 genes <- read.csv('expression_results.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
```
## 1. Density distribution of expression values within genes$KITA_02 on a log scale
```{r}
d <- density(genes$KITA_02) # returns the density data
plot(density(log2(genes$KITA_02[(genes$KITA_02>0)]))) # plots the results
```
##  2. Density distribution of PF_BASES
```{r}
d <- density(samples$PF_BASES) # returns the density data
 plot(d) # plots the results
```
## 3. Scatterplot of two columns against one another, KITA_01 vs. KITA_03 on a log2 scale
```{r}
plot(log2(genes$KITA_01[(genes$KITA_01>10 |genes$KITA_03>10 )]),log2(genes$KITA_03[(genes$KITA_01>10 |genes$KITA_03>10 )]))
```
## 4. Heatmap showing how well all the samples correlate
```{r}
library(reshape)
library(ggplot2)
samples$uid=rownames(samples)
genes_summary<-data.frame(
 UID=rownames(samples),
 min=minBySample <- sapply(genes, function(x) min(x[x > 0])),
 max=maxBySample <- sapply(genes, function(x) max(x))
 )
corr<-cor(genes)
melted_corr <- melt(corr)
 head(melted_corr)
```
```{r}
ggplot(melted_corr , aes(x = X1, y = X2)) +
 geom_raster(aes(fill = value)) +
 scale_fill_gradient2(low="green", mid="white", high="red", midpoint=0.5) + theme_classic()
```
## 5. Hierarchal clustering dendrogram
```{r}
genes_transsample <- t(genes[c(rep(FALSE,19),TRUE), ])
clusters <- hclust(dist(genes_transsample))
library('dendextend')
dend <- as.dendrogram(clusters)
dend <- rotate(dend, 1:93)
dend <- color_branches(dend, k=4)
par(cex=0.5) # reduces font
plot(dend)
```
## 6. PCA
```{r}
library('plotly')
setwd("/Users/lundai/Desktop/RNASeqExample/")
samples <- read.csv('sample_info.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
genes <- read.csv('expression_results.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
min(genes[genes>0])
# [1] 8.05e-12
genes.log <-log2(genes+8.05e-12)
genes.log.small <- genes.log[seq(1, nrow(genes.log), 20), ]
pca <- prcomp(genes.log.small,center = TRUE,scale. = TRUE)

std_dev <- pca$sdev
  pr_var <- std_dev^2
  pr_var[1:10]
  prop_varex <- pr_var/sum(pr_var)
pcadf<-data.frame(pca$rotation)
samples$uid<-rownames(samples)
pcadf$uid<-rownames(pcadf)
samples<-inner_join(samples,pcadf,by="uid")
plot_ly(samples, x = ~PC2, y = ~PC3, z = ~PC4, size= ~reads,marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 25), color = ~Kit, colors = c('#F25385', '#0099E2')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC2'),
  yaxis = list(title = 'PC3'),
  zaxis = list(title = 'PC4')))
```
```

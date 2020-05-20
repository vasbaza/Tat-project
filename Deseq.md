```{r}
```
###Скачиваем библиотеки
````
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(AnnotationDbi)
library("org.Hs.eg.db", character.only = TRUE)

````
###Загружаем образцы 
````

control.1 <- read.delim("c1.counts.tab", header = TRUE)


control.2 <- read.delim("c2.counts.tab", header = TRUE)


control.3 <- read.delim("c3.counts.tab", header = TRUE)


Tat.0h.1 <- read.delim("0h1.counts.tab", header = TRUE)


Tat.0h.2 <- read.delim("0h2.counts.tab", header = TRUE)


Tat.0h.3 <- read.delim("0h3.counts.tab", header = TRUE)


Tat.16h.1 <- read.delim("16h1.counts.tab", header = TRUE)


Tat.16h.2 <- read.delim("16h2.counts.tab", header = TRUE)



Tat.16h.3 <- read.delim("16h3.counts.tab", header = TRUE)


Tat.s.1 <- read.delim("s1.counts.tab", header = TRUE)


Tat.s.2 <- read.delim("s2.counts.tab", header = TRUE)


Tat.s.3 <- read.delim("s3.counts.tab", header = TRUE)

raw_counts <- cbind(control.1, control.2, control.3, Tat.0h.1, Tat.0h.2, Tat.0h.3, Tat.16h.1, Tat.16h.2, Tat.16h.3, Tat.s.1, Tat.s.2, Tat.s.3)

````
###Deseq2
````

#Stable vs Control and O hours vs Control

condition = c('control', 'control', 'control', 'Tat.0h', 'Tat.0h', 'Tat.0h', 'Tat.16h', 'Tat.16h', 'Tat.16h', 'Tat.s', 'Tat.s', 'Tat.s')
col_matrix <- data.frame(sample = c('control.1', 'control.2', 'control.3', 'Tat.0h.1', 'Tat.0h.2', 'Tat.0h.3', 'Tat.16h.1', 'Tat.16h.2', 'Tat.16h.3', 'Tat.s.1', 'Tat.s.2', 'Tat.s.3'),
condition = c('control', 'control', 'control', 'Tat.0h', 'Tat.0h', 'Tat.0h', 'Tat.16h', 'Tat.16h', 'Tat.16h', 'Tat.s', 'Tat.s', 'Tat.s'))
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
colData = col_matrix,
design = ~ condition)

dds <- DESeq(dds)

#
#normalized_counts <-  counts(dds, normalized=TRUE)
#normalized_counts <- as.data.frame(normalized_counts)
#rownames(normalized_counts) <- rownames(res_s_c)
#
resultsNames(dds)

res_T0_c <- results(dds, name = "condition_Tat.0h_vs_control")
res_s_c <- results(dds, name = "condition_Tat.s_vs_control")
res_T0_c <- as.data.frame(res_T0_c)
res_s_c <- as.data.frame(res_s_c)

#16 hours vs 0 hours
raw_counts_16_0 <- cbind(Tat.0h.1, Tat.0h.2, Tat.0h.3, Tat.16h.1, Tat.16h.2, Tat.16h.3)
condition = c('Tat.0h', 'Tat.0h', 'Tat.0h', 'Tat.16h', 'Tat.16h', 'Tat.16h')
col_matrix <- data.frame(sample = c('Tat.0h.1', 'Tat.0h.2', 'Tat.0h.3', 'Tat.16h.1', 'Tat.16h.2', 'Tat.16h.3'),
condition = c('Tat.0h', 'Tat.0h', 'Tat.0h', 'Tat.16h', 'Tat.16h', 'Tat.16h'))
dds_16_0 <- DESeqDataSetFromMatrix(countData = raw_counts_16_0,
colData = col_matrix,
design = ~ condition)

dds_16_0 <- DESeq(dds_16_0)
resultsNames(dds_16_0)

res_T16_T0 <- results(dds_16_0, name = "condition_Tat.16h_vs_Tat.0h")
res_T16_T0 <- as.data.frame(res_T16_T0)

````
###Добавляем названия генов в кодировке Ensembl
````

FeatureID <- read.delim('FeatureID.tab', header = TRUE) # здесь названия взяты из таблицы с каунтами
rownames(res_T0_c) <- FeatureID$FeatureID
rownames(res_s_c) <- FeatureID$FeatureID
rownames(res_T16_T0) <- FeatureID$FeatureID

````
###Добавляем названия генов в кодировках Symbol, EntrezID, Genename
````
res_T0_c$symbol = mapIds(org.Hs.eg.db,
keys = row.names(res_T0_c),
column = "SYMBOL",
keytype = "ENSEMBL",
multiVals = "first")
res_T0_c$entrez = mapIds(org.Hs.eg.db,
keys = row.names(res_T0_c),
column = "ENTREZID",
keytype = "ENSEMBL",
multiVals = "first")
res_T0_c$name = mapIds(org.Hs.eg.db,
keys = row.names(res_T0_c),
column = "GENENAME",
keytype = "ENSEMBL",
multiVals = "first")
res_T0_c$foldChange <- 2^(res_T0_c$log2FoldChange)

res_s_c$symbol = mapIds(org.Hs.eg.db,
keys = row.names(res_s_c),
column = "SYMBOL",
keytype = "ENSEMBL",
multiVals = "first")
res_s_c$entrez = mapIds(org.Hs.eg.db,
keys = row.names(res_s_c),
column = "ENTREZID",
keytype = "ENSEMBL",
multiVals = "first")
res_s_c$name = mapIds(org.Hs.eg.db,
keys = row.names(res_s_c),
column = "GENENAME",
keytype = "ENSEMBL",
multiVals = "first")
res_s_c$foldChange <- 2**(res_s_c$log2FoldChange)

res_T16_T0$symbol = mapIds(org.Hs.eg.db,
keys = row.names(res_T16_T0),
column = "SYMBOL",
keytype = "ENSEMBL",
multiVals = "first")
res_T16_T0$entrez = mapIds(org.Hs.eg.db,
keys = row.names(res_T16_T0),
column = "ENTREZID",
keytype = "ENSEMBL",
multiVals = "first")
res_T16_T0$name = mapIds(org.Hs.eg.db,
keys = row.names(res_T16_T0),
column = "GENENAME",
keytype = "ENSEMBL",
multiVals = "first")
res_T16_T0$foldChange <- 2**(res_T16_T0$log2FoldChange)

````
###Сортируем по padj
````
res_s_c_padj <- res_s_c[res_s_c$padj <= 0.05,]
res_s_c_padj <- na.omit(res_s_c_padj)

res_T0_c_padj <- res_T0_c[res_T0_c$padj <= 0.05,]
res_T0_c_padj <- na.omit(res_T0_c_padj)

res_T16_T0_padj <- res_T16_T0[res_T16_T0$padj <= 0.05,]
res_T16_T0_padj <- na.omit(res_T16_T0_padj)

````
###Разбиваем на группы (повышение и понижение экспрессии)
````

s_c_up <- res_s_c_padj[res_s_c_padj$log2FoldChange > 0,]
s_c_up <- s_c_up[order(- (as.vector(s_c_up$stat))),]

s_c_down <- res_s_c_padj[res_s_c_padj$log2FoldChange < 0,]
s_c_down <- s_c_down[order(- (as.vector(s_c_down$stat))),]


T0_c_up <- res_T0_c_padj[res_T0_c_padj$log2FoldChange > 0,]
T0_c_up <- T0_c_up[order(- (as.vector(T0_c_up$stat))),]

T0_c_down <- res_T0_c_padj[res_T0_c_padj$log2FoldChange < 0,]
T0_c_down <- T0_c_down[order(- (as.vector(T0_c_down$stat))),]


T16_T0_up <- res_T16_T0_padj[res_T16_T0_padj$log2FoldChange > 0,]
T16_T0_up <- T16_T0_up[order(- (as.vector(T16_T0_up$stat))),]

T16_T0_down <- res_T16_T0_padj[res_T16_T0_padj$log2FoldChange < 0,]
T16_T0_down <- T16_T0_down[order(- (as.vector(T16_T0_down$stat))),]
```
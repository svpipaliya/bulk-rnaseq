### (c) Shweta V. Pipaliya 
### 2022-04-09 
## This script performs enrichment analyses using clusterprofiler, gseGO, and enrichR 
### run from command-line: Rscript --vanilla enrich_anal.svp.R
### Bioconductor packages in this script require R v.4.2.0

## Set seed 
set.seed(777)

## Set current working directory to output DGE summary results and plots
wd <- getwd() 
if (!is.null(wd)) setwd(wd) # set current directory as working directory; wd should have the two input files and this rscript
dir.create("dge_salmon.results.svp") # create a new results folder 
list.files(wd)

## Install and load required libraries (addition of tximport, readr, and biomaRt)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("tibble")
library(tibble)

#BiocManager::install("BiocParallel")
library(BiocParallel)
register(MulticoreParam(4)) # register the number of nodes on multicoreParam for DESeq2 run

#if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)

#if (!require(RColorBrewer)) install.packages('RColorBrewer')
library(RColorBrewer)

#BiocManager::install("DOSE")
library(DOSE)

#BiocManager::install("enrichplot")
library(enrichplot)

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE, force = TRUE)
library(organism, character.only = TRUE)

#install.packages("ggnewscale")
library(ggnewscale)

#install.packages("ggupset")
library(ggupset)

#BiocManager::install("pathview")
library(pathview)

#BiocManager::install("clusterProfiler")
library(clusterProfiler)

#BiocManager::install("AnnotationHub")
library(AnnotationHub)

#BiocManager::install("ensembldb")
library(ensembldb)

install.packages("devtools")
devtools::install_github("stephenturner/annotables")

## Perform OverExpression Analysis (ORA) using Clusterprofiler and Gene Ontology

## clusterProfiler does not work as easily using gene names, so we will turn gene names into Ensembl IDs using 
## clusterProfiler::bitr and merge the IDs back with the DE results
keytypes(org.Hs.eg.db)

## do this for the AM comparision
ids_AM <- bitr(rownames(dge.AM.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_AM <- which(duplicated(ids_AM$SYMBOL) == FALSE)

ids_AM <- ids_AM[non_duplicates_AM, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_AM <- merge(x=dge.AM.vp, y=ids_AM, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
AM_sigOE <- subset(merged_gene_ids_AM, padj < 0.05)

AM_sigOE_genes <- as.character(AM_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
AM_allOE_genes <- as.character(merged_gene_ids_AM$ENSEMBL)

## Run GO enrichment analysis 
ego_AM <- enrichGO(gene = AM_sigOE_genes, 
                   universe = AM_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_AM <- data.frame(ego_AM)

write.csv(cluster_summary_AM, "/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/enrichment_analyses/clusterProfiler_AM.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_AM, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
AM_OE_foldchanges <- AM_sigOE$log2FoldChange

names(AM_OE_foldchanges) <- AM_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_AM, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AM_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
AM_OE_foldchanges <- ifelse(AM_OE_foldchanges > 2, 2, AM_OE_foldchanges)
AM_OE_foldchanges <- ifelse(AM_OE_foldchanges < -2, -2, AM_OE_foldchanges)

cnetplot(ego_AM, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AM_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the AS comparision
ids_AS <- bitr(rownames(dge.AS.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_AS <- which(duplicated(ids_AS$SYMBOL) == FALSE)

ids_AS <- ids_AS[non_duplicates_AS, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_AS <- merge(x=dge.AS.vp, y=ids_AS, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
AS_sigOE <- subset(merged_gene_ids_AS, padj < 0.05)

AS_sigOE_genes <- as.character(AS_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
AS_allOE_genes <- as.character(merged_gene_ids_AS$ENSEMBL)

## Run GO enrichment analysis 
ego_AS <- enrichGO(gene = AS_sigOE_genes, 
                   universe = AS_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_AS <- data.frame(ego_AS)

write.csv(cluster_summary_AS, "/enrichment_analyses/clusterProfiler_AS.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_AS, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
AS_OE_foldchanges <- AS_sigOE$log2FoldChange

names(AS_OE_foldchanges) <- AS_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_AS, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AS_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
AS_OE_foldchanges <- ifelse(AS_OE_foldchanges > 2, 2, AS_OE_foldchanges)
AS_OE_foldchanges <- ifelse(AS_OE_foldchanges < -2, -2, AS_OE_foldchanges)

cnetplot(ego_AS, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AS_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the AC comparision
ids_AC <- bitr(rownames(dge.AC.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_AC <- which(duplicated(ids_AC$SYMBOL) == FALSE)

ids_AC <- ids_AC[non_duplicates_AC, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_AC <- merge(x=dge.AC.vp, y=ids_AC, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
AC_sigOE <- subset(merged_gene_ids_AC, padj < 0.05)

AC_sigOE_genes <- as.character(AC_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
AC_allOE_genes <- as.character(merged_gene_ids_AC$ENSEMBL)

## Run GO enrichment analysis 
ego_AC <- enrichGO(gene = AC_sigOE_genes, 
                   universe = AC_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_AC <- data.frame(ego_AC)

write.csv(cluster_summary_AC, "/enrichment_analyses/clusterProfiler_AC.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_AC, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
AC_OE_foldchanges <- AC_sigOE$log2FoldChange

names(AC_OE_foldchanges) <- AC_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_AC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AC_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
AC_OE_foldchanges <- ifelse(AC_OE_foldchanges > 2, 2, AC_OE_foldchanges)
AC_OE_foldchanges <- ifelse(AC_OE_foldchanges < -2, -2, AC_OE_foldchanges)

cnetplot(ego_AC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=AC_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the MS comparision
ids_MS <- bitr(rownames(dge.MS.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_MS <- which(duplicated(ids_MS$SYMBOL) == FALSE)

ids_MS <- ids_MS[non_duplicates_MS, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_MS <- merge(x=dge.MS.vp, y=ids_MS, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
MS_sigOE <- subset(merged_gene_ids_MS, padj < 0.05)

MS_sigOE_genes <- as.character(MS_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
MS_allOE_genes <- as.character(merged_gene_ids_MS$ENSEMBL)

## Run GO enrichment analysis 
ego_MS <- enrichGO(gene = MS_sigOE_genes, 
                   universe = MS_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_MS <- data.frame(ego_MS)

write.csv(cluster_summary_MS, "/enrichment_analyses/clusterProfiler_MS.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_MS, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
MS_OE_foldchanges <- MS_sigOE$log2FoldChange

names(MS_OE_foldchanges) <- MS_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_MS, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=MS_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
MS_OE_foldchanges <- ifelse(MS_OE_foldchanges > 2, 2, MS_OE_foldchanges)
MS_OE_foldchanges <- ifelse(MS_OE_foldchanges < -2, -2, MS_OE_foldchanges)

cnetplot(ego_MS, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=MS_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the MC comparision
ids_MC <- bitr(rownames(dge.MC.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_MC <- which(duplicated(ids_MC$SYMBOL) == FALSE)

ids_MC <- ids_MC[non_duplicates_MC, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_MC <- merge(x=dge.MC.vp, y=ids_MC, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
MC_sigOE <- subset(merged_gene_ids_MC, padj < 0.05)

MC_sigOE_genes <- as.character(MC_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
MC_allOE_genes <- as.character(merged_gene_ids_MC$ENSEMBL)

## Run GO enrichment analysis 
ego_MC <- enrichGO(gene = MC_sigOE_genes, 
                   universe = MC_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_MC <- data.frame(ego_MC)

write.csv(cluster_summary_MC, "/enrichment_analyses/clusterProfiler_MC.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_MC, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
MC_OE_foldchanges <- MC_sigOE$log2FoldChange

names(MC_OE_foldchanges) <- MC_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_MC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=MC_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
MC_OE_foldchanges <- ifelse(MC_OE_foldchanges > 2, 2, MC_OE_foldchanges)
MC_OE_foldchanges <- ifelse(MC_OE_foldchanges < -2, -2, MC_OE_foldchanges)

cnetplot(ego_MC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=MC_OE_foldchanges, 
         vertex.label.font=6)

## Repeat for the SC comparison
ids_SC <- bitr(rownames(dge.SC.vp), 
               fromType = "SYMBOL", 
               toType = c("ENSEMBL", "ENTREZID"), 
               OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_SC <- which(duplicated(ids_SC$SYMBOL) == FALSE)

ids_SC <- ids_SC[non_duplicates_SC, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_SC <- merge(x=dge.SC.vp, y=ids_SC, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
SC_sigOE <- subset(merged_gene_ids_SC, padj < 0.05)

SC_sigOE_genes <- as.character(SC_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
SC_allOE_genes <- as.character(merged_gene_ids_SC$ENSEMBL)

## Run GO enrichment analysis 
ego_SC <- enrichGO(gene = SC_sigOE_genes, 
                   universe = SC_allOE_genes, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_SC <- data.frame(ego_SC)

write.csv(cluster_summary_SC, "/enrichment_analyses/clusterProfiler_SC.csv")

## Dotplot gives top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not padj.
dotplot(ego_SC, showCategory=50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
SC_OE_foldchanges <- SC_sigOE$log2FoldChange

names(SC_OE_foldchanges) <- SC_sigOE$Row.names

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_SC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=SC_OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a minimum and maximum fold change value
SC_OE_foldchanges <- ifelse(SC_OE_foldchanges > 2, 2, SC_OE_foldchanges)
SC_OE_foldchanges <- ifelse(SC_OE_foldchanges < -2, -2, SC_OE_foldchanges)

cnetplot(ego_SC, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=SC_OE_foldchanges, 
         vertex.label.font=6)

## Perform Gene Set Enrichment Analysis using ClusterProfiler

## Start with the asymptomatic v mild analysis - prepare input
df.am <- read.csv("asymptomatic_v_mild_dge.csv")
colnames(df.am)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.am <- df.am[-which(df.am$padj > 0.05), ]

# convert column 1 into rownames
#df.am.rownames <- df.am %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.am.gsea <- df.am[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
df.am.gsea <- df.am.gsea[-3614,]

# feature 1: numeric vector
geneList = df.am.gsea[,2]

# feature 2: named vector
names(geneList) = as.character(df.am.gsea[,1])

# omit NA values
geneList_am <- na.omit(geneList_am)

# sort the list in decreasing order (required for clusterProfiler)
sort.am_gene_list = sort(geneList_am, decreasing = TRUE)
#sort.am_gene_list = order(gene_list$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_am, decreasing = TRUE), "geneList_GSEAsorted_AM.csv")

# perform gene set enrichment analysis using gseGOvi - This step can take upto 2-3 hours
gse_am <- gseGO(geneList= sort.am_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_am
# visualize the results using a dotplot
dotplot(gse_am, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_am, showCategory = 10)

## Repeat for asymptomatic v severe analysis - prepare input
df.as <- read.csv("asymptomatic_v_severe_dge.csv")
colnames(df.as)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.as <- df.as[-which(df.as$padj > 0.05), ]

# convert column 1 into rownames
#df.am.rownames <- df.am %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.as.gsea <- df.as[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.as.gsea <- df.as.gsea[-3614,]

# feature 1: numeric vector
geneList_as = df.as.gsea[,2]

# feature 2: named vector
names(geneList_as) = as.character(df.as.gsea[,1])

# omit NA values
geneList_as <- na.omit(geneList_as)

# sort the list in decreasing order (required for clusterProfiler)
sort.as_gene_list = sort(geneList_as, decreasing = TRUE)
#sort.as_gene_list = order(gene_list_as$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_as, decreasing = TRUE), "geneList_GSEAsorted_AS.csv")

# perform gene set enrichment analysis using gseGOvi
gse_as <- gseGO(geneList= sort.as_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_as
# visualize the results using a dotplot
dotplot(gse_as, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_as, showCategory = 10)

## Repeat for asymptomatic v critical analysis - prepare input
df.ac <- read.csv("asymptomatic_v_critical_dge.csv")
colnames(df.ac)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.ac <- df.ac[-which(df.ac$padj > 0.05), ]

# convert column 1 into rownames
#df.ac.rownames <- df.ac %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.ac.gsea <- df.ac[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.ac.gsea <- df.ac.gsea[-3614,]

# feature 1: numeric vector
geneList_ac = df.ac.gsea[,2]

# feature 2: named vector
names(geneList_ac) = as.character(df.ac.gsea[,1])

# omit NA values
geneList_ac <- na.omit(geneList_ac)

# sort the list in decreasing order (required for clusterProfiler)
sort.ac_gene_list = sort(geneList_ac, decreasing = TRUE)
#sort.ac_gene_list = order(gene_list_ac$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_ac, decreasing = TRUE), "geneList_GSEAsorted_AC.csv")

# perform gene set enrichment analysis using gseGOvi
gse_ac <- gseGO(geneList= sort.ac_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_ac
# visualize the results using a dotplot
dotplot(gse_ac, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_ac, showCategory = 10)

## Repeat for mild v critical analysis - prepare input
df.mc <- read.csv("mild_v_critical_dge.csv")
colnames(df.mc)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.mc <- df.mc[-which(df.mc$padj > 0.05), ]

# convert column 1 into rownames
#df.mc.rownames <- df.mc %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.mc.gsea <- df.mc[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.mc.gsea <- df.mc.gsea[-3614,]

# feature 1: numeric vector
geneList_mc = df.mc.gsea[,2]

# feature 2: named vector
names(geneList_mc) = as.character(df.mc.gsea[,1])

# omit NA values
geneList_mc <- na.omit(geneList_mc)

# sort the list in decreasing order (required for clusterProfiler)
sort.mc_gene_list = sort(geneList_mc, decreasing = TRUE)
#sort.mc_gene_list = order(gene_list_mc$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_mc, decreasing = TRUE), "geneList_GSEAsorted_MC.csv")

# perform gene set enrichment analysis using gseGOvi
gse_mc <- gseGO(geneList= sort.mc_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_mc
# visualize the results using a dotplot
dotplot(gse_mc, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_mc, showCategory = 10)

## Repeat for mild v severe analysis - prepare input
df.ms <- read.csv("mild_v_severe_dge.csv")
colnames(df.ms)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.ms <- df.ms[-which(df.ms$padj > 0.05), ]

# convert column 1 into rownames
#df.ms.rownames <- df.ms %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.ms.gsea <- df.ms[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.ms.gsea <- df.ms.gsea[-3614,]

# feature 1: numeric vector
geneList_ms = df.ms.gsea[,2]

# feature 2: named vector
names(geneList_ms) = as.character(df.ms.gsea[,1])

# omit NA values
geneList_ms <- na.omit(geneList_ms)

# sort the list in decreasing order (required for clusterProfiler)
sort.ms_gene_list = sort(geneList_ms, decreasing = TRUE)
#sort.ms_gene_list = order(gene_list_ms$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_ms, decreasing = TRUE), "geneList_GSEAsorted_MS.csv")

# perform gene set enrichment analysis using gseGOvi
gse_ms <- gseGO(geneList= sort.ms_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_ms
# visualize the results using a dotplot
dotplot(gse_ms, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

## Repeat for severe v critical analysis - prepare input
df.sc <- read.csv("severe_v_critical_dge.csv")
colnames(df.sc)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.sc <- df.sc[-which(df.sc$padj > 0.05), ]

# convert column 1 into rownames
#df.sc.rownames <- df.sc %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

# keep only the symbol and the L2FC column
df.sc.gsea <- df.sc[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
#df.sc.gsea <- df.sc.gsea[-3614,]

# feature 1: numeric vector
geneList_sc = df.sc.gsea[,2]

# feature 2: named vector
names(geneList_sc) = as.character(df.sc.gsea[,1])

# omit NA values
geneList_sc <- na.omit(geneList_sc)

# sort the list in decreasing order (required for clusterProfiler)
sort.sc_gene_list = sort(geneList_sc, decreasing = TRUE)
#sort.sc_gene_list = order(gene_list_sc$log2FoldChange, decreasing = TRUE)

write.csv(sort(geneList_sc, decreasing = TRUE), "geneList_GSEAsorted_SC.csv")

# perform gene set enrichment analysis using gseGOvi
gse_sc <- gseGO(geneList= sort.sc_gene_list, 
                ont = "ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE,
                pAdjustMethod = "none", 
                OrgDb = "org.Hs.eg.db", 
                nPermSimple = 1000)
gse_sc
# visualize the results using a dotplot
dotplot(gse_sc, showCategory = 15, split=".sign") +
  facet_grid(.~.sign)

# visualize using the encrichment map
emapplot(gse_sc, showCategory = 10)

## METHOD USING ENTREZIDS
# convert geneIDs to ENTREZIDs
#hs <- org.Hs.eg.db
#my.symbols.am <- c(names(am_gene_list))
#gene_data_am <- AnnotationDbi::select(hs, 
#keys = my.symbols,
#columns = c("ENTREZID", "SYMBOL"), 
#keytype = "SYMBOL")

# merge dataframes
#merge_df.am <- merge(df.am.gsea, gene_data_am, by="SYMBOL")


# USING ENSEMBL
#ids_AM <- bitr(rownames(dge.AM.vp), 
# fromType = "SYMBOL", 
#toType = c("ENSEMBL", "ENTREZID"), 
#OrgDb = "org.Hs.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
#non_duplicates_AM <- which(duplicated(ids_AM$SYMBOL) == FALSE)

#ids_AM <- ids_AM[non_duplicates_AM, ] 

## Merge the Ensembl IDs with the results     
#merged_gene_ids_AM <- merge(x=dge.AM.vp, y=ids_AM, by.x="row.names", by.y="SYMBOL")             

## Extract significant results
#AM_sigOE <- subset(merged_gene_ids_AM, padj < 0.05)

#AM_sigOE_genes <- as.character(AM_sigOE$ENSEMBL)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
#AM_allOE_genes <- as.character(merged_gene_ids_AM$ENSEMBL)

# extract the log 2 fold change column
#am_gene_list <- AM_sigOE$log2FoldChange

# omit NA values
#am_gene_list <- na.omit(am_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
#sort.am_gene_list = sort(am_gene_list, decreasing = TRUE)

# perform gene set enrichment analysis
#gse_am <- gseGO(geneList= sort.am_gene_list, 
#ont = "ALL", 
#keyType = "ENSEMBL", 
#nPerm = 10000,
#minGSSize = 3, 
#maxGSSize = 800, 
#pvalueCutoff = 0.05, 
#pAdjustMethod = "BH", 
#OrgDb = org.Hs.eg.db)

# visualize the results using a dotplot

## Output results from GO analysis to a table
#cluster_summary_AM <- data.frame(ego_AM)

#write.csv(cluster_summary_AM, "/Users/ShwetaPipaliya/Documents/PostDoc/projects/viral_transcriptomics/analyses/sars_cov_2/PRJEB43380/enrichment_analyses/clusterProfiler_AM.csv")


#### CODE BELOW IS NULL = REQUIRES ANNOTATION HUB WHICH I HAVENT FIGURED OUT JUST YET ####
## Visualize enriched terms from the GSE list of DGE genes using DOSE (Disease Ontology Semantic and Enrichment Analysis)
## Explore the grch38 table loaded by the annotables library
#ah = AnnotationHub()
#ah

## Return the IDs for the gene symbols in the DE results
#idx <- grch38$symbol %in% rownames(AM)

#ids <- grch38[idx, ]

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
#non_duplicates <- which(duplicated(ids$symbol) == FALSE)

#ids <- ids[non_duplicates, ] 

## Merge the IDs with the results 
#res_ids <- inner_join(res_tableOE_tb, ids, by=c("gene"="symbol")) 

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
#allOE_genes <- as.character(res_ids$ensgene)

## Extract significant results
#sigOE <- filter(res_ids, padj < 2.9e-06) # use genes that meet the genome wide significance threshold

#sigOE_genes <- as.character(sigOE$ensgene)

### END ###
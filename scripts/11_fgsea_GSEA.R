### (c) Shweta V. Pipaliya 
### 2022-04-09 
## This script performs gene set enrichment analysis using fgsea package
## Workfllow based on Stephen Turner's tutorial: https://stephenturner.github.io/deseq-to-fgsea/
### How to run this script from command-line: Rscript --vanilla 11_fgsea_GSEA.R
### Bioconductor packages in this script require atleast R v.4.2.0

## Set seed 
set.seed(777)

## Set current working directory to output DGE summary results and plots
wd <- getwd() 
if (!is.null(wd)) setwd(wd) # set current directory as working directory; wd should have the two input files and this rscript
dir.create("fgsea.results.svp") # create a new results folder 
list.files(wd)

## Install and load required libraries 
#if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

#if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)

#organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE, force = TRUE)
#library(organism, character.only = TRUE)
library(org.Hs.eg.db)

#if (!require(tidyverse)) install.packages('fgsea')
library(fgsea)

#install.packages("msigdbr")
library(msigdbr)

## Start with the asymptomatic v severe analysis - prepare input
# load in csv
df.as <- read.csv("asymptomatic_v_severe_dge.csv")
colnames(df.as)[1] = "SYMBOL" #rename column 1 as SYMBOL

# remove all values that do not have a padj of 0.05
df.as <- df.as[-which(df.as$padj > 0.05), ]

# get the gene symbol and test statistic. Also remove the NAs and deal with multiple test statistics
#df.as.1 <- df.as %>% 
  #dplyr::select(SYMBOL, stat) %>% 
  #na.omit() %>% 
  #distinct() %>% 
  #group_by(SYMBOL) %>%
  #summarize(stat=mean(stat))
#df.as.1

# instead of using LFC, we will use the stat values. So here we will just extract the SYMBOl and the stat valyes and remove the NA
df.as <- df.as[c("SYMBOL","log2FoldChange")]

# remove row with the empty name 
df.as.gsea <- df.am.gsea[-2202,]

# Create a named vector of the test statistic using the deframe function
ranks <- deframe(df.as)
head(ranks, 20)

# Use the Hallmark gene Set from MSigDBB to summarize and represent specific well-defined biological states and processes
# the gmtPathways() function will take a GMT  file downloaded from MSigDB and turn it into a list

#pathways.hallmark <- gmtPathways(system.file("extdata", "")

# retrieve human H (hallmark gene sets)
msigdbr_df <- msigdbr(species = "Homo sapiens", category = "H")
head(msigdbr_df)

# fixing format to work with fgsea
pathwaysH = split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name, )

# run fgsea enrichment
fgseaRes <- fgsea(pathways=pathwaysH, stats=ranks)

# Tidy the results
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

# plot the normalized enrichment scores and colour based on if or not the pathway was significant
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

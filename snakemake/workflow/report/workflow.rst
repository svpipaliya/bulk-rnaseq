This workflow performs differential expression analysis on single- or paired-end RNA-seq data. After adapter removal with Cutadapt, reads were mapped and gene counts were generated with STAR. Gene counts of replicated were summed up. Integrated normalization and differential expression analysis was conducted with DESeq2 following standard procedure as outlined in the manual.


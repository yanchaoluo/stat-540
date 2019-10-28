Paper Review by Yanchao Luo
================
20 February, 2019

[Gene Expression by Mouse Inner Ear Hair Cells during Development](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4405555/)

Goals and findings of the paper
-------------------------------

The study is on the gene-expression during mouse development of inner ear hair cells. They developed an enzymatic treatment to dissociate cells of the sensory epithelia and then use FACS method to purify hair cells (HCs). Surrounding cells (SCs) refer to the other cells in the sensory epithelium. Then, they compared gene-expression between hair cells (HCs) and surrounding cells (SCs) to study genes that cause deafness. They find hereditary deafness genes highly expressed in hair cells than surrounding cells. In addition, they conclude that understanding the difference gene-expression between hair cells (HCs) and surrounding cells (SCs) could help find the cause of deafness.

Data description
----------------

RNA-Seq was used to generate mRNAs data. This data could be found in the Shared Harvard Inner Ear Laboratory Database. The data contains information including cell types (HCs and SCs), tissues sources (cochlea and utricle), developmental stages (E16, P0, P4, and P7). In addition, information such as RT-PCR, Quantitative PCR, in situ hybridization, Immunocytochemistry were also collected. They used a transgenic mouse strain expressing GFP and also used strain expressing tdTomato to check results. To purify HCs, cells from cochleae and utricles were separated using enzymatic dissociation followed by FACS. Transcriptomes of the purified HCs and SCs were then studied by HTS. Each of 16 samples was sequenced to a depth of 40–100 million reads.

Analysis/methods
----------------

**Principal Component Analysis (PCA)**

Normalized gene expression values were assembled into a matrix with rows of different sample types and columns of genes. Principal component analysis (PCA) was performed across the 16 samples using the statistics package in R to understand the primary determinants of difference in genes-expression. The result showed that the first component (PC1; accounting for 22% of the variance) is highly associated with the cell type. The second component (PC2; accounting for 17% of the variance) is related to the developmental stage approximately (P4, P7 at the top and P0 at the bottom). PCA also used to validate the reproducibility of results across samples. Samples from the same cell type tended to cluster with each other in the PCA plot.

**Enrichment of genes**

Fold change (FC) is a ratio of GFP+ to GFP- across organs and ages. They were used FC to identify gene enrichment. The result showed that among the 20,207 annotated mouse genes, only 2008 genes were not expressed. 5430 were found to be preferentially expressed in HCs, 3230 were found in SCs and 9539 were expressed in both.

Further analysis was to identify the enrichments in cochlear HCs and utricular HCs. The result showed that some of the genes in functions unique to OHCs are associated with cochlear HCs. Components of microtubule-based cilia are correlated with utricular HCs.

Moreover, they used a histogram to show the enrichment of different developmental stages. They also applied DAVID to visualize how enrichment changes with age.

**Hierarchical clustering**

The dChip V2010.01 was used to visualize the expression data by hierarchically clustering genes. The hierarchical clustering heat maps showed subsets of genes enriched in SCs, in all HCs, in cochlear HCs, or in utricular HCs. Firstly, they computed log transformations of the normalized expression values. These expression values were standardized for each gene (scaled to have mean 0 and SD 1). Secondly, they defined the distance between two genes as a correlation of the expression of the 3′-tags of two genes across samples. Then, they computed the distance between a cluster and a cluster (centroid linkage method). The result showed that genes critical for mechanotransduction are expressed in HCs but not SCs.

**Deafness genes**

They ranked the 18,199 genes by enrichment level in HCs and showed that Genes expressed in HCs are a rich pool of deafness-gene candidates.

Conclusion
----------

This is a study focused on gene expression during mouse development. They conclude that understanding the difference gene-expression between hair cells (HCs) and surrounding cells (SCs) could help find the cause of deafness. The main methods of this study include multiple testing, PCA and Hierarchical clustering. One limitation is that there is no much information about multiple testing. Furthermore, they could construct a multiple linear regression of deafness based on types, tissues sources, developmental stages to confirm the result.

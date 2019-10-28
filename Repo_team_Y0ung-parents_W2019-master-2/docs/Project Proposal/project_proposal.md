# Project proposal - Y0ung-parents

## Motivation and background work

The most common neurodevelopmental disorders, Intellectual Disability (ID) and Autsim Spectrum Disorder (ASD), arise from three causative factors: genes, epigenetics and the environment[1]. Large-scale whole exome sequencing studies have found there is no single gene, but a collection of rare variants distributed across many genes that confer the manifestation of ASD[2],[3]. The variants occur in regions of multiple general transcription factors, which lend to alterations in global levels of gene expression regulation. Expression profiles have been determined with the use of RNA sequencing data from blood samples to develop predictive risk assessments of ASD in children. Previous studies have targeted peripheral blood lymphocytes[4] and whole blood[5] to generate transcriptome signatures, yet there has been no investigation as to whether one model may successfully be applied to a second, independent dataset. We are interested in comparing the global gene expression profiles of children with ASD versus a control group from these two studies, by recreating the predictive model using the data of the whole blood study and applying machine learning to evaluate the robustness of its true predictive power. To further substantiate the veracity of the models, we will use the metadata to evaluate global expression signatures across age. Accounting for this variable is necessary as gene regulation in individuals with ASD has been shown to fluctuate over maturation[6]. The validation of predictive risk models is fundamental to the refinement of ASD diagnostic tools and may lend to the development of gene-targeted treatments.


## Division of labour 

| Name | Background | Affiliation | Responsibilites | 
| ------------- | ------------- | ------------- | ------------- |
| Benson Chang | Biotechnology, Biological Engineering | Genome Science and Technology |  Analysis literature validation, Data visualization, Statistical analysis, Data interpretation |
| Virginia Pichler |  Genomic Epidemiology| Microbiology and Immunology | Literature search, Project planning, Data inspection/Quality control |
| Yanchao Luo | Statistical methods and Empirical Likelihood | Department of Statistics | Data cleaning & processing, Statistical analysis (hypothesis testing), Data visualization | 
| Anthony Jakob | Biological Engineering | Biomedical Engineering | Machine learning model |

## Dataset

### General description and characteristics of the data

We are working with two different datasets, both data from microarrays for expression profiling.
In both datasets, RNA was amplified, labeled and hybridized, the source of the RNA being peripheral blood lymphocytes for the first dataset and peripheral blood cells for the second dataset.
In both datasets, the subjects are children suffering from Autism Spectrum Disorder (ASD) or not.

In both datasets, the samples assign a measured expression value to the corresponding gene reference ID.
In the first sample, the measured value is the Robust Microarray Analysis (RMA) signal intensity.
In the second sample, the measured value is the Probe Logarithmic Intensity Error (PLIER) signal intensity.
The samples from both datasets each contain 54613 rows and the same gene reference ids, as the same kits were used.


### Technology used to generate the data

#### Dataset 1 [7]: 
Expression profiling by array
- Use of Qiagen Qiaquick kit on blood draw
- Double round amplification, followed by biotin-labelling using Affymetrix's GeneCHip Two-Cyle Target Labeling kit
- Checks for evenly distributed range of transcript size, verification of fragmentation
- Hybridization cocktails hybridized on Affymetrix Human Genome U133 Plus 2.0 Array (in situ oligonucleotide)
- Washing and scanning on Affymetrix's GeneChip Scanner 3000 7G
- Extraction of raw signal intensity from scanned images of the array in two batches (on two machines).
- No adjustment made for possible batch effects as the algorithms alter the gene expression distribution.
- Analysis of covariance of batch numbers.
- Analysis of gene expression distribution with MAS 5.0 (no alteration by algorithm)
- Analysis of gene expression with RMA (uses quantile normalization and my remove group level differences in gene expression distribution), looking for specific gene expression differences between groups.

#### Dataset 2 [8]:
Expression profiling by array
- Trizol extraction of total RNA according to manufacturer
- Generation of biotin-labeled cRNA according to Affymetrix protocols
- Quantification (A260) and fragmentation of cRNA
- Hybridization of fragmented cRNA on GeneChips (Affymetrix Human Genome U133 Plus 2.0 Array)
- Scanning using Affymetrix GeneChip scanner 3000 at 2.5 microm resultion
- Recording of PLIER of the samples as signal intensity

### Description of data

| Data | Dataset 1 			  | Dataset 2 			   |
| ---- | -------------------- | ---------------------- |
| Size | 663.4 MB 			  | 1.2 GB				   |
| Number of samples | 146	  | 285					   |
| Number of genes / samples | 54'613 | 54'612		   |
| Data in sample | 2 columns (ID_REF and VALUE) | 2 columns (ID_REF and VALUE) |
| Age of child | Yes		  | Yes 				   |
| Age of parents | Yes		  | No 					   |
| Gender | No		  		  | Yes 				   |
| Ethnicity | No			  | Yes					   |
| Specific diagnosis | No	  | Yes 				   |
| Other diseases | No		  | Yes					   |
| Scan batch number | Yes     | Not needed			   |

## Aims and methodology

We are interested in the differential gene expression of patients with/without autism. We are also interested in how the gene expression profiles of patients with/without autism differ at certain age groups. we have found two datasets from two different studies that explore highly similar variables (gene expression of 54613 genes, age, control vs. autism).

The first dataset is limited in information compared to the second dataset, so the bulk of our primary analysis will be performed on the second dataset, while the first dataset can be used to verify the predictive modelling generated from the second dataset. We will use normalization methods to ensure our comparative analysis of both datasets is consistent. The primary study from the second dataset has suggested that certain genes are significantly correlated with disorders associated with autism, we will analyze and verify whether this set of genes is also observed in the sample patients with autism from the first dataset. We will conduct principal component analysis (PCA) on gene expression profiles in the second dataset to identify variability as well as to explore possible batch effect of the data since second dataset contains information on sample batches. We will apply clustering in the gene expression samples to see if the identified groups correspond to control patients and patients with autism.

Both datasets did not look at the significance of sample patient age affecting gene expression profiles, even though the second dataset acknowledged that patient age at blood draw significantly influenced different gene expression levels of certain genes. We will rebuild a prediction model based on the second dataset, but also factor in subgrouping samples into different age groups. Statistical methods involved in building our new prediction model will include linear regression and logistic regression (because our response will only have two values - control vs autism). 

Building our new prediction model will also involve machine learning tools to help us verify the results of our prior clustering analysis. The first aim in this domain will be to replicate the model and findings of the second study, and further aims include the utilization of different machine learning models (random forests, PCA, non-negative matrix factorization, sparse models, deep neural networks, clustering, ...) and feature selection to find groups of genes correlated with ASD or not, potentially revealing unknown interactions between genes associated with ASD, and comparing the findings with the previous statistical analyses. Our prediction model will also be verified by predicting control vs autism of the samples in the first dataset to test the sensitivity and accuracy of our model.

## References

[1] [Srivastava, Anand K., and Charles E. Schwartz. "Intellectual disability and autism spectrum disorders: causal genes and molecular mechanisms." Neuroscience & Biobehavioral Reviews 46 (2014): 161-174.](https://www.sciencedirect.com/science/article/pii/S0149763414000773)

[2] [Anney, Richard, Lambertus Klei, Dalila Pinto, Joana Almeida, Elena Bacchelli, Gillian Baird, Nadia Bolshakova et al. "Individual common variants exert weak effects on the risk for autism spectrum disorders." Human molecular genetics 21, no. 21 (2012): 4781-4792.](https://academic.oup.com/hmg/article/21/21/4781/781997)

[3] [Liu, Li, Aniko Sabo, Benjamin M. Neale, Uma Nagaswamy, Christine Stevens, Elaine Lim, Corneliu A. Bodea et al. "Analysis of rare, exonic variation amongst subjects with autism spectrum disorders and population controls." PLoS genetics 9, no. 4 (2013): e1003443.](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003443)

[4] [Alter, Mark D., Rutwik Kharkar, Keri E. Ramsey, David W. Craig, Raun D. Melmed, Theresa A. Grebe, R. Curtis Bay et al. "Autism and increased paternal age related changes in global levels of gene expression regulation." PloS one 6, no. 2 (2011): e16715.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3040743/)

[5] [Kong, Sek Won, Christin D. Collins, Yuko Shimizu-Motohashi, Ingrid A. Holm, Malcolm G. Campbell, In-Hee Lee, Stephanie J. Brewster et al. "Characteristics and predictive value of blood transcriptome signature in males with autism spectrum disorders." PLoS One 7, no. 12 (2012): e49475.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0049475)

[6] [Chow, Maggie L., Tiziano Pramparo, Mary E. Winn, Cynthia Carter Barnes, Hai-Ri Li, Lauren Weiss, Jian-Bing Fan et al. "Age-dependent brain gene expression and copy number anomalies in autism suggest distinct pathological processes at young versus mature ages." PLoS genetics 8, no. 3 (2012): e1002592.​](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002592​)

[7] [Dataset 1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25507)

[8] [Dataset 2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18123)

[9] [Cluster Analysis for Gene Expression Data](https://cse.buffalo.edu/DBGROUP/bioinformatics/papers/survey.pdf)

[10] [Principal Component Analysis for clustering gene expression data](http://faculty.washington.edu/kayee/pca/bioinformaticspca.pdf)

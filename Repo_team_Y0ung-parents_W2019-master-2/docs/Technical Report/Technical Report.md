Technical Report
================
Introduction
---------------------------
Autism Spectrum Disorder (ASD) large-scale whole exome sequencing studies have found there is no single gene, but a collection of rare variants distributed across many genes that confer its manifestation [1],[2]. The variants occur in regions of multiple general transcription factors, which lend to alterations in global levels of gene expression regulation. Expression profiles have been generated with RNA sequencing data from blood samples to develop predictive risk assessments of ASD in children. Previous studies have targeted peripheral blood lymphocytes (leukocyte type) [3] and whole peripheral blood (erythrocytes, thrombocytes, and leukocytes) [4] to generate transcriptome signatures, yet there has been little investigation as to whether a predictive model designed on one dataset may successfully be applied to a second, independent dataset. We are interested in comparing the global gene expression profiles of children with ASD versus a control group from these two studies, by generating our own predictive model using the combined data of both studies through machine learning methods to evaluate the robustness of its predictive power. To further substantiate the veracity of the model, we will use the metadata to evaluate global expression signatures across age, batch and diagnosis (autism vs control) and compare it to linear regression results of these variables.  Age was selected as a variable as gene regulation in individuals with ASD has been shown to fluctuate over development and maturation [5]. The validation of predictive risk models is fundamental to the refinement of ASD diagnostic tools and may lend to the development of gene-targeted treatments.


Data Cleaning
---------------------------

First, we combined data from two publicly available datasets from the Gene Expression Omnibus (GEO) database (GEO accession numbers GSE18123 and GSE25507). The second dataset was subdivided into two separate datasets, which were generated using different sequencing platforms in the study. No rubric was provided to merge these datasets by IDs successfully. For this reason, we selected the first of the two datasets (GSE25507), which uses the same ID references as our first dataset (GSE18123). We choose the interested variables such as age, batch, diagnosis from the two metadata.

We convert diagnosis to a categorical variable which has two levels(autism and control).

GSE25507 has batch information(batch 1 and batch 2); however, GSE18123 has no information about batch, so we assign "none" for the GSE18123 metadata. We want to use PCA or cluster analysis to identify whether there is a batch effect of our combined data.

There are some missing values (15 samples) for our age variable. We consider these missing points might correlate with other variables(such as gender or batch). This satisfies the condition that the missing value is missing at random. Then we perform the multiple imputations with five repetitions for our age variable. The two pre-processed GEO datasets were merged after these the step.

Normalization: We first took log2 transformation for both gene expression datasets. Since the scales of the two datasets are different, we decided that quantile normalization was appropriate for our data after considering a few normalization methods. The density plot of average gene expression value can be found in our analysis. For further details of this step, please refer to [preprocessing](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/code/preprocessing/Pre_processing.md).


Statistical Methods
-------------------

Afterward, we performed our exploratory analysis on the data and model the data with Limma. Several analytical approaches we applied include PCA, agglomerative hierarchical clustering, modeling with limma.

Principal Component Analysis (PCA): PCA was conducted on gene expression to identify the variability explained by factors available in the combined dataset. The first three principal components capture around 60% of the variability in our data. We also plotted the PC1 vs. PC2 for age, batch, and diagnosis separately. The results of these plots show that age, batch, and diagnosis are all not related to PC1 and PC2.

Hierarchical Clustering: The clustering method we applied was Average method, and our distance metric was “Euclidean”. We did so to determine if our variables such as age, batch, and diagnosis are well-defined clusters.

Linear regression: Linear regression was performed by using the “limma” package in R to identify the top differentially gene expression between control and autism cases. We identified 15 different genes, given a p-value cutoff of 0.1. Then we plotted the top 10 genes and used Hierarchical Clustering algorithms to cluster the top 15 genes that showed differential expression between control and autism cases.

Normality check: We need to check whether the residuals are normally distributed when we use linear regression. We use the Anderson-Darling and the Shapiro-Wilks test to check the residuals. We also us a Bonferroni correction on the significance threshold to control the number of false positives.

More details for our analysis can be found in [exploratory and limma analyses](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/code/exploratory_and_limma_analysis/exploratory_and_Limma_analysis.md

Statistical Result
------------------

We used Limma to identify the top differentially genes To do so we fit each gene with a linear regression model as follows:
**Y**<sub>gene</sub> = *β*<sub>0</sub> + *β*<sub>1</sub>**a****g****e** + *β*<sub>2</sub>**b****a****t****c****h** + *β*<sub>3</sub>**d****i****a****g****n****o****s****i****s** + *ϵ*

where,

**Y**<sub>gene</sub> represents the vector of expression levels for each gene.

*β*<sub>0</sub> is the intercept vector.

**a****g****e** is a continuous variable that indicates people age between 1 and 17.5. *β*<sub>1</sub> is the coefficient of the **a****g****e** .

**b****a****t****c****h** is a binary variable with 3 levels(batch 1, batch 2 and none). *β*<sub>2</sub> is the coefficient of the **b****a****t****c****h**.

**d****i****a****g****n****o****s****i****s** is a binary variable with two levels (autism and control). *β*<sub>3</sub> is the coefficient of the **d****i****a****g****n****o****s****i****s**.

*ϵ* is a variable of residuals for a certain gene.

We identified 15 different genes between control and autism cases (p-value cutoff = 0.01) by using the multiple linear regression. One of our hypotheses was that different gene expression would be detectable in comparing control and autism cases, concerning age and batch. Additionally, our PCA results show that age, diagnosis (control and autism) and batch are all not related to first two of the variability we are observing in our data.


Geneset Enrichment Analysis
------------------

A geneset enrichment analysis was employed to evaluate the expression values of the genes for GO term multifunctionality. Input data for this analysis was generated in the [preprocessing](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/code/preprocessing/Pre_processing.md) step, where the expression data is combined from both datasets (GEO accession numbers GSE18123 and GSE25507 (P1)), cleaned and normalised. To further prepare the data, a differential expression analysis was run by creating a design matrix with the metadata, lmfit() to fit a linear model, and ebayes() to compute logFC values and the probes IDs were annotated with their corresponding gene symbols. The input data was subsequently subsetted to gene symbol and logFC. 

Gene multifunctionality scores were inputted from the STAT540 Github repository and used to check for the degree of multifunctional bias. The method employed was a Spearman’s correlation, which outputted a coefficient value close to zero, indicating little to no multifunctional bias. 

Finally, a Precision-Recall method was utilised to calculate the p-value of GO terms after a multifunctionality correction. This enrichment method put emphasis on highly ranked genes, identifying enrichment that could potentially impact the weight of the predictor genes generated from the machine learning analysis. GO.xml file was accessed using ErmineR to place genes into sets.  As expected from the Spearman’s correlation coefficient, GO terms were only minimally corrected in the enrichment analysis, with the largest log p-value adjustment at 0.4714262. 

These results overall indicate a low level of multifunctional bias. Furthermore, the top ranked genesets for bias did not contain any of the predictor genes generated in the machine learning analysis. 

[Full geneset enrichment analysis](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/code/geneset_enrichment_analysis/Geneset_Enrichment_Analysis.md).


Machine Learning
------------------

Using Python and the library scikit-learn, three binary classification models to predict diagnosis based on gene expression value were fitted to the data. Input for this analysis was generated in the [preprocessing](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/code/preprocessing/Pre_processing.md) step, where the expression data is combined from both datasets (GEO accession numbers GSE18123 and GSE25507). The data was cleaned, normalized, log2-transformed and the gene expression values were set to mean 0 and standard deviation 1 for each gene. 70% of the total combined data was used for training and 30% was used for testing. The motivation behind this approach was to automatically let the model select the genes most relevant for prediction using an appropriate regularizer for sparse solutions. Age was included to the data as an additional feature.

Model 1: Logistic regression (`sklearn.linear_model.LogisticRegression`) with L1-penalty, inverse regularization strength C 0.1, maximum 100 iterations to train the model. The value for C was chosen using 4-fold cross-validation. The best model found uses 50 genes and achieves an average classification accuracy of 84% on the combined dataset. 

Model 2: Linear Support Vector Machine (SVM) trained with Stochastic Gradient Descent (`sklearn.linear_model.SGDClassifier`) with hinge loss, Elastic Net penalty with a L1-to-L2 ratio of 0.9, inverse regularization strength C 0.1, maximum 100 iterations of Stochastic Gradient Descent to train the model. The value for C was chosen using 4-fold cross-validation. The best model found uses 39 genes and achieves an average classification accuracy of 75% on the combined dataset.

Model 3: Random Forest (`sklearn.ensemble.RandomForestClassifier`) with 100 trees per forest and a maximum depth of 20 for every tree. The best model found uses 3949 genes and achieves and average classification accuracy of 85% on the combined dataset. A prior gene selection is necessary for random forests to use a reasonable amount of genes. Since the goal of this step was to come up with said gene selection, this model was not further investigated.

The models' selected genes and the top differentially expressed genes between control and autism cases identified by the linear regression were compared. 4 genes were found to overlap between model 1 and 2 (*GATA2, CXCR3, STATH, MMP27*), 3 genes between model 1 and the linear regression (*AC092718.4, PCM1, ORMDL1*) and 1 gene between model 2 and the linear regression (*HIST1H2BG*). Some of the identified genes have been associated with ASD in the literature, however most of them have a very broad range of functions in the human body. Age was not selected in any model.

[Full Machine Learning Analysis](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/code/ml/ML_results.md).


References
------------------

[1] Anney, Richard, Lambertus Klei, Dalila Pinto, Joana Almeida, Elena Bacchelli, Gillian Baird, Nadia Bolshakova et al. "Individual common variants exert weak effects on the risk for autism spectrum disorders." Human molecular genetics 21, no. 21 (2012): 4781-4792.

[2] Liu, Li, Aniko Sabo, Benjamin M. Neale, Uma Nagaswamy, Christine Stevens, Elaine Lim, Corneliu A. Bodea et al. "Analysis of rare, exonic variation amongst subjects with autism spectrum disorders and population controls." PLoS genetics 9, no. 4 (2013): e1003443.

[3] Alter, Mark D., Rutwik Kharkar, Keri E. Ramsey, David W. Craig, Raun D. Melmed, Theresa A. Grebe, R. Curtis Bay et al. "Autism and increased paternal age related changes in global levels of gene expression regulation." PloS one 6, no. 2 (2011): e16715.

[4] Kong, Sek Won, Christin D. Collins, Yuko Shimizu-Motohashi, Ingrid A. Holm, Malcolm G. Campbell, In-Hee Lee, Stephanie J. Brewster et al. "Characteristics and predictive value of blood transcriptome signature in males with autism spectrum disorders." PLoS One 7, no. 12 (2012): e49475.

[5] Chow, Maggie L., Tiziano Pramparo, Mary E. Winn, Cynthia Carter Barnes, Hai-Ri Li, Lauren Weiss, Jian-Bing Fan et al. "Age-dependent brain gene expression and copy number anomalies in autism suggest distinct pathological processes at young versus mature ages." PLoS genetics 8, no. 3 (2012): e1002592.


Datasets
------------------

GSE25507 Dataset

GSE18123 Dataset


Programmes and Packages
------------------

[cluster](https://cran.r-project.org/web/packages/cluster/cluster.pdf)   
[pvclust](https://cran.r-project.org/web/packages/pvclust/pvclust.pdf)   
[xtable](https://cran.r-project.org/web/packages/xtable/xtable.pdf)   
[limma](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf)   
[GEOquery](https://github.com/seandavi/GEOquery)   
[knitr](https://cran.r-project.org/web/packages/knitr/knitr.pdf)   
[pheatmap](https://cran.r-project.org/web/packages/pheatmap/pheatmap.pdf)   
[stringr](https://cran.r-project.org/web/packages/stringr/stringr.pdf)     
[ggplot2](https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf)   
[reshape2](https://cran.r-project.org/web/packages/reshape2/reshape2.pdf)   
[tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)    
[ermineR](https://github.com/PavlidisLab/ermineR)     
[devtools](https://cran.r-project.org/web/packages/devtools/index.html)   
[Python](https://www.python.org/)    
[NumPy](http://www.numpy.org/)    
[scikit-learn](https://scikit-learn.org/stable/index.html) 

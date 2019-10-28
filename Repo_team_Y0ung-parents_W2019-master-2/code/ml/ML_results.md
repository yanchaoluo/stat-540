

# Machine Learning Methods

The motivation behind using machine learning methods is to build a binary classifier that predicts PSD or not based on the gene expression profile of a patient. Training the model with an L1 or LASSO regularizer (penalize the use of more genes) induces that only "relevant" genes be included in the final model.

In all following models, 70% of the samples were used for training and the remaining 30% were used to validate the model.

The second study contains two datasets: one was used to train the model, and one to validate the model. Unfortunately, different platforms were used to acquire the gene expression profiles, leading to a different nomenclature of the gene names. As such, the second part of the second dataset could not be used in any analysis.

## Procedure
The expression data was pre-processed in R to center and normalize it.
The samples were used as examples in the model. The expression value of every gene was considered to be a feature of every example fed to the model.

## Logistic Regression (scikit-learn LogisticRegression)

### Implementation details
- 4-fold cross-validation using bootstrapping led to an optimal regularization constant C of 0.1 for L1-regularization.
- 100 iterations to train each model

### Results
Training on only the second dataset led to a model averaging 89% classification accuracy using only 40 genes. This model however performed poorly on the first dataset, averaging only 56% classification accuracy, and 66% on the combined datasets. This is a clear example of overfitting and thus the model is of little value for further investigations.

Training on the combined datasets led to a model averaging 84% classification accuracy using only 50 genes. Training and validating on the first dataset only led to a 73% classification accuracy, and 82% when only using the second dataset.

These 50 genes are:

| Name | 
| ------------- |
|FILNC1|
|PTGFR|
|RTN2|
|STXBP5-AS1|
|LINC02366|
|POLD3|
|MMP27|
|AC009063.1|
|UBR4|
|NUDT6|
|ITIH5|
|PIF1|
|AL162171.1|
|CXCR3|
|COX8C|
|MYO1G|
|OR2F1|
|ERAP2|
|PCM1|
|MTMR3|
|LINC00901|
|ORMDL1|
|LYZL1|
|ENGASE|
|SCG3|
|RBCK1|
|SLC25A53|
|DTNA|
|COX8C|
|ZFHX3|
|YES1|
|STATH|
|SYCE1|
|WDR48|
|LYZL2|
|GATA2|
|TTC26|
|AC092718.4|
|OR2F1|
|STATH|


## Linear SVM trained with Stochastic Gradient Descent (scikit-learn SGDClassifier)

### Implementation details
- Hinge loss (leads to linear SVM) and elastic net penalty, with a L1-to-L2-ratio of 0.9
- 4-fold cross-validation using bootstrapping led to an optimal regularization constant C of 0.1 for elastic net regularization
- Maximum 100 iterations of stochastic gradient descent

### Results
Training on only one dataset proved difficult, as the model selected huge numbers of relevant genes (around 30'000) to compensate for the limited amount of samples.

Training on the combined datasets however led to a model averaging 75% classification accuracy using only 39 genes. Training and validating on the first dataset only led to a 75% classification accuracy, and 77% when only using the second dataset.

Training a logistic regression model using only these 39 selected genes improved the average classification accuracy to 78%.

These 39 genes are:

| Name | 
| ------------- |
|GOLGB1|
|AC068620.2|
|SLC8A1-AS1|
|HLA-DPB2|
|TEK|
|WWOX|
|CCNL1|
|AL606752.1|
|HLA-DPB2|
|CMIP|
|HLA-DPB2|
|SCAF11|
|MIR4645|
|MMP27|
|ALB|
|SPSB2|
|JCAD|
|HLA-DPB2|
|AC136621.1|
|SLC6A2|
|NOP14-AS1|
|TREML4|
|HIST1H2BG|
|CXCR3|
|LOXL1-AS1|
|SERPINB9P1|
|SRSF2|
|STATH|
|D21S2088E|
|HLA-DPB2|
|ASPN|
|GATA2|
|NR2F2|
|HLA-DPB2|
|STATH|

## Random Forest (scikit-learn RandomForestClassifier)

### Results
Random forests do not automatically select features in the same sense as the previous models. If the trees used are deep enough, a random forest will eventually use all features. 

Restricting the depth of the trees, an average of 85% classification accuracy was achieved on the combined datasets using 3949 genes. The number of used genes is too large and thus further analyses using random forests were abandoned.

### Implementation details
- Maximal depth of 20 nodes for each tree
- 100 estimators (trees) in each forest

## Technicalities
- For any approach, 100 models of the same type were trained, and the best one was saved if it used less than 100 genes. This process was repeated multiple times.
- On single dataset models, only genes with expression values of more than 150 in all samples were used to build the model
- On combined dataset models, all genes were used to build the model

# Analysis

## Accuracy

The achieved classification accuracy of 75 - 80% improves on the paper's model's 68% accuracy, using even less genes (39 or 50 instead of 55). 

## Common genes

We compare the genes used in the models with the genes found to be statistically relevant:

| Log. reg. vs SGD | Log. reg. vs Stat | SGD vs Stat | All | 
| --- | ----------- | --------- | ----- |
|GATA2| AC092718.4	| HIST1H2BG	|  		|
|CXCR3| PCM1		| 			|		|
|STATH| ORMDL1		| 			|		|	
|MMP27| 			|			|		|


The list of the 13 statistically relevant genes is provided below for completeness:

| Statistical analysis |
| ------ |
|PNN|
|LUC7L|
|LUC7L3|
|TIGD1|
|MALAT1|
|MLLT6|
|TLL2|
|HIST1H2BG|
|PCM1|
|HNRNPA2B1|
|ORMDL1|
|GSAP|
|AC092718.4|

It is to note that in the statistical logisitic regression, we try to predict the gene expression profile based on the presence of the disease, while in the models, we try to predict the presence of the disease given the gene expression profile.

There is little overlap between the genes used by the different models, and no single gene was found in all three cases.

The selected genes seem uncorrelated to ASD at first glance: GATA2 is a transcription factor which can lead to a wide range of diseases if mutated, PCM1 is essential for the localization of centrosomal proteins and HIST1H2BG is a histone-coding gene.

## Age

The age of the patient was added to the data as another feature to train the models (i.e. it acts as if it were the expression value of an additional gene). Age was not selected as a relevant feature for prediction of ASD in any model.

# References
- [Python](https://www.python.org/)
- [numpy](http://www.numpy.org/), used to simplify numerical operations
- [scikit-learn](https://scikit-learn.org/stable/index.html), used as machine learning library

# Individual report of Yanchao Luo
## Summary of contributions of each group member
### Describes the tasks and contributions of each of your group members.

**Propose the project** 

We found the dataset together. Virginia proposed the initial idea that we want to use RNA sequencing data from blood samples to develop predictive risk assessments of ASD in children. Then, after exploring the data, I came up with two statistical questions. One was different gene expression would be detectable in comparing control and autism cases, concerning age and batch. Another one was that age is correlated with autism using logistical regression. Additionally, Virginia and Benson did a lot of research on the biological reason behind many parts of our analysis. Anthony considered utilizing a machine learning method to construct predict model.

**Clean the data**

I did this part for our project. We kept the useful information from two datasets. I have completed the following steps to process the data. 1. Convert the raw data to log2 transformation 2. Use quantile normalization to make two gene dataset have a similar distribution. 3. Scale the data with mean zero and variance 1 for Principal component analysis and Limma model. 4. Convert diagnosis to a categorical variable only with two levels(autism or control), and convert batch to a categorical variable with three levels. 5. Use multiple imputations to handle the missing value in the age variable.

**Statistical Analysis**

All of us handled this part. I did some statistical analysis including the PCA, Hierarchical Clustering, Linear regression (Limma) and normality check. Benson tried to do the pathway analysis, but the result showed no significant differentially expressed genes. Virginia did Geneset Enrichment Analysis to evaluate the expression values of the genes for GO term multifunctionality and compared the result with machine learning analysis. Anthony separated the data to 70% of the samples were used for training, and the remaining 30% were used to validate the model for the prediction part. He also used 4-fold cross-validation to optimal regularization for L1-regularization and then constructed a logistic regression.

**Machine Learning Analysis**

 Anthony used three binary classification models to predict diagnosis based on gene expression value were fitted to the data. Except for the logistic regression,  Linear Support Vector Machine and Random Forest were constructed in this part.
 
**Interpretation of Results:** 

We all worked on this together. I focus on the statistical result; Anthony focuses on the machine learning result. Virginia worked on the gene Enrichment result. Benson compared all the methods and found some meaningful reasons from a biological perspective, and he also explained the result from the pathway analysis.

**Github Repository organization**

Anthony, Virginia and I worked to ensure each of the folders had a README with a description of the contents. All of us to make sure links are available in the GitHub repo.

**Presentation**

Virginia conducted a literature review for the basis of motivation in analyzing global gene expression in individuals with autism, and also presented the Geneset Enrichment Analysis. I showed how we cleaned and merged the data and PCA result. Benson presented the clustering analysis result, the linear regression, the KEGG-pathway classification and then discussed the limitation for our project. Anthony gave the machine learning result.

### Do you think the task assignments were fairly assigned and appropriate given each member’s background and skills? If no, how would you change it?

Yes, I believe the task assignments were fair. I focus on the data process and statistics part; Virginia and Benson concentrate on the biological perspective; Anthony concentrates on the machine learning part.

## Your specific contributions and comments

### Explain what are your contributions to the project?

I played a role in the statistical analysis for this project. I cleaned the whole dataset, and the detail could be found on the above clean the data part.
I did the following steps for statistical analysis:
1. Use PCA analysis to identify the variable effect.
2. Perform Hierarchical Clustering to determine if our variables such as age, batch, and diagnosis are well-defined.  
3. Linear regression to identify the top differentially gene expression between control and autism cases. 
4. Normality check for residual.

I also participate in interpreting the result and presentation.
For the proposal, progress and final report,  I wrote the data cleaning, statistical method, and result of the statistical part.
For the Github Repository organization, I wrote some descriptions in the README file.
For the presentation, I presented the part of statistical analysis.

### What worked well and what did not? What was the most challenging or rewarding moment during your group project?

the challenging for this project list as follows,
1. We fail to use all the dataset, and we only work on the data which has the same ID references.
2. We know the cell type heterogeneity can cause the differential expression to the differences in the cell proportions between test and control. We look upon the cell type heterogeneity. However, we cannot control the proportion from peripheral blood and all peripheral blood.
3. Initially, I set the false discovery rate to be 0.05, and we found no differential gene expression. We use 0.1 instead of 0.05. Therefore the result may not highly statistical significance. 
4. We build the Principal component analysis on the data set after quantile normalization with mean zero and variance 1 to identify the batch effect. However, the result shows that there is no relation between the batch and the principal components.

### How did you find having members of different backgrounds for a scientific project?

It was great to have members from a different background. I learned a lot of biological knowledge from my teammate. I also really enjoyed learning some machine learning methods using Python.
I am looking forward to using this new knowledge in my future. 

### What have you learned from this group project?

Except learning the biological knowledge and some machine learning methods, this project also helps me to review the statistical techniques I learned from this course. It also provides a reliable way for me to handle the high dimensional data in the future. For this group project, It also tells me how to separate a complex job by part.

### Any other comments on how the group project could have been structured differently (i.e. delivery requirements, assessment)

In general, I think the group project is excellent for us to learn and practice the material in this course. I just have a few comments for the delivery of the group project.
We have a month to do the project and then write a process report. However, after getting some suggestions from TAs and instructor, we only have around two weeks to complete this project. It would be a benefit for us to give more time between the process report and the final report. Another suggestion is maybe we could provide a guest talk before the initial proposal, and then we could have a rough idea about how the group project will look like.
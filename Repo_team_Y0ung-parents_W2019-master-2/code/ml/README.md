# Machine Learning code

The results can be found [here](ML_results.md).

## Python files
- [extract.py](extract.py) and [extractR.py](extractR.py) are used to extract data from either the raw .txt or R into a processable form
- [logReg.py](logReg.py) fits a logistic regression model to the data
- [sgd.py](sgd.py) fits a model to the data with stochastic gradient descent
- [randomForest.py](randomForest.py) fits a random forest model to the data
- [analysis.py](analysis.py) evaluates a given model

## .pkl files
- [bestLogReg.pkl](bestLogReg.pkl) contains the names of the genes (in Affymetrix code name) that are used in the best achieved logistic regression model
- [bestSGD.pkl](bestSGD.pkl) contains the names of the genes (in Affymetrix code name) that are used in the best achieved stochastic gradient descent model

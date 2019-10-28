import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import scale
import argparse
import multiprocessing as mp

"""
Anthony Jakob, 11.02.2019

Perform analysis
"""

def load_dataset_np(filename):
    return np.load(filename)

def load_dataset_pickle(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

def classification_error(y, yhat):
    return np.mean(y!=yhat)

def classification_accuracy(y, yhat):
    return 1 - np.mean(y!=yhat)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--dataX', required=True)
    parser.add_argument('-y', '--dataY', required=True)

    io_args = parser.parse_args()
    filename = io_args.dataX  # data_X.npy
    filenameMeta = io_args.dataY  # data_y.pkl

    X = load_dataset_np(filename)
    dataset = load_dataset_pickle(filenameMeta)
    y, metadata, gene_names = dataset["y"], np.array(dataset["metadata"]), dataset["gene_names"]

    verbose = True
    logData = True
    isScaled = True

    if not logData:
        # Remove lowly expressed genes
        X_new = []
        genes_new = []
        for i in range(X.shape[1]):
            if all(gene >= 150 for gene in X[:,i]):
                X_new.append(X[:,i])
                genes_new.append(genes[i])
        
        X = np.array(X_new).transpose()
        genes = genes_new

    if not isScaled:
        # Center and scale data
        X = scale(X)

    N = len(y)
    split = int(N*7/10)

    nModels = 100
    best_acc = 0
    best_acc_nnz = 10000
    best_nnz = 10000
    best_nnz_acc = 0

    for i in range(nModels):
        if verbose: print("Iteration %d" % i)

        # shuffe data
        randomize = np.arange(N)
        np.random.shuffle(randomize)
        X, y, metadata = X[randomize], y[randomize], metadata[randomize]

        # split data
        Xtrain, ytrain = X[:split], y[:split]
        Xtest, ytest = X[split:], y[split:]

        # Fit model
        model = RandomForestClassifier(n_estimators=100, max_depth=20)
        model.fit(Xtrain, ytrain)

        # Test model
        acc = classification_accuracy(model.predict(Xtest), ytest)
        nnz = model.n_features_
        if acc > best_acc:
            best_acc = acc
            best_acc_nnz = nnz
        if best_nnz > nnz:
            best_nnz = nnz
            best_nnz_acc = acc

    # Print data
    print("Best validation accuracy: %.3f with %d non-zeros" % (best_acc, best_acc_nnz))
    print("Best non-zeros: %d with validation accuracy %.3f" % (best_nnz, best_nnz_acc))

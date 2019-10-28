import numpy as np
import pickle
from sklearn.linear_model import LogisticRegression, ElasticNet, SGDClassifier
from sklearn import svm
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

def cross_validate_async(X, y, C):
    n = len(y)
    n_blocks = 4
    len_blocks = int(n * 1.0 / n_blocks)

    # do crossvalidation leaving out one block at the time
    cumulative_error = 0
    for i in range(0, n, len_blocks):
        # select training examples
        indices_val = np.array([j >= i and j <= i + len_blocks - 1 for j in range(n)])
        indices_train = np.invert(indices_val)

        model = SGDClassifier(loss='hinge', penalty='elasticnet', alpha=C, l1_ratio=0.9, max_iter=100)
        model.fit(X[indices_train, :], y[indices_train])
        cumulative_error += classification_error(model.predict(X[indices_val,:]), y[indices_val])

    cumulative_error /= 1.0 * n_blocks
    return { 'error': cumulative_error, 'C': C }

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--dataX', required=True)
    parser.add_argument('-y', '--dataY', required=True)

    io_args = parser.parse_args()
    filename = io_args.dataX  # data_X.npy
    filenameMeta = io_args.dataY  # data_y.pkl
    
    verbose = True
    logData = True
    scaled = True
    doCrossVal = False

    # load data
    X = load_dataset_np(filename)
    dataset = load_dataset_pickle(filenameMeta)
    y = dataset["y"]
    metadata = np.array(dataset["metadata"])
    genes = dataset["gene_names"]  

    if not logData:
        # Remove lowly expressed genes
        X_new = []
        genes_new = []
        for i in range(X.shape[1]):
            # min threshold expression is 150 for all samples
            if all(gene >= 150 for gene in X[:,i]):
                X_new.append(X[:,i])
                genes_new.append(genes[i])
        
        X = np.array(X_new).transpose()
        genes = genes_new

    if not scaled:
        # Center and scale data
        X = scale(X)

    # Cross-validation to select optimal C
    if doCrossVal:
        # Shuffle data
        randomize = np.arange(len(y))
        np.random.shuffle(randomize)
        X, y, metadata = X[randomize], y[randomize], metadata[randomize]

        # Cross validate in parallel
        cores = max(mp.cpu_count() - 1, 1)  # ensure at least one process
        if verbose: print("Running on %i cores" % cores)

        with mp.Pool(processes=cores) as pool:
            Cs = [pool.apply_async(cross_validate_async, (X, y, 10 ** c,)) for c in range(-4, 4)]
            results = [result.get() for result in Cs]

            best_err = np.inf
            for result in results:
                if result['error'] < best_err:
                    best_err = result['error']
                    C = result['C']
    else:
        C = 0.1

    if verbose: print("Running with C = %.3f" % C)

    # Use 70% of the data to train and the rest to test
    N = len(y)
    split = int(N*7/10)

    nModels = 100
    best_acc = 0
    best_acc_nnz = 10000
    best_acc_genes = []
    best_nnz = 10000
    best_nnz_acc = 0
    best_nnz_genes = []

    for i in range(nModels):
        if verbose: print("Iteration %d" % i)

        # Shuffle data
        randomize = np.arange(len(y))
        np.random.shuffle(randomize)
        X, y, metadata = X[randomize], y[randomize], metadata[randomize]

        Xtrain, ytrain = X[:split], y[:split]
        Xtest, ytest = X[split:], y[split:]
        
        # Fit model
        model = SGDClassifier(loss='hinge', penalty='elasticnet', alpha=C, l1_ratio=0.9, max_iter=100)
        model.fit(Xtrain, ytrain)

        # Test model
        acc = classification_accuracy(model.predict(Xtest), ytest)
        nnz = (model.coef_ != 0).sum()
        if acc > best_acc:
            best_acc = acc
            best_acc_nnz = nnz
            best_acc_genes = (model.coef_ != 0)[0]
        if best_nnz > nnz:
            best_nnz = nnz
            best_nnz_acc = acc
            best_nnz_genes = (model.coef_ != 0)[0]

    # Display some stats
    print("Best validation accuracy: %.3f with %d non-zeros" % (best_acc, best_acc_nnz))
    print("Best non-zeros: %d with validation accuracy %.3f" % (best_nnz, best_nnz_acc))

    # Save model if it uses less than 100 genes
    if best_acc_nnz < 100:
        dt = np.array(genes)[best_acc_genes[:len(best_acc_genes)-1]]
        if best_acc_genes[len(best_acc_genes)-1]:
            dt.append('Age')
        with open('best_genes_sgd_' + str(best_acc_nnz) + '_' + str(best_acc) + '.pkl', 'wb') as f:
            pickle.dump(dt, f)

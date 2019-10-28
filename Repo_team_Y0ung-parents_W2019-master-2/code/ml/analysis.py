import numpy as np
import pickle
import argparse
from sklearn.linear_model import LogisticRegression, ElasticNet, SGDClassifier

def load_dataset_np(filename):
    return np.load(filename)

def load_dataset_pickle(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

def classification_error(y, yhat):
    return np.mean(y!=yhat)

def classification_accuracy(y, yhat):
    return 1 - np.mean(y!=yhat)

def intersect(a, b):
    matching = []
    for e in a:
        if e in b:
            matching.append(e)
    return e


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--dataX', required=True)
    parser.add_argument('-y', '--dataY', required=True)

    io_args = parser.parse_args()
    filename = io_args.dataX  # data_X.npy
    filenameMeta = io_args.dataY  # data_y.pkl
        
    # Load data
    X = load_dataset_np(filename)
    dataset = load_dataset_pickle(filenameMeta)
    y = dataset["y"]
    metadata = np.array(dataset["metadata"])
    genes = dataset["gene_names"]

    # Model to test
    modelGenes = load_dataset_pickle("bestLogReg.pkl")
    indices = []

    # Best genes from different analyses
    genesStat = ['209997_x_at', '223187_s_at', '225932_s_at', '223546_x_at', '226675_s_at', '212036_s_at', '1553099_at', '241792_x_at', '215843_s_at', '226334_s_at', '210387_at', '224567_x_at', '209710_at', '213142_x_at', '225628_s_at', '228559_at', '202933_s_at', '1557446_x_at', '219758_at', '208835_s_at', '234789_at', '213931_at', '220544_at', '225390_s_at', '225594_at', '239964_at', '1569482_at', '232931_at', '219594_at', '220183_s_at', '213364_s_at', '242389_at', '1566974_at', '208720_s_at', '46665_at', '220577_at', '1560684_x_at', '1555475_x_at', '214802_at', '225767_at']
    bestLogReg = load_dataset_pickle("bestLogReg.pkl")
    bestSGD = load_dataset_pickle("bestSGD.pkl")

    # Gene names (not in order)
    geneNamesStat = ['PNN', 'AHSA2P', 'CREBZF', 'RF02193', 'RBM39', 'TCL6', 'LUC7L', 'TTLL3', 'AL034370.1', 'LUC7L3', 'LUC7L3', 'LUC7L3', 'TIGD1', 'TREML3P', 'NUDT6', 'MALAT1', 'MALAT1', 'MLLT6', 'TSKS', 'SNX1', 'TLL2', 'NBEAP1', 'HIST1H2BG', 'MLLT6', 'EXOC7', 'PCM1', 'HNRNPA2B1', 'KLF13', 'ORMDL1', 'YES1', 'ID2', 'SEMA4C', 'GVINP1', 'NINJ2', 'GSAP', 'GATA2', 'TTC26', 'AC092718.4', 'SNRNP200', 'KLF13']
    geneNamesLogReg = ['FILNC1', 'PTGFR', 'RTN2', 'STXBP5-AS1', 'LINC02366', 'POLD3', 'MMP27', 'AC009063.1', 'UBR4', 'NUDT6', 'ITIH5', 'PIF1', 'AL162171.1', 'CXCR3', 'COX8C', 'MYO1G', 'OR2F1', 'ERAP2', 'PCM1', 'MTMR3', 'LINC00901', 'ORMDL1', 'LYZL1', 'ENGASE', 'SCG3', 'RBCK1', 'SLC25A53', 'DTNA', 'COX8C', 'ZFHX3', 'YES1', 'STATH', 'SYCE1', 'WDR48', 'LYZL2', 'GATA2', 'TTC26', 'AC092718.4', 'OR2F1', 'STATH']
    geneNamesSGD = ['GOLGB1', 'AC068620.2', 'SLC8A1-AS1', 'HLA-DPB2', 'TEK', 'WWOX', 'CCNL1', 'AL606752.1', 'HLA-DPB2', 'CMIP', 'HLA-DPB2', 'SCAF11', 'MIR4645', 'MMP27', 'ALB', 'SPSB2', 'JCAD', 'HLA-DPB2', 'AC136621.1', 'SLC6A2', 'NOP14-AS1', 'TREML4', 'HIST1H2BG', 'CXCR3', 'LOXL1-AS1', 'SERPINB9P1', 'SRSF2', 'STATH', 'D21S2088E', 'HLA-DPB2', 'ASPN', 'GATA2', 'NR2F2', 'HLA-DPB2', 'STATH']

    best = intersect(geneNamesLogReg, geneNamesSGD)
    matchingLogReg = intersect(geneNamesLogReg, geneNamesStat)
    matchingSGD = intersect(geneNamesSGD, geneNamesStat)
    matchingAll = intersect(intersect(geneNamesLogReg, geneNamesStat), geneNamesSGD)

    index = 0
    for gene in genes:
        if gene in modelGenes:
            indices.append(index)
        index += 1

    # Reduce X
    X = np.array([[example[idx] for idx in indices] for example in X])
    split = int(len(y)*7/10)
    n = 1000
    acc = 0.0
    nnz = 0.0

    bestAcc = 0
    bestCoef = []
    avgCoef = []

    # format for markdown
    # for g in modelGenes: print('|' + g.strip("'") + '|')

    for i in range(n):
        # Shuffle data
        randomize = np.arange(len(y))
        np.random.shuffle(randomize)
        X, y, metadata = X[randomize], y[randomize], metadata[randomize]

        # Split data
        Xtrain, ytrain = X[:split], y[:split]
        Xtest, ytest = X[split:], y[split:]
                
        # Train model
        model = LogisticRegression()
        model.fit(Xtrain, ytrain)

        # Test model
        acc_i = classification_accuracy(model.predict(Xtest), ytest) 
        acc += acc_i
        nnz += (model.coef_ != 0).sum()
        if len(avgCoef) < 1:
            avgCoef = model.coef_
        else:
            avgCoef += model.coef_

        # Save best model
        if acc_i > bestAcc:
            bestAcc = acc_i
            bestCoef = model.coef_

    avgCoef = avgCoef / n

    # Print stats
    print("validation accuracy of avg model: %.3f with 50 non-zeros" % classification_accuracy(((np.sign(Xtest@avgCoef.ravel().transpose())+1)/2).astype(int), ytest))
    print("avg validation accuracy: %.3f with %d non-zeros" % (acc/n, nnz/n))

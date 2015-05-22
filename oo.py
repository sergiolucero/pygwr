from openopt import NLP
import numpy as np
numAssets = 25
assets = ['Asset_' + str(currentAsset) for currentAsset in range(0, numAssets)]
def objectiveFunction(x, mu, sigma, aversion = 1):
    x_ = x.reshape(-1,1)
    return float(np.dot(x_.T,mu) - aversion/2 * np.dot(x_.T, np.dot(sigma, x_)))
def objectiveGradient(x, mu, sigma, aversion = 1):
    x_ = x.reshape(-1,1)
    return mu - aversion * np.dot(sigma, x_)
def nlConstraint(x, sigma, conVal = 0.02):
    x_ = x.reshape(-1, 1)
    return np.dot(x_.T, np.dot(sigma, x_)) - pow(conVal,2)/21
def nlConstraintGradient(x, sigma, conVal = 0.02):
    x_ = x.reshape(-1, 1)
    return 2 * np.dot(sigma, x_)
A = np.identity(numAssets)
b = np.tile(0.2, (numAssets,1))
Aeq = np.tile(1.0, (1, len(assets)))
beq = np.array([0.]).reshape(-1,1)
numOpt = 1000
x_init = np.tile(0., (len(assets), 1))
for currentOpt in range(1, numOpt):
    print 'Running opt ' + str(currentOpt)
    assetReturns = np.random.randn(1000, numAssets) * 0.5*0.05/pow(252,0.5)
    currentMu = assetReturns.mean(axis = 0).reshape(-1, 1)*21
    currentSigma = np.cov(assetReturns.T)*21
    problem = NLP(f = objectiveFunction, x0 = x_init, A = A, b = b, Aeq  = Aeq, beq = beq, \
    goal = 'max', df = objectiveGradient, ftol = 1e-5, scale = 1, maxIter = 1e5, iprint = 0,\
    c = nlConstraint, dc = nlConstraintGradient , contol = 1e-6, maxFunEvals = 1e6)
    problem.args.f = (currentMu, currentSigma)
    problem.args.c = (currentSigma,)
    problem.solve('ralg')
    if problem.stopcase == 1:
        x_init = problem.xf.reshape(-1, 1)
    else:
        print "Not feasible"
        break

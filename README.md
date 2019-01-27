# M-estimator-stable-tail-dependence-function-
Please see 'MEstmatorForTailDependenceInArbitraryDimension -EINMAHL et al 2012' for details, you can download this paper from 'https://projecteuclid.org/euclid.aos/1349196391'. All code devote to replicating the result in their paper.

1. Generate data with specific dependence structure (logistic distribution with Frechet marginal distribution, 2-factor model)

2. Define M-estimator which trys to minimize the distance between empirical stable tail dependence function and the theretical one by introducing an intermediary function. 

3. compare bias and variance of this estimator with Censured MLE, and find the best threshold and intemediary function.

/* 
* Calibration Factor Estimation
* ---------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description 
* ---------------------------------------------------
* This model samples the posterior probability distribution
* for the calibration factor between observed fluorescence 
* and fluorophore copy number. This model assumes that 
* the noise in measurement is negligible compared to the 
* observed fluorescence, making individual fluorescence 
* measurements behave as delta functions. 
*/
functions{
    /** 
    * Approximate the Binomial distirubution for continuous variables 
    * as a ratio of Gamma functions 
    * 
    * @param I1: Observed fluorescence of daughter cell 1. 
    * @param I2: Observed fluorescence of daughter cell 2.
    * @param alpha: Fluorescenc calibration factor in units of a.u. / molecule
    * @param N: Total number of measurements 
    **/
    real GammaApproxBinom_lpdf(vector I1, vector I2, real alpha, int N)
    {
        vector[N] n1 = I1 ./ alpha;
        vector[N] n2 = I2 ./ alpha;
        vector[N] ntot = n1 + n2;
        return -N * log(alpha) + sum(lgamma(ntot + 1) - 
                    lgamma(n1 + 1) - lgamma(n2 + 1) - 
                    ntot * log(2));
    } 
}
     
data {
    int<lower=0> N; // Number of data points
    vector<lower=0>[N] I1; // Observed fluorescence of daughter cell 1
    vector<lower=0>[N] I2; // Observed fluorescence of daughter cell 2
}


parameters {
    real<lower=1> alpha;
}

model {   
    alpha ~ normal(0, 500);
    I1 ~ GammaApproxBinom(I2, alpha, N);  
}

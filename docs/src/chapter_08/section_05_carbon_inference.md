## Parameter Estimation of DNA Binding Energies and Comparison Across Carbon Source 

In Chapter 4, we conclude that the biophysical parameters defining
the fold-change input-output function are unperturbed between different
carbon sources. This conclusion is reached primarily by comparing how
well the fold-change and the free energy shift $\Delta F$ is predicted
using the parameter values determined in glucose supported medium at
37$^\circ$ C. In this section, we redetermine the DNA binding energy for
each carbon source condition and test its ability to predict the
fold-change of the other conditions.

### Reparameterizing the Fold-Change Input-Output Function

As described previously, the fold-change in gene expression is defined by the
total repressor copy number $R$, the energetic difference between the active
and inactive states of the repressor $\Delta\varepsilon_{AI}$, and the
binding energy of the repressor to the DNA $\Delta\varepsilon_{RA}$. Using
fluorescence microscopy, we can directly measure the average repressor copy
number per cell, reducing the number of variable parameters to only the
energetic terms.

Estimating both $\Delta\varepsilon_{RA}$ and $\Delta\varepsilon_{AI}$
simultaneously is fraught with difficulty as the parameters are highly
degenerate [@razo-mejia2018]. We can avoid this degeneracy by
reparameterizing the fold-change input-output function as
$$
\text{fold-change} =\left( 1 +
    \frac{R}{N_{NS}}e^{-\beta\epsilon}\right)^{-1},
$${#eq:fc_reparam}
where $\epsilon$ is the effective energetic parameter 
$$
\epsilon = \Delta\varepsilon_{R} - k_BT \log\left(1 +
    e^{-\beta\Delta\varepsilon_{AI}}\right).
$${#eq:ep_reparam}
Thus, to further elucidate any changes to the parameter values due to
changing the carbon source, we can infer the best-fit value of $\epsilon$ for
each condition and explore how well it predicts the fold-change in other
conditions.

### Statistical Inference of $\epsilon$

We are interested in the probability distribution of the parameter
$\epsilon$ given knowledge of the repressor copy number $R$ and a
collection of fold-change measurements $\mathbf{fc}$. This quantity can be calculated
via by Bayes' theorem as
$$
g(\epsilon\,\vert\,R, \mathbf{fc}) = \frac{f(\mathbf{fc}\,\vert\, R,
    \epsilon)g(\epsilon)}{f(\mathbf{fc})},
$${#eq:epsilon_bayes} 
where $g$ and $f$ represent probability densities of parameters and data,
respectively. In this context, the denominator term $f(\mathbf{fc})$ serves
only as normalization constant and can be neglected. As is described in detail in Refs.
[@razo-mejia2018; @chure2019], we assume that a given set of fold-change
measurements are normally distributed about the theoretical value $\mu$
defined by the fold-change function. The likelihood can be mathematically
defined as 
$$
f(\mathbf{fc}\,\vert\,\epsilon, R) =
\frac{1}{\left(2\pi\sigma^2\right)^{N/2}}\prod\limits_i^N
\exp\left[\frac{-\left(\text{fc}_i - \mu(\epsilon,
R)\right)^2}{2\sigma^2}\right],
$${#eq:lin_fc_like}
where $N$ is the total number of fold-change measurements, and $\sigma$ is
the standard deviation of the observations about the true mean and is another
parameter that must be included in the estimation. As the fold-change in gene
expression in this work covers several orders of magnitude (from $\approx
10^{-3} - 10^{0}$), it is better to condition the parameters on the log
fold-change rather than linear scaling, translating @Eq:lin_fc_like to 
$$
f(\mathbf{fc}^*\,\vert\, \epsilon, R) =
\frac{1}{\left(2\pi\sigma\right)^{N/2}}\prod\limits_i^N\exp\left[-\frac{\left(\text{fc}^*_i
- \mu^*(\epsilon, R)\right)^2}{2\sigma^2}\right],
$${#eq:log_fc_like}
where $\mathbf{fc}^*$ and $\mu^*(\epsilon, R)$
are the transformations 
$$
\mathbf{fc}^* = \log(\mathbf{fc})
    \label{eq:log_trans_fc}$$ and $$\mu^*(\epsilon, R) = -\log\left(1 +
    \frac{R}{N_{NS}}e^{-\beta\epsilon}\right),
$${#eq:log_trans_mu}
respectively.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With a likelihood in hand, we can now turn
towards defining functional forms for the prior distributions $g(\sigma)$ and
$g(\epsilon)$. For these definitions, we can turn to those used in
@chure2019,
$$
g(\epsilon) \sim \text{Normal}(\mu=-12, \sigma=6)
$${#eq:chure2019_ep_prior}
and
$$
g(\sigma) \sim \text{HalfNormal}(\phi=0.1),
$${#eq:chure2019_sig_prior} 
where we introduce the shorthand notation of "Normal" and "HalfNormal."
Combining @Eq:log_fc_like and @Eq:chure2019_sig_prior -
@Eq:chure2019_ep_prior yields the complete posterior distribution for
estimating the DNA binding energy for each carbon source medium. The complete
posterior distribution was sampled using Markov chain Monte Carlo in the Stan
probabilistic programming language [@carpenter2017].

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The sampled posterior distributions for
$\epsilon$ and $\sigma$ for each carbon source condition are shown in
@Fig:carbon_cornerplot and are summarized in Table 8.1. The posterior
distributions of $\epsilon$ across the conditions are approximately equal
with highly overlapping 95\% credible regions. The predictive capacity of
each estimate of $\epsilon$ is shown in @Fig:pairwise_carbon where all
fold-change measurements fall upon the theoretical prediction regardless of
which carbon source that the parameter value was conditioned upon. With this
analysis, we can say with quantitative confidence that the biophysical
parameters are indifferent to the physiological changes resulting from
variation in carbon quality.

| **Growth Condition** | **Parameter** | **Value** |
|:-- | :--: | --: |
| Glucose, 37$^\circ$ C | $\epsilon$ | $-14.5^{+0.2}_{-0.3}\, k_BT$|
| | $\sigma$ | $0.3^{+0.1}_{-0.1}$ | 
| Glycerol, 37$^\circ$ C | $\epsilon$ | $-14.6_{-0.1}^{+0.1}\, k_BT$|
| | $\sigma$ | $0.39^{+0.10}_{-0.08}$ | 
| Acetate, 37$^\circ$ C | $\epsilon$ | $-14.1_{-0.3}^{+0.2}\, k_BT$ |
| | $\sigma$ | $0.35_{-0.1}^{+0.1}$ | 
: Summarized parameter estimates of $\epsilon$ and $\sigma$ given a single
growth condition. Reported as median and upper/lower bounds of 95\% credible
region.

![**Posterior probability distributions of effective DNA binding energy
$\epsilon$ and standard deviation $\sigma$.** Top plot shows the marginal
probability distribution of $\epsilon$ conditioned only on the purple (glucose),
glycerol (green), or acetate (brown) measurements. Bottom left plot shows the
joint probability distribution between the effective DNA binding energy
$\epsilon$ and the standard deviation $\sigma$. Bottom right shows the marginal
posterior distribution over the standard deviation
$\sigma$. The [Python code (`ch8_figS9.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS9.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch8_figS9){#fig:carbon_cornerplot short-caption="Posterior
probability distributions of effective DNA binding and standard deviation for
different carbon sources."}

![**Pairwise estimation and prediction of DNA binding energies.** Rows indicate
the strain to which the effective DNA binding energy $\epsilon$ was estimated and
columns are the strains whose fold-change is predicted. Shaded lines represent
the 95\% credible region of the prediction given the estimated value of
$\epsilon$. Points and error correspond to the median and standard error of five
to eight biological replicates. The [Python code (`ch8_figS10.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS10.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch8_figS10){#fig:pairwise_carbon
short-caption="Pairwise estimation and prediction of DNA binding energies
estimated from different carbon sources."}

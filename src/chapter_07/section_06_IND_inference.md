## Bayesian Parameter Estimation for Inducer Binding Domain Mutants 

In Chapter 3, we put forward two naïve hypotheses for which parameters of our
fold-change equation are affected by mutations in the inducer binding domain
of the repressor. The first hypothesis was that only the inducer dissociation
constants, $K_A$ and $K_I$, were perturbed from their wild-type values.
Another hypothesis was that the inducer dissociation constants were affected
in addition to the energetic difference between the active and inactive
states of the repressor, $\Delta\varepsilon_{AI}$.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In this section, we first derive the
statistical model for each hypothesis and then perform a series of diagnostic
tests that expose the inferential limitations of each model. With well
calibrated statistical models, we then apply each to an induction profile of
the inducer binding mutant Q291K and assess the validity of each hypothesis.
To understand the statistical models for each hypothesis, only the subsection
*Building A Generative Statistical Model* is necessary.

### Building a Generative Statistical Model 

For both hypotheses, we assume that the underlying physical model is the same
while a subset of the parameters are modified. As the fold-change
measurements for each biological replicate are statistically independent, we
can assume that they are normally distributed about the theoretical
fold-change value. Thus, for each model, we must include a parameter $\sigma$
which is the standard deviation of the distribution of fold-change
measurements. For the first hypothesis, in which only $K_A$ and $K_I$ are
changed, we are interested in sampling the posterior distribution 
$$
g(K_A,K_I, \sigma\,\vert\, y) \propto f(y\,\vert\, K_A, K_I,
\sigma)g(K_A)g(K_I)g(\sigma),
$${#eq:hyp1_generic_bayes}
where $y$ corresponds to the set of fold-change measurements. In the above
model, we have assumed that the priors for $K_A$ and $K_I$ are independent.
It is possible that it is more appropriate to assume that they are dependent
and that a single prior distribution captures both parameters, $g(K_A, K_I)$.
However, assigning this prior is more difficult and requires strong knowledge
*a priori* about the relationship between them. Therefore, we continue under
the assumption that the priors are independent.

The generic posterior given in @Eq:hyp1_generic_bayes can be extended to 
evaluate the second hypothesis in which $\Delta\varepsilon_{AI}$ is also modified,
$$
g(K_A, K_I, \Delta\varepsilon_{AI}, \sigma\,\vert\,y) \propto f(y\,\vert\, K_A, K_I, \Delta\varepsilon_{AI}, \sigma)g(K_A)g(K_I)g(\Delta\varepsilon_{AI})g(\sigma)
$${#eq:hyp2_generic_bayes}
where we have included $\Delta\varepsilon_{AI}$ as an estimated parameter and
assigned a prior distribution.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As we have assumed that the fold-change
measurements across replicates are independent and normally distributed, the
likelihoods for each hypothesis can be written as
$$
f(y\,\vert\, K_A, K_I, \sigma) \sim \text{Normal}\{\mu(K_A, K_I), \sigma\},
$${#eq:hyp1_likelihood} for the first hypothesis and
$$
f(y\,\vert\, K_A, K_I, \Delta\varepsilon_{AI}, \sigma) \sim \text{Normal}\{\mu(K_A, K_I, \Delta\varepsilon_{AI}), \sigma\},
$${#eq:hyp2_likelihood} for the second. Here, we have assigned
$\mu(\dots)$ as the mean of the normal distribution as a function of the
parameters defined by our fold-change equation.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With a likelihood distribution in hand, we now
turn toward assigning functional forms to each prior distribution. As we have
used in the previous sections of this chapter, we can assign a half-normal prior for
$\sigma$ with a standard deviation of 0.1, namely,
$$
g(\sigma)\sim \text{HalfNormal}\{0, 0.1\}.
$${#eq:IND_halfnorm_sigma} 
It is important to note that the inducer dissociation constants $K_A$ and $K_I$ are scale
invariant, meaning that a change from $0.1\,\mu$M to $1\,\mu$M yields a
decrease in affinity equal to a change from $10\,\mu$M to $100\,\mu$M.
As such, it is better to sample the dissociation constants on a
logarithmic scale. We can assign a log normal prior for each
dissociation constant as
$$
g(K_A) = {1 \over K_A\sqrt{2\pi\phi^2}}\exp\left[-{(\log {K_A \over 1\,\mu\text{M}} - \mu_{K_A})^2 \over 2\phi^2}\right],
$${#eq:lognormal_explicit} 
or with the short-hand notion of
$$
g(K_A) \sim \text{LogNormal}\{\mu_{K_A}, \phi\}
$${#eq:lognormal_shorthand}
For $K_A$, we assigned a mean $\mu_{K_A} = 2$ and a standard deviation
$\phi=2$. For $K_I$, we chose a mean of $\mu_{K_I} = 0$ and $\phi = 2$,
capturing our prior knowledge that $K_A > K_I$ for the wild-type LacI. While
the prior distributions are centered differently, they both show extensive
overlap, permitting mutations in which $K_A < K_I$. For
$\Delta\varepsilon_{AI}$, we assign a normal distribution of the prior
centered at 0 with a standard deviation of $5\, k_BT$,
$$
g(\Delta\varepsilon_{AI}) \sim \text{Normal}\{0, 5\}.
$${#eq:epAI_prior} 
This permits values of $\Delta\varepsilon_{AI}$ that are above or below zero, meaning that the
inactive state of the repressor can be either more or less energetically
favorable to the active state. A standard deviation of $5\,k_BT$ permits
a wide range of energies with $+5\,k_BT$ and $-5\,k_BT$ corresponding to
$\approx 99.5\%$ and $\approx 0.5\%$ of the repressors being active in
the absence of inducer, respectively.

### Prior predictive checks 

To ensure that these choices of prior distributions are appropriate, we
performed prior predictive checks for each hypothesis as previously
described in the second section of this chapter. We drew 1000 values from the prior
distributions shown in @Fig:ind_prior_predictive (A) for $K_A$, $K_I$, and
$\Delta\varepsilon_{AI}$. Using the draws from the $K_A$, and $K_I$
priors alone, we generated data sets of $\approx 70$ measurements. The
percentiles of the fold-change values drawn for the 1000 simulations is
shown in the top panel of @Fig:ind_prior_predictive (B).

It can be seen that in the absence of inducer, the fold-change values
are close to zero and are with distributed about the leakiness value due
to $\sigma$. This is in contrast to the data sets generated when
$\Delta\varepsilon_{AI}$ is permitted to vary along with $K_A$ and
$K_I$. In the bottom panel of @Fig:ind_prior_predictive (B), the fold-change when $c = 0$
can extend above 1.0 which is possible only when
$\Delta\varepsilon_{AI}$ is included, which sets what fraction of the
repressors is active. Under both hypotheses, the 99$^\text{th}$
percentile of the fold-change extends to just above 1 or just below 0,
which matches our intuition of how the data should behave. Given these
results, we are satisfied with these choices of priors and continue onto
the next level of calibration of our model.

![**Prior predictive checks for two hypotheses of inducer binding domain
mutants.** (A) Probability density functions for $K_A$, $K_I$,
$\Delta\varepsilon_{AI}$, and $\sigma$. Black points correspond to draws
from the distributions used for prior predictive checks. (B) Percentiles
of the simulated data sets using draws from the $K_A$ and $K_I$
distributions only (top, red bands) and using draws from $K_A$, $K_I$,
and $\Delta\varepsilon_{AI}$ (bottom, blue
bands).](ch7_figS13){#fig:ind_prior_predictive short-caption="Prior predictive
checks for two hypotheses of inducer binding domain mutants."}

### Simulation Based Calibration 

With an appropriate choice of priors, we turn to simulation based
calibration to root out any pathologies lurking in the model itself or
the implementation through MCMC. For each parameter under each model, we
compute the $z$-score and shrinkage of each inference, shown in @Fig:ind_sbc.
Under the first hypothesis in which $K_A$ and
$K_I$ are the only perturbed parameters [@Fig:ind_sbc(A)], we see all parameters have $z-$scores
clustered around 0, indicating that the value of the ground-truth is
being accurately estimated through the inference. While the shrinkage
for $\sigma$ is close to 1 (indicating the prior is being informed by
the data), the shrinkage for $K_A$ and $K_I$ is heavily tailed with some
values approaching zero. This is true for both statistical models,
indicating that for some values of $K_A$ and $K_I$, the parameters are
difficult to pin down with high certainty. In the application of these
models to data, this will be revealed as large credible regions in the
reported parameters. Under the second hypothesis in which all allosteric
parameters are allowed to change, we see moderate shrinkage for
$\Delta\varepsilon_{AI}$ [purple points in @Fig:ind_sbc(B)] with the minimum shrinkage being around
0.5. The samples resulting in low shrinkage correspond to values of
$\Delta\varepsilon_{AI}$ that are highly positive or highly negative, in
which small changes in the active fraction of repressors cannot be
accurately measured through our model. However, the median shrinkage for
$\Delta\varepsilon_{AI}$ is approximately 0.92, meaning that the the
data highly informed the prior distributions for the majority of the
inferences. The rank distributions for all parameters under each model
appear to be highly uniform, indicating that both statistical models are
computationally tractable.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With knowledge of the caveats of
estimating $K_A$ and $K_I$ for both models, we proceed with our analysis and
examine how accurately these models can capture the phenomenology of the data.

![**Simulation based calibration of statistical models for inducer binding
domain mutants.** (A) Sensitivity statistics and rank distribution for a
statistical model in which $K_A$ and $K_I$ are the only parameters
permitted to vary. (B) Sensitivity statistics and rank distribution for
a model in which all allosteric parameters $K_A$, $K_I$, and
$\Delta\varepsilon_{AI}$ are allowed to be modified by the mutation.
Gray envelope in the bottom plots correspond to the 99$^\text{th}$
percentile of variation expected from a true uniform
distribution.](ch7_figS14){#fig:ind_sbc short-caption="Simulation based
calibration of statistical models for inducer binding domain mutants."}

### Posterior Predictive Checks 

With a properly calibrated statistical model for each hypothesis, we now
apply it to a representative dataset. While each model was applied to
each inducer binding domain mutant, we only show the application to the
mutant Q291K with 260 repressors per cell paired with the native *lac*
operator O2.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The results from applying the statistical model in which only $K_A$ and
$K_I$ can change is shown in @Fig:kaki_post_predictive. The joint and marginal
distributions for each parameter [@Fig:kaki_post_predictive (A)] reveal a strong correlation
between $K_A$ and $K_I$ whereas all other parameters are symmetric and
independent. While the joint and marginal distributions look well
behaved, the percentiles of the posterior predictive checks
[@Fig:kaki_post_predictive (B)] are more suspect. While all
data falls within the 95$^\text{th}$ percentile, the overall trend of
the data is not well predicted. Furthermore, the percentiles expand far
below zero, indicating that the sampling of $\sigma$ is compensating for
the leakiness in the data being larger than it should be if only $K_A$
and $K_I$ were the changing parameters.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We see significant improvement when $\Delta\varepsilon_{AI}$ is
permitted to vary in addition to $K_A$ and $K_I$. @Fig:kaki_epai_post_predictive
(A) shows the joint and marginal distributions between all parameters from the MCMC sampling. We
still see correlation between $K_A$ and $K_I$, although it is not as
strong as in the case where they are the only parameters allowed to
change due to the mutation. We also see that the marginal distribution
for $\sigma$ has shrunk significantly compared to the marginal
distribution in @Fig:kaki_post_predictive (A). The percentiles of the
posterior predictive checks, shown in @Fig:kaki_epai_post_predictive (B)
are much more in line with the experimental measurements, with the 5$\text{th}$ percentile
following the data for the entire induction profile.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In this section we have presented two hypotheses for the minimal
parameter set needed to describe the inducer binding mutations, derived
a statistical model for each, thoroughly calibrated its behavior, and
applied it to a representative data set. The posterior predictive checks
(@Fig:kaki_post_predictive and @Fig:kaki_epai_post_predictive) help us understand which
hypothesis is more appropriate for that particular mutant. The
incredibly wide percentiles and significant change in the leakiness that
result from a model in which only $K_A$ and $K_I$ are perturbed suggests
that more than those two parameters should be changing. We see
significant improvement in the description of the data when
$\Delta\varepsilon_{AI}$ is altered, indicating that it is the more
appropriate hypothesis of the two.

![**Posterior predictive checks for inducer binding domain mutants where
only $K_A$ and $K_I$ are changed.** (A) MCMC sampling output for each
parameter. Joint distributions are colored by the value of the log
posterior with increasing probability corresponding to transition from
blue to yellow. (B) Percentiles of the data generated from the
likelihood distribution for each sample of $K_A$, $K_I$, and $\sigma$.
Overlaid points are the experimentally observed
measurements.](ch7_figS15){#fig:kaki_post_predictive short-caption="Posterior
predictive checks for inducer binding domain mutants where only $K_A$ and $K_I$
are changed."}


![**Posterior predictive checks for inducer binding domain mutants where
all allosteric parameters can change.** (A) MCMC sampling output for all
parameters. Joint distributions are colored by the value of the log
posterior with increasing probability corresponding to the transition
from blue to yellow. Marginal distributions are shown adjacent to each
joint distribution. (B) Percentiles of the data generated from the
likelihood for each sample of $K_A$, $K_I$, $\Delta\varepsilon_{AI}$,
and $\sigma$. The corresponding experimental data for Q291K are shown as
black open-faced circles.](ch7_figS16){#fig:kaki_epai_post_predictive
short-caption="Posterior predictive checks for inducer binding domain mutants
where all allosteric parameter can change."}

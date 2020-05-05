## Bayesian Parameter Estimation For DNA Binding Mutants 

In this section, we outline the statistical model used in this work to
estimate the DNA binding energy for a given mutation in the DNA binding
domain. The methodology presented here is similar to that performed in
Chapter 2 and outlined in accompanying Chapter 6. In the following text, we
take a very detailed approach to vetting the robustness of our statistical
inference machinery as determination of parameter values is critical to
assessing the effects of mutations. Similarly to what is presented in Chapter
6, we begin with a derivation of our statistical model using Bayes' theorem
and then perform a series of principled steps to validate our choices of
priors, ensure computational feasibility, and assess the validity of the
model given the collected data. This work follows the analysis pipeline
outlined by Michael Betancourt in his case-study entitled ["Towards A
Principled Bayesian Workflow."
](https://betanalpha.github.io/assets/case_studies/principled_bayesian_workflow.html)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The second subsection *Building a Generative
Statistical Model* lays out the statistical model used in this work to
estimate the DNA binding energy and the error term $\sigma$. The subsequent
subsections -- *Prior Predictive Checks*, *Simulation Based Calibration*, and
*Posterior Predictive Checks* -- define and summarize a series of
tests that ensure that the parameters of the statistical model can be
identified and are computationally tractable. To understand how we defined
our statistical model, only the second subsection is needed.

### Calculation of the Fold-Change in Gene Expression 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We appreciate the subtleties of the efficiency
of photon detection in the flow cytometer, fluorophore maturation and
folding, and autofluorescence correction, and we understand the importance in
modeling the effects that these processes have on the reported value of the
fold-change. However, in order to be consistent with the methods used in the
literature, we took a more simplistic approach to calculate the fold-change.
Given a set of fluorescence measurements of the constitutive expression
control ($R = 0$), an autofluorescence control (no YFP), and the experimental
strain ($R > 0$), we calculate the fold-change as
$$
\text{fold-change} = {\langle I_\text{cell}(R > 0)\rangle - \langle I_\text{autofluorescence}\rangle \over \langle I_\text{cell}(R = 0) \rangle - \langle I_\text{autofluorescence}\rangle}.
$${#eq:experimental_foldchange} 
It is important to note here that for a given biological replicate, we
consider only a point estimate of the mean fluorescence for each sample and
perform a simple subtraction to adjust for background fluorescence. For the
analysis going forward, all mentions of measured fold-change are determined
by this calculation.

### Building a Generative Statistical Model 

To identify the minimal parameter set affected by a mutation, we assume
that mutations in the DNA binding domain of the repressor alters only
the DNA binding energy $\Delta\varepsilon_{RA}$, while the other
parameters of the repressor are left unperturbed from their wild-type
values. As a first approach, we can assume that all of the other
parameters are known without error and can be taken as constants in our
physical model. Ultimately, we want to know how probable a particular
value of $\Delta\varepsilon_{RA}$ is given a set of experimental
measurements $y$. Bayes' theorem computes this distribution, termed the
*posterior distribution* as
$$
g(\Delta\varepsilon_{RA}\,\vert\, y) = {f(y\,\vert\,
 \Delta\varepsilon_{RA}) g(\Delta\varepsilon_{RA}) \over f(y)}
$${#eq:simple_bayes_DNA}
where we have used $g$ and $f$ to represent probability densities over
parameters and data, respectively. The expression $f(y\,\vert
\Delta\varepsilon_{RA})$ captures the likelihood of observing our data set
$y$ given a value for the DNA binding energy under our physical model. All
knowledge we have of what the DNA binding energy *could* be, while remaining
completely ignorant of the experimental measurements, is defined in
$g(\Delta\varepsilon_{RA})$, referred to as the *prior distribution*.
Finally, the likelihood that we would observe the data set $y$ while being
ignorant of our physical model is defined by the denominator $f(y)$. In this
work, this term serves only as a normalization factor and as a result will be
treated as a constant. We can therefore say that the posterior distribution
of $\Delta\varepsilon_{RA}$ is proportional to the joint distribution between
the likelihood and the prior,
$$
g(\Delta\varepsilon_{RA}\,\vert\, y) \propto f(y\,\vert\,\Delta\varepsilon_{RA})g(\Delta\varepsilon_{RA}).
$${#eq:epRA_post_prop}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We are now tasked with translating this
generic notation into a concrete functional form. Our physical model derived
in Chapter 2 given by computes the average fold-change in gene expression.
Speaking practically, we make several replicate measurements of the
fold-change to reduce the effects of random errors. As each replicate is
independent of the others, it is reasonable to expect that these measurements
will be normally distributed about the theoretical value of the fold-change
$\mu$, computed for a given $\Delta\varepsilon_{RA}$. We can write this
mathematically for each measurement as 
$$
f(y\,\vert\, \Delta\varepsilon_{RA}) = {1 \over
    (2\pi\sigma^2)^{N/2}}\prod\limits_{i}^N \exp\left[-(y_i -
    \mu(\Delta\varepsilon_{RA}))^2 \over 2\sigma^2\right],
$${#eq:likelihood_DNA_explicit} 
where $N$ is the number of measurements in $y$ and $y_i$ is the $i^\text{th}$
experimental fold-change measurement. We can write this likelihood in
shorthand as
$$
f(y\,\vert\,\Delta\varepsilon_{RA}) = \text{Normal}\{\mu(\Delta\varepsilon_{RA}),
\sigma\}\, 
$${#eq:likelihood_DNA_short}
which we will use for the remainder of this section.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Using a normal distribution for our likelihood
has introduced a new parameter $\sigma$ which describes the spread of our
measurements about the true value. We must therefore include it in our
parameter estimation and assign an appropriate prior distribution such that
the posterior distribution becomes
$$
g(\Delta\varepsilon_{RA}, \sigma\,\vert y) \propto f(y\,\vert\,
\Delta\varepsilon_{RA}, \sigma)g(\Delta\varepsilon_{RA})g(\sigma).
$${#eq:full_post_sig_prior} 
We are now tasked with assigning functional forms to the priors $g(\Delta\varepsilon_{RA})$ and
$g(\sigma)$. Though one hopes that the result of the inference is not
too dependent on the choice of prior, it is important to choose one that
is in agreement with our physical and physiological intuition of the
system.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We can impose physically reasonable bounds on the possible values of the
DNA binding energy $\Delta\varepsilon_{RA}$. We can say that it is
unlikely that any given mutation in the DNA binding domain will result
in an affinity greater than that of biotin to streptavidin ($1\,
\mathrm{fM} \approx -35\, \mathrm{k_BT}$, BNID 107139 [@milo2010]), one
of the strongest known non-covalent bonds. Similarly, it's unlikely that
a given mutation will result in a large, positive binding energy,
indicating non-specific binding is preferable to specific binding
($\sim 1$ to $10\,\mathrm{k_BT}$). While it is unlikely for the DNA
binding energy to exceed these bounds, it's not impossible, meaning we
should not impose these limits as hard boundaries. Rather, we can define
a weakly informative prior as a normal distribution with a mean and
standard deviation as the average of these bounds,
$$
g(\Delta\varepsilon_{RA}) \sim \text{Normal}\{-12, 12\}
$${#eq:epRA_prior}
whose probability density function in shown in @Fig:epRA_prior_pred (A).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;By definition, fold-change is restricted to the
bounds $[0, 1]$. Measurement noise and fluctuations in autofluorescence
background subtraction means that experimental measurements of fold-change
can extend beyond these bounds, though not substantially. By definition, the
scale parameter $\sigma$ must be positive and greater than zero. We also know
that for the measurements to be of any use, the error should be less than the
available range of fold-change, $1.0$. We can choose such a prior as a half
normal distribution
$$
g(\sigma) = {1 \over \phi}\sqrt{2 \over \pi}\exp\left[-{\sigma^2 \over
2\phi^2}\right];\, \forall\, \sigma \geq 0
$${#eq:sig_hn_explicit}
where $\phi$ is the standard deviation. By choosing $\phi = 0.1$, it is
unlikely that $\sigma \geq 1$ yet not impossible, permitting the occasional
measurement significantly outside of the theoretical bounds. The probability
density function for this prior is shown in @Fig:epRA_prior_pred (B).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;While these choices for the priors seem
reasonable, we can check their appropriateness by using them to simulate a
data set and checking that the hypothetical fold-change measurements obey our
physical and physiological intuition.

![**Prior distributions and prior predictive check for estimation of the DNA
binding energy.** (A) Prior probability density function for DNA binding energy
$\Delta\varepsilon_{RA}$ as $\sim$ Normal($-12, 12$). (B) Prior probability
density function for the standard deviation in measurement noise $\sigma \sim$
HalfNormal($0, 0.1$). (C) Percentiles of values drawn from the likelihood
distribution given draws from prior distributions given $R=260$, $K_A=139
\,\mu$M, $K_I=0.53\,\mu$M, and $\Delta\varepsilon_{AI} = 4.5\,k_BT$, which match
the parameters used in @razo-mejia2018. Black points at top of (A) and (B)
represent draws used to generate fold-change measurements from the likelihood
distribution. Percentiles in (C) generated from 800 draws from the prior
distributions. For each draw from the prior distributions, a data set of 70
measurements over 12 IPTG concentrations (ranging from $0$ to $5000\,\mu$M) were
generated from the likelihood. The [Python code (`ch7_figS2.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS2.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS2){#fig:epRA_prior_pred
short-caption="Prior distributions and prior predictive check for estimation of
the DNA binding energy."}

### Prior Predictive Checks 

If our choice of prior distribution for each parameter is appropriate,
we should be able to simulate data sets using these priors that match
our expectations. In essence, we would hope that these prior choices
would generate some data sets with fold-change measurements above 1 or
below zero, but they should be infrequent. If we end up getting
primarily negative values for fold-change, for example, then we can
surmise that there is something wrong in our definition of the prior
distribution. This method, coined a *prior predictive check*, was first
put forward in @good1950 and has received newfound attention in computational statistics.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We perform the simulation in the following
manner. We first draw a random value for $\Delta\varepsilon_{RA}$ out of its
prior distribution stated in @Eq:epRA_prior and calculate what the mean
fold-change should be given our physical model. With this in hand, we draw a
random value for $\sigma$ from its prior distribution, specified in
@Eq:sig_hn_explicit. We then generate a simulated data set by drawing
$\approx$ 70 fold-change values across twelve inducer concentrations from the
likelihood distribution which we defined in @Eq:likelihood_DNA_short. This
roughly matches the number of measurements made for each mutant in this work.
We repeat this procedure for 800 draws from the prior distributions, which is
enough to observe the occasional extreme fold-change value from the
likelihood. As the DNA binding energy is the only parameter of our physical
model that we are estimating, we had to choose values for the others. We kept
the values of the inducer binding constants $K_A$ and $K_I$ the same as the
wild-type repressor ($139\,\mu$M and $0.53\,\mu$M, respectively). We chose to
use $R=260$ repressors per cell as this is the repressor copy number we used
in the main text to estimate the DNA binding energies of the three mutants.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The draws from the priors are shown in
@Fig:epRA_prior_pred (A) and (B) as black points above the corresponding
distribution. To display the results, we computed the percentiles of the
simulated data sets at each inducer concentration. These percentiles are
shown as red shaded regions in @Fig:epRA_prior_pred (C). The 5$\text{th}$
percentile (dark purple band) has the characteristic profile of an induction
curve. Given that the prior distribution for $\Delta\varepsilon_{RA}$ is
centered at $-12\, k_BT$ and we chose $R = 260$, we expect the generated data
sets to cluster about the induction profile defined by these values. More
importantly, approximately 95\% of the generated data sets fall between
fold-change values of -0.1 and 1.1, which is within the realm of possibility
given the systematic and biological noise in our experiments. The
99$^\text{th}$ percentile maximum is approximately $1.3$ and the minimum
approximately $-0.3$. While we could tune our choice of prior further to
minimize draws this far from the theoretical bounds, we err on the side of
caution and accept these values as it is possible that fold-change
measurements this high or low can be observed, albeit rarely.
Through these prior predictive checks, we feel confident that these
choices of priors are appropriate for the parameters we wish to
estimate. We can now move forward and make sure that the statistical
model as a whole is valid and computationally tractable.

### Sensitivity Analysis and Simulation Based Calibration 

Satisfied with our choice of prior distributions, we can proceed to
check other properties of the statistical model and root out any
pathologies lurking in our model assumptions.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To build trust in our model, we could generate
a data set $\tilde{y}$ with a *known* value for $\sigma$ and
$\Delta\varepsilon_{RA}$, estimate the posterior distribution
$g(\Delta\varepsilon_{RA}, \sigma\,\vert\,\tilde{y})$, and determine how well
we were able to retrieve the true value of the parameters. However, running
this once or twice for handpicked values of $\sigma$ and
$\Delta\varepsilon_{RA}$ won't reveal edge-cases in which the inference
fails, some of which may exist in our data. Rather than performing this
operation once, we can run this process over a variety of data sets where the
ground truth parameter value is drawn from the prior distribution (as we did
for the prior predictive checks). For an arbitrary parameter $\theta$, the
joint distribution between the ground truth value $\tilde\theta$, the
inferred value $\theta$, and the simulated data set $\tilde y$ can be written
as
$$
\pi(\theta, \tilde{y}, \tilde{\theta}) =
g(\theta\,\vert\,\tilde{y})f(\tilde{y}\,\vert\,\tilde{\theta})g(\tilde{\theta}).
$${#eq:bayesian_joint_dist}
If this process is run for a large number of simulations,
@Eq:bayesian_joint_dist can be marginalized over all data
sets $\tilde{y}$ and all ground truth values $\tilde{\theta}$ to yield
the original prior distribution,
$$
\int d\tilde{\theta}\int d\tilde{y} \pi(\theta,\tilde{y}, \tilde{\theta}) =
g(\theta).
$${#eq:self_consistency_relation}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This result, described by @talts2018, holds true for *any* statistical model
and is a natural self consistency property of Bayesian inference. Any
deviation between the distribution of our inferred values for $\theta$ and
the original prior distribution $g(\theta)$ indicates that either our
statistical model is malformed or the computational method is not behaving as
expected. There are a variety of ways we can ensure that this condition is
satisfied, which we outline below.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Using the data set generated for the prior predictive checks (shown in
@Fig:epRA_prior_pred (C)), we sampled the posterior
distribution and compute $\Delta\varepsilon_{RA}$ and $\sigma$ for each
simulation and checked that they matched the original prior
distribution. To perform the inference, we use Markov chain Monte Carlo
(MCMC) to sample the posterior distribution. Specifically, we use the
Hamiltonian Monte Carlo algorithm implemented in the Stan probabilistic
programming language [@carpenter2017]. The specific code files can be
accessed through the [paper website](http://rpgroup.caltech.edu/mwc_mutants) or
the associated [GitHub repository](https://github.com/rpgroup-pboc/mwc_mutants).
The original prior distribution and the distribution of inferred parameter
values can be seen in @Fig:epRA_sbc (A) and (B). For both $\Delta\varepsilon_{RA}$
and $\sigma$, we can accurately recover the ground truth distribution
(purple) via sampling with MCMC (orange). For $\Delta\varepsilon_{RA}$, there
appears to be an upper and lower limit past which we are unable to
accurately infer the binding energy. This can be seen in both the
histogram [@Fig:epRA_sbc(A)] and the empirical cumulative distribution
[@Fig:epRA_sbc (B)] as deviations from the ground truth when
DNA binding is below $\approx -25\,k_BT$ or above $\approx -5\,k_BT$.
These limits hinder our ability to comment on exceptionally strong or
weak binding affinities. However, as all mutants queried in this work
exhibited binding energies between these limits, we surmise that the
inferential scheme permits us to draw conclusions about the inferred DNA
binding strengths.

![**Comparison of averaged posterior and prior distributions for
$\Delta\varepsilon_{RA}$ and $\sigma$.** (A) Distribution of the average
values for the DNA binding energy $\Delta\varepsilon_{RA}$ (orange)
overlaid with the ground truth distribution (purple). (B) Data averaged
posterior (orange) for the standard deviation of fold-change measurements
overlaid with the ground truth distribution (purple). Top and bottom show
the same data with different
visualizations. The [Python code (`ch7_figS3.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS3.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS3){#fig:epRA_sbc short-caption="Comparison of averaged
posterior and prior distributions for DNA binding energy and homoscedastic error."}


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Rather than examining the agreement of the
data-averaged posterior and the ground truth prior distribution solely by
eye, we can compute summary statistics using the mean $\mu$ and standard
deviation $\sigma$ of the posterior and prior distributions which permit
easier identification of pathologies in the inference. One such quantity is
the posterior $z$-score, which is defined as
$$
z = {\mu_\text{posterior} - \tilde\theta \over \sigma_\text{posterior}}.
$${#eq:zscore}
This statistic summarizes how accurately the posterior recovers the ground
truth value beyond simply reporting the mean, median, or mode of the
posterior distribution. $Z$-scores around $0$ indicate that the posterior is
concentrating tightly about the true value of the parameter whereas large
values (either positive or negative) indicate that the posterior is
concentrating elsewhere. A useful feature of this metric is that the width of
the posterior is also considered. It is possible that the posterior could
have a mean very close to the ground truth value, but have an incredibly
narrow distribution/spread such that it does not overlap with the
ground-truth. Only comparing the mean value to the ground truth would suggest
that the inference "worked". However with a small standard deviation
generates a very large $z$-score, telling us that something has gone awry.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If our inferential model is behaving properly, the width of the
posterior distribution should be significantly smaller than the width of
the prior, meaning that the posterior is being informed by the data. The
level to which the posterior is being informed by the data can be easily
calculated given knowledge of both the prior and posterior distribution.
This quantity, aptly named the shrinkage $s$, can be computed as
$$
s = 1 - {\sigma^2_\text{posterior} \over \sigma^2_\text{prior}}.
$${#eq:shrinkage}
When the shrinkage is close to zero, the variance of the posterior is
approximately the same as the variance of the prior, model is not being
properly informed by the data. When $s\approx 1$, the variance of the
posterior is much smaller than the variance of the prior, indicating that the
it is being highly informed by the data. A shrinkage less than 0 indicates
that the posterior is wider than the prior distribution, revealing a severe
pathology in either the model itself or the implementation.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In @Fig:epRA_sbc_sensitivity, we compute these summary
statistics for each parameter. For both $\Delta\varepsilon_{RA}$ and
$\sigma$, we see clustering of the $z$-score about 0 with the extrema
reaching $\approx \pm 3$. This suggests that for the vast majority of
our simulated data sets, the posterior distribution concentrated about
the ground truth value. We also see that for both parameters, the
posterior shrinkage $s$ is $\approx 1$, indicating that the posterior is
being highly informed by the data. There is a second distribution
centered $\approx$ 0.8 for $\Delta\varepsilon_{RA}$, indicating that for
a subset of the data sets, the posterior is only $\approx 80$\% narrower
than the prior distribution. These samples are those that were drawn
outside of the limits of $\approx -25$ to $-5\, k_BT$ where the
inferential power is limited. Nevertheless, the posterior still
significantly shrank, indicating that the data strongly informs the
posterior.

![**Inferential sensitivity for estimation of $\Delta\varepsilon_{RA}$ and
$\sigma$**. The posterior $z$-score for each posterior distribution
inferred from a simulated data set is plotted against the shrinkage for
(A) the DNA binding energy $\Delta\varepsilon_{RA}$ and (B) the standard
deviation of fold-change measurements
$\sigma$. The [Python code (`ch7_figS4.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS4.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch7_figS4){#fig:epRA_sbc_sensitivity short-caption="Inferential
sensitivity for estimation of the DNA binding energy and homoscedastic error."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The general self-consistency condition given by
@Eq:self_consistency_relation provides another route to
ensure that the model is computationally tractable. Say that we draw a
value for the DNA binding energy from the prior distribution, simulate a
data set, and sample the posterior using MCMC. The result of this
sampling is a collection of $N$ values of the parameter which may be
above, below, or equal to the ground-truth value. From this set of
values, we select $L$ of them and rank order them by their value. @talts2018 
derived a general theorem which states that the number of samples less than the ground truth value of the parameter
(termed the rank statistic) is uniformly distributed over the interval
$[0, L]$. As @Eq:self_consistency_relation *must* hold true for any
statistical model, deviations from uniformity signal that there is a
problem in the implementation of the statistical model. How the
distribution deviates is also informative as different types of failures
result in different distributions. The nature of these deviations, along
with a more formal proof of the uniform distribution of rank statistics
can be found in @talts2018 where it was originally derived.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Given the sampling statistics for each of the
simulated data sets, we took 800 of the MCMC samples of the posterior
distribution for each of the 800 simulated data sets and computed the rank
statistic. The distributions are shown in 
@Fig:epRA_sbc_rank as both histograms and ECDFs for the DNA
binding energy and standard deviation. The distribution of rank
statistics for both parameters appears to be uniform. The purple band
overlaying the histograms (top row) as well as the purple envelopes
overlaying the ECDFs (bottom row) represent the $99^\text{th}$
percentile expected from a true uniform distribution. The uniformity of
this distribution, along with the well-behaved $z$-scores and shrinkage
for each parameter, tells us that there are no underlying pathologies in
our statistical model and that it is computationally tractable. However,
this does not mean that it is correct. Whether this model is valid for
the actual observed data is the topic of the next section.

![**Rank distribution of the posterior samples from simulated data.** Top
row shows a histogram of the rank distribution with $n=20$ bins. Bottom
row is the cumulative distribution for the same data. Purple bands
correspond to the 99th percentile of expected variation from a uniform
distribution. (A) Distribution for the DNA binding energy
$\Delta\varepsilon_{RA}$ and (B) for the standard deviation
$\sigma$. The [Python code (`ch7_figS5.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS5.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch7_figS5){#fig:epRA_sbc_rank short-caption="Rank distribution of the
posterior samples from simulated data for the DNA binding energy and
homoscedastic error."}


### Parameter Estimation and Posterior Predictive Checks 

We now turn to applying our vetted statistical model to experimental
measurements. While the same statistical model was applied to all three
DNA binding mutants, here we only focus on the mutant Q18M for brevity.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Using a single induction profile, we sampled the posterior distribution
over both the DNA binding energy $\Delta\varepsilon_{RA}$ and the
standard deviation $\sigma$ using MCMC implemented in the Stan
programming language. The output of this process is a set of 4000
samples of both parameters along with the value of their log posterior
probabilities, which serves as an approximate measure of the probability
of each value. The individual samples are shown in @Fig:epRA_ppc. The joint distribution between
$\Delta\varepsilon_{RA}$ and $\sigma$ is shown in the lower left hand
corner, and the marginal distributions for each parameter are shown
above and to the right of the joint distribution, respectively. The
joint distribution is color coded by the value of the log posterior,
with bright orange and dark purple corresponding to high and low probability,
respectively. The symmetric shape of the joint distribution is a telling
sign that there is no correlation between two parameters. The marginal
distributions for each parameter are also relatively narrow, with the
DNA binding energy covering a range of $\approx 0.6\, k_BT$ and $\sigma$
spanning $\approx 0.02$. To more precisely quantify the uncertainty, we
computed the shortest interval of the marginal distribution for each
parameter contains 95\% of the probability. The bounds of this interval,
coined the Bayesian credible region, can accommodate asymmetry in the
marginal distribution since the upper and lower bounds of the estimate
are reported. In the main text, we reported the DNA binding energy
estimated from these data to be $15.43^{+0.06}_{-0.06}\, k_BT$, where
the first value is the median of the distribution and the super- and
subscripts correspond to the upper and lower bounds of the credible
region, respectively.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;While looking at the shape of the posterior
distribution can be illuminating, it is not enough to tell us if the
parameter values extracted make sense or accurately describe the data on
which they were conditioned. To assess the validity of the statistical model
in describing actual data, we again turn to simulation, this time using the
posterior distributions for each parameter rather than the prior
distributions. The likelihood of our statistical model assumes that across
the entire induction profile, the observed fold-change is normally
distributed about the theoretical prediction with a standard deviation
$\sigma$. If this is an accurate depiction of the generative process, we
should be able to draw values from the likelihood using the sampled values
for $\Delta\varepsilon_{RA}$ and $\sigma$ that are indistinguishable from the
actual experimental measurements. This process is known as a *posterior
predictive check* and is a Bayesian method of assessing goodness-of-fit.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For each sample from the posterior, we computed the theoretical mean
fold-change given the sampled value for $\Delta\varepsilon_{RA}$. With this mean in hand, we used the
corresponding sample for $\sigma$ and drew a data set from the
likelihood distribution the same size as the real data set used for the
inference. As we did this for every sample of our MCMC output (a total
of $\approx 4000$), it is more instructive to compute the percentiles of
the generated data than to show the entire output. In 
@Fig:epRA_ppc (B), the percentiles of the generated data sets
are shown overlaid with the data used for the inference. We see that all
of the data points fall within the 99$^\text{th}$ percentile of
simulated data sets with the $5^\text{th}$ percentile tracking the mean
of the data at each inducer concentration. As there are no systematic
deviations or experimental observations that fall far outside those
generated from the statistical model, we can safely say that the
statistical model derived here accurately describes the observed data.

![**Markov Chain Monte Carlo (MCMC) samples and posterior predictive check
for DNA binding mutant Q18M.** (A) Marginal and joint sampling
distributions for DNA binding energy $\Delta\varepsilon_{RA}$ and
$\sigma$. Each point in the joint distribution is a single sample.
Marginal distributions for each parameter are shown adjacent to joint
distribution. Color in the joint distribution corresponds to the value
of the log posterior with the progression of dark purple to bright orange
corresponding to increasing probability. (B) The posterior predictive
check of model. The measurements of the fold-change in gene expression
are shown as black open-faced circles. The percentiles are shown as
colored bands and indicate the fraction of simulated data drawn from the
likelihood that fall within the shaded region. The [Python code (`ch7_figS6.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS6.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS6){#fig:epRA_ppc short-caption="Markov chain Monte Carlo samples and
posterior predictive check for DNA binding mutant Q18M."}


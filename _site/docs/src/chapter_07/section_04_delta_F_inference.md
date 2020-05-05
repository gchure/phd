## Inferring the Free Energy From Fold-Change Measurements 

In this section, we describe the statistical model to infer the free
energy $F$ from a set of fold-change measurements. We follow the same
principled workflow as described previously for the DNA binding
estimation, including declaration of the generative model, prior
predictive checks, simulation based calibration, and posterior
predictive checks. Finally, we determine an empirical limit in our
ability to infer the free energy and define a heuristic which can be
used to identify measurements that are likely inaccurate. To understand
the statistical model and the empirical limits of detection, only the
subsections *Building A Generative Model* and *Sensitivity Limits and
Systematic Errors in Inference* are necessary.

### Building A Generative Model 

In Chapter 2, we showed that the fold-change equation can be rewritten in the
form of a Fermi function, 
$$
\text{fold-change} = {1 \over 1 + e^{-F / k_BT}},
$${#eq:mut_si_fermi}
where $F$ corresponds to the free energy
difference between the repressor bound and unbound states of the
promoter. While the theory prescribes a way for us to calculate the free
energy based on our knowledge of the biophysical parameters, we can
directly calculate the free energy of a measurement of fold-change by
simply rearranging @Eq:mut_si_fermi
as
$$
F = -k_BT\log\left({1 \over \text{fold-change}} - 1\right).
$${#eq:empirical_F}
With perfect measurement of the fold-change in
gene expression (assuming no experimental or measurement noise), the
free energy can be directly calculated. However, actual measurements of
the fold-change in gene expression can extend beyond the theoretical
bounds of $0$ and $1$, for which the free energy is mathematically
undefined.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As the fold-change measurements between
biological replicates are independent, it is reasonable to assume that they
are normally distributed about a mean value $\mu$ with a standard deviation
$\sigma$. While the mean value is restricted to the bounds of \[$0, 1$\],
fold-change measurements outside of these bounds are still possible given
that they are distributed about the mean with a scale of $\sigma$. Thus, if
we have knowledge of the mean fold-change in gene expression about which the
observed fold-change is distributed, we can calculate the mean free energy as
$$
F = -k_BT\log\left({1 \over \mu} - 1\right).
$${#eq:empirical_F_mu}
For a given set of fold-change measurements $y$, we wish to infer the
posterior probability distribution for $\mu$ and $\sigma$, given by
Bayes' theorem as
$$
g(\mu, \sigma\, \vert\, y) \propto f(y\,\vert\, \mu, \sigma)g(\mu)g(\sigma),
$${#eq:fc_mu_bayes}
where we have dropped the normalization constant $f(y)$ and assigned a proportionality between the posterior and
joint probability distribution. Given that the measurements are
independent, we define the likelihood $f(y\,\vert\, \mu, \sigma)$ as a
normal distribution,
$$
f(y\,\vert\,\mu\,\sigma) \sim \text{Normal}\{\mu, \sigma\}.
$${#eq:fc_mu_likelihood}
While the mean $\mu$ is restricted to
the interval \[$0, 1$\], there is no reason *a priori* to think that it
is more likely to be closer to either bound. To remain uninformative and
be as permissive as possible, we define a prior distribution for $\mu$
as a Uniform distribution between 0 and 1,
$$
g(\mu)=\begin{cases}{1\over \mu_\text{max} - \mu_\text{min}} & \mu_\text{min} < \mu < \mu_\text{max}\\
0 & \text{otherwise}\end{cases}.
$${#eq:fc_mu_prior}
Here, $\mu_\text{min} = 0$ and $\mu_\text{max} = 1$, reducing $g(\mu)$ to 1. For
$\sigma$, we can again assume a half-normal distribution with a standard
deviation of 0.1 as was used for estimating the DNA binding energy
@Eq:sig_hn_explicit,
$$
g(\sigma) = \text{HalfNormal}\{0, 0.1\}.
$${#eq:fc_sigma_prior}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With a full generative model defined, we can
now use prior predictive checks to ensure that our choices of prior are
appropriate for the inference.

### Prior Predictive Checks 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To check the validity of the chosen priors, we pulled 1000 combinations
of $\mu$ and $\sigma$ from their respective distributions
(@Fig:empirical_F_prior_pred (A)) and subsequently drew a set
of 10 fold-change values (a number comparable to the number of
biological replicates used in this work) from a normal distribution
defined by $\mu$ and $\sigma$. To visualize the range of values
generated from these checks, we computed the percentiles of the
empirical cumulative distributions of the fold-change values, as can be
seen in @Fig:empirical_F_prior_pred (C). Approximately 95\% of the the
generated fold-change measurements were between the theoretical bounds
of \[$0, 1$\] whereas 5\% of the data sets fell outside with the maximum
and minimum values extending to $\approx
1.2$ and $-0.2$, respectively. Given our familiarity with these
experimental strains and the detection sensitivity of the flow
cytometer, these excursions beyond the theoretical bounds agree with our
intuition. Satisfied with our choice of prior distributions, we can
proceed to check the sensitivity and computational tractability of our
model through simulation based calibration.

![**Prior predictive checks for inference of the mean fold-change.** (A) The
prior distributions for $\mu$ (left) and $\sigma$ (right). The vertical
axis is proportional to the probability of the value. Black points above
distributions correspond to the values used to perform the prior
predictive checks. (B) Percentiles of the data generated for each draw
from the prior distributions shown as a cumulative distribution.
Percentiles were calculated for 1000 generated data sets, each with 10
fold-change measurements drawn from the likelihood given the drawn
values of $\mu$ and
$\sigma$. The [Python code                                                
(`ch7_figS7.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS7.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS7){#fig:empirical_F_prior_pred short-caption="Prior predictive checks
for inference of the mean fold-change"}




### Simulation Based Calibration 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To ensure that the parameters can be estimated with confidence, we
sampled the posterior distribution of $\mu$ and $\sigma$ for each data
set generated from the prior predictive checks. For each inference, we
computed the $z$-score and shrinkage for each parameter, shown in @Fig:empirical_F_sensitivity(A). For both parameters, the
$z$-scores are approximately centered about zero, indicating that the
posteriors concentrate about the ground truth value of the parameter.
The $z$-scores for $\sigma$ [green points in Fig.
@Fig:empirical_F_sensitivity) (A)] appear to be slightly off
centered with more negative values than positive. This suggests that
$\sigma$ is more likely to be slightly overestimated in some cases. The
shrinkage parameter for $\mu$ (red points) is very tightly distributed
about $1.0$, indicating that the prior is being strongly informed by the
data. The shrinkage is more broadly distributed for for $\sigma$ with a
minimum value of $\approx 0.5$. However, the median shrinkage for
$\sigma$ is $\approx 0.9$, indicating that half of the inferences shrank
the prior distribution by at least 90\%. While we could revisit the model
to try and improve the shrinkage values, we are more concerned with
$\mu$ which shows high shrinkage and zero-centered $z$-scores.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To ensure that the model is computationally tractable, we computed the rank
statistic of each parameter for each inference. The empirical cumulative
distributions for $\mu$ (black) and $\sigma$ (red) can be seen in
@Fig:empirical_F_sensitivity (B). Both distributions appear to be uniform,
falling within the 99$^\text{th}$ percentile of the variation expected from a
true uniform distribution. This indicates that the self-consistency relation
defined by @Eq:self_consistency_relation. holds for this statistical model.
With a computationally tractable model in hand, we can now apply the
statistical model to our data and verify that data sets drawn from the
data-conditioned posterior are indistinguishable from the experimental
measurements.

![**Sensitivity measurements and rank statistic distribution of the
statistical model estimating $\mu$ and $\sigma$.** (A) Posterior $z$-score
of each inference plotted against the posterior shrinkage factor for the
parameters $\mu$ (blue points) and $\sigma$ (green points). (B)
Distribution of rank statistics for $\mu$ (red) and $\sigma$ (black).
Purple envelope represents the 99$^\text{th}$ percentile of a true uniform
distribution. The [Python code                                                
(`ch7_figS8.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS8.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS8){#fig:empirical_F_sensitivity
short-caption="Sensitivity measurements and rank statistic distribution of the
statistical model estimating $\mu$ and $\sigma$."}


### Posterior Predictive Checks 

The same statistical model was applied to every unique set of fold-change
measurements used in this work. Here, we focus only on the set of fold-change
measurements for the double mutant Y17I-Q291V at 50 $\mu$M IPTG. The samples
from the posterior distribution conditioned on this dataset can be seen in
@Fig:empirical_F_post_pred (A). The joint distribution, shown in the lower
left-hand corner, appears fairly symmetric, indicating that $\mu$ and
$\sigma$ are independent. There is a slight asymmetry in the sampling of
$\sigma$, which can be more clearly seen in the corresponding marginal
distribution to the right of the joint distribution.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For each MCMC sample of $\mu$ and $\sigma$, we drew 10 samples from a
normal distribution defined by these parameters. From this collection of
data sets, we computed the percentiles of the empirical cumulative
distribution and plotted them over the data, as can be seen in 
@Fig:empirical_F_post_pred (B). We find that the observed
data falls within the 99$^\text{th}$ percentile of the generated data
sets. This illustrates that the model can produce data which is
identically distributed to the actual experimental measurements,
validating our choice of statistical model.

![**MCMC sampling output and posterior predictive checks of the
statistical model for the mean fold-change $\mu$ and standard deviation
$\sigma$.** (A) Corner plot of sampling output. The joint distribution
between $\sigma$ and $\mu$ is shown in the lower left hand corner. Each
point is an individual sample. Points are colored by the value of the
log posterior with increasing probability corresponding to transitions
from purple to orange. Marginal distributions for each parameter are shown
adjacent to the joint distribution. (B) Percentiles of the cumulative
distributions from the posterior predictive checks are shown as shaded
bars. Data on which the posterior was conditioned are shown as white
orange circles and lines. The [Python code                                                
(`ch7_figS9.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS9.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS9){#fig:empirical_F_post_pred
short-caption="MCMC sampling output and posterior predictive checks of the
statistical model for the mean fold-change and standard deviation."}


### Sensitivity Limits and Systematic Errors in Inference 


Considering the results from the prior predictive checks, simulation
based calibration, and posterior predictive checks, we can say that the
statistical model for inferring $\mu$ and $\sigma$ fold-change from a
collection of noisy fold-change measurements is valid and
computationally tractable. Upon applying this model to the experimental
data of the wild-type strain (where the free energy is theoretically
known), we observed that systematic errors arise when the fold-change is
exceptionally high or low, making the resulting inference of the free
energy inaccurate.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To elucidate the source of this systematic error, we return to a
simulation based approach in which the true free energy is known (black
points in @Fig:empirical_F_error (A)). For a range of free energies,
we computed the theoretical fold-change prescribed by
@Eq:mut_si_fermi. For each free energy value, we pulled a value for
$\sigma$ from the prior distribution defined in
@Eq:sig_hn_explicit and generated a data set of 10
measurements by drawing values from a normal distribution defined by the
true fold-change and the drawn value of $\sigma$ (purple points in
@Fig:empirical_F_error (A)). We then sampled the statistical
model over these data and inferred the mean fold-change $\mu$ (orange points in
@Fig:empirical_F_error (A)). By eye, the inferred points
appear to collapse onto the master curve, in many cases overlapping the
true values. However, the points with a free energy less than
$\approx-2\, k_BT$ and greater than $\approx 2\, k_BT$ are slightly
above or below the master curve, respectively. This becomes more obvious
when the inferred free energy is plotted as a function of the true free
energy, shown in @Fig:empirical_F_error (B). Points in which the difference
between $\mu$ and the nearest boundary ($0$ or $1$) is less than the
value of $\sigma$ are shown as purple or green. When this condition is
met, the inferred mean free energy strays from the true value,
introducing a systematic error. This suggests that the spread of the
fold-change measurements sets the detection limit of fold-change close
to either boundary. Thus, the narrower the spread in the fold-change the
better the estimate of the fold-change near the boundaries.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;These systematic errors can be seen in
experimental measurements of the wild-type repressor. Data from
@razo-mejia2018 in which the IPTG titration profiles of seventeen different
bacterial strains were measured is shown collapsed onto the master curve in
@Fig:empirical_F_error (C) as red points. Here, each point corresponds
to a single biological replicate. The inferred mean fold-change $\mu$ and
95\% credible regions are shown as purple, blue, or green points. The color
of these points correspond to the relative value of $\mu$ or $1 - \mu$ to
$\sigma$. The discrepancy between the predicted and inferred free energy of
each measurement set can be seen in @Fig:empirical_F_error (D). The
significant deviation from the predicted and inferred free energy occurs past
the detection limit set by $\sigma$. In this work, we therefore opted to not
display inferred free energies at the extrema where the inferred fold-change
was closer to the boundaries than the corresponding standard deviation, as it
reflects limitations in our measurement rather than a deviation from the
theoretical predictions.

![**Identification of systematic error in simulated and real data when
considering the free energy.** (A) The true fold-change (black points),
simulated fold-change distribution (purple points), and inferred
mean fold-change (orange) is plotted as a function of the true free
energy. Error bars on inferred fold-change correspond to the 95\%
credible region of the mean fold-change $\mu$. (B) Inferred free energy
plotted as a function of the true free energy. Black line indicates
perfect agreement between the ground truth free energy and inferred free
energy. Blue points correspond to the inferred free energy where the
median values of the parameters satisfy the condition $\mu > \sigma$ and
$1 - \mu > \sigma$. Purple points correspond to the inferred mean
fold-change $\mu < \sigma$. Green points correspond to those where the
inferred mean fold-change $1 - \mu < \sigma$. Error bars correspond to
the bounds of the 95\% credible region. (C) Biological replicate data
from @razo-mejia2018 (red points) plotted as a
function of the theoretical free energy. Inferred mean fold-change $\mu$
and the 95\% credible region are shown as blue points. Purple and green
points are colored by the same conditions as in (B). (D) Inferred free
energy as a function of the predicted free energy colored by the
satisfied condition. Error bars are the bounds of the 95\% credible
region. All inferred values in (A -- D) are the median values of the
posterior distribution. The [Python code                                                
(`ch7_figS10.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS10.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS10){#fig:empirical_F_error
short-caption="Identification of systematic error in simulated and real data
when considering the free energy."}

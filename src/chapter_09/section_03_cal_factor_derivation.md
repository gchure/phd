## Calibration of a Standard Candle 

&nbsp;&nbsp;&nbsp;&nbsp;To estimate the single-cell MscL abundance via
microscopy, we needed to determine a calibration factor that could translate
arbitrary fluorescence units to protein copy number. To compute this
calibration factor, we relied on *a priori* knowledge of the mean copy number
of MscL-sfGFP for a particular bacterial strain in specific growth
conditions. In @bialecka-fornal2012, the
average MscL copy number for a population of cells expressing an MscL-sfGFP
fusion (*E. coli* K-12 MG1655 $\phi$(*mscL-sfGFP*)) cells was measured using
quantitative Western blotting and single-molecule photobleaching assays. By
growing this strain in identical growth and imaging conditions, we can make
an approximate measure of this calibration factor. In this section, we derive
a statistical model for estimating the most-likely value of this calibration
factor and its associated error.

### Definition of a Calibration Factor 
&nbsp;&nbsp;&nbsp;&nbsp;We assume that all detected fluorescence signal from
a particular cell is derived from the MscL-sfGFP protein, after background
subtraction and correction for autofluorescence. The arbitrary units of
fluorescence can be directly related to the protein copy number via a
calibration factor $\alpha$,
$$
I_\text{tot} = \alpha N_\text{tot},
$${#eq:mscl_ian}
where $I_\text{tot}$ is the total cell fluorescence and $N_\text{tot}$ is the
total number of MscL proteins per cell. @bialecka-fornal2012 report the
average cell MscL copy number for the population rather than the
distribution. Knowing only the mean, we can rewrite [@Eq:mscl_ian] as
$$
\langle I_\text{tot}\rangle = \alpha \langle N_\text{tot} \rangle,
$${#eq:mscl_avg_ian}
assuming that $\alpha$ is a constant value that does not change from cell to
cell or fluorophore to fluorophore.

&nbsp;&nbsp;&nbsp;&nbsp;The experiments presented in this work were performed
using non-synchronously growing cultures. As there is a uniform distribution
of growth phases in the culture, the cell size distribution is broad. As
described in the main text, the cell size distribution of a population is
broadened further by modulating the MscL copy number with low copy numbers
resulting in aberrant cell morphology. To speak in the terms of an effective
channel copy number, we relate the average areal intensity of the population
to the average cell size,
$$
\langle I_\text{tot} \rangle = \langle I_A \rangle \langle A \rangle = \alpha \langle N_\text{tot} \rangle,
$${#eq:area_conversion}
where $\langle I_A\rangle$ is the average areal intensity in arbitrary units
per pixel of the population and $\langle A \rangle$ is the average area of a
segmented cell. As only one focal plane was imaged in these experiments, we
could not compute an appropriate volume for each cell given the highly
aberrant morphology. We therefore opted to use the projected two-dimensional
area of each cell as a proxy for cell size. Given this set of measurements,
the calibration factor can be computed as
$$
\alpha = {\langle I_A \rangle\langle A \rangle \over \langle N_\text{tot} \rangle}.
$${#eq:mscl_simple_cal_factor}
&nbsp;&nbsp;&nbsp;&nbsp;While it is tempting to use @Eq:mscl_simple_cal_factor
directly, there are multiple sources of error that are important to propagate
through the final calculation. The most obvious error to include is the
measurement error reported in @bialecka-fornal2012 for the average
MscL channel count . There are also slight variations
in expression across biological replicates that arise from a myriad of
day-to-day differences. Rather than abstracting all sources of error away
into a systematic error budget, we used an inferential model derived from
Bayes' theorem that allows for the computation of the probability
distribution of $\alpha$.

### Estimation of $\alpha$ for a Single Biological Replicate 

&nbsp;&nbsp;&nbsp;&nbsp;A single data set consists of several hundred
single-cell measurements of intensity, area of the segmentation mask, and
other morphological quantities. The areal density $I_A$ is computed by
dividing the total cell fluorescence by the cell area $A$. We are interested
in computing the probability distributions for the calibration factor
$\alpha$, the average cell area $\langle A \rangle$, and the mean number of
channels per cell $\langle N_\text{tot} \rangle$ for the data set as a whole
given only $I_A$ and $A$. Using Bayes' theorem, the probability distribution
for these parameters given a single cell measurement, hereafter called the
posterior distribution, can be written as
$$
g(\alpha, \langle A \rangle, \langle N_\text{tot} \rangle\,\vert\, A, I_A) = {f(A, I_A\,\vert
\, \alpha,\langle A \rangle, \langle N_\text{tot} \rangle) g(\alpha,\langle A \rangle, \langle N_\text{tot} \rangle) \over f(\alpha, I_A)},
$${#eq:simple_bayes}
where $g$ and $f$ represent probability density functions over parameters and
data, respectively. The term $f(A, I_A\,\vert\, \alpha, \langle A \rangle,
\langle N_\text{tot} \rangle)$ in the numerator represents the likelihood of
observing the areal intensity $I_A$ and area $A$ of a cell for a given value
of $\alpha$, $\langle A \rangle$, and $\langle N_\text{tot} \rangle$. The
second term in the numerator $g(\alpha,\langle A \rangle, \langle
N_\text{tot} \rangle)$ captures all prior knowledge we have regarding the
possible values of these parameters knowing nothing about the measured data.
The denominator, $f(I_A, A)$ captures the probability of observing the data
knowing nothing about the parameter values. This term, in our case, serves
simply as a normalization constant and is neglected for the remainder of this
section.

&nbsp;&nbsp;&nbsp;&nbsp;To determine the appropriate functional form for the
likelihood and prior, we must make some assumptions regarding the biological
processes that generate them. As there are many independent processes that
regulate the timing of cell division and cell growth, such as DNA replication
and peptidoglycan synthesis, it is reasonable to assume that for a given
culture, the distribution of cell size would be normally distributed with a
mean of $\langle A \rangle$ and a variance $\sigma_{\langle A \rangle}$.
Mathematically, we can write this as
$$
f(A\,\vert\,\langle A \rangle, \sigma_{\langle A \rangle}) \propto {1 \over \sigma_{\langle A \rangle} }\exp\left[-{(A - \langle A \rangle)^2 \over 2\sigma_{\langle A \rangle}^2}\right],
$${#eq:area_likelihood}
where the proportionality results from dropping normalization constants for
notational simplicity.

&nbsp;&nbsp;&nbsp;&nbsp; While total cell intensity is intrinsically
dependent on the cell area, the areal intensity $I_A$ is independent of cell
size. The myriad processes leading to the detected fluorescence, such as
translation and proper protein folding, are largely independent, allowing us
to assume a normal distribution for $I_A$ as well with a mean $\langle I_A
\rangle$ and a variance $\sigma_{I_A}^2$. However, we do not have knowledge
of the average areal intensity for the standard candle strain *a priori*.
This can be calculated knowing the calibration factor, total MscL channel
copy number, and the average cell area as
$$
I_A =  {\alpha\langle N_\text{tot} \rangle \over \langle A \rangle}.
$${#eq:avg_rho}
Using @Eq:avg_rho to calculate the expected areal intensity for the population, we can write the likelihood as a Gaussian distribution, 
$$
f(I_A\,\vert\,\alpha,\langle A \rangle,\langle N_\text{tot} \rangle, \sigma_{I_A}) \propto {1 \over \sigma_{I_A} }\exp\left[-{\left(I_A - {\alpha \langle N_\text{tot} \rangle\over \langle A \rangle}\right)^2 \over 2 \sigma_{I_A}^2}\right].
$${#eq:rho_likelihood}

&nbsp;&nbsp;&nbsp;&nbsp;With these two likelihoods in hand, we are tasked with determining the appropriate priors. As we have assumed normal distributions for the likelihoods of $\langle A \rangle$ and $I_A$, we have included two additional parameters, $\sigma_{\langle A \rangle}$ and $\sigma_{I_A}$, each requiring their own prior probability distribution. It is common practice to assume maximum ignorance for these variances and use a Jeffreys prior [@sivia2006],
$$
g(\sigma_{\langle A \rangle}, \sigma_{I_A}) = {1 \over \sigma_{\langle A \rangle}\sigma_{I_A} }.
$${#eq:jeffreys}

&nbsp;&nbsp;&nbsp;&nbsp;The next obvious prior to consider is for the average channel copy number $\langle N_\text{tot} \rangle$, which comes from @bialecka-fornal2012. In this work, they report a mean $\mu_N$  and variance $\sigma_N^2$, allowing us to assume a normal distribution for the prior,
$$
g(\langle N_\text{tot}\rangle\,\vert\, \mu_N,\sigma_N) \propto {1 \over \sigma_N}\exp\left[-{(\langle N_\text{tot} \rangle - \mu_N)^2 \over 2 \sigma_N^2}\right].
$${#eq:informative_prior}
For $\alpha$ and $\langle A \rangle$, we have some knowledge of what these parameters can and cannot be. For example, we know that neither of these parameters can be negative. As we have been careful to not overexpose the microscopy images, we can say that the maximum value of $\alpha$ would be the bit-depth of our camera. Similarly, it is impossible to segment a single cell with an area larger than our camera's field of view (although there are biological limitations to size below this extreme). To remain maximally uninformative, we can assume that the parameter values are uniformly distributed between these bounds, allowing us to state 
$$
g(\alpha) = \begin{cases} {1 \over \alpha_\text{max} - \alpha_\text{min} } & \alpha_\text{min} \leq \alpha \leq \alpha_\text{max} \\
0 & \text{otherwise}
\end{cases},
$${#eq:alpha_uniform_prior}
for $\alpha$ and 
$$
g(\langle A \rangle) = \begin{cases} {1 \over \langle A \rangle_\text{max} - \langle A \rangle_\text{min} } & \langle A \rangle_\text{min} \leq \langle A \rangle \leq \langle A \rangle_\text{max}\\
0 & \text{otherwise}
\end{cases}
$${#eq:area_uniform_prior}
for $\langle A \rangle$.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Piecing @Eq:area_likelihood through @Eq:area_uniform_prior together generates
a complete posterior probability distribution for the parameters given a
single cell measurement. This can be generalized to a set of $k$ single cell
measurements as
$$
\begin{aligned}
g(\alpha,&\langle A \rangle, \langle N_\text{tot} \rangle, \sigma_{I_A},
\sigma_{\langle A \rangle}\,\vert\, [I_A, A], \mu_N, \sigma_N) \propto {1 \over
(\alpha_\text{max} - \alpha_\text{min})(\langle A \rangle_\text{max} - \langle A
\rangle_\text{min})}\times \\ &{1 \over (\sigma_{I_A}\sigma_{\langle A \rangle})^{k+1}
}{1 \over \sigma_N}\exp\left[- {(\langle N_\text{tot}\rangle - \mu_N)^2 \over 2\sigma_N^2}\right]
\\ &\prod\limits_i^k\exp\left[-{(A^{(i)} - \langle A \rangle)^2 \over 2\sigma_{\langle A \rangle}^2} - {\left(I_A^{(i)} - {\alpha \langle N_\text{tot}\rangle \over \langle A \rangle}\right)^2 \over 2\sigma_{I_A}^2}\right] \end{aligned},
$${#eq:single_rep_post}
where $[I_A, A]$ represents the set of $k$ single-cell measurements.

&nbsp;&nbsp;&nbsp;&nbsp;As small variations in the day-to-day details of cell
growth and sample preparation can alter the final channel count of the
standard candle strain, it is imperative to perform more than a single
biological replicate. However, properly propagating the error across
replicates is not trivial. One option would be to pool together all
measurements of $n$ biological replicates and evaluate the posterior given in
@Eq:single_rep_post. However, this by definition assumes that there is no
difference between replicates. Another option would be to perform this
analysis on each biological replicate individually, and then compute a mean
and standard deviation of the resulting most-likely parameter estimates for
$\alpha$ and $\langle A \rangle$. While this is a better approach than simply
pooling all data together, it suffers a bias from giving each replicate equal
weight, skewing the estimate of the most-likely parameter value if one
replicate is markedly brighter or dimmer than the others. Given this type of
data and a limited number of biological replicates ($n = 6$ in this work), we
chose to extend the Bayesian analysis presented in this section to model the
posterior probability distribution for $\alpha$ and $\langle A \rangle$ as a
hierarchical process in which $\alpha$ and $\langle A \rangle$ for each
replicate is drawn from the same distribution.

### A Hierarchical Model for Estimating $\alpha$

&nbsp;&nbsp;&nbsp;&nbsp;In the previous section, we assumed maximally
uninformative priors for the most-likely values of $\alpha$ and $\langle A
\rangle$. While this is a fair approach to take, we are not completely
ignorant with regard to how these values are distributed across biological
replicates. A major assumption of our model is that the most-likely value of
$\alpha$ and $\langle A \rangle$ for each biological replicate are
comparable, so long as the experimental error between them is minimized. In
other words, we are assuming that the most-likely value for each parameter
for each replicate is drawn from the same distribution. While each replicate
may have a unique value, they are all related to one another. Unfortunately,
proper sampling of this distribution requires an extensive amount of
experimental work, making inferential approaches more attractive.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This approach, often called a multi-level or hierarchical model, is
schematized in @Fig:hierarchical_model. Here, we use an informative prior
for $\alpha$ and $\langle A \rangle$ for each biological replicate. This
informative prior probability distribution can be described by summary
statistics, often called hyper-parameters, which are then treated as the
"true" value and are used to calculate the channel copy number. As this
approach allows us to get a picture of the probability distribution of the
hyper-parameters, we are able to report a point estimate for the most-likely
value along with an error estimate that captures all known sources of
variation.

![**Schematic of hierarchical model structure.** The hyper-parameter probability
distributions (top panel) are used as an informative prior for the
most-likely parameter values for each biological replicate (middle panel).
The single-cell measurements of cell area and areal intensity (bottom panel)
are used as data in the evaluation of the
likelihood.](ch9_figS3){#fig:hierarchical_model short-caption="Schematic of
hierarchical model structure for estimating MscL-sfGFP fluorescence calibration factor."}

&nbsp;&nbsp;&nbsp;&nbsp;The choice for the functional form for the
informative prior is often not obvious and can require other experimental
approaches or back-of-the-envelope estimates to approximate. Each experiment
in this work was carefully constructed to minimize the day-to-day variation.
This involved adhering to well-controlled growth temperatures and media
composition, harvesting of cells at comparable optical densities, and
ensuring identical imaging parameters. As the experimental variation is
minimized, we can use our knowledge of the underlying biological processes to
guess at the approximate functional form. For similar reasons presented in
the previous section, cell size is controlled by a myriad of independent
processes. As each replicate is independent of another, it is reasonable to
assume a normal distribution for the average cell area for each replicate.
This normal distribution is described by a mean $\tilde{\langle A \rangle}$
and variance $\tilde{\sigma}_{\langle A \rangle}$. Therefore, the prior for
$\langle A \rangle$ for $n$ biological replicates can be written as
$$
g(\langle A \rangle\, \vert\, \tilde{\langle A \rangle}, \tilde{\sigma}_{\langle A \rangle}) \propto {1 \over \tilde{\sigma}_{\langle A \rangle}^n}\prod\limits_{j=1}^{n}\exp\left[-{(\langle A \rangle_j - \tilde{\langle A \rangle})^2 \over 2 \tilde{\sigma}_{\langle A \rangle}^2}\right].
$${#eq:hyperprior_area}
In a similar manner, we can assume that the calibration factor for each
replicate is normally distributed with a mean $\tilde{\alpha}$ and variance
$\tilde{\sigma}_\alpha$,
$$
g(\alpha\,\vert\,\tilde{\alpha}, \tilde{\sigma}_\alpha) \propto {1 \over \tilde{\sigma}_\alpha^n}\prod\limits_{j=1}^n \exp\left[-{(\alpha_j - \tilde{\alpha})^2 \over 2\tilde{\sigma}_\alpha^2}\right].
$${#eq:hyperprior_alpha}

&nbsp;&nbsp;&nbsp;&nbsp;With the inclusion of two more normal distributions,
we have introduced four new parameters, each of which need their own
prior. However, our knowledge of the reasonable values for the
hyper-parameters has not changed from those described for a single replicate.
We can therefore use the same maximally uninformative Jeffreys priors given
in @Eq:jeffreys for the variances and the uniform distributions given in
@Eq:alpha_uniform_prior and @Eq:area_uniform_prior for the means.
Stitching all of this work together generates the full posterior probability
distribution for the best-estimate of $\tilde{\alpha}$ and $\tilde{\langle A
\rangle}$ shown in @Eq:mscl_avg_ian given $n$ replicates of $k$ single cell
measurements,
$$
\begin{aligned}
g(\tilde{\alpha}, \tilde{\sigma}_\alpha, \tilde{\langle A \rangle}, \tilde{\sigma}_{\langle A \rangle}, &\{\langle N_\text{tot} \rangle, \langle A \rangle, \alpha, \sigma_{I_A}\}\,\vert\, [I_A, A], \mu_N, \sigma_N) \propto\\
&{1 \over (\tilde{\alpha}_\text{max} - \tilde{\alpha}_\text{min})(\tilde{\langle A \rangle}_\text{max} - \tilde{\langle A \rangle}_\text{min})\sigma_N^n(\tilde{\sigma}_\alpha\tilde{\sigma}_{\langle A \rangle})^{n + 1} }\,\times\\
&\prod\limits_{j=1}^n\exp\left[-{(\langle N \rangle_j^{(i)} - \mu_N)^2 \over 2\sigma_N^2} - {(\alpha_j - \tilde{\alpha})^2 \over 2\tilde{\sigma}_\alpha^2} - {(\langle A \rangle_j - \tilde{\langle A \rangle})^2 \over 2\tilde{\sigma}_{\langle A \rangle}^2}\right]\,\times\,\\ 
&{1 \over (\sigma_{ {I_A}_j}\sigma_{\langle A \rangle_j})^{k_j + 1} }\prod\limits_{i=1}^{k_j}\exp\left[-{(A_j^{(i)} - \langle A \rangle_j)^2 \over 2\sigma^{(i)2}_{\langle A \rangle_j} } - {\left({I_A}_{j}^{(i)} - {\alpha_j \langle N_\text{tot}\rangle_j \over \langle A \rangle_j}\right)\over 2\sigma_{ {I_A}_j}^{(i)2} }\right]
\end{aligned}
$${#eq:cal_factor_posterior}
where the braces $\{\dots\}$ represent the set of parameters for biological
replicates and the brackets $[\dots]$ correspond to the set of single-cell
measurements for a given replicate.

&nbsp;&nbsp;&nbsp;&nbsp;While @Eq:cal_factor_posterior is not analytically
solvable, it can be easily sampled using Markov chain Monte Carlo (MCMC). The
results of the MCMC sampling for $\tilde{\alpha}$ and $\tilde{\langle A
\rangle}$ can be seen in @Fig:mscl_posterior_samples. From this approach, we
found the most-likely parameter values of $3300^{+700}_{-700}$ a.u. per MscL
channel and $5.4^{+0.4}_{-0.5}$ $\mu$m$^2$ for $\tilde{\alpha}$ and
$\tilde{\langle A \rangle}$, respectively. Here, we have reported the median
value of the posterior distribution for each parameter with the upper and
lower bound of the 95\% credible region as superscript and subscript,
respectively. These values and associated errors were used in the calculation
of channel copy number.

![**Posterior distributions for hyper-parameters and replicate parameters.**
(A) The posterior probability distribution for $\tilde{\alpha}$ and
$\tilde{\langle A \rangle}$. Probability increases from light to dark red.
The replicate parameter (purple) and hyper-parameter (orange) are marginalized
posterior probability distributions for $\alpha$ (B) and $\langle A \rangle$
(C). The [Python code (`ch9_figS4.py`)](https://github.com/gchure/phd/blob/master/src/chapter_09/code/ch9_figS4.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).
](ch9_figS4){#fig:mscl_posterior_samples short-caption="Posterior
distributions for hyper-parameters and replicate parameters."}

### Effect of Correction

&nbsp;&nbsp;&nbsp;&nbsp; The posterior distributions for $\alpha$ and $\langle A
\rangle$ shown in @Fig:mscl_posterior_samples were used directly to compute the
most-likely channel copy number for each measurement of the Shine-Dalgarno
mutant strains, as is described in the coming section.
The importance of this correction can be seen in @Fig:area_correction. Cells
with low abundance of MscL channels exhibit notable morphological defects, as
illustrated in @Fig:area_correction (A). While these would all be considered
single cells, the two-dimensional area of each may be comparable to two or three
wild-type cells. For all of the Shine-Dalgarno mutants, the distribution of
projected cell area has a long tail, with the extremes reaching 35
$\mu\text{m}^2$ per cell (@Fig:area_correction (B)). Calculating the total number
of channels per cell  does nothing to decouple this correlation between cell
area and measured cell intensity. @Fig:area_correction (C) shows the correlation
between cell area and the total number of channels without normalizing to an
average cell size $\langle A \rangle$ differentiated by their survival after an
osmotic downshock. This correlation is removed by calculating an effective
channel copy number shown in @Fig:area_correction (D). 

![**Influence of area correction for Shine-Dalgarno mutants.** (A)
Representative images of aberrant cell morphologies found in low-expressing
Shine-Dalgarno mutants. (B) Empirical cumulative distribution of
two-dimensional projected cell area for the standard candle strain MLG910
(gray line) and for all Shine-Dalgarno mutants (red line). (C) The
correlation between channel copy number and cell area without the area
correction. (D) The correlation between effective channel copy number and
cell area with the area correction applied. The [Python code (`ch9_figS5.py`)](https://github.com/gchure/phd/blob/master/src/chapter_09/code/ch9_figS5.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).
 ](ch9_figS5){#fig:area_correction
short-caption="Influence of area correction for Shine-Dalgarno mutants."}

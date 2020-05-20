## Statistical Inference of Entropic Costs

In the main text, we describe how a simple rescaling of the energetic
parameters $\Delta\varepsilon_R$ and $\Delta\varepsilon_{AI}$ is not
sufficient to describe the fold-change in gene expression when the
growth temperature is changed from 37$^\circ$ C. In this section, we
describe the inference of hidden entropic parameters to
phenomenologically describe the temperature dependence of the
fold-change in gene expression.

### Definition of Hidden Entropic Costs

The values of the energetic parameters $\Delta\varepsilon_R$ and
$\Delta\varepsilon_{AI}$ were determined in a glucose supplemented medium
held at 37$^\circ$ C which we denote as $T_{ref}$. A null model to describe
temperature dependence of these parameters is to rescale them to the changed
temperature $T_{exp}$ as 
$$
\Delta\varepsilon^* =
\frac{T_{ref}}{T_{exp}}\Delta\varepsilon,
$${#eq:epstar} 
where $\Delta\varepsilon^*$ is either $\Delta\varepsilon_R$ or
$\Delta\varepsilon_{AI}$. However, we found that this null model was not
sufficient to describe the fold-change in gene expression, prompting the
formulation of a new phenomenological description.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Our thermodynamic model for the fold-change in
gene expression coarse-grains the regulatory architecture to
a two-state model, meaning many of the rich features of regulation such as
vibrational entropy, the material properties of DNA, and the occupancy of the
repressor to the DNA are swept into the effective energetic parameters. As
temperature was never perturbed when this model was developed, modeling these
features was not necessary. However, we must now return to these features to
consider what may be affected.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Without assigning a specific mechanism, we can
say that there is a temperature-dependent entropic parameter that was
neglected in the estimation of the energetic parameters in
@garcia2011 and @razo-mejia2018.
In this case, the inferred energetic parameter $\Delta\varepsilon^*$ is
composed of enthalpic ($\Delta H$) and entropic ($\Delta S$) parameters,
$$\Delta\varepsilon^* = \Delta H - T\Delta S. \label{eq:utds}$$ For a set of
fold-change measurements at a temperature $T_{exp}$, we are interested in
estimating values for $\Delta H$ and $\Delta S$ for each energetic parameter.
Given measurements from Refs. [@garcia2011; @razo-mejia2018], we know at
37$^\circ C$ what $\Delta\varepsilon_{RA}$ and $\Delta\varepsilon_{AI}$ are
inferred to be, placing a constraint on the possible values of $\Delta H$ and
$\Delta S$,
$$
\Delta H_R = \Delta\varepsilon_{R} + T_{ref}\Delta S_R = T_{ref}\Delta S_R -
13.9\, k_BT,
$${#eq:epra_constraint}
and
$$
\Delta H_{AI} = \Delta\varepsilon_{AI} + T_{ref}\Delta S_{AI} = T_{ref}\Delta S_{AI} + 4.5\,k_BT
$${#eq:epai_constraint}
for the DNA binding energy and allosteric state energy difference, respectively.

### Statistical Inference of $\Delta S_R$ and $\Delta S_{AI}$ 

Given the constraints from @Eq:epra_constraint and @Eq:epai_constraint, we
are interested in inferring the entropic parameters $\Delta S_R$ and $\Delta
S_{AI}$ given literature values for $\Delta\varepsilon_{RA}$ and
$\Delta\varepsilon_{AI}$ and the set of fold-change measurements
$\mathbf{fc}$ at a given temperature $T_{exp}$. The posterior probability
distribution for the entropic parameters can be enumerated via Bayes' theorem
as
$$
g(\Delta S_R, \Delta S_{AI}\,\vert\, \mathbf{fc}) =
\frac{f(\mathbf{fc}\,\vert\,\Delta S_R, \Delta S_{AI})g(\Delta S_R, \Delta
S_{AI})}{f(\mathbf{fc})},
$${#eq:temp_generic_bayes}
where $g$ and $f$ are used to denote probability densities over parameters
and data, respectively. As we have done elsewhere in this SI text, we treat
$f(\mathbf{fc})$ as a normalization constant and neglect it in our estimation
of $g(\Delta S_R, \Delta S_{AI}\,\vert\, \mathbf{fc})$. Additionally, as is
discussed in detail earlier in this chapter, we consider the $\log$
fold-change measurements to be normally distributed about a mean
$\text{fc}^*$ defined by the fold-change input-output function and a standard
deviation $\sigma$. Thus, the likelihood for the fold-change in gene
expression is
$$
f(\mathbf{fc}\,\vert\,\Delta S_R, \Delta S_{AI}, \sigma, T_{exp}) = 
\frac{1}{\left(2\pi\sigma^2\right)^{N/2}}\prod\limits_i^N
\exp\left[-\frac{\left[\log \text{fc}_i - \log \text{fc}^*(\Delta S_R, \Delta
S_{AI}, T_{exp})\right]^2}{2\sigma^2}\right].
$${#eq:delS_like}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In calculating the mean $\text{fc}^*$, the effective energetic
parameters $\Delta\varepsilon_R^*$ and $\Delta\varepsilon_{AI}^*$ can be
defined and constrained using @Eq:epra_constraint and @Eq:epai_constraint as
$$
\Delta\varepsilon_{R}^* = \Delta H_R - T_{exp}\Delta S_R = \Delta S_R(T_{ref} -
T_{exp}) - 13.9\, k_BT_{ref},
$${#eq:eprastar_constraint}
and
$$
\Delta\varepsilon_{AI}^* = \Delta H_{AI} - T_{exp}\Delta S_{AI}  = \Delta
S_{AI}(T_{ref} - T_{exp}) + 4.5\, k_BT_{ref}.
$${#eq:epaistar_constraint}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As the enthalpic parameters are calculated
directly from the constraints of @Eq:epra_constraint and @Eq:epai_constraint,
we must only estimate three parameters, $\Delta S_R$, $\Delta S_{AI}$, and
$\sigma$, each of which need a functional form for the prior distribution.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*A priori*, we know that both $\Delta S_R$ and
$\Delta S_{AI}$ must be small because $T_{ref}$ and $T_{exp}$ are defined in
K. As these entropic parameters can be either positive or negative, we can
define the prior distributions $g(\Delta S_R)$ and $g(\Delta S_{AI})$ as a
normal distribution centered at zero with a small standard deviation,
$$
g(\Delta S_R) \sim \text{\itshape Normal}(\mu=0, \sigma=0.1),
$${#eq:delSR_prior}
and
$$
g(\Delta S_{AI}) \sim \text{\itshape Normal}(\mu=0, \sigma=0.1).
$${#eq:delSAI_prior}
The standard deviation $\sigma$ can be defined as a half-normal distribution centered at 0 with
a small standard deviation $\phi$,
$$
g(\sigma) \sim \text{\itshape HalfNormal}(\phi=0.1).
$${#eq:unlabeled_sigma}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With the priors and likelihood functions in
hand, we sampled the posterior distribution using Markov chain Monte Carlo as
implemented in the Stan probabilistic programming language [@carpenter2017].
We performed three different estimations -- one inferring the parameters
using only data at 32$^\circ$ C, one using only data from 42$^\circ$ C, and
one using data sets from both temperatures pooled together.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The sampling results can be seen in
@Fig:entropy_corner. The estimation of $\Delta S_R$ is distinct for each
condition whereas the sampling for $\Delta S_{AI}$ is the same for all
conditions and is centered at about 0. The latter suggests that the value of
$\Delta\varepsilon_{AI}$ determined at 37$^\circ$ C is not dependent on
temperature within the resolution of our experiments. The difference between
the estimated value of $\Delta S_{R}$ between temperatures suggests that
there is another component of the temperature dependence that is not captured
by the inclusion of a single entropic parameter. @Fig:entropy_pairwise shows
that estimating $\Delta S$ from one temperature is not sufficient to predict
the fold-change in gene expression at another temperature. The addition of
the entropic parameter leads to better fit of the 32$^\circ$ C condition than
the simple rescaling of the energy as described by @Eq:epstar
(@Fig:entropy_pairwise, dashed line), but poorly predicts the behavior at
42$^\circ C$. Performing the inference on the combined 32$^\circ$ C and
42$^\circ$ C data strikes a middle ground between the predictions
resulting from the two temperatures alone [@Fig:entropy_pairwise, grey shaded
region].

### Entropy as a Function of Temperature

In the previous section, we made an approximation of the energetic
parameters $\Delta\varepsilon_{R}$ and $\Delta\varepsilon_{AI}$ to be
defined by an enthalpic and entropic term, both of which being
independent of temperature. However, entropy can be (and in many cases
is) dependent on the system temperature, often in a non-trivial manner.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To explore the effects of a temperature-dependent entropy in our
prediction of the fold-change, we perform a Taylor expansion of @Eq:utds about
the entropic parameter with respect to temperature keeping only the
first order term such that 
$$
\Delta S = \Delta S_0 + \Delta S_1 T,
$${#eq:taylor_exp}
where $\Delta S_0$ is a constant, temperature-independent entropic term and
$\Delta S_1$ is the entropic contribution per degree Kelvin. With this simple
relationship enumerated, we can now define the temperature-dependent
effective free energy parameters $\Delta\varepsilon_R^*$ and $\Delta
\varepsilon_{AI}^*$ as
$$
\Delta\varepsilon_{R}^* = ({S_0}_R + {S_1}_{R} T_{exp})( T_{ref} - T_{exp}) -
13.9\, k_BT,  
$${#eq:epRA_star_taylor}
and
$$
\Delta\varepsilon_{AI}^* = ({S_0}_{AI} + {S_1}_{AI} T_{exp})( T_{ref} - T_{exp})
+ 4.5\, k_BT,  
$${#eq:epAI_star_taylor}
respectively, again relying on the constraints defined by @Eq:eprastar_constraint
and @Eq:epaistar_constraint.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Using a similar inferential approach as
described in the previous section, we sample the posterior distribution of
these parameters using Markov chain Monte Carlo and compute the fold-change
and shift in free energy for each temperature. As seen in @Fig:t_dependent,
there is a negligible improvement in the description of the data by including
this temperature dependent entropic parameter.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;These results together suggest that our
understanding of temperature dependence in this regulatory architecture is
incomplete and requires further research from both theoretical and
experimental standpoints.

![**Sampled posterior probability distributions of entropic penalty parameter
inference.** Marginal and joint distributions conditioned only on data
collected at 32$^\circ$ C, only on 42 $^\circ$ C, or on both temperatures are
shown in blue, red, and black, respectively. The values of $\Delta S_R$ and
$\Delta S_{AI}$ are given in $k_BT / K$ where $K$ is 1 degree
Kelvin. The [Python code (`ch8_figS11.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS11.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch8_figS11){#fig:entropy_corner short-caption="Sampled posterior
probability distributions of entropic penalty parameter inference."}


![**Pairwise predictions of fold-change in gene expression at different
temperatures.** Each row represents the condition used to infer the entropic
parameters, and columns are the conditions that are being predicted. Color
shaded regions represent the 95\% credible region of the predicted
fold-change given the inferred values of $\Delta S_R$ and $\Delta S_{AI}$ for
each condition. The black dashed line represents the predicted fold-change in
gene expression by a simple rescaling of the binding energy determined at
37$^circ$ C. The grey shaded region is the 95\% credible region of the
fold-change given estimation of $\Delta S_R$ and $\Delta S_{AI}$ conditioned
on both temperatures pooled together. The [Python code (`ch8_figS12.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS12.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch8_figS12){#fig:entropy_pairwise
short-caption="Pairwise predictions of fold-change in gene expression at
different temperatures."}

![**Fold-change and shift in free energy including a temperature-dependent
entropic contribution.** Top row illustrates the estimated fold-change in
gene expression at 32$^\circ$ C (left) and 42$^\circ$ C (right). Bottom row
shows the estimated shift in free energy for each temperature. The grey
transparent line in all plots shows the relevant quantity including only a
temperature-independent entropic parameter. Purple hashed line illustrates
the relevant quantity including a temperature-dependent entropic
contribution. The width of each curve indicates the 95\% credible region of
the relevant quantity. Points in each correspond to the experimental
measurements. The points in the top row represent the mean and standard error
of five to eight biological replicates. Points in the bottom row correspond
to the median value of the inferred free energy shift and the mean of five to
eight biological replicates for the repressor copy number. Vertical error
represents the upper and lower bounds of the 95\% credible region of the
parameter. Horizontal error bars correspond to the standard error of five to
eight biological replicates. The [Python code (`ch8_figS13.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS13.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch8_figS13){#fig:t_dependent
short-caption="Fold-change and shift in free energy including a
temperature-dependent entropic contribution."}


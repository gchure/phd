## Global Fit of All Parameters

In Chapter 2, we used the repressor copy numbers $R$ and
repressor-DNA binding energies $\Delta\varepsilon_{RA}$ as reported by
@garcia2011. However, any error in these previous measurements of $R$ and
$\Delta\varepsilon_{RA}$ will necessarily propagate into our own
fold-change predictions. In this section, we take an alternative approach
to fitting the physical parameters of the system to that used in the
main text. First, rather than fitting only a single strain, we fit the
entire data set in along with microscopy data for the synthetic operator
Oid. In addition, we also
simultaneously fit the parameters $R$ and $\Delta\varepsilon_{RA}$
using the prior information given by the previous measurements. By using
the entire data set and fitting all of the parameters, we obtain the
best possible characterization of the statistical mechanical parameters
of the system given our current state of knowledge. 

&nbsp;&nbsp;&nbsp;&nbsp;To fit all of the parameters simultaneously, we
follow a similar approach to the one detailed in the Materials \& Methods of
Chapter 2. Briefly, we perform a Bayesian parameter estimation of the
dissociation constants $K_A$ and $K_I$, the six different repressor copy
numbers $R$ corresponding to the six *lacI* ribosomal binding sites used in
our work, and the four different binding energies $\Delta \varepsilon_{RA}$
characterizing the four distinct operators used to make the experimental
strains. As in the main text, we fit the logarithms $\tilde{k}_A = -\log
\frac{K_A}{1\,\text{M}}$ and $\tilde{k}_I = -\log \frac{K_I}{1\,\text{M}}$ of
the dissociation constants which grants better numerical stability.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We assume that deviations of the experimental fold-change from
the theoretical predictions are normally distributed with mean zero and
standard deviation $\sigma$. We begin by writing Bayes’ theorem,
$$
g(\tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
  \mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma \mid D) =
  \frac{f(D \mid \tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
  \mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma) g(\tilde{k}_A,
  \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
  \mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}},
  \sigma)}{f(D)},
$${#eq:bayes_global}
where $\textbf{\textit{R}}$ is an array
containing the six different repressor copy numbers to be fit,
$\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}$ is an
array containing the four binding energies to be fit, and $D$ is the
experimental fold-change data. The term $P(\tilde{k}_A, \tilde{k}_I,
\mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma \mid D)$
gives the probability distributions of all of the parameters given the
data. The prefixes $g$ and $f$ denote probability densities of parameters and
data, respectively. The term
$f(D \mid \tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma)$
represents the likelihood of having observed our experimental data given
some value for each parameter. $g(\tilde{k}_A, \tilde{k}_I,
\mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma)$
contains all the prior information on the values of these parameters.
Lastly, $f(D)$ serves as a normalization constant and is neglected.


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Given $n$ independent measurements of the
fold-change, the first term in @Eq:bayes_global can be written as
$$
\begin{aligned}
    f(D \mid \tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
    \mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma) =
    \frac{1}{(2\pi\sigma^2)^{\frac{n}{2}}}\prod\limits_{i=1}^n \exp
    \left[-\frac{(\text{fc}^{(i)}_{\exp} - \text{fc}(\tilde{k}_A, \tilde{k}_I,
    R^{(i)}, \Delta\varepsilon_{RA}^{(i)}, c^{(i)}))^2}{2\sigma^2}\right],
\end{aligned}
$${#eq:global_likelihood} 
where $\text{fc}^{(i)}_{\text{exp}}$ is the $i^{\text{th}}$ experimental
fold-change and $\text{fc}(\cdot\cdot\cdot)$ is the theoretical
prediction. Note that the standard deviation $\sigma$ of this
distribution is not known and hence needs to be included as a parameter
to be fit.

&nbsp;&nbsp;&nbsp;&nbsp;The second term in @Eq:bayes_global represents the prior information of the parameter
values. We assume that all parameters are independent of each other, such
that $g(\tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma)$ =
$g(\tilde{k}_A ) \cdot P(\tilde{k}_I ) \cdot \prod_i P(R^{(i)}) \cdot \prod_j
g(\Delta\varepsilon_{RA}^{(j)}) \cdot g(\sigma),$ where the superscript
$(i)$ indicates the repressor copy number of index $i$ and the
superscript $(j)$ denotes the binding energy of index $j$. As above,
we note that a prior must also be included for the unknown parameter
$\sigma$.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Because we know nothing about the values of $\tilde{k}_A$,
$\tilde{k}_I$, and $\sigma$ before performing the experiment, we
assign maximally uninformative priors to each of these parameters. More
specifically, we assign uniform priors to $\tilde{k}_A$ and
$\tilde{k}_I$ and a Jeffreys prior to $\sigma$, indicating that
$K_A$, $K_I$, and $\sigma$ are scale parameters [@sivia2006]. We do, however,
have prior information for the repressor copy numbers and the
repressor-DNA binding energies from @garcia2011. This prior knowledge is included
within our model using an informative prior for these two parameters,
which we assume to be Gaussian. Hence each of the $R^{(i)}$ repressor
copy numbers to be fit satisfies 
$$
g(R^{(i)}) = \frac{1}{\sqrt{2\pi\sigma_{R_i}^2}} \exp \left(-\frac{(R^{(i)} -
    \bar{R}^{(i)})^2}{2 \sigma_{R_i}^2} \right),
$${#eq:global_prior_r}
where $\bar{R}^{(i)}$ is the mean repressor copy number and $\sigma_{R_i}$
is the variability associated with this parameter as reported in @garcia2011. Note
that we use the given value of $\sigma_{R_i}$ from previous
measurements rather than leaving this as a free parameter.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Similarly, the binding energies $\Delta\varepsilon_{RA}^{(j)}$ are
also assumed to have a Gaussian informative prior of the same form. We
write it as
$$
g(\Delta\varepsilon_{RA}^{(j)}) = \frac{1}{\sqrt{2\pi\sigma_{\varepsilon_j}^2}}
\exp \left(- \frac{(\Delta\varepsilon_{RA}^{(j)} -
    \Delta\bar{\varepsilon}_{RA}^{(j)})^2}{2 \sigma_{\varepsilon_j}^2} \right),
$${#eq:global_prior_epra}
where $\Delta\bar{\varepsilon}_{RA}^{(j)}$ is the binding energy and
$\sigma_{\varepsilon_j}$ is the variability associated with that
parameter around the mean value as reported in @garcia2011.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The $\sigma_{R_i}$ and $\sigma_{\varepsilon_j}$ parameters will
constrain the range of values for $R^{(i)}$ and
$\Delta\varepsilon_{RA}^{(j)}$ found from the fitting. For example, if
for some $i$ the standard deviation $\sigma_{R_i}$ is very small, it
implies a strong confidence in the previously reported value.
Mathematically, the exponential in @Eq:global_prior_r will ensure that the best-fit
$R^{(i)}$ lies within a few standard deviations of $\bar{R}^{(i)}$.
Since we are interested in exploring which values could give the best
fit, the errors are taken to be wide enough to allow the parameter
estimation to freely explore parameter space in the vicinity of the best
estimates. Putting all these terms together, we use Markov chain Monte
Carlo to sample the posterior distribution
$P(\tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta \boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma \mid
D)$, enabling us to determine both the most likely value for each
physical parameter as well as its associated credible region. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;@Fig:induction_global_fit shows the result of
this global fit. When compared with the results of Chapter 2, we can see that fitting for the binding
energies and the repressor copy numbers improves the agreement between the
theory and the data. Table 6.2 summarizes the values of the parameters as
obtained with this MCMC parameter inference. We note that even though we
allowed the repressor copy numbers and repressor-DNA binding energies to
vary, the resulting fit values were very close to the previously reported
values. The fit values of the repressor copy numbers were all within one
standard deviation of the previous reported values provided in @garcia2011.
And although some of the repressor-DNA binding energies differed by a few
standard deviations from the reported values, the differences were always
less than $1~k_BT$, which represents a small change in the biological scales
we are considering. The biggest discrepancy between our fit values and the
previous measurements arose for the synthetic Oid operator, which we discuss
in more detail in the following section of this chapter.

| **Parameter** | **Reported Values** [@garcia2011] | **Global Fit**|
|:--|:--:|--:|
| $K_A$ | - | $205^{+11}_{-12}\,\mu$M |
| $K_I$ | - | $0.73^{0.04}_{-0.04}\,\mu$M |
| $R_{22}$ | $22 \pm 4$ per cell | $20^{+1}_{-1}$ per cell|
| $R_{60}$ | $60 \pm 20$ per cell | $74^{+4}_{-3}$ per cell|
| $R_{124}$ | $124 \pm 30$ per cell | $130^{+6}_{-6}$ per cell|
| $R_{260}$ | $260 \pm 40$ per cell | $257^{+9}_{-11}$ per cell|
| $R_{1220}$ | $1220 \pm 160$ per cell | $1191^{+32}_{-55}$ per cell|
| $R_{1740}$ | $1740 \pm 340$ per cell | $1599^{+75}_{-87}$ per cell|
| O1 $\Delta\varepsilon_{RA}$ | $-15.3 \pm 0.2\,k_BT$ | $-15.22^{+0.1}_{-0.1}\, k_BT$ | 
| O2 $\Delta\varepsilon_{RA}$ | $-13.9\pm 0.2 \, k_BT$| $-13.06^{+0.1}_{-0.1}\, k_BT$ |
| O3 $\Delta\varepsilon_{RA}$ | $-9.7\pm 0.1\, k_BT$ | $-9.4^{+0.1}_{-0.1}\, k_BT$ |
| Oid $\Delta\varepsilon_{RA}$| $-17.0 \pm 0.2\, k_BT$| $-17.7^{+0.2}_{-0.1}\, k_BT$ |
: Global parameter estimates and comparison to previously reported values. 

![**Global fit of dissociation constants, repressor copy numbers, and binding
energies.** Theoretical prediction resulting from simultaneous estimation of
the dissociation constants $K_A$ and $K_I$, the six repressor copy numbers $R$,
and the four repressor-DNA binding energies $\Delta\varepsilon_{RA}$ using the
entire dataset. Points and errors represent the mean and standard error of ~10
biological replicates for O1, O2, O3, and 3 biological replicates for
Oid. The [Python code (`ch6_fig12.py`)](https://github.com/gchure/phd/blob/master/src/chapter_06/code/ch6_figS12.py) used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch6_figS12){#fig:induction_global_fit short-caption="Induction curves
using global parameter estimates."}


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;@Fig:properties_global_fit shows the same key properties as Fig. 2.7 , but
uses the parameters obtained from this global fitting approach. We note that
even by increasing the number of degrees of freedom in our fit, the result
does not change substantially. or the O3 operator data, again,
agreement between the predicted $[EC_{50}]$ and the effective Hill
coefficient remains poor due to the theory being unable to capture the steepness
of the response curves.

![**Key properties of induction profiles as predicted with a global fit using all
data.** Data for (A) leakiness, (B) saturation, and (C) dynamic range are
computed directly from measured fold-change. Points and errors correspond to the
mean and standard error of 10 - 11 biological replicates. Points in (D) and (E)
for the [EC$_{50}$] and the effective Hill coefficient, respectively, represent
the estimated value using parameter estimates of $K_A$ and $K_I$ for that
particular strain. Errors represent the width of the 95\% credible region. In
all plots, curves represent the theoretical predictions given the parameter
estimates conditioned on all data sets. The [Python code
(`ch6_figS13.py`)](https://github.com/gchure/phd/blob/master/src/chapter_06/code/ch6_figS13.py)
used to generate this figure can be found on the thesis
[GitHub repository](https://github.com/gchure/phd).](ch6_figS13){#fig:properties_global_fit
short-caption="Key properties of induction profiles as predicted with a global
fit using all data."}

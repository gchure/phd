## Inferring Allosteric Parameters from Previous Data {#sec:ch2_si_epAI}

The fold-change profile described by features three unknown parameters
$K_A$, $K_I$, and $\Delta\varepsilon_{AI}$. In this section, we
explore different conceptual approaches to determining these parameters.
We first discuss how the induction titration profile of the simple
repression constructs used in this paper are not sufficient to determine
all three MWC parameters simultaneously, since multiple degenerate sets
of parameters can produce the same fold-change response. We then utilize
an additional data set from Brewster et al. (2014) to determine the
parameter $\Delta\varepsilon_{AI} = 4.5~k_BT$, after which the
remaining parameters $K_A$ and $K_I$ can be extracted from any
induction profile with no further degeneracy.

### Degenerate Parameter Values

In this section, we discuss how multiple sets of parameters may yield
identical fold-change profiles. More precisely, we shall show that if we
try to fit the data in to the fold-change and extract the three unknown
parameters ($K_A$, $K_I$, and $\Delta\varepsilon_{AI}$), then
multiple degenerate parameter sets would yield equally good fits. In
other words, this data set alone is insufficient to uniquely determine
the actual physical parameter values of the system. This problem
persists even when fitting multiple data sets simultaneously as in
Section “”.

In , we fit the $R=260$ data by fixing $\Delta\varepsilon_{AI}$ to
the value shown on the $x$-axis and determine the parameters $K_A$
and $K_I$ given this constraint. We use the fold-change function but
with $\beta
\Delta\varepsilon_{RA}$ modified to the form $\beta
\Delta\tilde{\varepsilon}_{RA}$ in to account for the underlying
assumptions used when fitting previous data (see Section “” for a full
explanation of why this modification is needed).

The best-fit curves for several different values of
$\Delta\varepsilon_{AI}$ are shown in . Note that these fold-change
curves are nearly overlapping, demonstrating that different sets of
parameters can yield nearly equivalent responses. Without more data, the
relationships between the parameter values shown in represent the
maximum information about the parameter values that can be extracted
from the data. Additional experiments which independently measure any of
these unknown parameters could resolve this degeneracy. For example, NMR
measurements could be used to directly measure the fraction
$(1 + e^{-\beta  \Delta\varepsilon_{AI}})^{-1}$ of active repressors in the
absence of IPTG [@gardino2003; @boulton2016]. 

### Computing $\Delta\varepsilon_{AI}$

As shown in the previous section, the fold-change response of a single
strain is not sufficient to determine the three MWC parameters ($K_A$,
$K_I$, and $\Delta\varepsilon_{AI}$), since degenerate sets of
parameters yield nearly identical fold-change responses. To circumvent
this degeneracy, we now turn to some previous data from the *lac* system
in order to determine the value of $\Delta\varepsilon_{AI}$ in for the
induction of the Lac repressor. Specifically, we consider two previous
sets of work from: (1) @garcia2011  and (2) @brewster2014, both of which measured fold-change with the same simple
repression system in the absence of inducer ($c=0$) but at various
repressor copy numbers $R$. The original analysis for both data sets
assumed that in the absence of inducer all of the Lac repressors were in
the active state. As a result, the effective binding energies they
extracted were a convolution of the DNA binding energy
$\Delta\varepsilon_{RA}$ and the allosteric energy difference
$\Delta\varepsilon_{AI}$ between the Lac repressor’s active and
inactive states. We refer to this convoluted energy value as
$\Delta \tilde{\varepsilon}_{RA}$. We first disentangle the
relationship between these parameters in Garcia and Phillips and then
use this relationship to extract the value of $\Delta\varepsilon_{AI}$
from @brewster2014.

Garcia and Phillips determined the total repressor copy numbers $R$ of
different strains using quantitative Western blots. Then they measured
the fold-change at these repressor copy numbers for simple repression
constructs carrying the O1, O2, O3, and Oid *lac* operators integrated
into the chromosome. These data were then fit to the following
thermodynamic model to determine the repressor-DNA binding energies
$\Delta\tilde{\varepsilon}_{RA}$ for each operator,

$$
\text{fold-change}(c=0) = \left(
1+\frac{R}{N_{NS}}e^{-\beta \Delta\tilde{\varepsilon}_{RA}} \right)^{-1}.
$${#eq:ind_si_simple_fc}

Note that this functional form does not exactly match our fold-change in
the limit $c=0$, 

$$
\text{fold-change}(c=0) = \left(
1+\frac{1}{1+e^{-\beta \Delta \varepsilon_{AI}}}\frac{R}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}} \right)^{-1},
$${#eq:ind_si_leakiness}

since it is missing the factor
$\frac{1}{1+e^{-\beta \Delta\varepsilon_{AI}}}$ which specifies what
fraction of repressors are in the active state in the absence of
inducer, 

$$
\frac{1}{1+e^{-\beta \Delta\varepsilon_{AI}}} = p_A(0). 
$${#eq:ind_si_r_active}

In other words, Garcia and Phillips assumed that in the absence of inducer, all
repressors were active. In terms of our notation, the convoluted energy
values $\Delta\tilde{\varepsilon}_{RA}$ extracted by Garcia and
Phillips (namely, $\Delta\tilde{\varepsilon}_{RA}=-15.3~k_B T$ for O1
and $\Delta\tilde{\varepsilon}_{RA}=-17.0~k_B T$ for Oid) represent

$$
\beta \Delta\tilde{\varepsilon}_{RA} = \beta \Delta\varepsilon_{RA} - \log
\left( \frac{1}{1 + e^{-\beta \Delta \varepsilon_{AI}}} \right).
$${#eq:ind_si_tilde_varepsilon}

Note that if $e^{-\beta \Delta \varepsilon_{AI}} \ll 1$, then nearly
all of the repressors are active in the absence of inducer so that
$\Delta\tilde{\varepsilon}_{RA} \approx \Delta\varepsilon_{RA}$. In
simple repression systems where we definitively know the value of
$\Delta
\varepsilon_{RA}$ and $R$, we can use to determine the value of
$\Delta \varepsilon_{AI}$ by comparing with experimentally determined
fold-change values. However, the binding energy values that we use from
@garcia2011 are effective parameters
$\Delta\tilde{\varepsilon}_{RA}$. In this case, we are faced with an
undetermined system in which we have more variables than equations, and
we are thus unable to determine the value of
$\Delta \varepsilon_{AI}$. In order to obtain this parameter, we must
turn to a more complex regulatory scenario which provides additional
constraints that allow us to fit for $\Delta \varepsilon_{AI}$.

A variation on simple repression in which multiple copies of the
promoter are available for repressor binding (for instance, when the
simple repression construct is on plasmid) can be used to circumvent the
problems that arise when using $\Delta \tilde{\varepsilon}_{RA}$. This
is because the behavior of the system is distinctly different when the
number of active repressors $p_A(0) R$ is less than or greater than
the number of available promoters $N$. Repression data for plasmids
with known copy number $N$ allows us to perform a fit for the value of
$\Delta\varepsilon_{AI}$.

To obtain an expression for a system with multiple promoters $N$, we
follow @weinert2014, writing the fold-change in terms of the
the grand canonical ensemble as 

$$
\text{fold-change} = \frac{1}{1 + \lambda_r e^{-\beta \Delta \varepsilon_{RA}}},
$${#eq:ind_si_fugacity}

where $\lambda_r = e^{\beta \mu}$ is the fugacity and $\mu$ is the
chemical potential of the repressor. The fugacity will enable us to
easily enumerate the possible states available to the repressor.

To determine the value of $\lambda_r$, we first consider that the
total number of repressors in the system, $R_{\text{tot}}$, is fixed
and given by 

$$
R_{\text{tot}} = R_S + R_{NS}, 
$${#eq:ind_si_fugacity_rtot}

where $R_S$ represents the number of repressors specifically bound to the
promoter and $R_{NS}$ represents the number of repressors nonspecifically bound
throughout the genome. The value of $R_S$ is given by 

$$
R_S = N \frac{\lambda_r e^{-\beta \Delta \varepsilon_{RA}}}{1 + \lambda_r
e^{-\beta \Delta \varepsilon_{RA}}},
$${#eq:ind_si_fugacity_rs}

where $N$ is the number of available promoters in the cell. Note that
in counting $N$, we do not distinguish between promoters that are on
plasmid or chromosomally integrated provided that they both have the
same repressor-operator binding energy (Weinert et al. 2014). The value
of $R_{NS}$ is similarly give by 

$$
R_{NS} = N_{NS} \frac{\lambda_r}{1 + \lambda_r},
$${#eq:ind_si_fugacity_rns}

where $N_{NS}$ is the number of non-specific sites in the cell (recall that we use
$N_{NS} = 4.6 \times 10^6$ for *E. coli*).

Substituting in into the modified yields the form

$$
p_A(0) R_{\text{tot}} = \frac{1}{1 + e^{-\beta \Delta \varepsilon_{AI}}} \left(
N \frac{\lambda_r e^{-\beta \Delta \varepsilon_{RA}}}{1 + \lambda_r e^{-\beta
\Delta \varepsilon_{RA}}} + N_{NS} \frac{\lambda_r}{1 + \lambda_r} \right),
$${#eq:ind_si_fugacity_pa0}

where we recall from that
$\beta \Delta \varepsilon_{RA} =   \beta \Delta \tilde\varepsilon_{RA} + \log{\left(\frac{1}{1 + e^{-\beta \Delta \varepsilon_{AI}}}\right)}.$
Numerically solving for $\lambda_r$ and plugging the value back into
yields a fold-change function in which the only unknown parameter is
$\Delta \varepsilon_{AI}$.

With these calculations in hand, we can now determine the value of the
$\Delta \varepsilon_{AI}$ parameter. shows how different values of
$\Delta\varepsilon_{AI}$ lead to significantly different fold-change
response curves. Thus, analyzing the specific fold-change response of
any strain with a known plasmid copy number $N$ will fix
$\Delta\varepsilon_{AI}$. Interestingly, the inflection point of
occurs near $p_A(0) R_{\text{tot}} = N$ (as shown by the triangles in
), so that merely knowing where the fold-change response transitions
from concave down to concave up is sufficient to obtain a rough value
for $\Delta\varepsilon_{AI}$. We note, however, that for
$\Delta\varepsilon_{AI} \geq 5\, k_BT$, increasing $\Delta\varepsilon_{AI}$ further does not affect
the fold-change because essentially every repressors will be in the
active state in this regime. Thus, if the $\Delta\varepsilon_{AI}$ is
in this regime, we can only bound it from below.

We now analyze experimental induction data for different strains with
known plasmid copy numbers to determine $\Delta\varepsilon_{AI}$.
shows experimental measurements of fold-change for two O1 promoters with
$N=64$ and $N=52$ copy numbers and one Oid promoter with $N=10$
from Brewster et al. (2014). By fitting these data to , we extracted the
parameter value $\Delta\varepsilon_{AI} = 4.5\,k_B T$. Substituting
this value into shows that 99\% of the repressors are in the active state
in the absence of inducer and $\Delta\tilde{\varepsilon}_{RA} \approx
\Delta\varepsilon_{RA}$, so that all of the previous energies and calculations
made by @garcia2011 and @brewster2014 were accurate.

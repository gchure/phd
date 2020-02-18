## Theoretical Model

This work considers the inducible simple repression regulatory motif
depicted in Fig.
@fig:induction_theory (A) from a
thermodynamic perspective which has been thoroughly dissected and tested
experimentally [@garcia2011; @brewster2014; @razo-mejia2018] and is described in
depth in Chapter 1. The result of this extensive theory-experiment
dialogue is a succinct input-output function schematized in Fig.
@fig:induction_theory (B) that computes the fold-change in gene expression relative to an unregulated promoter.
This function is of the form 

$$
\text{fold-change} = \left(1 + {R_A \over
N_{NS}}e^{-\beta\Delta\varepsilon_{RA}}\right)^{-1},
$${#eq:foldchange}

where $R_A$ is the number of active repressors
per cell, $N_{NS}$ is the number of non-specific binding sites for the
repressor, $\Delta\varepsilon_{RA}$ is the binding energy of the
repressor to its specific binding site relative to the non-specific
background, and $\beta$ is defined as ${1 \over k_B T}$ where
$k_B$ is the Boltzmann constant and $T$ is the temperature. While
this theory requires knowledge of the number of *active* repressors, we
often only know the total number $R$ which is the sum total of active
and inactive repressors. We can define a prefactor $p_\text{act}(c)$
which captures the allosteric nature of the repressor and encodes the
probability a repressor is in the active (repressive) state rather than
the inactive state for a given inducer concentration $c$, namely,
$$
p_\text{act}(c) = {\left(1 + {c \over K_A}\right)^n \over \left(1 + {c \over
K_A}\right)^n + e^{-\beta\Delta\varepsilon_{AI}}\left(1 + {c \over
K_I}\right)^n}.
$${#eq:pact}

Here, $K_A$ and $K_I$ are the dissociation
constants of the inducer to the active and inactive repressor,
$\Delta\varepsilon_{AI}$ is the energetic difference between the
repressor active and inactive states, and $n$ is the number of
allosteric binding sites per repressor molecule ($n=2$ for LacI). With
this in hand, we can define $R_A$ in Eq. @eq:foldchange as $R_A = p_\text{act}(c)
R$.

![](Fig1.pdf)

<span id="fig:induction_theory" label="fig:induction_theory">fig:induction_theory</span>

A key feature of [eq:foldchange](#eq:foldchange) and
[eq:pact](#eq:pact) is that the diverse phenomenology of the gene
expression induction profile can be collapsed onto a single master curve
by rewriting the input-output function in terms of the free energy $F$
also called the Bohr parameter ,
\text{fold-change} = \left(1 + e^{-\beta F}\right)^{-1},
\label{eq:collapse} where
F = -k_BT \log p_\text{act}(c) - k_BT\log\left({R \over N_{NS}}\right) + \Delta\varepsilon_{RA}.
\label{eq:Bohr} Hence, if different combinations of parameters yield
the same free energy, they will give rise to the same fold-change in
gene expression, enabling us to collapse multiple regulatory scenarios
onto a single curve. This can be seen in
Fig. [fig:induction_theory](#fig:induction_theory)(C) where
eighteen unique inducer titration profiles of a LacI simple repression
architecture collected and analyzed in Razo-Mejia *et al.* 2018 
collapse onto a single master curve. The tight distribution about this
curve reveals that the fold-change across a variety of genetically
distinct individuals can be adequately described by a small number of
parameters. Beyond predicting the induction profiles of different
strains, the method of data collapse inspired by
[eq:collapse](#eq:collapse) and [eq:Bohr](#eq:Bohr) can be used
as a tool to identify mechanistic changes in the regulatory architecture
. Similar data collapse approaches have been used previously in such a
manner and have proved vital for distinguishing between changes in
parameter values and changes in the fundamental behavior of the system .

Assuming that a given mutation does not result in a non-functional
protein, it is reasonable to say that any or all of the parameters in
[eq:foldchange](#eq:foldchange) can be affected by the mutation,
changing the observed induction profile and therefore the free energy.
To examine how the free energy of a mutant $F^\mathrm{(mut)}$ differs
from that of the wild-type $F^\mathrm{(wt)}$, we define
$\Delta F = F^\mathrm{(mut)}-
F^\mathrm{(wt)}$, which has the form \begin{aligned}
\Delta F = -k_BT\log\left({p_\text{act}^\mathrm{(mut)}(c) \over  p_\text{act}^\mathrm{(wt)}(c)}\right) &- k_BT \log\left({R^\mathrm{(mut)}\over R^\mathrm{(wt)}}\right)\\
 &+ (\Delta\varepsilon_{RA}^\mathrm{(mut)}- \Delta\varepsilon_{RA}^\mathrm{(wt)}).
    \end{aligned}
    \label{eq:delF}

$\Delta F$ describes how a mutation translates a point across the
master curve shown in Fig.
[fig:induction_theory](#fig:induction_theory)(C). As we will show
in the coming paragraphs (illustrated in Fig.
[fig:deltaF_theory](#fig:deltaF_theory)), this formulation coarse
grains the myriad parameters shown in
[eq:foldchange](#eq:foldchange) and [eq:pact](#eq:pact) into
three distinct quantities, each with different sensitivities to
parametric changes. By examining how a mutation changes the $\Delta F$
as a function of the inducer concentration, one can draw conclusions as
to which parameters have been modified based solely on the shape of the
curve. To help the reader understand how various perturbations to the
parameters tune the free energy, we have hosted an interactive figure on
the dedicated [paper
website](http://www.rpgroup.caltech.edu/mwc_mutants) which makes
exploration of parameter space a simpler task.

The first term in [eq:delF](#eq:delF) is the log ratio of the
probability of a mutant repressor being active relative to the wild type
at a given inducer concentration $c$. This quantity defines how
changes to any of the allosteric parameters – such as inducer binding
constants $K_A$ and $K_I$ or active/inactive state energetic
difference $\Delta\varepsilon_{AI}$ – alter the free energy $F$,
which can be interpreted as the free energy difference between the
repressor bound and unbound states of the promoter. Fig.
[fig:deltaF_theory](#fig:deltaF_theory) (A) illustrates how changes
to the inducer binding constants $K_A$ and $K_I$ alone alter the
induction profiles and resulting free energy as a function of the
inducer concentration. In the limit where $c
= 0$, the values of $K_A$ and $K_I$ do not factor into the
calculation of $p_\text{act}(c)$ given by [eq:pact](#eq:pact),
meaning that $\Delta\varepsilon_{AI}$ is the lone parameter setting
the residual activity of the repressor. Thus, if only $K_A$ and
$K_I$ are altered by a mutation, then $\Delta F$ should be
$0\, k_BT$ when $c = 0$, illustrated by the overlapping red, purple,
and grey curves in the right-hand plot of Fig.
[fig:deltaF_theory](#fig:deltaF_theory)(A). However, if
$\Delta\varepsilon_{AI}$ is influenced by the mutation (either alone
or in conjunction with $K_A$ and $K_I$), the leakiness will change,
resulting in a non-zero $\Delta F$ when $c=0$. This is illustrated
in Fig. [fig:deltaF_theory](#fig:deltaF_theory) (B) where
$\Delta\varepsilon_{AI}$ is the only parameter affected by the
mutation.

It is important to note that for a mutation which perturbs only the
inducer binding constants, the dependence of $\Delta F$ on the inducer
concentration can be non-monotonic. While the precise values of $K_A$
and $K_I$ control the sensitivity of the repressor to inducer
concentration, it is the ratio $K_A / K_I$ that defines whether this
non-monotonic behavior is observed. This can be seen more clearly when
we consider the limit of saturating inducer concentration,
\lim\limits_{c \rightarrow \infty} \log\left({p_\text{act}^\mathrm{(mut)}\over p_\text{act}^\mathrm{(wt)}}\right) \approx \log\left[{1 + e^{-\beta\Delta\varepsilon_{AI}^\mathrm{(wt)}} \left({K_A^\mathrm{(wt)}\over K_I^\mathrm{(wt)}}\right)^n \over 1 + e^{-\beta\Delta\varepsilon_{AI}^\mathrm{(wt)}} \left({K_A^\mathrm{(mut)}\over K_I^\mathrm{(mut)}}\right)^n}\right]
    \label{eq:kaki_sat_c}, which illustrates that $\Delta F$ returns
to zero at saturating inducer concentration when the ratio $K_A / K_I$
is the same for both the mutant and wild-type repressors, so long as
$\Delta\varepsilon_{AI}$ is unperturbed. Non-monotonicity can *only*
be achieved by changing $K_A$ and $K_I$ and therefore serves as a
diagnostic for classifying mutational effects reliant solely on
measuring the change in free energy. A rigorous proof of this
non-monotonic behavior given changing $K_A$ and $K_I$ can be found
in the SI text.

The second term in [eq:delF](#eq:delF) captures how changes in the
repressor copy number contributes to changes in free energy. It is
important to note that this contribution to the free energy change
depends on the total number of repressors in the cell, not just those in
the active state. This emphasizes that changes in the expression of the
repressor are energetically divorced from changes to the allosteric
nature of the repressor. As a consequence, the change in free energy is
constant for all inducer concentrations, as is schematized in Fig.
[fig:deltaF_theory](#fig:deltaF_theory)(C). Because the magnitude
of the change in free energy scales logarithmically with changing
repressor copy number, a mutation which increases expression from 1 to
10 repressors per cell is more impactful from an energetic standpoint
($k_BT \log(10) \approx 2.3\,  k_BT$) than an increase from 90 to 100
($k_BT \log(100/90) \approx 0.1\, k_BT$). Appreciable changes in the
free energy only arise when variations in the repressor copy number are
larger than or comparable to an order of magnitude. Changes of this
magnitude are certainly possible from a single point mutation, as it has
been shown that even synonymous substitutions can drastically change
translation efficiency .

The third and final term in [eq:delF](#eq:delF) is the difference in
the DNA binding energy between the mutant and wild-type repressors. All
else being equal, if the mutated state binds more tightly to the DNA
than the wild type
($\Delta\varepsilon_{RA}^\mathrm{(wt)}> \Delta\varepsilon_{RA}^\mathrm{(mut)}$),
the net change in the free energy is negative, indicating that the
repressor bound states become more energetically favorable due to the
mutation. Much like in the case of changing repressor copy number, this
quantity is independent of inducer concentration and is therefore also
constant Fig. [fig:deltaF_theory](#fig:deltaF_theory)(D).
However, the magnitude of the change in free energy is linear with DNA
binding affinity while it is logarithmic with respect to changes in the
repressor copy number. Thus, to change the free energy by $1\, k_BT$,
the repressor copy number must change by a factor of $\approx 2.3$
whereas the DNA binding energy must change by $1\, k_BT$.

The unique behavior of each quantity in [eq:delF](#eq:delF) and its
sensitivity with respect to the parameters makes $\Delta F$ useful as
a diagnostic tool to classify mutations. Given a set of fold-change
measurements, a simple rearrangement of [eq:collapse](#eq:collapse)
permits the direct calculation of the free energy, assuming that the
underlying physics of the regulatory architecture has not changed. Thus,
it becomes possible to experimentally test the general assertions made
in Fig. [fig:deltaF_theory](#fig:deltaF_theory).

![image](Fig2.pdf)

<span id="fig:deltaF_theory" label="fig:deltaF_theory">fig:deltaF_theory</span>


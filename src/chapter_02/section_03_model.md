---
bibliography: ../references.bib
---

## Theoretical Model 

### Inducible Transcriptional Repression Via The MWC Model of Allostery
We begin by considering a simple repression genetic architecture in which the
binding of an allosteric repressor occludes the binding of RNA polymerase
(RNAP) to the DNA [@ackers1982; @buchler2003]. When an effector (hereafter
referred to as an “inducer" for the case of induction) binds to the
repressor, it shifts the repressor’s allosteric equilibrium towards the
inactive state as specified by the MWC model [@monod1965]. This causes the
repressor to bind more weakly to the operator, which increases gene
expression. Simple repression motifs in the absence of inducer have been
previously characterized by an equilibrium model where the probability of
each state of repressor and RNAP promoter occupancy is dictated by the
Boltzmann distribution [@ackers1982; @buchler2003; @vilar2003; @bintu2005a;
@garcia2011; @brewster2014] [we note that non-equilibrium models of simple
repression have been shown to have the same functional form that we derive
below [@phillips2015]]. We extend these models to consider allostery by
accounting for the equilibrium state of the repressor through the MWC model.

Thermodynamic models of gene expression begin by enumerating all possible
states of the promoter and their corresponding statistical weights. As shown
in Fig. @fig:states_weights (A), the promoter can either be empty, occupied
by RNAP, or occupied by either an active or inactive repressor. The
probability of binding to the promoter will be affected by the protein copy
number, which we denote as $P$ for RNAP, $R_{A}$ for active repressor, and
$R_{I}$ for inactive repressor. We note that repressors fluctuate between the
active and inactive conformation in thermodynamic equilibrium, such that
$R_{A}$ and $R_{I}$ will remain constant for a given inducer concentration
[@monod1965]. We assign the repressor a different DNA binding affinity in the
active and inactive state. In addition to the specific binding sites at the
promoter, we assume that there are $N_{NS}$ non-specific binding sites
elsewhere (i.e. on parts of the genome outside the simple repression
architecture) where the RNAP or the repressor can bind. All specific binding
energies are measured relative to the average non-specific binding energy.
Thus, $\Delta\varepsilon_{P}$ represents the energy difference between the
specific and non-specific binding for RNAP to the DNA. Likewise,
$\Delta\varepsilon_{RA}$ and $\Delta\varepsilon_{RI}$ represent the
difference in specific and non-specific binding energies for repressor in the
active or inactive state, respectively.

![States and weights for the simple repression motif. (A) Occupancy states of
the promoter. RNAP (light blue) and a
repressor compete for binding to a promoter of interest. There are $R_A$
repressors in the active state (red) and $R_I$ repressors in the inactive
state (purple). The difference in energy between a repressor bound to the
promoter of interest versus another non-specific site elsewhere on the DNA
equals $\Delta\varepsilon_{RA}$ in the active state and
$\Delta\varepsilon_{RI}$ in the inactive state; the $P$ RNAP have a
corresponding energy difference $\Delta\varepsilon_{P}$ relative to
non-specific binding on the DNA. $N_{NS}$ represents the number of
non-specific binding sites for both RNAP and repressor. (B) Allosteric states of
the repressor. A repressor has an
active conformation (red, left column) and an inactive conformation (purple,
right column), with the energy difference between these two states given by
$\Delta\varepsilon_{AI}$. The inducer (orange circle) at concentration $c$ is
capable of binding to the repressor with dissociation constants $K_A$ in the
active state and $K_I$ in the inactive state. The eight states for a dimer
with $n=2$ inducer binding sites are shown along with the sums of the active
and inactive states.](ch2_fig2){#fig:states_weights short-caption="States and
statistical weights for the simple repression motif."}

Thermodynamic models of transcription [@ackers1982; @buchler2003; @vilar2003;
@bintu2005; @bintu2005a; @kuhlman2007; @daber2011a; @garcia2011;
@brewster2014; @weinert2014] posit that gene expression is proportional to
the probability that the RNAP is bound to the promoter $p_\text{bound}$,
which is given by

$$
p_\text{bound} = \frac{\frac{P}{N_{NS}}e^{-\beta
\Delta\varepsilon_{P}}}{1+\frac{R_A}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}}+\frac{R_I}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RI}}+\frac{P}{N_{NS}}e^{-\beta\Delta\varepsilon_{P}}},
$${#eq:pbound_def}

with $\beta = 1/k_BT$ where $k_B$ is the Boltzmann constant
and $T$ is the temperature of the system. As $k_BT$ is the natural
unit of energy at the molecular length scale, we treat the products
$\beta \Delta\varepsilon_{j}$ as single parameters within our model.
Measuring $p_\text{bound}$ directly is fraught with experimental
difficulties, as determining the exact proportionality between
expression and $p_\text{bound}$ is not straightforward. Instead, we
measure the fold-change in gene expression due to the presence of the
repressor. We define fold-change as the ratio of gene expression in the
presence of repressor relative to expression in the absence of repressor
(i.e. constitutive expression), namely,

$$
\text{fold-change} \equiv \frac{p_\text{bound}(R > 0)}{p_\text{bound}(R = 0)}.
$${#eq:fold_change_definition}

We can simplify this expression using two well-justified approximations:
firstly, $(P / N_{NS})e^{-\beta\Delta\varepsilon_{P}}\ll$ 1 implying that the
RNAP binds weakly to the promoter ($N_{NS} = 4.6 \times 10^6$, $P \approx
10^3$ [@klumpp2008], $\Delta\varepsilon_{P} \approx -2\,\, \text{to} \, -5\,
k_BT$ [@brewster2012], so that $(P/N_{NS})e^{-\beta\Delta\varepsilon_{P}}
\approx 0.01$) and (2) $(R_I/N_{NS})e^{-\beta \Delta\varepsilon_{RI}} \ll 1 +
(R_A /N_{NS}) e^{-\beta\Delta\varepsilon_{RA}}$ which reflects our assumption
that the inactive repressor binds weakly to the promoter of interest. Using
these approximations, the fold-change reduces to the form

$$
\text{fold-change} \approx \left(1+\frac{R_A}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}}\right)^{-1} \equiv \left( 1+p_\text{act}(c)
\frac{R}{N_{NS}}e^{-\beta\Delta\varepsilon_{RA}} \right)^{-1},
$${#eq:fold_change_RA}

 where in the last step we
have introduced the fraction $p_\text{act}(c)$ of repressors in the active
state given a concentration $c$ of inducer, such that
$R_A(c)=p_\text{act}(c) R$. Since inducer binding shifts the repressors from
the active to the inactive state, $p_\text{act}(c)$ grows smaller as $c$
increases .

We use the MWC model to compute the probability $p_\text{act}(c)$ that a
repressor with $n$ inducer binding sites will be active. The value of
$p_\text{act}(c)$ is given by the sum of the weights of the active repressor
states divided by the sum of the weights of all possible repressor
states [see Fig. @fig:states_weights (B)], namely, 

$$
p_\text{act}(c)=\frac{\left(1+\frac{c}{K_A}\right)^n}{\left(1+\frac{c}{K_A}\right)^n+e^{-\beta
\Delta \varepsilon_{AI} }\left(1+\frac{c}{K_I}\right)^n},
$${#eq:eq_pactive}

where $K_A$ and $K_I$ represent the dissociation constant between the inducer
and repressor in the active and inactive states, respectively, and $\Delta
\varepsilon_{AI} = \varepsilon_{I} - \varepsilon_{A}$ is the free energy
difference between a repressor in the inactive and active state [the quantity
$e^{-\Delta \varepsilon_{AI}}$ is sometimes denoted by $L$ [@monod1965;
@marzen2013] or $K_{\text{RR}*}$ [@daber2011a]]. In this equation,
$c/K_A$ and $c/K_I$ represent the change in free energy when
an inducer binds to a repressor in the active or inactive state,
respectively, while $e^{-\beta \Delta \varepsilon_{AI}}$ represents the
change in free energy when the repressor changes from the active to inactive
state in the absence of inducer. Thus, a repressor which favors the active
state in the absence of inducer ($\Delta \varepsilon_{AI} > 0$) will be
driven towards the inactive state upon inducer binding when $K_I < K_A$. The
specific case of a repressor dimer with $n=2$ inducer binding sites is shown
in Fig. @fig:states_weights (B).

Substituting $p_\text{act}(c)$ from into yields the general formula for
induction of a simple repression regulatory architecture [@phillips2015],
namely,
$$
\text{fold-change} =
\left(1+\frac{\left(1+\frac{c}{K_A}\right)^n}{\left(1+\frac{c}{K_A}\right)^n+e^{-\beta
\Delta \varepsilon_{AI}
}\left(1+\frac{c}{K_I}\right)^n}\frac{R}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}} \right)^{-1}.
$${#eq:fold_change_full}

While we have used the specific case of simple repression with induction
to craft this model, the same mathematics describe the case of
corepression in which binding of an allosteric effector stabilizes the
active state of the repressor and decreases gene expression (see Fig. @fig:inducible_types).
Interestingly, we shift from induction (governed by $K_I < K_A$) to
corepression ($K_I > K_A$) as the ligand transitions from
preferentially binding to the inactive repressor state to stabilizing
the active state. Furthermore, this general approach can be used to
describe a variety of other motifs such as activation, multiple
repressor binding sites, and combinations of activator and repressor
binding sites [@bintu2005a; @brewster2014; @weinert2014].

The formula presented in Eq. @eq:fold_change_full enables us to make precise
quantitative statements about induction profiles. Motivated by the broad
range of predictions implied by Eq. @eq:fold_change_full, we designed a
series of experiments using the *lac* system in *E. coli* to tune the control
parameters for a simple repression genetic circuit. As discussed in Fig.
@fig:inducible_types (C), previous studies from our lab have provided
well-characterized values for many of the parameters in our experimental
system, leaving only the values of the the MWC parameters ($K_A$, $K_I$, and
$\Delta \varepsilon_{AI}$) to be determined. We note that while previous
studies have obtained values for $K_A$, $K_I$, and $L=e^{-\beta \Delta
\varepsilon_{AI}}$ [@ogorman1980; @daber2011a], they were either based upon
biochemical experiments or *in vivo* conditions involving poorly
characterized transcription factor copy numbers and gene copy numbers. These
differences relative to our experimental conditions and fitting techniques
led us to believe that it was important to perform our own analysis of these
parameters. After inferring these three MWC parameters (see , Section “” for
details regarding the inference of $\Delta \varepsilon_{AI}$, which was
fitted separately from $K_A$ and $K_I$), we were able to predict the
input/output response of the system under a broad range of experimental
conditions. For example, this framework can predict the response of the
system at different repressor copy numbers $R$, repressor-operator affinities
$\Delta\varepsilon_{RA}$, inducer concentrations $c$, and gene copy numbers
(see Chapter @sec:ch2_si, Sec. @sec:induction_fugacity).

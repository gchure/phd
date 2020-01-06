---
pagetitle: none
---

## Results

### Characterizing Transcription Factor Induction using the Monod-Wyman-Changeux (MWC) Model

We begin by considering a simple repression genetic architecture in
which the binding of an allosteric repressor occludes the binding of RNA
polymerase (RNAP) to the DNA . When an effector (hereafter referred to
as an “inducer" for the case of induction) binds to the repressor, it
shifts the repressor’s allosteric equilibrium towards the inactive state
as specified by the MWC model . This causes the repressor to bind more
weakly to the operator, which increases gene expression. Simple
repression motifs in the absence of inducer have been previously
characterized by an equilibrium model where the probability of each
state of repressor and RNAP promoter occupancy is dictated by the
Boltzmann distribution  (we note that non-equilibrium models of simple
repression have been shown to have the same functional form that we
derive below ). We extend these models to consider allostery by
accounting for the equilibrium state of the repressor through the MWC
model.

Thermodynamic models of gene expression begin by enumerating all
possible states of the promoter and their corresponding statistical
weights. As shown in [@fig:ind_states_weights], the promoter can either be empty, occupied by
RNAP, or occupied by either an active or inactive repressor. The
probability of binding to the promoter will be affected by the protein
copy number, which we denote as $P$ for RNAP, $R_{A}$ for active
repressor, and $R_{I}$ for inactive repressor. We note that repressors
fluctuate between the active and inactive conformation in thermodynamic
equilibrium, such that $R_{A}$ and $R_{I}$ will remain constant for
a given inducer concentration . We assign the repressor a different DNA
binding affinity in the active and inactive state. In addition to the
specific binding sites at the promoter, we assume that there are
$N_{NS}$ non-specific binding sites elsewhere (i.e. on parts of the
genome outside the simple repression architecture) where the RNAP or the
repressor can bind. All specific binding energies are measured relative
to the average non-specific binding energy. Thus,
$\Delta\varepsilon_{P}$ represents the energy difference between the
specific and non-specific binding for RNAP to the DNA. Likewise,
$\Delta\varepsilon_{RA}$ and $\Delta\varepsilon_{RI}$ represent the
difference in specific and non-specific binding energies for repressor
in the active or inactive state, respectively.

![States and weights for the simple repression motif.** RNAP (light
blue) and a repressor compete for binding to a promoter of interest.
There are $R_A$ repressors in the active state (red) and $R_I$
repressors in the inactive state (purple). The difference in energy
between a repressor bound to the promoter of interest versus another
non-specific site elsewhere on the DNA equals $\Delta\varepsilon_{RA}$
in the active state and $\Delta\varepsilon_{RI}$ in the inactive
state; the $P$ RNAP have a corresponding energy difference
$\Delta\varepsilon_{P}$ relative to non-specific binding on the DNA.
$N_{NS}$ represents the number of non-specific binding sites for both
RNAP and repressor. A repressor has an active conformation (red, left
column) and an inactive conformation (purple, right column), with the
energy difference between these two states given by $\Delta
\varepsilon_{AI}$. The inducer (blue circle) at concentration $c$
is capable of binding to the repressor with dissociation constants
$K_A$ in the active state and $K_I$ in the inactive state. The eight
states for a dimer with $n=2$ inducer binding sites are shown along
with the sums of the active and inactive
states.](ch2_fig2){#fig:ind_states_weights short-caption="States and weights for the simple repression motif."}


Thermodynamic models of transcription  posit that gene expression is
proportional to the probability that the RNAP is bound to the promoter
$p_{\text{bound}}$, which is given by 

$$
p_\text{bound}=\frac{\frac{P}{N_{NS}}e^{-\beta
\Delta\varepsilon_{P}}}{1+\frac{R_A}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}}+\frac{R_I}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RI}}+\frac{P}{N_{NS}}e^{-\beta\Delta\varepsilon_{P}}},
$$ {#eq:pbound_def}

with $\beta = \frac{1}{k_BT}$ where $k_B$ is the Boltzmann constant
and $T$ is the temperature of the system. As $k_BT$ is the natural
unit of energy at the molecular length scale, we treat the products
$\beta \Delta\varepsilon_{j}$ as single parameters within our model.
Measuring $p_{\text{bound}}$ directly is fraught with experimental
difficulties, as determining the exact proportionality between
expression and $p_{\text{bound}}$ is not straightforward. Instead, we
measure the fold-change in gene expression due to the presence of the
repressor. We define fold-change as the ratio of gene expression in the
presence of repressor relative to expression in the absence of repressor
(i.e. constitutive expression), namely,

$$
\text{fold-change} \equiv \frac{p_\text{bound}(R > 0)}{p_\text{bound}(R = 0)}.
$${#eq:ind_fc_def}

We can simplify this expression using two well-justified approximations:
(1) $\frac{P}{N_{NS}}e^{-\beta\Delta\varepsilon_{P}}\ll 1$ implying
that the RNAP binds weakly to the promoter
($N_{NS} = 4.6 \times 10^6$, $P \approx 10^3$ ,
$\Delta\varepsilon_{P} \approx -2 \,\, \text{to} \, -5~k_B
T$ , so that $\frac{P}{N_{NS}}e^{-\beta\Delta\varepsilon_{P}}
\approx 0.01$) and (2)
$\frac{R_I}{N_{NS}}e^{-\beta \Delta\varepsilon_{RI}} \ll
1 + \frac{R_A}{N_{NS}} e^{-\beta\Delta\varepsilon_{RA}}$ which reflects
our assumption that the inactive repressor binds weakly to the promoter
of interest. Using these approximations, the fold-change reduces to the
form 

$$
\text{fold-change} \approx \left(1+\frac{R_A}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}}\right)^{-1} \equiv \left( 1+p_A(c)
\frac{R}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}} \right)^{-1},
$$ {#eq:ind_fc_appx}

where in the last step we have introduced the fraction $p_\text{act}(c)$ of
repressors in the active state given a concentration $c$ of inducer, such
that $R_A(c)=p_\text{act}(c) R$. Since inducer binding shifts the repressors
from the active to the inactive state, $p_\text{act}(c)$ grows smaller as $c$
increases.


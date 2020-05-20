## Generalizability of Data Collapse to Other Regulatory Architectures

In Chapters 2, 3, and 4, we have stated that the input-output function for
the fold-change in gene expression can be rewritten in the form of a Fermi
function. However, this result can be derived directly by coarse graining the
transcription factor bound and unbound states of the systems' partition
function. In this section, we show how coarse graining of promoter occupancy
states results in a general Fermi function so long as the transcription
factor bound states and transcriptionally active states do not overlap.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As shown in Chapter 2, the partition function $\mathcal{Z}$ for the simple
repression motif (ignoring allosteric control) can be enumerated as
$$
\mathcal{Z} = \overbrace{1 + {P \over
N_{NS}}e^{-\beta\Delta\varepsilon_P}}^\text{repressor unbound states} +
\underbrace{{R \over N_{NS}}e^{-\beta\Delta\varepsilon_{R}}}_\text{repressor
bound state} = \mathcal{Z}_{\neg r} + \mathcal{Z}_{r},
$${#eq:z_simple_rep} 
where $R$ is the number of repressors, $P$ is
the number of polymerases, $N_{NS}$ is the number of nonspecific binding
sites, and $\Delta\varepsilon_P$ and $\Delta\varepsilon_R$ are the
binding energies of the polymerase and repressor to the DNA,
respectively. The states can be grouped into either repressor bound
states (right-hand terms) or repressor unbound states (left hand terms),
denoted as $\mathcal{Z}_{r}$ and $\mathcal{Z}_{\neg r}$, respectively.
The probability of the repressor *not* being bound to the promoter can
be computed as
$$
P(\neg r) = {\mathcal{Z}_{\neg r} \over \mathcal{Z}_{\neg r} + \mathcal{Z}_r}.
$${#eq:p_notr}

&nbsp;&nbsp;&nbsp;&nbsp;In this coarse-grained description, the
transcriptionally active states are separate from the repressor bound
states. Thus, to calculate the fold-change in gene expression, we can
compute the ratio of $P(\neg r \vert R > 0)$ when repressor is present
to $P(\neg r\,\vert\, R = 0)$ when repressor is absent. As the latter
term is equal to 1, the fold-change in gene expression is equivalent to
@Eq:p_notr. We can compute an *effective* free energy
$\tilde{F}$ of each coarse-grained state as
$$
\tilde{F}_{\neg r} = -k_BT \log \mathcal{Z}_{\neg r}\, ;\, \tilde{F}_{r} = -k_BT \log \mathcal{Z}_r.
$${#eq:eff_energy} 
With an effective energy in hand, we can write @Eq:p_notr (which is equivalent to the fold-change) in terms of the Boltzmann weights of these two coarse-grained states as
$$
\text{fold-change} = {e^{-\beta \tilde{F}_{\neg r}} \over e^{-\beta \tilde{F}_{\neg r}} + e^{-\beta\tilde{F}_r}}.
$${#eq:fc_boltz} 
A simple rearrangement of this result produces a Fermi function
$$
\text{fold-change} = {1 \over 1 + e^{-\beta \tilde{F}}},
$${#eq:fermi_fn}
where $\tilde{F}$ is the effective free energy of the system defined as
$$
\tilde{F} = -k_BT\left( \log \mathcal{Z}_r - \log \mathcal{Z}_{\neg r}\right).
$${#eq:f_def}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In the case of the simple repression motif, we apply the weak promoter
approximation which is based on the fact that
$(P/N_{NS})e^{-\beta\Delta\varepsilon_P} \ll 1$, reducing
$\mathcal{Z}_{\neg r} = 1$. With this approximation, the effective free
energy $\tilde{F}$ defined in @Eq:f_def can be interpreted as the free energy between the repressed and transcriptionally active states. In general, this coarse-graining approach can be applied to any regulatory architecture in which the transcription factor bound states and transcriptionally active states
share no overlapping substates. @Fig:twostate_archs shows various repression based
architectures along with the coarse graining of states and the effective
free energy.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This approach cannot be applied to
architectures in which the transcription factor bound and transcriptionally
active states cannot be separated. One such example is the case of simple
activation in which the binding of an activator increases gene expression.
For this architecture, the total partition function [as derived in @bintu2005] is
$$
\mathcal{Z} = \overbrace{1 + {P \over N_{NS}}e^{-\beta\Delta\varepsilon_{P}}}^\text{activator unbound states} +
              \underbrace{{A \over N_{NS}}e^{-\beta\Delta\varepsilon_A} + {A \over N_{NS}}{P \over N_{NS}}e^{-\beta(\Delta\varepsilon_{A} + 
              \Delta\varepsilon_{P} + \varepsilon_\text{interaction})}}_\text{activator bound states} = \mathcal{Z}_{\neg a} + \mathcal{Z}_a.
$${#eq:act_zed}
Here, we have used $A$ to denote the number of activators, $\Delta\varepsilon_{A}$ is the binding energy of the activator to its specific binding site, and
$\varepsilon_\text{interaction}$ is the interaction energy between the
activator and polymerase which makes the activator-polymerase-DNA state
more energetically favorable to the polymerase-DNA state. Unlike in the
case of repression, the transcriptionally active states are present in
the $\mathcal{Z}_a$. While we can compute the probability of the
activator not being bound $P(\neg a)$, this is not equivalent to the
fold-change in gene expression as was the case for simple repression.
Rather, the fold-change in gene expression can be written as
$$
\text{fold-change} = {\mathcal{Z}_a -1 + \mathcal{Z}_{\neg a} - {A \over N_{NS}}e^{-\beta\Delta\varepsilon_{A}} \over \mathcal{Z}_{a} + \mathcal{Z}_{\neg a}}{\mathcal{Z}_a \over \mathcal{Z}_a - 1}.
$${#eq:unlabbed}
This is a more cumbersome expression than in the case for simple
repression, and cannot be massaged into a one-parameter description of
the fold-change. Thus, the analysis presented here is only applicable to
cases in which the occupancy of the promoter by either the polymerase or
the transcription factor are mutually exclusive.

![**Various repression-based regulatory architectures and their
coarse-grained states.** The left column shows a schematic of the
regulatory architecture. The middle column shows the equilibrium states
of the system coarse grained into repressor bound (orange boxes) and
repressor unbound (purple boxes) states. States in which the polymerase is
bound is shaded in white and are neglected by the weak promoter
approximation. The right column illustrates that the formulation of the
effective free energy is the same for all architectures, although the
formulation of $\mathcal{Z}_{\neg r}$ and $\mathcal{Z}_r$ are different
between the examples. The bottom row illustrates the coarse graining of
a simple repression motif in which the same repressor regulates multiple
copies of the same gene. Rather than considering all states, the gene
expression can be calculated by considering one promoter and an effective
copy number, described by the repressor fugacity $\lambda_r$, and is
derived in @weinert2014. Grayed-out states illustrate this further
level of coarse graining.](ch7_figS23){#fig:twostate_archs short-caption="Various repression-based regulatory architectures and their coarse-grained states."}

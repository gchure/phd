## Discussion

Allosteric regulation is often couched as “biological action at a
distance". Despite extensive knowledge of protein structure and
function, it remains difficult to translate the coordinates of the
atomic constituents of a protein to the precise parameter values which
define the functional response, making each mutant its own intellectual
adventure. Bioinformatic approaches to understanding the
sequence-structure relationship have permitted us to examine how the
residues of allosteric proteins evolve, revealing conserved regions
which hint to their function. Co-evolving residues reveal sectors of
conserved interactions which traverse the protein that act as the
allosteric communication channel between domains [@suel2002; @mclaughlin2012; @reynolds2011]. Elucidating these
sectors has advanced our understanding of how distinct domains "talk" to
one another and has permitted direct engineering of allosteric responses
into non-allosteric enzymes [@raman2014; @poelwijk2011a; @raman2016]. Even so, we are left without a
quantitative understanding of how these admittedly complex networks set
the energetic difference between active and inactive states or how a
given mutation influences binding affinity. In this context, a
biophysical model in which the various parameters are intimately
connected to the molecular details can be of use and can lead to
quantitative predictions of the interplay between amino-acid identity
and system-level response.

By considering how each parameter contributes to the observed change in
free energy, we are able to tease out different classes of parameter
perturbations which result in stereotyped responses to changing inducer
concentration. These characteristic changes to the free energy can be
used as a diagnostic tool to classify mutational effects. For example,
we show in Fig. @fig:deltaF_theory that
modulating the inducer binding constants $K_A$ and $K_I$ results in
non-monotonic free energy changes that are dependent on the inducer
concentration, a feature observed in the inducer binding mutants
examined in this work. Simply looking at the inferred $\Delta F$ as a
function of inducer concentration, which requires no fitting of the
biophysical parameters, indicates that $K_A$ and $K_I$ must be
modified considering those are the only parameters which can generate
such a response.

Another key observation is that a perturbation to only $K_A$ and
$K_I$ requires that the $\Delta F = 0\, k_BT$ at $c = 0$. Deviations from
this condition imply that more than the inducer binding constants must
have changed. If this shift in $\Delta F$ off of $0\, k_BT$ at $c = 0$ is
not constant across all inducer concentrations, we can surmise that the
energy difference between the allosteric states
$\Delta\varepsilon_{AI}$ must also be modified. We again see this
effect for all of our inducer mutants. By examining the inferred
$\Delta F$, we can immediately say that in addition to $K_A$ and
$K_I$, $\Delta\varepsilon_{AI}$ must decrease relative to the
wild-type value as $\Delta F > 0$ at $c = 0$. When the allosteric
parameters are fit to the induction profiles, we indeed see that this is
the case, with all four mutations decreasing the energy gap between the
active and inactive states. Two of these mutations, Q291R and Q291K,
make the inactive state of the repressor *more* stable than the active
state, which is not the case for the wild-type repressor [@razo-mejia2018].

Our formulation of $\Delta F$ indicates that shifts away from $0\, k_BT$
that are independent of the inducer concentration can only arise from
changes to the repressor copy number and/or DNA binding specificity,
indicating that the allosteric parameters are untouched. We see that for
three mutations in the DNA binding domain, $\Delta F$ is the same
irrespective of the inducer concentration. Measurements of $\Delta F$
for these mutants with repressor copy numbers across three orders of
magnitude yield approximately the same value, revealing that
$\Delta\varepsilon_{RA}$ is the sole parameter altered via the
mutations.

We note that the conclusions stated above can be qualitatively drawn
without resorting to fitting various parameters and measuring the
goodness-of-fit. Rather, the distinct behavior of $\Delta F$ is
sufficient to determine which parameters are changing. Here, these
conclusions are quantitatively confirmed by fitting these parameters to
the induction profile, which results in accurate predictions of the
fold-change and $\Delta F$ for nearly every strain across different
mutations, repressor copy numbers, and operator sequence, all at
different inducer concentrations. With a collection of evidence as to
what parameters are changing for single mutations, we put our model to
the test and drew predictions of how double mutants would behave both in
terms of the titration curve and free energy profile.

A hypothesis that arises from our formulation of $\Delta F$ is that a
simple summation of the energetic contribution of each mutation should
be sufficient to predict the double mutants (so long as they are in
separate domains). We find that such a calculation permits precise and
accurate predictions of the double mutant phenotypes, indicating that
there are no epistatic interactions between the mutations examined in
this work. With an expectation of what the free energy differences
should be, epistatic interactions could be understood by looking at how
the measurements deviate from the prediction. For example, if epistatic
interactions exist which appear as a systematic shift from the predicted
$\Delta F$ independent of inducer concentration, one could conclude
that DNA binding energy is not equal to that of the single mutation in
the DNA binding domain alone. Similarly, systematic shifts that are
dependent on the inducer concentration (i.e. not constant) indicate that
the allosteric parameters must be influenced. If the expected difference
in free energy is equal to $0\, k_BT$ when $c=0$, one could surmise that
the modified parameter must not be $\Delta\varepsilon_{AI}$ nor
$\Delta\varepsilon_{RA}$ as these would both result in a shift in
leakiness, indicating that $K_A$ and $K_I$ are further modified.

Ultimately, we present this work as a proof-of-principle for using
biophysical models to investigate how mutations influence the response
of allosteric systems. We emphasize that such a treatment allows one to
boil down the complex phenotypic responses of these systems to a
single-parameter description which is easily interpretable as a free
energy. The general utility of this approach is illustrated in Fig.
@fig:all_data_collapse where gene
expression data from previous work  along with all of the measurements
presented in this work collapse onto the master curve defined by
Eq. @eq:collapse. While our model coarse grains many of
the intricate details of transcriptional regulation into two states (one
in which the repressor is bound to the promoter and one where it is
not), it is sufficient to describe a swath of regulatory scenarios. As
discussed in the supplemental Chapter 7, any architecture in which the
transcription-factor bound and transcriptionally active states of the
promoter can be separated into two distinct coarse-grained states can be
subjected to such an analysis.

Given enough parametric knowledge of the system, it becomes possible to
examine how modifications to the parameters move the physiological
response along this reduced one-dimensional parameter space. This
approach offers a glimpse at how mutational effects can be described in
terms of energy rather than Hill coefficients and arbitrary prefactors.
While we have explored a very small region of sequence space in this
work, coupling of this approach with high-throughput sequencing-based
methods to query a library of mutations within the protein will shed
light on the phenotypic landscape centered at the wild-type sequence.
Furthermore, pairing libraries of protein and operator sequence mutants
will provide insight as to how the protein and regulatory sequence
coevolve, a topic rich with opportunity for a dialogue between theory
and experiment.

![**Data collapse of the simple repression regulatory architecture. All
data are means of biological replicates.** Where present, error bars
correspond to the standard error of the mean of five to fifteen
biological replicates. Red triangles indicate data from Garcia and
Phillips  obtained by colorimetric assays. Blue squares are data from
Brewster et al. acquired from video microscopy. Green circles are data
from Razo-Mejia et al.  obtained via flow cytometry. All other symbols
correspond to the work presented here. An interactive version of this
figure can be found on the [paper
website](https://www.rpgroup.caltech.edu/mwc_mutants) where the
different data sets can be viewed in more
detail.](ch3_fig6){#fig:all_data_collapse short-caption="Data collapse of the
simple repression regulatory architecture. All data are means of biological
replicates."}
## The Janus Face of Molecules

Monod is perhaps most famous for his discovery of allostery, to which he
famously referred  to as "the second secret of
life" [@ullmann2011; @monod1965]. It is fair to say that this
"secret" has been now been declassified. Allosteric regulation can be found
in all domains of life across varied types of biological processes. Allostery
can be found governing the behavior of ion channels [@einav2017,
@auerbach2012a], enzymatic reactions [@einav2016], chemotaxis [@keymer2006],
G-protein coupled receptors [@canals2012], quorum sensing [@swem2008], and
transcriptional regulation [@huang2018, @lindsley2006a], to name a few of
many examples. Despite the objective complexity in the molecular structures
of all of these allosteric molecules, they can be frequently be reduced to
simple cartoons where the details of conformational changes, substrate
binding affinities, and more can be massaged into a small set of key details.
@Fig:allostery (A) shows the molecular structures of a variety of allosteric
transcriptional repressors (top). While each has their own fascinating
structure and continuum of conformational states, all can be coarse grained
into a simple cartoon representation (bottom) with an active (red) and
inactive (purple) state, both of which possess binding pockets
(semicircular notches) for an inducer molecule (orange).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Much as we can reduce the complexity of allosteric molecules schematically,
we can enumerate simple mathematical models that describe their behavior.
Thermodynamic models built on an assumption of quasi-equilibrium are
routinely used to describe complex biological phenomena despite the reality
that being in thermodynamic equilibrium is synonymous with being dead.
Even with this glaring assumption, such models have been shown to be
exceptionally predictive for a variety of complex systems, especially in
modeling molecular binding reactions [@dill2010] and allostery writ large
[@swem2008; @keymer2006; @einav2017; @einav2016; @phillips2015]. As the timescales of
binding and unbinding reactions are orders of magnitude smaller than that of
many other processes in the cell, it is fair to make the approximation that
molecular binding is in equilibrium. Under this assumption, we are granted
the powerful mathematical privilege to say that the probability of a given state
of the system $P_\text{state}$ follows a Boltzmann distribution,
$$
P_\text{state} = \frac{e^{-\frac{\epsilon_\text{state}}{k_BT}}}{\mathcal{Z}},
$${#eq:boltzmann}
where $\epsilon_\text{state}$ is the energy of that state, $k_B$ is the Boltzmann constant,
and $T$ is the system temperature. The denominator $\mathcal{Z}$ is the
partition function of the system and is the sum 
$$
\mathcal{Z} = \sum\limits_{i \in \text{states}} e^{\frac{\epsilon_i}{k_BT}},
$${#eq:partition_function_def}
ensuring that the distribution is normalized. Therefore, if we are interested in computing
the probability of a given allosteric protein being in the active state, we
merely  have to enumerate all of the Boltzmann weights (given by the numerator in @Eq:boltzmann)
and compute
$$
P_\text{active}  = \frac{\text{sum over all possible active states}}{\text{sum
over all possible states}}.
$${#eq:verbal_pact}
This probability, defined as a function of the inducer concentration, is shown
schematically in @Fig:allostery (B). While we have passed over some of the more
subtle details of this calculation, the plot in @Fig:allostery (B) presents a
*quantitative* prediction of how the activity of an allosteric molecule should
scale as a function of the inducer,in this case becoming less active as more
inducer is present).

![**A Coarse grained representation of an allosteric molecule.** (A) Crystal
structures of a variety of allosteric transcription factors are shown at the
top. In this thesis, we coarse grain away many of the details to a minimal model
(bottom) where the protein can be represent as being either active (red) or 
inactive (purple), both of which can bind an inducer molecule (orange).
(B) By making an assumption of quasi-equilibrium, we can trivially compute a
mathematical description of the active probability of an allosteric protein as a
function of the inducer concentration (top). In this particular case, the inactive
state becomes more probable relative to the active state at higher concentrations
of inducer molecule.](ch1_fig3){#fig:allostery short-caption="A Coarse grained
representation of an allosteric molecule."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In **Chapter 2** and the associated
supplementary
**Chapter 6** of this dissertation, we use the Monod-Wyman-Changeux model of
allostery [@monod1965] to build a predictive model of transcriptional regulation
where the level of gene expression changes
in response to changing activity of an allosteric transcriptional repressor.
Using the same tricks given by @Eq:boltzmann and @Eq:verbal_pact, we expand 
upon a previously characterized thermodynamic model of the simple repression
motif. This motif, schematized in @Fig:induction_intro (A) is not just a
convenient abstraction of a regulatory architecture. Rather, this motif is
the most ubiquitous regulatory scheme in *E. coli* [@gama-castro2016;
@ireland2020] and has been the target of much theoretical and
experimental dissection [@vilar2003; @buchler2003; @bintu2005; @garcia2011;
@brewster2014; @phillips2019]. However inclusion of allostery in a mathematical
sense had yet to be experimentally dissected.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;At the beginning of 2016, Manuel Razo-Mejia, Stephanie L. Barnes, Nathan M.
Belliveau, Tal Einav, and I joined forces and set out to build a complete
theoretical model for allosteric transcriptional regulation coupled with a
thorough experimental dissection. This was no small task and would have likely
taken a full Ph.D.'s worth of effort for a single person to do. Yet, within a
year of project inception we had submitted a manuscript to preprint servers
where all of us were annotated as equal contributors. This experience defined
how I view collaboration in scientific research and serves
as a shining example of scientific socialism. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Together, we enumerated a complete thermodynamic model for the inducible
simple repression motif and defined a succinct input-output function for the
fold-change in gene expression (schematized in @Fig:induction_intro (B)).
This model, which is explored in depth in Chapter 2, is defined by a
minimal set of biophysical parameters, many of which can be directly measured
using standard tricks of molecular biology and biochemistry. With a model in
hand, we turned to a collection of 17 unique *E. coli* strains, each with
different copy numbers of the repressor protein and different regulatory DNA
sequences. Using our theoretical model, we inferred the lone two biophysical
parameters which we did not know *a priori* from a single experimental strain
(white points in middle panel of @Fig:induction_intro (C)), and tested our
predictions on all other experimental strains. We found the model to be
remarkably predictive, suggesting that our "toy" model of an allosteric
repressor captured the underlying physics of the system.

![**Experimental dissection of the inducible simple repression input-output
function.** (A) Schematic diagram of the inducible simple repression motif. (B)
Schematic diagram of the input-output function as is derived in Chapter 2. (C)
Experimental measurements of the fold-change in gene expression using the
*lac* repressor from *E coli*. Different rows correspond to different operator
sequences and therefore different values for $\Delta\varepsilon_{RA}$. Different
colors correspond to different values for the average repressor copy number $R$.
While filled points in the middle panel represent the experimental strain used to
infer the values of the inducer dissociation constants. All points correspond to
the mean of at least 10 biological replicates.](ch1_fig4){#fig:induction_intro
short-caption="Experimental dissection of the inducible simple repression
input-output function."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A key feature of this work is a derivation of
thermodynamic state variable of this regulatory architecture which we term
the *free energy*. This parameter provides an intuition for the effective
free energy difference between states of the promoter in which the repressor
is bound relative those states in which the repressor is not bound to the
promoter. This parameter accounts for all of the ways in which one can tune
the parameter values and still achieve the same fold-change in gene expression, as is diagrammed
in @Fig:collapse_intro (A). While we leave the details of this derivation to
Chapter 2, we emphasize
that this formalism provides a means by which all of the experimental
measurements plotted in @Fig:induction_intro (C) can be collapsed onto a
master curve defined *only* by the free energy, which is illustrated in
@Fig:collapse_intro (B). This scaling, often referred to as "data collapse"
in physics, concretely shows that one has identified the *natural variable*
of the system. With this scaling function in hand, we are able to make a
measurement of the biophysical parameters, compute the free energy, and make
a concrete prediction of what the fold-change in gene expression will be. Or,
as we will see in the following section, details of the biophysical
parameters can be determined directly from an empirical measurement of the
free energy.

![**Collapse of individual induction profiles onto a simple scaling function.**
(A) Many different combinations of parameter values yield the same value of the
fold-change in gene expression, shown as red and blue horizontal planes. Any
point on those planes corresponds to a single value of the free energy (middle)
and will appear on the master curve. (B) Data presented in @Fig:induction_intro
(C) collapsed onto the master curve defined by the predicted value of the free
energy.](ch1_fig5){#fig:collapse_intro short-caption="Collapse of individual
induction profiles onto a simple scaling function."}

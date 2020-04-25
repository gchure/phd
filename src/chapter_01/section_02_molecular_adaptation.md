## Topic I: The Janus Face of Molecules

Monod famously referred to allostery as "the second secret of life"
[@monod1965]. This secret has been now been declassified as allosteric
regulation can be found across all domains of life across myriad types of
bioological processes. Allostery can be found governing the behavior of ion
channels [@einav2017, @auerbach2012a], enzymatic reactions [@einav2016],
chemotaxis [@keymer2006], G-protein coupled receptors [@canals2012], quorum
sensing [@swem2008], and transcriptional regulation [@huang2018,
@lindsley2006a], to name a few of many examples. Despite the objective
complexity of the structures of all of these allosteric molecules, they can
be frequently be reduced to simple cartoons where the details of
conformational changes, substrate binding affinities, and more can be
massaged into a small set of key details. @Fig:allostery (A) shows the
molecular structures of a variety of allosteric transcriptional regulators
(top). While each has their own fascinating structure and continuum of
conformational states, all can be coarse grained into a simple cartoon of an
active (red) or inactive (purple) molecule, both of which states possess
binding pockets (semicircular notches) for an inducer molecule (orange).

Much as we can reduce the complexity of allosteric molecules schematically, we
can enumerate simple mathematical models that describe their behavior. Thermodynamic models built on an assumption of quasi-equilibrium are routinely used to
describe complex biological phenomena despite the reality that being in
thermodynamic equilibrium is synonymous with being dead. Despite this glaring
assumption, such models have been shown to be exceptionally predictive for a
variety of complex systems, especially in modeling molecular binding reactions 
[@dill2010] and allostery writ large [@swem2008; @keymer2006,  @einav2017,
@einav2016]. As the timescales of binding and unbinding reactions are orders of
magnitude smaller than that of many other processes in the cel, it is fair to
make the approximation that molecular binding is in equilibrium. Under this
assumption, we are granted the powerful mathematical privilege to say that the
probability of any state of the system $P_\text{state}$ follows a Boltzmann distribution,
$$
P_\text{state} \propto e^{-\frac{\epsilon_\text{state}}{k_BT}},
$${#eq:boltzmann}
where $\epsilon_\text{state}$ is the energy of that state, $k_B$ is the Boltzmann constant,
and $T$ is the system temperature. Therefore, if we are interested in computing
the probability of a given allosteric protein being in the active state, we
merely  have to enumerate all of the Boltzmann weights (given by @Eq:boltzmann)
and compute
$$
P_\text{active}  = \frac{\text{sum over all possible active states}}{\text{sum
over all possible states}}
$${#eq:verbal_pact}.
This probability, defined as a function of the inducer concentration is shown
schematically in @Fig:allostery(B). While we have passed over some of the more
subtle details of this calculation, the plot in @Fig:allostery (B) presents a
*quantitative* prediction of how the activity of an allosteric molecule should
scale as a function of the inducer (in this case, becoming less active as more
inducer is present).

![**A Coarse grained representation of an allosteric molecule.** (A) Crystal
structures of a variety of allosteric transcription factors are shown at the
top. In this thesis, we coarse grain away many of the details to a minimal model
(bottom) where the protein can be represent as being either active (red),
inactive (purple), both of which states can bind an inducer molecule (orange).
(B) By making an assumption of quasi-equilibrium, we can trivially compute a
mathematical description of the active probability of an allosteric protein as a
function of the inducer concentration. In this particular case, the inactive
state becomes stabilized relative to the active state at higher concentrations
of inducer molecule](ch1_fig3){#fig:allostery short-caption="A Coarse grained
representation of an allosteric molecule."}

In **Chapter 2** of this dissertation, we use the Monod-Wyman-Changeux model of
allosery [@monod1965] to build a predictive model of gene expression changes in
response to changing activity of an allosteric transcriptional repressor.
Using the same tricks given by @Eq:boltzmann and @Eq:verbal_pact, we extend upon
a previously characterized thermodynamic model  of the simple repression motif.
This motif, schematized in @Fig:induction_intro (A) is not a convenient
abstraction of a regulatory architecture. Rather, this motif is the most
ubiquitous regulatory schemes in *E. coli* [@gama-castro2016, @ireland2020].
This architecture has been the target of much theoretical and experimental
dissection [@vilar2003; @buchler2003; @bintu2005; @garcia2011; @brewster2014; 
@phillips2019], though previous inclusion of allostery into a thermodynamic model of
of this motif had left much to be desired. 


At the beginning of 2016, Manuel Razo-Mejia, Stephanie L. Barnes, Nathan M.
Belliveau, Tal Einav, and I joined forces and set out to build a complete
theoretical model for allosteric transcriptional regulation coupled with a
thorough experimental dissection. This was no small task and would have likely
taken a full Ph.D.'s worth of effort for a single person to do. Yet, within a
year of project inception we had submitted a version to the preprint servers
where all of us were annotated as equal contributors. This experience defined
how I view collaboration in scientific research and serves, at least in my mind,
as a shining example of scientific socialism. 

Together, we enumerated a complete thermodynamic model for the inducible
simple repression motif and defined a succinct input-output function for the
fold-change in gene expression [schematized in @Fig:induction_intro (B)].
This model, which is explored in depth in **Chapter 2** is defined by a
minimal set of biophysical parameters, many of which can be directly measured
using standard tricks of molecular biology and biochemistry. With a model in
hand, we turned to a collection of 17 unique *E. coli* strains, each with
different copy numbers of the repressor protein and different regulatory DNA
sequences. Using our theoretical model, we inferred the lone two biophysical
parameters which we did not know *a priori* from a single experimental strain
[white points in middle panel of @Fig:induction_intro (C)], and tested our
predictions on all other experimental strains.

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

A key feature of this work is a derivation of thermodynamic state variable of
this regulatory architecture which we term the "free energy". This parameter
provides an intuition for the effective free energy difference between states of
the promoter in which thee repressor is bound relative those states in which the
repressor is not bound to the  promoter. This parameter accounts for all of the
ways in which one can tune the parameter values and still achieve the same
fold-change, as is diagrammed in @Fig:collapse_intro (A). While we leave the details of this derivation
to **Chapter 2** and the corresponding supplement in **Chapter 6**, we emphasize
that this formalism provides a means by which all of the experimental
measurements plotted in @Fig:induction_intro (C) can be collapsed onto a master
curve defined *only* by the free energy, which is illustrated in
@Fig:collapse_intro (B). This scaling, often referred to as "data collapse" in
physics, concretely shows that one has identified the "natural variable" of the
system. With this scaling function in hand, we are able to make a measurement of
the biophysical parameters, compute the free energy, and make a concrete
prediction of what the fold-change in gene expression will be. Or, as we will
see in the following section, details of the biophysical parameters can be
determined directly from an empirical measurement of the free energy. 

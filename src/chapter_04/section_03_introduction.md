## Introduction

Cellular physiology is inextricably tied to the extracellular environment.
Fluctuations in nutrient availability and variations in temperature, for
example, can drastically modulate the cell’s growth rate, which is often used
as a measure of the evolutionary fitness [@schaechter1958]. In response to such
environmental insults, cells have evolved myriad clever mechanisms by which
they can adapt to their changing surroundings, many of which involve
restructuring their proteome such that critical processes (i.e. protein
translation) are allocated the necessary resources. Recent work exploring this
level of adaptation using mass spectrometry, ribosomal profiling, and RNA
sequencing have revealed that various classes of genes (termed "sectors") are
tuned such that the protein mass fraction of the translational machinery is
prioritized over the metabolic and catabolic machinery in nutrient replete
environments [@scott2014; @klumpp2014; @hui2015; @schmidt2016; @li2014a]. This
drastic reorganization is mediated by the regulation of gene expression,
relying on the concerted action of myriad transcription factors. Notably, each
gene in isolation is regulated by only one or a few components
[@gama-castro2016]. The most common regulatory architecture in *Escherichia
coli* is the simple repression motif in which a transcriptional repressor binds
to a single site in the promoter region, occluding binding of an RNA polymerase
[@rydenfelt2014; @phillips2019]. The simple activation architecture, in which
the simultaneous binding of an activator and an RNA polymerase amplifies gene
expression, is another common mode of regulation. Combinatorial regulation such
as dual repression, dual activation, or combined activation and repression can
also be found throughout the genome, albeit with lower
frequency [@phillips2019]. The ubiquity of the simple repression and simple
activation motifs illustrate that, for many genes, the complex systems-level
response to a physiological perturbation boils down the binding and unbinding
of a single regulator to its cognate binding sites.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Despite our knowledge of these modes of
regulation, there remains a large disconnect between concrete, physical
models of their behavior and experimental validation. The simple repression
motif is perhaps the most thoroughly explored theoretically and
experimentally [@phillips2019] where equilibrium thermodynamic [@garcia2011;
@garcia2012; @brewster2014; @razo-mejia2018; @barnes2019] and kinetic
[@jones2014; @ko1991; @kepler2001; @michel2010] models have been shown to
accurately predict the level of gene expression in a variety of contexts.
While these experiments involved variations of repressor copy number,
operator sequence, concentration of an external inducer, and amino acid
substitutions, none have explored how the physiological state of the cell as
governed by external factors influences gene expression. This is arguably one
of the most critical variables one can experimentally tune to understand the
roles these regulatory architectures play in cellular physiology writ
large.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In this work, we interrogate the adaptability
of a simple genetic circuit to various physiological stressors, namely carbon
source quality and growth temperature. Following the aforementioned
thermodynamic models, we build upon this theory-experiment dialogue by using
environmental conditions as an experimentally tunable variable and determine
their influence on various biophysical parameters. Specifically, we use
physiological stressors to tune the growth rate. One mechanism by which we
modulate the growth rate is by exchanging glucose in the growth medium for
the poorer carbon sources glycerol and acetate, which decrease the growth
rate by a factor of $\approx$ 1.5 and $\approx$ 4 compared to glucose,
respectively. We hypothesize that different carbon sources should, if
anything, only modulate the repressor copy number seeing as the relationship
between growth rate and total protein content has been rigorously quantified
[@schaechter1958; @schmidt2016; @li2014a; @jun2018]. Using single-cell
time-lapse fluorescence microscopy, we directly measure the copy number of
the repressor in each condition. Under a simple hypothesis, all other
parameters should be unperturbed, and we can thus rely on previously
determined values to make parameter-free predictions of the fold-change in
gene expression.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Despite the decrease in growth rate, both the
fold-change in gene expression and the repressor copy number remains largely
unaffected. We confirm this is the case by examining how the effective free
energy of the system changes between carbon sources, a method we have used
previously to elucidate parametric changes due to mutations within a
transcription factor [@chure2019] and has been extensively discussed in
Chapter 3. This illustrates that the energetic parameters defining the
fraction of active repressors and their affinity for the DNA are ignorant of
the carbon-dependent physiological states of the cell. Thus, in this context,
the values of the biophysical parameters determined in one condition can be
used to draw predictions in others.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We then examine how variations in temperature
influence the transcriptional output. Unlike in the case of carbon source
variation, temperature dependence is explicit in our model: the repressor-DNA
binding energy and the energetic difference between the active and inactive
states of the repressor are scaled to the thermal energy of the system at
37$^\circ$ C. This is defined via the Boltzmann distribution which states
that the probability of a state $p_{state}$ is related to the energy of that
state $\varepsilon_{state}$ as
$$
p_{state} \propto e^{-\varepsilon_{state} / k_BT},
$${#eq:boltz}
where $k_B$ is the Boltzmann constant and $T$ is the temperature of the system.
Given knowledge of $T$ for a particular experiment, we can easily draw
predictions of the fold-change in gene expression. However, we find that the
fold-change in gene expression is inconsistent with this simple model,
revealing an incomplete description of the energetics. We then examine how
entropic effects neglected in the initial estimation of the energetic
parameters may play an important role; a hypothesis that is supported when we
examine the change in the effective free energy.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The results presented here are, to our
knowledge, the first attempts to systematically characterize the
growth-dependent effects on biophysical parameters in thermodynamic models of
transcription. While some parameters of our model are affected by changing
the growth rate, they change in ways that are expected or fall close within
our *a priori* predictions, suggesting that such modeling can still be
powerful in understanding how adaptive processes influence physiology at the
level of molecular interactions.

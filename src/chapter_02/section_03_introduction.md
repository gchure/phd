## Introduction

Understanding how organisms sense and respond to changes in their
environment has long been a central theme of biological inquiry. At the
cellular level, this interaction is mediated by a diverse collection of
molecular signaling pathways. A pervasive mechanism of signaling in these
pathways is allosteric regulation, in which the binding of a ligand
induces a conformational change in some target molecule, triggering
a signaling cascade [@lindsley2006a]. One of the most important examples
of such signaling is offered by transcriptional regulation, where
a transcription factors' propensity to bind to DNA will be altered upon
binding to an allosteric effector.

Despite the ubiquity of allostery, we largely lack a formal, rigorous, and
generalizable framework for studying its effects across the broad variety of
contexts in which it appears. A key example of this is transcriptional
regulation, in which allosteric transcription factors can be induced or
corepressed by binding to a ligand. An allosteric transcription factor can
adopt multiple conformational states, each of which has its own affinity for
the ligand and for its DNA target site. *In vitro* studies have rigorously
quantified the equilibria of different conformational states for allosteric
transcription factors and measured the affinities of these states to the
ligand [@harman2001; @lanfranco2017]. In spite of these experimental
observations, the lack of a coherent quantitative model for allosteric
transcriptional regulation has made it impossible to predict the behavior of
even a simple genetic circuit across a range of regulatory parameters,
physiological states of the organism, and evolutionary isoforms of the
regulatory sequences.

The ability to predict circuit behavior robustly— that is, across both
broad ranges of parameters and regulatory architectures —is important for
multiple reasons. First, in the context of a specific gene, accurate
prediction demonstrates that all components relevant to the genes'
behavior have been identified and characterized to sufficient
quantitative precision. Second, in the context of genetic circuits in
general, robust prediction validates the model that generated the
prediction. Possessing a validated model also has implications for future
work. For example, when we have sufficient confidence in the model,
a single data set can be used to accurately extrapolate a system’s
behavior in other conditions. Moreover, there is an essential distinction
between a predictive model, which is used to predict a system’s behavior
given a set of input variables, and a retroactive model, which is used to
describe the behavior of data that has already been obtained. We note
that even some of the most careful and rigorous analysis of
transcriptional regulation often entails only a retroactive reflection on
a single experiment. This raises the fear that each regulatory
architecture may require a unique analysis that cannot carry over to
other systems, a worry that is exacerbated by the prevalent use of
phenomenological functions (e.g. Hill functions or ratios of polynomials)
that can analyze a single data set but cannot be used to extrapolate
a system’s behavior in other conditions [@setty2003; @poelwijk2011a;
@vilar2013; @rogers2015; @rohlhill2017].

This work explores what happens when theory takes center stage, namely, we
first write down the equations governing a system and describe its expected
behavior across a wide array of experimental conditions, and only then do we
set out to experimentally confirm these results. Building upon previous work
[@garcia2011; @brewster2014; @weinert2014] and the work of Monod, Wyman, and
Changeux [@monod1965], we present a statistical mechanical rendering of
allostery in the context of induction and corepression (shown schematically
in @Fig:inducible_types and henceforth referred to as the MWC model) and use
it as the basis of parameter-free predictions which we then test
experimentally. More specifically, we study the simple repression motif – a
widespread bacterial genetic regulatory architecture in which binding of a
transcription factor occludes binding of an RNA polymerase, thereby
inhibiting transcription initiation. The MWC model stipulates that an
allosteric protein fluctuates between two distinct conformations – an active
and inactive state – in thermodynamic equilibrium [@monod1965]. During
induction, for example, effector binding increases the probability that a
repressor will be in the inactive state, weakening its ability to bind to the
promoter and resulting in increased expression. To test the predictions of
our model across a wide range of operator binding strengths and repressor
copy numbers, we design a genetic construct in *Escherichia coli* in which
the binding probability of a repressor regulates gene expression of a
fluorescent reporter.

In total, the work presented here demonstrates that one extremely compact set
of parameters can be applied self-consistently and predictively to different
regulatory situations including simple repression on the chromosome, cases in
which decoy binding sites for repressor are put on plasmids, cases in which
multiple genes compete for the same regulatory machinery, cases involving
multiple binding sites for repressor leading to DNA looping, and induction by
signaling [@garcia2011; @garcia2011b; @brewster2012; @brewster2014;
@boedicker2013a; @boedicker2013]. Thus, rather than viewing the behavior of
each circuit as giving rise to its own unique input-output response, the MWC
model provides a means to characterize these seemingly diverse behaviors
using a single unified framework governed by a small set of parameters.

![ **Transcriptional regulatory motifs involving an allosteric repressor.**
(A) We consider a promoter regulated solely by an allosteric repressor in
which the active (repressive, red blobs) state of the repressor is
energetically favorable in the absence (induction, left panel) or presence
(corepression, right panel) of an allosteric effector. Both inducible
repression and corepression are ubiquitous regulatory strategies in *E.
coli*, several examples of which are given in the tables below each panel.
(B) A representative regulatory response (fold-change in gene expression) of
the two architectures shown in Panel (A) as a function of the corresponding
allosteric effector concentration. Properties of interest to this work are
shown schematically upon the regulatory response. (C) Historical progression
of thermodynamic modeling of the inducible simple-repression regulatory
architecture. @garcia2011 used colorimetric assays and quantitative Western
blots to investigate how single-site repression is modified by the repressor
copy number and repressor-DNA binding energy. @brewster2014 used video
microscopy to probe how the copy number of the promoter and presence of
competing repressor binding sites affect gene expression. Building upon these
works, we use flow cytometry to determine the inducer-repressor dissociation
constants and demonstrate that with these parameters we can predict *a
priori* the behavior of the system for any repressor copy number, DNA binding
energy, gene copy number, and inducer
concentration.](ch2_fig1){#fig:inducible_types short-caption="Transcriptional
regulatory architectures involving an allosteric repressor."}

## Introduction

Understanding how organisms sense and respond to changes in their environment
has long been a central theme of biological inquiry. At the cellular level,
this interaction is mediated by a diverse collection of molecular signaling
pathways. A pervasive mechanism of signaling in these pathways is allosteric
regulation, in which the binding of a ligand induces a conformational change
in some target molecule, triggering a signaling cascade [@lindsley2006a]. One
of the most important examples of such signaling is offered by transcriptional
regulation, where a transcription factor’s propensity to bind to DNA will be
altered upon binding to an allosteric effector.

Despite allostery’s ubiquity, we lack a formal, rigorous, and generalizable
framework for studying its effects across the broad variety of contexts in
which it appears. A key example of this is transcriptional regulation, in
which allosteric transcription factors can be induced or corepressed by
binding to a ligand. An allosteric transcription factor can adopt multiple
conformational states, each of which has its own affinity for the ligand and
for its DNA target site. *In vitro* studies have rigorously quantified the
equilibria of different conformational states for allosteric transcription
factors and measured the affinities of these states to the ligand
[@harman2001; @lanfranco2017]. In spite of these experimental observations,
the lack of a coherent quantitative model for allosteric transcriptional
regulation has made it impossible to predict the behavior of even a simple
genetic circuit across a range of regulatory parameters.

The ability to predict circuit behavior robustly— that is, across both
broad ranges of parameters and regulatory architectures —is important
for multiple reasons. First, in the context of a specific gene, accurate
prediction demonstrates that all components relevant to the gene’s
behavior have been identified and characterized to sufficient
quantitative precision. Second, in the context of genetic circuits in
general, robust prediction validates the model that generated the
prediction. Possessing a validated model also has implications for
future work. For example, when we have sufficient confidence in the
model, a single data set can be used to accurately extrapolate a
system’s behavior in other conditions. Moreover, there is an essential
distinction between a predictive model, which is used to predict a
system’s behavior given a set of input variables, and a retroactive
model, which is used to describe the behavior of data that has already
been obtained. We note that even some of the most careful and rigorous
analysis of transcriptional regulation often entails only a retroactive
reflection on a single experiment. This raises the fear that each
regulatory architecture may require a unique analysis that cannot carry
over to other systems, a worry that is exacerbated by the prevalent use
of phenomenological functions (e.g. Hill functions or ratios of
polynomials) that can analyze a single data set but cannot be used to
extrapolate a system’s behavior in other conditions [@setty2003; @poelwijk2011a; @vilar2013; @rogers2015; @rohlhill2017].

This work explores what happens when theory takes center stage, namely, we first write down the equations governing a system and describe its expected behavior across a wide array of experimental conditions, and only then do we set out to experimentally confirm these results. Building upon previous work [@garcia2011; @brewster2014; @weinert2014] and the work of Monod, Wyman, and Changeux [@monod1965], we present a statistical mechanical rendering of allostery in the context of induction and corepression (shown schematically in 
@Fig:inducible_types and henceforth referred to as the MWC model) and use it as
the basis of parameter-free predictions which we then test experimentally. More
specifically, we study the simple repression motif – a widespread bacterial
genetic regulatory architecture in which binding of a transcription factor
occludes binding of an RNA polymerase, thereby inhibiting transcription
initiation. The MWC model stipulates that an allosteric protein fluctuates
between two distinct conformations – an active and inactive state – in
thermodynamic equilibrium [@monod1965]. During induction, for example, effector
binding increases the probability that a repressor will be in the inactive
state, weakening its ability to bind to the promoter and resulting in increased
expression. To test the predictions of our model across a wide range of operator binding strengths and repressor copy numbers, we design an *E. coli* genetic
construct in which the binding probability of a repressor regulates gene
expression of a fluorescent reporter.

In total, the work presented here demonstrates that one extremely
compact set of parameters can be applied self-consistently and
predictively to different regulatory situations including simple
repression on the chromosome, cases in which decoy binding sites for
repressor are put on plasmids, cases in which multiple genes compete for
the same regulatory machinery, cases involving multiple binding sites
for repressor leading to DNA looping, and induction by signaling [@garcia2011;
@garcia2011b; @brewster2012; @brewster2014; @boedicker2013a; @boedicker2013] . Thus,
rather than viewing the behavior of each circuit as giving rise to its
own unique input-output response, the MWC model provides a means to
characterize these seemingly diverse behaviors using a single unified
framework governed by a small set of parameters.

![ We consider a promoter regulated solely by an allosteric repressor.
When bound, the repressor prevents RNAP from binding and initiating
transcription. Induction is characterized by the addition of an effector
which binds to the repressor and stabilizes the inactive state (defined
as the state which has a low affinity for DNA), thereby increasing gene
expression. In corepression, the effector stabilizes the repressor’s
active state and thus further reduces gene expression. We list several
characterized examples of induction and corepression that support
different physiological roles in <span>*E. coli*</span> . A schematic
regulatory response of the two architectures shown in Panel plotting the
fold-change in gene expression as a function of effector concentration,
where fold-change is defined as the ratio of gene expression in the
presence versus the absence of repressor. We consider the following key
phenotypic properties that describe each response curve: the minimum
response (leakiness), the maximum response (saturation), the difference
between the maximum and minimum response (dynamic range), the
concentration of ligand which generates a fold-change halfway between
the minimal and maximal response ($[EC_{50}]$), and the log-log slope
at the midpoint of the response (effective Hill coefficient). (C) Over
time we have refined our understanding of simple repression
architectures. A first round of experiments used colorimetric assays and
quantitative Western blots to investigate how single-site repression is
modified by the repressor copy number and repressor-DNA binding energy .
A second round of experiments used video microscopy to probe how the
copy number of the promoter and presence of competing repressor binding
sites affect gene expression, and we use this data set to determine the
free energy difference between the repressor’s inactive and active
conformations . Here we used flow cytometry to determine the
inducer-repressor dissociation constants and demonstrate that with these
parameters we can predict *a priori* the behavior of the system for any
repressor copy number, DNA binding energy, gene copy number, and inducer
concentration.](ch2_fig1){#fig:inducible_types short-caption="Transcriptional
regulatory architectures involving an allosteric repressor."}

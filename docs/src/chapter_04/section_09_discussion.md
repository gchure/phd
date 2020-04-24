## Discussion


The past century of work in bacterial physiology has revealed a rich
molecular complexity that drives cellular growth rate [@jun2018]. A key
finding of this body of work is that the composition of the proteome is
highly dependent on the specific growth condition, with entire classes
of genes being up- or down-regulated to ensure that enough resources are
allocated towards maintaining a pool of active ribosomes
[@hui2015; @scott2014]. These studies have led to a coarse-grained view
of global gene expression where physiological perturbations
substantially change the molecular composition of the cell, obfuscating
the utility of using thermodynamic models of individual regulatory
elements across physiological states. In this work, we rigorously
examine how robust the values of the various biophysical parameters are
to changes in cellular physiology.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We first examined how nutrient fluctuations
dictate the output of this architecture. We took three carbon sources with
distinct metabolic pathways and varying quality and measured the level of
gene expression, hypothesizing that the values of the biophysical parameters
to be independent of the growth medium. We found that even when the growth
rate is varied across a wide range (220 minute doubling time in acetate to 60
minute doubling time in glucose supported medium), there is no significant
change to the fold-change in gene expression or in the expression of the
transcription factor itself, within the resolution of our experiments. Given
numerous quantitative studies of the proteomic composition reveal a
dependence on protein content with growth rate [@hui2015; @schmidt2016;
@li2014], we find this robustness to be striking.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;@schmidt2016 found that the native expression of
LacI has a weak positive correlation with the growth rate. The native
LacI promoter region is solely regulated by activation via the cAMP
Receptor Protein (CRP), a broadly acting dual regulator in *E. coli*
[@gama-castro2016]. This is in contrast to the LacI expression system
used in the present work where the promoter is negatively regulated by
the TetR repressor, itself expressed from a low-copy number plasmid.
Furthermore, the expression of LacI in this work is tuned by the
addition of the allosteric effector of TetR, ATC, adding yet another
layer of allosteric regulation on LacI expression. The significant
difference in the regulatory mechanisms between the native and synthetic
circuit used in this work makes the two findings difficult to directly
compare. Regardless, our finding that the fold-change in gene expression
is unaltered from one carbon source to another illustrates that the
values of the biophysical parameters $\Delta\varepsilon_ {R}$ and
$\Delta\varepsilon_{AI}$ remain unperturbed, permitting quantitative
prediction of gene expression across numerous physiological states.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;However, in varying the temperature, we find that the predictive utility
of the biophysical parameters values determined at 37$^\circ$ C is
diminished, indicating that there are hidden effects not explicitly
accounted for in our thermodynamic model. The measurements of the
fold-change in gene expression are under- or over-estimated when the
temperature is increased or decreased, respectively, when one simply
rescales the energetic terms by the relative change in temperature.
There are many features of transcriptional regulation that are not
explicitly considered in our coarse-graining of the architecture into a
two-state model. Recently, it has been suggested that the phenomenon of
allostery writ large should be framed in the language of an ensemble of
states rather than a simple active/inactive distinction [@motlagh2014].
While our recent work illustrates that a two-state rendering of an
allosteric repressor is highly predictive in a variety of situations
[@razo-mejia2018; @chure2019], we must now consider details which are
dependent on the temperature of the system. In @Fig:deltaF_temp, we demonstrate that incorporating a
temperature-dependent entropic cost to the energetic terms significantly
improves the description of the experimental data. This is not to say,
however, that this is now an open-and-closed case for what precisely
defines this entropic cost. Rather, we conclude that the phenomenology
of the temperature dependence can be better described by the inclusion
of a correction factor that is linearly dependent on the system
temperature. Biology is replete with phenomena which can introduce such
an effect, including changes to the material properties of the repressor
and DNA [@goethe2015; @mondal2011], excluded volume effects
[@driessen2014], and solubilities
[@elf2007; @kao-huang1977; @yakovchuk2006]. Understanding the
mechanistic underpinnings of temperature dependence in elasticity theory
was borne out of similar phenomenological characterization
[@friedel1974] and required a significant dialogue between theory and
experiment [@phillips2001]. Further work is now needed to develop a
theory of temperature effects in the regulation of gene expression.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The effective free energy $F$, as defined in @Eq:bohr_parameter, is a state variable of the simple
repression regulatory architecture. This is illustrated in @Fig:collapse
 where fold-change measurements from a wide array of conditions (and measurement techniques) can be collapsed onto
the same theoretical description. Evolutionary perturbations (such as
mutations in the operator or repressor sequence), physiological changes
(such as modulations of the growth rate), or changes in the level of
activity of the repressor (due to changes in inducer concentration) do
not change the fundamental physics of the system and can all be
described by changes in the free energy relative to one another. While
such a statement is not "surprising", we can now say it with
quantitative confidence and use this principle to probe the degree to
which physiological perturbations influence the biophysical parameters
writ large. With such a framework in hand, we are in the auspicious
position to take a predictive approach towards understanding how this
regulatory architecture evolves in experimental settings, shedding light
on the interplay between biophysical parameters, organismal fitness, and
the fundamental forces of evolution.

![**A singular theoretical description for the molecular biophysics of
physiological and evolutionary adaptation in the simple repression regulatory
architecture.** Measurements of the fold-change in gene expression varying the
sequence of the operator site [orange pentagon, @garcia2011], concentration of
the extracellular inducer [green squares, @razo-mejia2018 and Chapter 2],
amino-acid sequence of the repressor [blue points, @chure2019 and Chapter 3],
and the various growth conditions queried in this work can be collapsed as a
function of the effective free energy. Error bars correspond to the standard
error of five to 10 biological replicates.](ch4_fig6){#fig:collapse
short-caption="A singular theoretical description for the molecular biophysics
of physiological and evolutionary adaptation in the simple repression motif."}

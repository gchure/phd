## Topic IV: The Physiological Adaptability of Transient Molecular Interactions

In **Chapers 4 and 5** (and the associated supplementary **Chapters 8 and 9**),
we explore the final level of adaptation in @Fig:adaptation_levels --
physiological adaptation. We do so in two distinctly different systems. The
first (Chapter 4) builds upon our discussion of transcriptional regulation, but
now examines how robust the biophysical parameters of the thermodynamic model
are to changes in physiology, either by changing the available carbon source or
by changing the temperature. Secondly, we examine
physiological adaptation in the context of osmoregulation -- a true matter of
life and death in the single-celled world.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Up to this point in our travels through scientific history, we have examined
Monod's growth curves  in various pairwise combinations of sugars. A feature of
note is that the presence of diauxic shifts can be seen in various organisms and
for many different types of sugars such as sucrose/arabinose, glucose/sorbitol, 
and glucose/lactose pairings [@monod1947]. These combinations reveal that cells
are able to juggle dual-input logic systems where the "decision" to digest one
carbon source or another relies on monitoring changes in concentrations of
either sugar. In his 1947 treatise, Monod showed that this phenomenon was not
limited to dual-carbon mixtures and presented a "triauxic" growth curve of *E.
coli* grown on a glucose/sorbitol/glycerol mixture, shown in
@Fig:triauxic_growth. This result illustrated to Monod that the mechanisms
underlying enzymatic adaptation "have the character of *competitive
interactions* between different specific enzyme forming systems" [@monod1947]. 

![**A metabolic hierarchy in a three-component growth mixture of glucose,
sorbitol, and glycerol.** A "triauxic" growth curve illustrating a hierarchy of
carbon source metabolism. An *E. coli* culture was grown in an medium with equal
parts glucose, sorbitol and glycerol with utilization in that order. Auxic
transitions are shown as black points and white vertical lines. Regions of the
growth curve where glucose, sorbitol, and glycrol are primarily consumed are
colored in blue, orange, and green, respectively. Data digitized from
@monod1947.](ch1_fig9){#fig:triauxic_growth short-caption="A metabolic hierarchy
in a growth medium containing glucose, sorbitol, and glycerol."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;These competing interactions must be resilient
to a variety of physiological states. Despite the fact that the carbon atoms
in glucose, sorbitol, and glycerol are all ultimately incorporated into the
same biomolecules, their pathways to utilization are all distinct and include
a variety of different metabolic intermediates. Furthermore, the exponential
growth phases in @Fig:triauxic_growth for each carbon source have different
growth rates which itself results in large changes in cell volume
[@taheri-araghi2015; @jun2018], genome copy number [@nordstrom2006], and
global gene expression patterns [@li2014a; @schmidt2016; @hui2015]. Despite
these changes in cellular physiology, the regulatory systems underlying
enzymatic adaptation still function with binding of transcription factors
being ignorant of the majority of possible metabolic states of the cell.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Despite this empirical observation, it has been
commonly assumed that the utility of thermodynamic models of gene expression are
limited and that the precise values of the biophysical parameters are
directly tied to the physiological state in which they were determined. It has even been said that thermodynamic models of gene
expression have been a "tactical success, yet strategic failure" in building
an understanding of how genomes operate [@phillips2019].

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In **Chapter 4** of this dissertation, we
quantitatively assess these assumptions in the context of gene expression by
considering the theoretical models built in Chapters 2 and 3 and directly
measuring the adaptability of the inducible simple repression regulatory
architecture across different physiological states. Namely, we explore how
predictive our thermodynamic model can be when modulating either the quality of the
carbon source (glucose, glycerol, or acetate, @Fig:growth_intro (A)) or by
changing the temperature of the growth medium (32$^\circ$ C, 37$^\circ$ C, or
42$^\circ$ C, @Fig:growth_intro (B)). The culture doubling time varies by
nearly a factor of four across the different conditions, illustrating the
diversity in physiological states.

![**Control of cellular physiology via carbon source and temperature
variation.** (A) Carbon sources used in work presented in Chapter 4 (left) with
"star rating" indicating quality of the carbon source. Growth curves for the
three carbon sources, all at 37$^\circ$ C are shown on right-hand panel. (B)
Growth temperatures explored in Chapter 4 (left) with "star rating" indicating
fastest growth rate. Growth curves (right) are shown for the three temperatures,
all of which use glucose as the sole carbon source. For right-hand panels in (A)
and (B), optical density is computed relative to the initial optical density of
the culture.](ch1_fig10){#fig:growth_intro short-caption="Methods of
physiological control used in Chapter 4."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; How could this variation in cellular physiology
be incorporated into our thermodynamic model? Up to this point in this thesis, a
experiments have been conducted in a growth medium supplemented with glucose
held at a balmy 37$^\circ$ C. However, nowhere in the thermodynamic model
schematized in @Fig:induction_intro (B) is it specified *which* carbon source
must be present whereas the temperature of the system is explicitly included as
a multiplicative factor $\beta = \left(k_BT\right)^{-1}$ in front of the
exponentiated terms. These features of the model allow us to make explicit
predictions of how these perturbations should influence the observed fold-change
in gene expression, if at all.   


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The parameter that we can say *a priori* is
very likely to change is the repressor copy number $R$. In Chapters 2 and 3,
we knew the total repressor copy number from previous work where the copy
numbers were directly measured via quantitative Western blotting in a
particular physiological state [@garcia2011]. However, it has been known for
nearly three-quarters of a century that the total protein content of the cell
scales linearly with the growth rate [@schaechter1958; @jun2018], a
phenomenon that has recently been queried at the single-protein level through
proteomic methods [@schmidt2016; @li2014a;@valgepea2013a;@peebo2015]. Thus,
we cannot assume that the protein copy number of the strains used in Chapters
2 and 3 will not be perturbed. To account for this fact, we used a
fluorescence-based method to directly count the number of LacI repressors per
cell in each growth condition, a method which is discussed in extensive
detail in Chapter 9. This experimental approach, while necessary, is
extremely laborious. I am indebted to the work of Zofii A. Kaczmarek for her
heroic efforts in conducting a large number of the experiments presented in
Chapter 4.


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Our work revealed two key features of this
thermodynamic model of gene expression. First, we found that the values of
the biophysical parameters inferred from a single physiological sate were
remarkably predictive when the quality of the carbon source was decreased
(@Fig:growth_results_intro (A, left)). This indicates that this genetic circuit is
largely insulated from the metabolic state of the cell. This is exemplified
in our ability to collapse the measurements across different carbon sources
as a function of the free energy, shown in @Fig:growth_results_intro (A, right). For
the carbon sources studied in this Chapter, we conclude that this simple
thermodynamic model can be considered both tactically and strategically
successful.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Yet when it comes to temperature, we find
that a simple rescaling of the thermal energy of the system is *not*
sufficient to predict the output of this genetic circuit when the temperature
is varied (dashed lines in left-hand side of @Fig:growth_results_intro, (B)). This is
not necessarily a surprising result as binding of transcription factors is not
strictly an enthalpic process. Temperature is known to have a strong influence
on many material properties of DNA ,such as persistence length and salt
release [@goethe2015], excluded volume effects [@driessen2014], and
repressor-DNA solubility [@elf2007], to name a few of many effects. To phenomenologically
characterize the influence of temperature on the fold-change in gene expression,
we considered that there was a constant entropic penalty (though inclusion of a
temperature-dependent entropic cost is discussed in  Chapter 9). We found that
inclusion of this parameter markedly improved the description of the data (solid
lines in left-hand side of @Fig:growth_results_intro (B)) and permitted data
collapse within experimental noise of data collected at 37$^\circ$ C (right-hand
side of @Fig:growth_results_intro (B)).

![**Performance of a simple thermodynamic model of simple repression in
diverse physiological states.** (A) Fold-change in gene expression
measurements in different carbon sources plotted against the average
repressor copy number (left) and free energy (right). Black line in the
left-hand panel is the predicted fold-change assuming no parameters are
modified. (B) Fold-change measurements at different temperatures plotted as a
function of the repressor copy number (left) and free energy (right).
Dashed-lines in left-hand plot show the predicted fold-change with a simple
rescaling of the thermal energy. Solid lines are predicted fold-change upon
inclusion of an entropic penalty. Points on right-hand plot were computed
using parameters with an entropic penalty. All measurements and errors
displayed are the mean and standard error of three to eight biological
replicates. Light-grey points in right hand panels are data from @garcia2011,
@brewster2014, @razo-mejia2018, and @chure2019, all of which were measured in
glucose-supplemented media at 37$^\circ$
C.](ch1_fig11){#fig:growth_results_intro short-caption="Performance of a
simple thermodynamic model of simple repression in diverse physiological
states."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The inclusion of a phenomenological entropic
parameter is not by any means meant to shut the book on temperature effects
in this model. Rather, it serves as a representation of what *may* help
explain these effects and demands more focused theoretical and experimental
work. To say that the current disagreement between theory and experiments
embodies the "strategic failure" of thermodynamic models is, in my view,
disingenuous. To say so would be to deem the initial failures of elasticity
theory to properly predict the influence of impurities and temperature on the
elastic constants of materials as a "strategic failure" in material science.
Initial phenomenological models of the effects of impurities and temperature on
elastic properties of solids [@friedel1974] led to several decades of focused theoretical and experimental
work that resulted in a complete predictive and mechanistic description
[@phillips2001]. Now is the time for a similar approach to biology in the context of
temperature and the regulation of gene expression. 


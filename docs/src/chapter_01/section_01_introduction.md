## Introduction

From archaea thriving in hydrothermal vents on the ocean floor to aspen trees
dominating a Coloradan mountainside, all forms of life are unified in their
ability to sense the state of their environment and, if they are lucky, adapt
to whatever the conditions may be. It is relatively easy for us
anthropormorphize these organisms and use terms such as "sense" and "adapt",
but can we do the same for things at the molecular scale? Adaptive phenomena
can be found traversing the spatial and temporal scales in biology, ranging
the nanosecond scale conformational switching of proteins
[@Fig:adaptation_levels (A)], to the minutes to hours long rewiring of
metabolic networks [@Fig:adaptation_levels (B)], to the 3.5 billion years of
speciation which has reared the biological diversity we see to day
[@Fig:adaptation_levels(C)]. Much as our archaeon and aspen trees at the
beginning of this paragraph, individual proteins are also capable of sensing
changes in the environment and adapting.

![**The spatial, temporal, and mechanistic scale of adaptation.** (A)
Molecular adaptation in this work is defined through the lens of allostery
where the activity of a protein complex is modulated the reversible binding
of a small molecule. These binding and unbinding events lead to rapid changes
in protein conformation whose behavior (both energetic and temporal) is
comparable to that of thermal motion. (B) Physiological adaptation here is
defined as the rewiring of biochemical reaction networks that lead to changes
in cellular behavior (such as chemotaxis) or metabolic capacity (such as
aerobic to fermentative metabolism). (C) Evolutionary adaptation is recorded
in the variation in the genetic sequence of regulatory molecules. Variations
in sequence influence the function of the proteins and RNAs they encode which
ultimately define the cellular fitness.](ch1_fig1){#fig:adaptation_levels
short-caption="The spatial, temporal, and mechanistic scale of adaptation."}


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This idea of molecular adaption is not novel by
any means and demands a brief foray into the history of bacterial growth and
the dawn of regulatory biology. In the late 1890's, Emilé Duclaux and his
graduate student Frédéric Diénert performed a series of experiments that
revealed the phenomenon of *enzymatic adaptation* where "the production of
diastases [enzymes] depends on the manner of nutrition" in which the cultures
were grown [@loison2013]. This is one of the the first observations of the
fact that, while an organism may be able to digest a certain sugar, it may
not *always* be able to do so. Rather, there seemed to be certain conditions
in which the production or formation of these enzymes could occur. In his
doctoral thesis in 1900, Diénert proposed two mechanisms for the origin of
enzymatic adaptation observed for galactozymase [@loison2013]. Either (a) the
presence of galactose directly transformed enzymes already present in the
cell (*S. cerevisiae*) into galactozymase or (b) that the galactose activated
the production of galactozymase *de novo* [@dienert1900, @loison2013].

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Nearly half a century later, Jacques Monod
would rediscover the phenomenon of enzymatic adaptation, this time in the
context of bacterial growth. In his 1941 work, *Sur un phénoméne nouveau de
croissance complexe dans les cultures bactériennes*, Monod for the first time
reported on the phenomenon of diauxic growth, shown in @Fig:diauxie_fig (A).
He noted that for some mixtures of carbon sources, the culture grew
"kinetically normal" meaning exponential growth to a saturating level [blue
points, @Fig:diauxie_fig (A)]. However, some mixtures (such as sucrose and
arabinose) led biphasic growth mode where the culture would grow
exponentially, undergo a period where growth had ceased, followed by again by
another round of exponential growth [blue points, @Fig:diauxie_fig (A)]. In
particular mixtures of three carbon sources, Monod even found that "triauxic"
growth could be observed [@Fig:diauxie_fig (B)].

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Monod immediately recognized a connection
between diauxic growth and enzymatic adaptation [@loison2013]. Despite his
work appearing 40 years after the pioneering work of Ducleaux and Diéner,
there had been little progress towards a mechanistic, needless to say
quantitative, explanation for the phenomenon of enzymatic adaptation. In
fact, Monod was particularly disappointed by the teleological "explanations"
for the phenomenon where the cells somehow controlled their behavior to
perform only the chemical reactions that were needed [@loison2013]. The
teleological approach to much of biology during this time period, especially
in the French scientific community, severely bothered Monod. To him, this
kind of approach lacked the "postulate of objectivity" that other fields of
science (particularly physics) had adopted and belonged to a pre-scientific
era of biology [@loison2013]. Near the middle of the 20th century, Monod
published a remarkable treatise on the the phenomena of enzymatic adaptation
with the level of quantitative rigor he thought it deserved [@monod1947]. In
this work, he set out to progressively deconstruct and invalidate a series of
hypotheses for the phenomenon of enzymatic adaptation. In doing so, he laid
the groundwork for what would become (in his opinion) his greatest
contribution to science, the nature of allosteric transitions [@loison2013,
@monod1963, @monod1965], a topic that will feature prominently in the
remainder of this thesis.

![**The phenomenon of enzymatic adaptation revealed in bacterial growth
curves.** (A) Optical density measurements of *Bacillus subtilis* cultures grown in a
mixture of sucrose and either glucose (blue points) or arabinose (green points).
Biphasic growth can be observed in the sucrose/arabinose mixture where the pause
in growth (white vertical line) corresponds to enzymatic adaptation. Data
digitized from @monod1941. (C) A triauxic growth curve of *Escherichia coli*
cells grown on a mixture of glucose, sorbitol, and glycerol. Periods of
enzymatic adaptation are highlighted by white vertical
lines.](ch1_fig2){#fig:diauxie_fig short-caption="The phenomenon of enzymatic
adaptation revealed in bacterial growth curves."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The di- and triauxic growth transitions shown
in @Fig:diauxie_fig, in my opinion, illustrate adaptive processes across the
biological scales, as were schematized in @Fig:adaptation_levels. The central
aim of this dissertation is to explore the biophysical mechanisms by which
these levels of adaptation -- molecular, physiological, and evolutionary,
schematized in @Fig:adaptation_levels -- are interconnected. Furthermore, in
the spirit of Monod, we seek to make our exploration quantitative and
leverage the tools of statistical physics to provide precise predictions from
pen-and-paper theory that can be rigorously tested through experiment. The
remaining sections of this chapter will summarize this journey and present
the key findings and will finish with a discussion of how these types of
models can be used to explore the predictability of evolution.









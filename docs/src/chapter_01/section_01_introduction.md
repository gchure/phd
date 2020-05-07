## Introduction

From archaea thriving in hydrothermal vents on the ocean floor to aspen trees
dominating a Coloradan mountainside, all forms of life are unified in their
obedience to (and influence on) their environment. Whether it be the
stochastic events that can trigger mass extinctions [@jablonski2001] or the
rapid "terraforming" of environments by their inhabitants (such as the
rusting of the oceans [@luo2018]), organisms and their environments are
continually engaged in a dialogue that demands the ability to adapt. Over the
past 3.5 billion years of evolution, life has evolved myriad clever ways to
combat (and exploit) environmental fluctuations to amplify reproductive
success. The mechanisms behind this adaptation are diverse and traverse the
biological scales, ranging from nanosecond-scale conformational switching of
proteins (@Fig:adaptation_levels (A)), to reconfiguration of metabolic
networks to consume different sugars (@Fig:adaptation_levels (B)), to
evolutionary trajectories that only become visible over many generations
(@Fig:adaptation_levels (C)). While "adaptation" is typically only associated
with organisms (at least colloquially), one can use the same language to
describe the microscopic operations of molecules.

![**The spatial, temporal, and mechanistic scales of adaptation.** (A)
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
ultimately define cellular fitness.](ch1_fig1){#fig:adaptation_levels
short-caption="The spatial, temporal, and mechanistic scale of adaptation."}


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The idea of molecular adaption is not novel and
demands a brief foray into the history of bacterial growth and the dawn of
regulatory biology. In the late 1890's, Emilé Duclaux and his graduate
student Frédéric Diénert performed a series of experiments illustrating that
the common yeast could only consume galactose after an incubation period with
the sugar. This led to a general conclusion that "the production of diastases
[enzymes] depends on the manner of nutrition" in which the cultures were
grown [@loison2013], a phenomenon later coined *enzymatic adaptation*. This
is one of the the first observations of the fact that, while an organism may
be able to digest a certain sugar, it may not *always* be able to do so.
Rather, there seemed to be certain conditions in which the production or
formation of these enzymes could occur. In his doctoral thesis in 1900,
Diénert proposed two mechanisms for the origin of enzymatic adaptation
observed for galactozymase in *S. cerevisiae* [@loison2013]. Either (a) the
presence of galactose *directly* transformed enzymes already present in the
cell into galactozymase or (b) that the galactose activated the production of
galactozymase *de novo* [@dienert1900; @loison2013].

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Nearly half a century later, Jacques Monod
would rediscover the phenomenon of enzymatic adaptation, this time in the
context of bacterial growth. In his 1941 work, Monod for the first time
reported on the phenomenon of diauxic or biphasic growth, shown in @Fig:diauxie_fig (A).
He noted that for some mixtures of carbon sources, the culture grew
"kinetically normal" meaning they grew exponentially to saturation (blue
points, @Fig:diauxie_fig (A)). However, some mixtures (such as sucrose and
arabinose) led to biphasic growth where the culture would grow exponentially,
undergo a period where growth had ceased (lasting typically 20 - 60 minutes),
followed again by another round of exponential growth (green points,
@Fig:diauxie_fig (A)). Additionally, Monod showed that the onset of this
diauxic shift could be tuned by varying the relative concentrations of the
carbon sources, revealing a controllable chemical basis for the adaptation
(@Fig:diauxie_fig (B)).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Monod immediately made the connection
between diauxic growth and enzymatic adaptation [@loison2013]. Despite his
work appearing 40 years after the pioneering work of Ducleaux and Diénért,
there had been little progress towards a mechanistic, needless to say
quantitative, explanation for the phenomenon. In
fact, Monod was particularly disappointed by the teleological explanations
where the cells simply changed  their behavior to
perform only the chemical reactions that were "needed" [@loison2013]. The
teleological approach to much of biology during this time period, especially
in the French scientific community, severely bothered Monod. To him, this
kind of approach belonged to a pre-scientific era and lacked the "postulate of
objectivity" that other fields of science (particularly physics) had adopted [@loison2013]. Near the middle of the 20th century, Monod
published a 60-page treatise on the the phenomenon of enzymatic adaptation
with the level of quantitative rigor he thought it deserved [@monod1947]. In
this work, he set out to progressively deconstruct and invalidate a series of
hypotheses for the phenomenon of enzymatic adaptation. In doing so, he laid
the groundwork for what would become (in his opinion) his greatest
contribution to science, the nature of allosteric transitions [@loison2013;
@monod1963; @monod1965], a topic that will feature prominently in the
remainder of this thesis.

![**The phenomenon of enzymatic adaptation revealed in bacterial growth
curves.** (A) Optical density measurements of *Bacillus subtilis* cultures
grown in a mixture of sucrose and either glucose (blue points) or arabinose
(green points). Biphasic growth can be observed in the sucrose/arabinose
mixture where the pause in growth (white shaded region) corresponds to
enzymatic adaptation. Data digitized from @monod1941. (B) Diauxic growth
curves of *Escherichia coli* cells grown on a mixture of glucose and sorbitol
in different proportions. Data digitized from @monod1947. Periods of
enzymatic adaptation are highlighted by white vertical lines. In both plots,
a univariate spline interpolation was used to draw lines to reflect data
presentation of original literature. The [Python code
(`ch1_fig2.py`)](https://github.com/gchure/phd/blob/master/src/chapter_01/code/ch1_fig2.py)
used to generate this figure can be accessed via the thesis [GitHub
repository](https://github.com/gchure/phd).](ch1_fig2){#fig:diauxie_fig
short-caption="The phenomenon of enzymatic adaptation revealed in bacterial
growth curves."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The diauxic growth transitions shown in
@Fig:diauxie_fig illustrate adaptive processes across the biological scales,
as were schematized in @Fig:adaptation_levels. While it was not known to
Monod at the time, we now know that many cases of enzymatic adaptation are
driven by the regulation of gene expression. As the bacterial culture
approaches the auxic shift, the presence or absence of the substrate is
sensed by regulatory molecules that control whether the genes encoding the
enzymes for metabolism of the substrate are expressed. This represents the
level of *molecular adaptation* where, given binding or unbinding of the
substrate molecule, the activity of the regulatory protein is modulated. The
amino acid sequence of these proteinaceous regulators are the product of
billions of years of evolution with regions of the protein (such as the
DNA-binding and inducer-binding domains) undergoing *evolutionary
adaptation* that define how the regulatory molecule senses and responds to these
signals. Finally, the precision with which these genes are regulated is
determined by their sensitivity to physiological states, capturing the level
of *physiological adaptation*.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The central aim of this dissertation is to
explore the biophysical mechanisms by which these levels of adaptation --
molecular, physiological, and evolutionary -- are interconnected. Furthermore, in the spirit of
Monod, we seek to make our exploration quantitative and leverage the tools of
statistical physics to provide precise predictions from pen-and-paper theory
that can be rigorously tested through experiment. The remaining sections of
this chapter will outline the major topics of this thesis and place them in a
historical context alongside the work of Monod. This thesis is structured such
that **Chapters 1 -- 5** present a self-contained summary of how quantitative
methods can be used to interrogate and understand the molecular biophysics of
adaptation. The remainder, **Chapters 6 -- 9**, are detailed supplements to
Chapters 1 -- 5 and are targeted to readers who want to dig into the weeds of
statistical inference, error estimation, and analytic properties of the
theoretical models.  

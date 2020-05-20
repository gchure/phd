## On Facing the Elements 

The first four chapters of this work encompass myriad perspectives of adaptive
processes at the level of transcription regulation. However, just as important
as the regulation is the action of the gene that is ultimately expressed. While
Monod's work described in the preceding sections was focused on the expression
of enzymes, we now turn to yet another level of physiological adaptation in
bacteria -- the regulation of turgor pressure. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In the wild, microbes are constantly faced with an array of environmental
insults ranging from changes in temperature, availability of oxygen for aerobic
respiration, and even chemical warfare from neighboring microbial  communities
[@czaran2002]. One such environmental challenge microbes often face is the
variation in the osmolarity of their surroundings. Changes in ion concentrations
can result in large volumes of water rushing through the cell membrane, leading
to rupture of the membrane and ultimately cell death (@Fig:mscl_intro (A)). Unsurprisingly, all domains of life have
evolved clever mechanisms to combat these osmotic shocks and regulate their
internal turgor pressure. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;One such mechanism for osmoregulation in *E. coli* is through the action of
mechanosensitive ion channels -- large, transmembrane structures which sense
tension in the cell membrane. Exposure to a hypo-osmotic shock (where water
rushes across the cell membrane *into* the cell), a change in membrane tension
is sensed by these mechanosensitive channels, triggering a conformational change
which opens a pore in the membrane without rupture (@Fig:mscl_intro (B)). This acts as a pressure
release valve, providing a means for turgor pressure to be relieved without a
potentially fatal burst. This phenomenon represents yet another system in which
adaptation can be found at the molecular, evolutionary, and physiological
levels. In **Chapters 5 and 9** of this dissertation, we explore a fundamental
question -- how many mechanosensitive channels does a cell need to have an
appreciable chance at surviving an osmotic shock?


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To approach this question, Heun Jin Lee and I collaborated on the experimental
and data analysis components, respectively. This is a project that had been in
preparation for several years before I had the privilege of joining the team in
the summer of 2017. While the experimental techniques used to probe
transcriptional regulation were far from simple, they pale in comparison to
those employed by Heun Jin. This project required an enormous amount of
molecular biology to generate the necessary strains in which the number of
mechanosensitive channels could be tuned across three orders of magnitude and
measured with precision. This process involved reworking classic techniques in
molecular biology to remove the presence of osmotic shocks which would prove
fatal for strains with few or no mechanosensitive channels. On top of
the complex biochemistry, Heun Jin developed a clever microfluidic system where
osmotic shocks could be imaged in real time at the single cell level. While the
majority of Chapter 5 focuses on the analysis and interpretation of the data,
none of it would have been possible without Heun Jin's Herculean efforts. 

![**The connection between mechanosensitive channel copy number and
probability of survival.** (A) In the absence of mechanosensitive channels,
water rushing across the membrane during a hypo-osmotic shock can lead to
membrane rupture and large-scale release of intracellular components into the
extracellular space, resulting in cell death. (B) In the presence of
mechanosensitive channels (specifically, the major *E. coli* mechanosensitive
channel MscL as shown in yellow), increased membrane tension results in a
conformational change of the channel, resulting in the expulsion of water and
some small constituents of the intracellular milieu. The inferred survival
probability curves for slow and fast shock exchange rates are shown in (C)
and (D), respectively. Different shaded purple regions correspond to
different credible regions of the estimates. The [Python code
(`ch1_fig12.py`)](https://github.com/gchure/phd/blob/master/src/chapter_01/code/ch1_fig12.py)
used to generate this figure can be accessed via the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch1_fig12){#fig:mscl_intro
short-caption="The connection between mechanosensitive channel number and
probability of survival."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;While we leave the details of the inference to Chapter 5 and the supplemental
Chapter 9, the survival probability as a function of the total mechanosensitive
channel number is given for "slow" and "fast" osmotic shocks in @Fig:mscl_intro
(C) and (D), respectively. The credible regions in this plot illustrate that for
an $\approx 80\%$ chance of surviving either a slow or a fast osmotic shock, at
least $\approx 500$ channels are needed. This number is in agreement with recent
proteomic measurements in *E. coli* [@li2014; @schmidt2016; @soufi2015], but
are at odds with current theoretical models. While it is difficult to
theoretically define a survival/death criterion, current physical models predict
that only a few mechanosensitive channels (specifically MscL) are needed to relieve
even large increases in membrane tension. These findings illustrate another
avenue in which the disagreement between theory and careful, quantitative
experiments reveal gaps in our understanding of fundamental biological
phenomena.

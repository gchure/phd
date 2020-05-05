## On Using Free Energy to Examine Evolutionary Adaptation

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Allow us to briefly return to Monod and his
biphasic growth curves in the mid 1940's. At this point in scientific
history, the French vision of biology had taken a strongly finalistic and
vitalistic turn [@loison2013]. In particular, a neo-Lamarckian view had been
employed to explain the phenomenon of enzymatic adaptation where the the
enzymes appropriate for digesting the substrate could be spontaneously formed
out the bacterial cytoplasm and inherited by the cell's descendents,
completely independent of genes. In general, this approach to biology deeply
frustrated Monod and strongly influenced his desire to "physicalize" the
science [@loison2013]. One tool he knew was critical to this mission was the
burgeoning field of genetics. In the mid 1930's Monod undertook a short
retreat to Thomas Hunt Morgan's lab at Caltech where he was introduced to
genetics which he later remarked to as "biology's first discipline"
[@loison2013]. This visit had a profound impact on Monod, who reflected upon it
some three decades later:

> *"Upon my return to France, I had again taken up the study of bacterial growth.
> But my mind remained full of the concepts of genetics and I was confident of
> its ability to analyze and convinced that one day these ideas would be applied
> to bacteria."* [@monod1966]

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Once he returned to his study of bacterial
growth and enzymatic adaptation, he was confronted with incorporating the
role of genetic inheritance into his mechanistic explanations. In the mid
1940's, Monod and his coworkers had begun experimenting with a strain of *E.
coli* which was unable to digest lactose, termed $L-$. When grown on a mixture of glucose and lactose, this strain
would not display a diauxic shift and would only be able to consume the
glucose in the medium (@Fig:lacneg, black). However, Monod and his coworker Alice Audureau
discovered a mutation in this strain which *enabled* the digestion of
lactose, termed $L+$ [@monod1947]. The growth curve of this strain had the
striking feature of diauxic growth. Rather than this mutation merely enabling
the digestion of lactose, it did so in a non-constitutive manner and
preserved the phenomenon of adaptation. This was an important step forward in
Monod's understanding of enzymatic adaptation [@loison2013], revealing that
it was a "truly genetic property" [@monod1966].

![**Growth curves of lactose-positive and lactose-negative *E. coli* strains
on a glucose/lactose mixture.** Black curve shows the growth curve of an *E. coli* strain unable
to digest lactose grown on a glucose/lactose mixed medium. Red curve shows a
mutant of the same *E. coli* strain which is able to consume lactose. The latter
displays a diauxic growth cycle with an adaptive period (highlighted in white),
illustrating that enzymatic adaptation is a truly genetic property. Figure
adapted from @monod1947. Lines correspond to univariate splines fit to the data
to retain the data presentation from the literature. The [Python code
(`ch1_fig6.py`)](https://github.com/gchure/phd/blob/master/src/chapter_01/code/ch1_fig6.py)
used to generate this figure can be accessed via the thesis [GitHub
repository](https://github.com/gchure/phd).](ch1_fig6){#fig:lacneg short-caption="Growth curves of
lactose-positive and lactose-negative *E. coli* strains."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This finding illustrates the level of
evolutionary adaptation operating at the level of molecules. While it is
difficult to find any literature dissecting this particular $L+$ mutation,
it is not difficult to imagine several different mechanisms by how it could
be manifest. One such explanation is that this $L+$ mutation is within a
transcriptional regulator itself where a deficiency in the ability to respond
to the presence of lactose (and decreasing glucose concentration) had been restored. Such
mutations are the crux of **Chapter 3** and the corresponding supplemental
**Chapter 7** of this dissertation.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As summarized in the previous section and
discussed in depth in Chapter 2, Chapter 3 and the associated supplemental
Chapter 7 of this dissertation focus on the influence of mutation within the
allosteric transcription factor. Furthermore, Chapter 3 presents a generic
mechanism by which shifts in the free energy can be mapped directly to
changes in values of the biophysical parameters. This chapter, much like
Chapter 2, was borne out of a wonderful collaboration with Manuel Razo-Mejia,
Stephanie L. Barnes, Nathan M. Belliveau, Tal Einav, and Zofii A. Kaczmarek.
Being able to launch another collaborative effort afforded us the opportunity
to both develop a new theoretical interpretation for how mutations influence
the free energy and acquire enough experimental data to thoroughly test it.

![**Mutations lead to shifts in free energy, permitting prediction of double
mutant phenotypes.** Consider a wild-type
bacterium which on, on average, exhibits a fold-change of $\approx$ 0.3
and a free energy of $-1\, k_BT$ (grey point in (B)). We can consider that a single mutation (either
orange or purple) changes the mean fold-change and therefore the free energy,
translating the measurement about the master curve (black line in (B)). Assuming
there are no interactions between the two single mutations, a null
hypothesis predicts that for the double mutant (blue bacterium in (A) and point in
(B), the net free energy is simply the sum of the individual free energy
shifts.The [Python code
(`ch1_fig7.py`)](https://github.com/gchure/phd/blob/master/src/chapter_01/code/ch1_fig7.py)
used to generate this figure can be accessed via the thesis [GitHub
repository](https://github.com/gchure/phd).](ch1_fig7){#fig:pedagogical_delF_intro short-caption="Mutations lead to
predictive shifts in free energy."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The primary conceptual development of Chapter 3 is illustrated in
@Fig:pedagogical_delF_intro. Theoretically, we consider a bacterial strain with
an allosteric repressor (which we term the "wild-type" repressor) that has been
sufficiently characterized. Given enough parametric knowledge of the system, we
can easily compute both its predicted average fold-change in gene expression
along with the corresponding free energy. However, once a mutation has been
introduced *into the repressor protein*  (resulting in a non-synonymous amino
acid change), we are once again ignorant *a priori* of what changes, if any,
that mutation may have imparted on the system. In @Fig:pedagogical_delF_intro,
we examine two separate hypothetical mutations, shown in purple and orange, which
significantly change the character of the system by either increasing or
decreasing the fold-change in gene expression, respectively. If we assume that
these mutations do not change the underlying physics of the system, we are
permitted to use the theoretical framework outlined in Chapter 2 and in
@Fig:induction_intro to characterize each mutation and determine what
biophysical parameters have been changed. This permits us to calculate the new
free energy of the system ($F_\text{mutation 1}$) as well as the shift in free
energy from the wild-type value, 
$$
\Delta F_\text{mutation 1} = F_\text{mutation 1} - F_\text{wt}.
$${#eq:intro_delF_definition}
As will be described in detail in Chapter 3 and the supplemental Chapter 7, the
precise value of this free energy shift $\Delta F$ can be directly computed
given sufficient parametric knowledge.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This formalism provides a mathematical
hypothesis for how double mutants may behave. Given known values for $\Delta
F$ of each mutation in isolation, can we compute the shift in free energy of
the pairwise double mutant $\Delta F_\text{mutations 1 \& 2}$?
@Eq:intro_delF_definition presents a mathematical null hypothesis that the
net shift in the free energy is simply the sum of the individual shifts in
free energy,
$$
\Delta F_\text{mutations 1 \& 2} = \Delta F_\text{mutation 1} +
\Delta F_\text{mutation 2},
$${#eq:delF_null_intro}

assuming there is no interaction between the mutations. Given the fact that we can compute the fold-change in gene expression given
knowledge of the free energy, we can therefore predict the double mutant
phenotype *a priori*, a prediction not possible prior to this work.


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Over the course of two years (while this theory was in the works), the
experimental cast of characters (Stephanie L. Barnes, Nathan M. Belliveau,
Manuel Razo-Mejia, and Zofii A. Kaczmarek) and I made a series of mutations
in the LacI repressor that we had characterized in the work presented in
Chapter 2. These mutations included three point mutations in the DNA binding
domain of the repressor, four mutations in the inducer binding domain, nine
double mutants (one inducer binding and one DNA binding each), across four
repressor copy numbers and three operator sequences. While this process of
strain generation and data collection is not the primary focus of the work,
it took $\approx 80\%$ of the effort. Without them, this work would have
remained an untested theoretical novelty. While we leave many of the rich
details of this prediction to the reader in Chapter 3, we showcase our 
experimental success in  @Fig:double_muts_intro (B) where the predicted induction profiles
of nine double mutants (light blue shaded regions) are overlaid with their
experimental measurements (points). The near 100\% agreement between theory and
experiment illustrates the utility of using free energy shifts as a means to
predict new phenotypes.

![**Theoretical prediction and experimental validation of double mutant
phenotypes.** (A) Cartoon representation of the LacI repressor with mutations in
the inducer binding domain and DNA binding domain represented by hats and socks,
respectively. While the mutations have known chemical features, we characterize
each mutation as potentially modifying four biophysical parameters,
the dissociation constants ($K_A$, $K_I$), the relative energy difference
between active and inactive states of the repressor ($\Delta\varepsilon_{AI}$) for inducer binding mutants, or
the DNA affinity ($\Delta\varepsilon_{RA}$) for DNA binding mutants. (B) Predicted induction
profiles for pairwise double mutants are shown as blue shaded regions
representing the uncertainty in our predictions. Experimental measurements are
shown as blue points (means of at least 10 biological replicates). Each row
corresponds to a single DNA binding domain mutation and each column to a single
inducer binding domain mutation. The [Python code
(`ch3_fig5.py`)](https://github.com/gchure/phd/blob/master/src/chapter_03/code/ch3_fig5.py)
used to generate this figure can be accessed via the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch1_fig8){#fig:double_muts_intro
short-caption="Theoretical prediction and experimental validation of double
mutant phenotypes."}

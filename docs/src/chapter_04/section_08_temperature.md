## Fold-change Dependence on Temperature Variation

Unlike the identity of the carbon source, the temperature of the system
is explicitly stated in @Eq:foldchange_simple and @Eq:growth_pact
where $\Delta\varepsilon_{RA}$ and $\Delta\varepsilon_{AI}$ are defined
relative to the thermal energy of the system in which they were
determined. This scaling is mathematically quantified as $k_BT$ dividing
the exponentiated terms in @Eq:foldchange_simple  and @Eq:growth_pact. As
all biophysical parameters were determined at a reference temperature of
37$^\circ$ C, any change in the growth temperature must be included as a
correction factor. The simplest approach is to rescale the energy by the
relative change in temperature. This is a simple multiplicative factor
of $\phi_T = k_BT_{ref} / k_BT_{exp}$ where $T_{ref}$ is the reference temperature of
37$^\circ$ C and $T_{exp}$ is the experimental temperature. This is an
intuitive result since an increase in temperature relative to the
reference results in $\phi_T < 1$, weakening the binding. Similarly,
decreasing the temperature scales $\phi_T > 1$, strengthening the
binding relative to that of the reference temperature.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;@Fig:deltaF_temp (A) shows the measured fold-change in gene
expression (points) plotted against the theoretical prediction with this
correction factor (orange line). It is immediately evident that a simple
rescaling of the energetic parameters is not sufficient for the
32$^\circ$ C condition and slightly underestimates the fold-change in
the 42$^\circ$ C condition. To identify the source of this disagreement,
we can again examine the free energy shift $\Delta F$. As both
$\Delta\varepsilon_{AI}$ and $\Delta\varepsilon_{R}$ are scaled to the
thermal energy, $\Delta F$ defined as $F_{T} - F_{ref}$ can be directly
calculated as
$$
\Delta F = k_BT_{exp} \left[ -\log{1 + e^{-{\Delta\varepsilon_{AI} \over
 k_BT_{ref}}} \over 1 + e^{-{\phi_T\Delta\varepsilon_{AI} \over
 k_BT_{exp}}}}- \log{R_{T} \over R_{ref}}\right] +
 \Delta\varepsilon_{R}\left(1 -\phi_T\right).
$${#eq:delF_temp}
This prediction along with the empirically determined $\Delta F$ is
shown in @Fig:deltaF_temp (B). Again, we see that this simple
correction factor significantly undershoots or overshoots the observed
$\Delta F$ for 32$^\circ$ C and 42$^\circ$ C, respectively, indicating
that there are temperature dependent effects that are not accounted for
in the simplest null model of temperature dependence.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The model described by @Eq:foldchange_simple
and @Eq:growth_pact subsumes the myriad rich dynamical processes
underlying protein binding and conformational changes into two effective
energies, $\Delta\varepsilon_{R}$ and $\Delta\varepsilon_{AI}$. By no means
is this done to undercut the importance of these details in transcriptional
regulation. Rather, it reduces the degrees of freedom in this objectively
complex system to the set of the details critical to particular conditions in
which we want to draw predictions. All prior dissections of this
thermodynamic model have been performed at a single temperature, abrogating
the need to consider temperature dependent effects. As we now vary
temperature, we must consider details that are swept into the effective
energies.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The model presented here only considers entropy by enumerating the
multiplicity of states in which the repressor can bind to the DNA
nonspecifically, resulting in terms of the form ${R / N_{NS}}$.
However, there are many other temperature-dependent entropic
contributions to the effective energies such as the fraction of
repressors bound to DNA versus in solution [@elf2007; @kao-huang1977],
the vibrational entropy of the repressor [@goethe2015], or
conformational entropy of the genome [@driessen2014; @mondal2011]. We
can consider the effective energies $\Delta\varepsilon_R$ and
$\Delta\varepsilon_{AI}$ as having generic temperature
dependent-entropic components $\Delta S_{R}$ and $\Delta S_{AI}$,
$$
\Delta\varepsilon_R = \Delta H_R - T\Delta S_{R}, 
$${#eq:epra_entropy}
and
$$
\Delta\varepsilon_{AI} = \Delta H_{AI} - T \Delta S_{AI},
$${#eq:epai_entropy}
where $\Delta H_R$ and $\Delta H_{AI}$ is the enthalpic contribution to the energies $\Delta\varepsilon_R$ and
$\Delta\varepsilon_{AI}$, respectively. Given the fold-change
measurements at 32$^\circ$ C and 42$^\circ$ C, we estimated the entropic
parameters $\Delta S_R$ and $\Delta S_ {AI}$ under the constraints that
at 37$^\circ$ C, $\Delta\varepsilon_R = -13.9\, k_BT$ and
$\Delta\varepsilon_{AI} = 4.5\,k_BT$ (see supplemental Chapter 8 for more
discussion on this parameter estimation). The grey shaded lines in 
@Fig:deltaF_temp show the result of this fit where the width
represents the 95\% credible region of the prediction given the estimated
values of $\Delta S_R$ and $\Delta S_{AI}$. Including this
phenomenological treatment of the entropy improves the prediction of the
fold-change in gene expression (@Fig:deltaF_temp (A)) as well as shift in free energy
(@Fig:deltaF_temp (B)). This phenomenological description
suggests that even small shifts in temperature can drastically alter the
expression of a genetic circuit simply by tuning hidden entropic effects
rather than scaling the difference in affinity between specific and
nonspecific binding.

![**Temperature effects on the fold-change in gene expression and free energy.**
(A) The fold-change in gene expression for growth in glucose supplemented medium
at 32$^\circ$ C (left) and 42$^\circ$ C (right). Points and errors correspond to
the mean and standard error of five biological replicates. Predictions of the
fold-change are shown without correcting for temperature (purple), with
multiplicative scaling (orange), and with an entropic penalty (grey). The width
of the prediction of the entropic penalty is the 95\% credible region. (B)
Predicted and observed shifts in the free energy  for growth in glucose medium
at 32$^\circ$ C (left) and 42$^\circ$ C (right). Points correspond to the median
of the inferred shift in free energy. Vertical error bars indicate the bounds of
the 95\% credible region. Horizontal position and error corresponds to the mean
and standard error for the repressor count over five to eight biological
replicates. The [Python code
(`ch4_fig5.py`)](https://github.com/gchure/phd/blob/master/src/chapter_04/code/ch4_fig5.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch4_fig5){#fig:deltaF_temp short-caption="Temperature effects on the
fold-change in gene expression and free energy."}

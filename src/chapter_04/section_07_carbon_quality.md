## Fold-change dependence on carbon quality

Given *a priori* knowledge of the biophysical parameter values
[@garcia2011; @razo-mejia2018] present in @Eq:foldchange_simple 
and @Eq:growth_pact and direct measurement of the repressor copy number, we made measurements of
the fold-change in gene expression for each growth medium to test the
prediction (@Fig:carbon_quality_delF (A)). We find that the
measurements fall upon the predicted theoretical curve within error,
suggesting that the values of the energetic terms in the model are
insensitive to changing carbon sources. This is notable as glucose,
glycerol, and acetate are metabolized via different pathways, changing
the metabolite and protein composition of the cytosol
[@martinez-gomez2012; @kim2007]. This result underscores the utility of
these thermodynamic parameters as quantitative traits in the study of
growth-condition dependent gene expression.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We have shown in the preceding chapters that
@Eq:foldchange_simple can be rewritten into a Fermi function of the form
$$
\text{fold-change} = \frac{1}{1 + e^{-F / k_BT}},
$${#eq:growth_fermi} 
where $F$ is the effective free energy difference between the repressor bound
and unbound states of the promoter as described in Chapters 2 and 3.
For the case of an allosteric simple repression architecture, and given knowledge of the
values of the biophysical parameters, $F$ can be directly calculated as
$$
F = k_BT \left[-\log\left(1 + e^{-\Delta\varepsilon_{AI} / k_BT}\right)^{-1} -
\log\frac{R}{N_{NS}} + \frac{\Delta\varepsilon_{RA}}{k_BT}\right].
$${#eq:growth_bohr_parameter}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We have recently interrogated how this
formalism can be used to map mutations within the repressor to
biophysical parameters by examining the difference in free energy
between a mutant and reference (wild-type) strain, $\Delta F = F_{mut} -
F_{ref}$[@chure2019 and Chapter 3]. This approach revealed that different parametric
changes yield characteristic response functions to changing inducer
concentrations. Rather than using wild-type and mutant variants of the
repressor, we can choose a reference condition and compare how the free
energy changes between different growth media. Here, we choose the
reference condition to be a sample grown at 37$^\circ$ C with glucose as
the available carbon source and a repressor copy number $R = 100$ per
cell. Under the hypothesis that the only variable parameter in these
growth conditions is the repressor copy number $R$, the shift in free
energy $\Delta F$ becomes
$$
\Delta F = F_{C} - F_{ref} = -\log\frac{R_{C}}{R_{ref}},
$${#eq:delF_carbon}
where $F_C$ and $R_C$ correspond to the free energy and repressor copy number
of the different growth conditions. This concise prediction serves as a
quantitative measure of how robust the energetic parameters
$\Delta\varepsilon_{RA}$ and $\Delta\varepsilon_{AI}$ are in units of $k_BT$
(@Fig:carbon_quality_delF (B)). In using free energy shifts as a diagnostic,
one can immediately determine the effect of the perturbation on the parameter
values by quantifying the disagreement between the observed and predicted
$\Delta F$ as the parameters and the free energy are both in the same natural
units.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We inferred the observed free energy for the fold-change measurements shown
in @Fig:carbon_quality_delF (A)(as described previously [@chure2019] and in
Chapter 7) and compared it to the theoretical prediction of @Eq:delF_carbon,
shown in @Fig:carbon_quality_delF (C). Again, we see the observed change in
free energy is in strong agreement with our theoretical predictions. This
agreement indicates that the free energy shift $\Delta F$ can be used in
multiple contexts to capture the energetic consequences of physiological and
evolutionary perturbations between different states of the system. The
insensitivity of the biophysical parameters to these distinctly different
physiological states demonstrates that $\Delta\varepsilon_{AI}$ and
$\Delta\varepsilon_{R}$ are material properties of the repressor defined by
the intricate hydrogen bonding networks of its constituent amino acids rather
than by the chemical constituency of its surroundings. Tuning temperature,
however, can change these material properties.

![**Fold-change in gene expression and free energy shifts in different growth
media.** (A) measurements of the fold-change in gene expression plotted against
the measured repressor copy number. Black line is the prediction and the width
corresponds to the uncertainty in the DNA binding energy as reported in
@garcia2011. Points and errors re the mean and standard error of five to eight
replicates for each ATC induction condition. (B) An example of how fold-change
is mapped onto changes in the free energy. Top panel shows the fold-change in
gene expression at specific repressor copy numbers (points). Arbitrarily choosing
100 repressors per cell as the reference state (white point) any measurement with
lager or smaller fold-change has a negative or positive shift in the free
energy, shown as red and green panels respectively. If $R$ is the only changing
variable, the free energy shift should be a linear function of $R$ (bottom pane).
(C) The inferred free energy shift from the fold-change measurements in (A).
Black line is the prediction. Points correspond to the mean and standard error
of the repressor copy number measurements over five to eight biological
replicates. Vertical error bars correspond to the 95\% credible region of the
inferred energy. The [Python code
(`ch4_fig4.py`)](https://github.com/gchure/phd/blob/master/src/chapter_04/code/ch4_fig4.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch4_fig4){#fig:carbon_quality_delF
short-caption="Fold-change in gene expression and free energy shifts in
different growth media."}
## Additional Characterization of DNA Binding Mutants 
In Chapter 3, we estimated the DNA binding energy o f each mutant
using the mutant strains that had approximately 260 repressors per cell.
In this section, we examine the effect of the choice of fit strain on
the predictions of both the induction profiles and $\Delta F$ for each
DNA binding domain mutant.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We applied the statistical model derived in
Section 2 of this chapter for each unique strain of the DNA binding
mutants and estimated the DNA binding energy. The median of the
posterior distribution along with the upper and lower bounds of the 95%
credible region are reported in Table 7.1. We found that the choice of fitting strain
did not strongly influence the estimate of the DNA binding energy. The
largest deviations appear for the weakest binding mutants paired with
the lowest repressor copy number. In these cases, such as for Q18A, the
difference in binding energy between the repressor copy numbers is
$\approx 1\, k_BT$ which is small compared to the overall DNA binding
energy. Using these energies, we computed the predicted induction
profiles of each mutant with different repressor copy numbers, shown in
@Fig:DNA_profile_pairwise_comparisons. In this plot, the
rows correspond to the repressor copy number of the strain used to
estimate the DNA binding energy. The columns correspond to the repressor
copy number of the predicted strains. The diagonals, shaded in grey,
show the induction profile of the fit strain along with the
corresponding data. In all cases, we find that the predicted profiles
are relatively accurate with the largest deviations resulting from using
the lowest repressor copy number as the fit strain.

|**Mutant**| **Repressors** | **DNA Binding Energy** [$\mathbf{k_BT}$] |
|:--:| :--:| :--:|
| Q18A | 60 | $-9.8_{-0.2}^{+0.2}$ |
| | 124|  $-10.3_{-0.1}^{+0.1}$ |
| | 260| $-11.0_{-0.1}^{+0.1}$ |
| | 1220| $-11.3_{-0.1}^{+0.1}$ |
| Q18M | 60 | $-15.83^{+0.08}_{-0.08}$ | 
| | 124| $-15.7_{-0.1}^{+0.1}$|   
| | 260 | $-15.43_{-0.06}^{+0.07}$ | 
| | 1220 | $-15.27_{-0.07}^{_+0.07}$|
| Y17I| 60 | $-9.4_{-0.3}^{+0.3}$|
| | 124 | $-9.5_{-0.1}^{+0.1}$|
| | 260 | $-9.9_{-0.1}^{+0.1}$ |
| | 1220| $-10.1_{-0.2}^{+0.2}$ |                                
: Estimated DNA binding energy for DNA binding domain mutants with different
repressor copy numbers. Reported values are the median of the posterior
distribution with the upper and lower bounds of the 95\% credible regions.


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The predicted change in free energy $\Delta F$
using each fit strain can be seen in @Fig:DNA_delF_pairwise_comparison. In
this figure, the rows represent the repressor copy number of the strain to
which the DNA binding energy was fit whereas the columns correspond to each
mutant. In each plot, we have shown the data for all repressor copy numbers
with the fit strain represented by white filled circles. Much as for the
induction profiles, we see little difference in the predicted $\Delta F$ for
each strain, all of which accurately describe the inferred free energies. The
ability to accurately predict the majority of the induction profiles of each
mutant with repressor copy numbers ranging over two orders of magnitude
strengthens our assessment that for these DNA binding domain mutations, only
the DNA binding energy is modified.

![**Pairwise comparisons of DNA binding mutant induction profiles.** Rows
correspond to the repressor copy number of the strain used to estimate
the DNA binding energy for each mutant. Columns correspond to the
repressor copy number of the strains that are predicted. Diagonals in
which the data used to estimate the DNA binding energy are shown with a
gray background.](ch7_figS12){#fig:DNA_profile_pairwise_comparisons
short-caption="Pairwise comparisons of DNA binding mutant estimated induction profiles."}

![**Dependence of fitting strain on $\Delta F$ predictions of DNA binding
domain mutants.** Rows correspond to the repressor copy number used to
estimate the DNA binding energy. Columns correspond to the particular
mutant. Colored lines are the bounds of the 95% credible region of the
predicted $\Delta F$. Open face points indicate the strain to which the
DNA binding energy was fit.](ch7_figS13){#fig:DNA_delF_pairwise_comparison
short-caption="Dependence of fitting strain on $\Delta F$ predictions of DNA
binding domain mutants.**}




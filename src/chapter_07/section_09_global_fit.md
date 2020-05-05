## Parameter Estimation Using All Induction Profiles
In Chapter 3 and Sec. 7.1 and 7.2 of this chapter, we have
laid out our strategy for inferring the the various parameters of our
model to a single induction profile and using the resulting values to
predict the free energy and induction profiles of other strains. In this
section, we estimate the parameters using all induction profiles of a
single mutant and using the estimated values to predict the free energy
profiles.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The inferred DNA binding energies considering induction profiles of all
repressor copy numbers for the three DNA binding mutants are reported in
Table 7.5. These parameters are close to those
reported in Table 7.1 for each repressor copy number with Q18A
showing the largest differences. The resulting induction profiles and
predicted change in free energy for these mutants can be seen in @Fig:global_DNA_profiles. 
Overall, the induction profiles match the data to an appreciable
agree. We acknowledge that even when using *all* repressor copy numbers,
the fit to Q18A remains imperfect. However we contend that this
disagreement is comparable to that observed in @razo-mejia2018 which
described the induction profile of the wild-type repressor. We find that
the predicted change in free energy (bottom row in @Fig:global_DNA_profiles (B))
narrows compared to that in @Fig:DNA_delF_pairwise_comparison
and Fig. 3.3 of Chapter 3, confirming that considering all induction profiles improves our
inference of the most-likely DNA binding energy. There appears to be a
very slight trend in the $\Delta F$ for Q18A at higher inducer
concentrations, though the overall change in free energy from 0 to 5000
$\mu$M IPTG is small.

![**Induction profiles and predicted change in free energy using parameters
estimated from the complete data sets.** Top row shows fold-change
measurements (points) as mean and standard error with ten to fifteen biological replicates. Shaded lines correspond to the 95\% credible regions of
the induction profiles using the estimated values of the DNA binding energies
reported in Table 7.5. Bottom row shows the 95\% credible regions of the
predicted change in free energy (shaded lines) along with the inferred free
energy of data shown in the top row. In all plots, the inducer concentration
is shown on a symmetric log scale with linear scaling between 0 and
$10^{-2}\,\mu$M and log scaling elsewhere. The [Python code (`ch7_figS21.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS21.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS21){#fig:global_DNA_profiles short-caption="Induction profiles and predicted change in free energy using parameters estimated from the complete data sets."}

| **Mutant** |  $\Delta\varepsilon_{RA}$ [$k_BT$] |
| :--: |:--:| 
| Y17I | $-9.81^{+0.04}_{-0.08}$|
| Q18A | $-10.60^{+0.07}_{-0.07}$|
| Q18M | $-15.61^{+0.05}_{-0.05}$|
  : Estimated DNA binding energies for each DNA binding domain mutant
  using all repressor copy numbers

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We also estimated the allosteric parameters ($K_A$, $K_I$, and
$\Delta\varepsilon_{AI}$) for all inducer binding domain mutations using
the induction profiles of all three operator sequences. The values,
reported in Table 7.6  are very similar to those estimated from a
single induction profile (Table 7.3). We note that for Q291R, it is
difficult to properly estimate the values for $K_A$ and $K_I$ as the
observed induction profile is approximately flat. The induction profiles
and predicted change in free energy for each inducer binding mutant is
shown in @Fig:global_IND_profiles. We see notable improvement in the
agreement between the induction profiles and the observed data,
indicating that considering all data significantly shrinks the
uncertainty of each parameter. The predicted change in free energy is
also improved compared to that shown in @Fig:IND_deltaF_comparison. We emphasize that the observed
free energy difference for each point assumes no knowledge of the
underlying parameters and comes directly from measurements. The
remarkable agreement between the predicted free energy and the
observations illustrates that redetermining the allosteric parameters is
sufficient to describe how the free energy changes as a result of the
mutation.

![**Induction profiles and predicted change in free energy using parameters
estimated from the complete data sets for inducer binding domain mutants.** Top
row shows fold-change measurements (points) as mean and standard error with
ten to fifteen biological replicates. Shaded lines correspond to the 95\%
credible regions of the induction profiles using the estimated values of the
allosteric parameters reported in Table 7.6. Bottom row shows the 95\%
credible regions of the predicted change in free energy (shaded lines) along
with the inferred free energy of data shown in the top row. In all plots, the
inducer concentration is shown on a symmetric log scale with linear scaling
between 0 and 10$^{-2}\, mu$M and log scaling elsewhere. The [Python code (`ch7_figS22.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS22.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS22){#fig:global_IND_profiles short-caption="Induction profiles and predicted change in free energy using parameters estimated from the complete data sets for inducer binding domain mutants."}

| **Mutant** |  $K_A$ \[$\mu$M\]  |    $K_I$ \[$\mu$M\]   |  $\Delta\varepsilon_{AI}$ \[$k_BT$\]  |
|:----------:| :-----------------:| :--------------------:| :-----------------------------------: |
|   F161T    |  $300_{-60}^{+60}$ |  $12.7_{-0.1}^{+0.1}$ |         $-0.9^{+0.3}_{-0.3}$ |
|   Q291K    |       $>1$ mM      |   $330_{-70}^{+60}$   |        $-3.17^{+0.07}_{-0.07}$ |
|   Q291R    |       $>1$ mM      |        $>1$ mM        |         $-2.4_{-0.2}^{+0.2}$ |
|   Q291V    |       $>1$ mM      |    $53^{+17}_{-13}$   |           $0_{-0.3}^{+0.3}$ |
                                                          
  : Estimated values for $K_A$, $K_I$, and $\Delta\varepsilon_{AI}$ for
  inducer binding domain mutations using induction profiles of all
  operator sequences.

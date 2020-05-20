## Applicability of Theory to the Oid Operator Sequence

In addition to the native operator sequences (O1, O2, and O3) considered in
the main text, we were also interested in testing our model predictions
against the synthetic Oid operator. In contrast to the other operators, Oid
is one base pair shorter in length (20 bp), is fully symmetric, and is known
to provide stronger repression than the native operator sequences considered
so far. While the theory should be similarly applicable, measuring the lower
fold-changes associated with this YFP construct was expected to be near the
sensitivity limit for our flow cytometer, due to the especially strong
binding energy of Oid ($\Delta \varepsilon_{RA}=-17.0\,k_BT$) [@garcia2011].
Accordingly, fluorescence data for Oid were obtained using microscopy, which
is more sensitive than flow cytometry.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We follow the approach of the main text and
make fold-change predictions based on the parameter estimates from our strain
with $R=260$ and an O2 operator. These predictions are shown in @Fig:oid_refit, where we
also plot data taken in triplicate for strains containing $R= 22$, 60, and 124,
obtained by single-cell microscopy. We find that the data are
systematically below the theoretical predictions. We also considered our
global fitting approach (see previous section) to see whether we
might find better agreement with the observed data. Interestingly, we
find that the majority of the parameters remain largely unchanged, but
our estimate for the Oid binding energy $\Delta \varepsilon_{RA}$ is
shifted to $-17.7~k_BT$ instead of the value $-17.0~k_BT$ found in @garcia2011.
In @Fig:oid_refit, we again plot the Oid fold-change data but with theoretical
predictions, using the new estimate for the Oid binding energy from our
global fit, and find substantially better agreement.

![**Predictions of fold-change for strains with an Oid binding sequence
versus experimental measurements with different repressor copy
numbers.** Experimental data is plotted against the parameter-free
predictions that are based on our fit to the O2 strain with $R=260$.
Here we use the previously measured binding energy
$\Delta\varepsilon_{RA}=-17.0\,k_BT$ [@garcia2011]. The same experimental data is
plotted against the best-fit parameters using the complete O1, O2, O3,
and Oid data sets to infer $K_A$, $K_I$, repressor copy numbers, and
the binding energies of all operators. Here the major
difference in the inferred parameters is a shift in the binding energy
for Oid from $\Delta\varepsilon_{RA}=-17.0~k_BT$ to
$\Delta\varepsilon_{RA}=-17.7~k_BT$, which now shows agreement between
the theoretical predictions and experimental data. Shaded regions from
the theoretical curves denote the 95\% credible
region. The [Python code (`ch6_figS14.py`)](https://github.com/gchure/phd/blob/master/src/chapter_06/code/ch6_figS14.py) used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch6_figS14){#fig:oid_refit short-caption="Predictions of fold-change
for strains with an Oid binding sequence versus experimental measurements with
different repressor copy numbers."}

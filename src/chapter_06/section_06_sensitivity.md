## Fold-Change Sensitivity Analysis
In we found that the width of the credible regions varied widely depending on
the repressor copy number $R$ and repressor operator binding energy $\Delta
\varepsilon_{RA}$. More precisely, the credible regions were much narrower
for low repressor copy numbers $R$ and weak binding energy
$\Delta\varepsilon_{RA}$. In this section, we explain how this behavior comes
about. We focus our attention on the maximum fold-change in the presence of
saturating inducer. While it is straightforward to consider the width of the
credible regions at any other inducer concentration, shows that the credible
region are widest at saturation.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The width of the credible regions corresponds to how sensitive the
fold-change is to the fit values of the dissociation constants $K_A$
and $K_I$. To be quantitative, we define
$$
    \Delta \text{fold-change}_{K_A} \equiv \text{fold-change}(K_A,K_I^\text{fit}) - \text{fold-change}(K_A^\text{fit},K_I^\text{fit}),
$${#eq:foldchange_difference_ka}
the difference between the fold-change at a particular $K_A$ value
relative to the best-fit dissociation constant
$K_A^\text{fit}=139 \, \mu\text{M}$. For simplicity, we keep the inactive state dissociation
constant fixed at its best-fit value $K_I^\text{fit}=0.53 \, \mu\text{M}$. A
larger difference $\Delta \text{fold-change}_{K_A}$ implies a wider credible region.
Similarly, we define the analogous quantity
$$
\Delta \text{fold-change}_{K_I} = \text{fold-change}(K_A^{\text{fit}},K_I) -
\text{fold-change}(K_A^{\text{fit}},K_I^{\text{fit}})
$${#eq:foldchange_difference_ki}
to measure the sensitivity of the fold-change to $K_I$ at a fixed
$K_A^{\text{fit}}$. @Fig:kaki_sensitivity shows both of these quantities in the limit
$c \to \infty$ for different repressor-DNA binding energies
$\Delta\varepsilon_{RA}$ and repressor copy numbers $R$. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To understand how the width of the credible region scales with
$\Delta\varepsilon_{RA}$ and $R$, we can Taylor expand the
difference in fold-change to first order,
$\Delta \text{fold-change}_{K_A} \approx (\partial
    \text{fold-change}/\partial K_A) \left( K_A - K_A^{\text{fit}} \right)$,
where the partial derivative has the form
$$
\frac{\partial \text{fold-change}}{\partial K_A} = \frac{e^{-\beta \Delta \varepsilon_{AI}} \frac{n}{K_I}\left(\frac{K_A}{K_I}\right)^{n-1}}{\left( 1 + e^{-\beta \Delta \varepsilon_{AI}} \left(\frac{K_A}{K_I}\right)^n \right)^2}\frac{R}{N_{NS}}e^{-\beta \Delta\varepsilon_{RA}} \left(
    1+\frac{1}{1+e^{-\beta \Delta \varepsilon_{AI} }
    \left(\frac{K_A}{K_I}\right)^n }\frac{R}{N_{NS}}e^{-\beta
    \Delta\varepsilon_{RA}} \right)^{-2}.
$${#eq:foldchange_partial_ka}
Similarly, the Taylor expansion
$\Delta \text{fold-change}_{K_I} \approx (\partial
    \text{fold-change}/\partial K_I) \left( K_I - K_I^{\text{fit}} \right)$
features the partial derivative 
$$
\frac{\partial \text{fold-change}}{\partial K_I} = -\frac{e^{-\beta \Delta \varepsilon_{AI}} \frac{n}{K_I}\left(\frac{K_A}{K_I}\right)^{n}}{\left( 1 + e^{-\beta \Delta \varepsilon_{AI}} \left(\frac{K_A}{K_I}\right)^n \right)^2}\frac{R}{N_{NS}}e^{-\beta \Delta\varepsilon_{RA}} \left(
1+\frac{1}{1+e^{-\beta \Delta \varepsilon_{AI} } \left(\frac{K_A}{K_I}\right)^n
}\frac{R}{N_{NS}}e^{-\beta \Delta\varepsilon_{RA}} \right)^{-2}.
$${#eq:foldchange_partial_ki}
From @Eq:foldchange_partial_ka and @Eq:foldchange_partial_ki, we find that both $\Delta \text{fold-change}_{K_A}$ and
$\Delta \text{fold-change}_{K_I}$ increase in magnitude with $R$ and
decrease in magnitude with $\Delta\varepsilon_{RA}$. Accordingly, we
expect that the O3 strains (with the least negative
$\Delta\varepsilon_{RA}$) and the strains with the smallest repressor
copy number will lead to partial derivatives with smaller magnitude and
hence to tighter credible regions. Indeed, this prediction is carried
out in @Fig:kaki_sensitivity.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Lastly, we note that @Eq:foldchange_partial_ka and @Eq:foldchange_partial_ki
enable us to quantify the scaling relationship between the width of the credible
region and the two quantities $R$
and $\Delta\varepsilon_{RA}$. For example, for the O3 strains, where
the fold-change at saturating inducer concentration is $\approx 1$,
the right-most term in both equations which equals the fold-change
squared is roughly 1. Therefore, we find that both $\frac{\partial
    \text{fold-change}}{\partial K_A}$ and
$\partial \text{fold-change}/\partial K_I$ scale linearly with $R$
and $e^{-\beta \Delta\varepsilon_{RA}}$. Thus the width of the
$R=22$ strain will be roughly 1/1000 as large as that of the
$R=1740$ strain; similarly, the width of the O3 curves will be roughly
1/1000 the width of the O1 curves.

![**Determining how sensitive the fold-change values are to the fit values of
the dissociation constants.** The difference $\Delta
\text{fold-change}_{K_A}$ in fold change when the dissociation constant $K_A$
is slightly offset from its best-fit value $K_A=139^{+29}_{-22} \mu\text{M}$,
as given by . Fold-change is computed in the limit of saturating inducer
concentration ($c \rightarrow \infty$, see ) where the credible regions in
are widest. The O3 strain ($\Delta\varepsilon_{RA} = -9.7~k_B T$) is about
1/1000 as sensitive as the O1 operator to perturbations in the parameter
values, and hence its credible region is roughly 1/1000 as wide. All curves
were made using $R = 260$. As in (A), but plotting the sensitivity of
fold-change to the $K_I$ parameter relative to the best-fit value
$K_I=0.53^{+0.04}_{-0.04} \mu\text{M}$. Note that only the magnitude, and not
the sign, of this difference describes the sensitivity of each parameter.
Hence, the O3 strain is again less sensitive than the O1 and O2 strains. As
in (A), but showing how the fold-change sensitivity for different
repressor copy numbers. The strains with lower repressor copy number are less
sensitive to changes in the dissociation constants, and hence their
corresponding curves in have tighter credible regions. All curves were made
using $\Delta\varepsilon_{RA} = -13.9~k_B T$. As in (C), the sensitivity
of fold-change with respect to $K_I$ is again smallest (in magnitude) for the
low repressor copy number strains.](ch6_figS11){#fig:kaki_sensitivity
short-caption="Sensitivity analysis of the fold-change function to $K_A$ and
$K_I$ estimates."}

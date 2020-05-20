## Properties of Induction Titration Curves

In this section, we expand on the phenotypic properties of the induction
response that were explored in Chapter 2. We begin by expanding on our
discussion of dynamic range and then show the analytic form of the
$[EC_{50}]$ for simple repression.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As stated in the main text, the dynamic range
is defined as the difference between the maximum and minimum system response,
or equivalently, as the difference between the saturation and leakiness of
the system. The dynamic range is therefore given by
$$
\begin{aligned}
\text{dynamic range} = &\left(
1+\frac{1}{1+e^{-\beta \Delta \varepsilon_{AI} } \left(\frac{K_A}{K_I}\right)^n
}\frac{R}{N_{NS}}e^{-\beta \Delta\varepsilon_{RA}} \right)^{-1} -\\
& \left(
1+\frac{1}{1+e^{-\beta \Delta \varepsilon_{AI} }}\frac{R}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}} \right)^{-1}.
\end{aligned}
$${#eq:induction_si_dyn_rng}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The dynamic range, along with saturation and
leakiness were plotted with our experimental data in Fig. 2.7(A-C) as a function
of repressor copy number. @Fig:properties_expanded shows how these properties
are expected to vary as a function of the repressor-operator binding energy.
Note that the resulting curves for all three properties have the same shape as
in Fig. 2.7 (A-C), since the dependence of the fold-change upon the repressor
copy number and repressor-operator binding energy are both contained in a single
multiplicative term, $R e^{-\beta \Delta\varepsilon_{RA}}$. Hence, increasing
$R$ on a logarithmic scale is equivalent to decreasing
$\Delta\varepsilon_{RA}$ on a linear scale.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;An interesting aspect of the dynamic range is
that it exhibits a peak as a function of either the repressor copy number (or
equivalently of the repressor-operator binding energy). Differentiating the
dynamic range (@Eq:induction_si_dyn_rng) and setting it equal to zero, we find that this peak occurs at
$$
\frac{R^*}{N_{NS}} = e^{-\beta (\Delta\varepsilon_{AI} -
\Delta\varepsilon_{RA})} \sqrt{e^{\Delta\varepsilon_{AI}} + 1}
\sqrt{e^{\Delta\varepsilon_{AI}} + \left( \frac{K_A}{K_I} \right)^n}.
$${#eq:loc_max_dyn_rng}
The magnitude of the peak is given by 
$$
\text{max dynamic range} = \frac{\left( \sqrt{e^{\Delta\varepsilon_{AI}} + 1} -
\sqrt{e^{\Delta\varepsilon_{AI}} + \left( \frac{K_A}{K_I} \right)^n}
\right)^2}{\left( \frac{K_A}{K_I} \right)^n - 1},
$${#eq:mag_max_dyn_rng}
which is independent of the repressor-operator binding energy
$\Delta\varepsilon_{RA}$ or $R$, and will only cause a shift in the
location of the peak, but not its magnitude.

![**Dependence of leakiness, saturation, and dynamic range on the operator
binding energy and repressor copy number.** Increasing repressor copy number or
decreasing the repressor-operator binding energy suppresses gene expression and
decreases both the (A) leakiness and (B) saturation. (C) The dynamic range
retains its shape, but shifts right as the repressor copy number increases. The
peak in the dynamic range can be understood by considering the two extremes for
$\Delta\varepsilon_{RA}$: for small repressor-operator binding energies, the
leakiness is small, but the saturation increases with $\Delta\varepsilon_{RA}$,
thereby decreasing the dynamic range. Repressor copy number does not affect the
maximum dynamic range. The [Python code (`ch6_figS18-S19.py`)](https://github.com/gchure/phd/blob/master/src/chapter_06/code/ch6_figS18-S19.py) used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch6_figS18){#fig:properties_expanded short-caption="Dependence of
leakiness, saturation, and dynamic range on the operator binding energy and
repressor copy number."}
s
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We now consider the two remaining properties, the
$[EC_{50}]$ and effective Hill coefficient, which determine the horizontal
properties of a system - that is, they determine the range of inducer
concentration in which the system’s response goes from its minimum to maximum
values. The $[EC_{50}]$ denotes the inducer concentration required to generate
fold-change halfway between its minimum and maximum value and was
defined implicitly in Eq. 2.9. For the simple repression system, the
$[EC_{50}]$ is given by
$$
\frac{[EC_{50}]}{K_A} = \frac{\frac{K_A}{K_I} - 1}{\frac{K_A}{K_I} - \left(
\frac{\left( 1 + \frac{R}{N_{NS}} e^{-\beta \Delta\varepsilon_{RA}} \right) +
\left( \frac{K_A}{K_I} \right)^n \left( 2 e^{-\beta \Delta \varepsilon_{AI}} +
\left( 1 + \frac{R}{N_{NS}} e^{-\beta \Delta\varepsilon_{RA}} \right) \right)
}{ 2 \left( 1 + \frac{R}{N_{NS}} e^{-\beta \Delta\varepsilon_{RA}} \right) +
e^{-\beta \Delta \varepsilon_{AI}} + \left( \frac{K_A}{K_I} \right)^n e^{-\beta
\Delta \varepsilon_{AI}} } \right)^{\frac{1}{n}}} - 1.
$${#eq:ec50_explicit}
Using this expression, we can then find the effective Hill coefficient $h$,
which equals twice the log-log slope of the normalized fold-change evaluated
at $c = [EC_{50}]$. In Fig. 2.7 (D-E), we show how these two properties vary
with repressor copy number, and in @Fig:ec50_ehill, we demonstrate
how they depend on the repressor-operator binding energy. Both the
$[EC_{50}]$ and $h$ vary significantly with repressor copy number for
sufficiently strong operator binding energies. Interestingly, for weak
operator binding energies on the order of the O3 operator, it is predicted
that the effective Hill coefficient should not vary with repressor copy
number. In addition, the maximum possible Hill coefficient is roughly 1.75,
which stresses the point that the effective Hill coefficient should not be
interpreted as the number of inducer binding sites, which is exactly 2.

![**[EC$_{50}$] and effective Hill coefficient depend strongly on repressor copy
number and operator binding energy.** (A) [EC$_{50}$] values from very small and
tightly clustered to relatively large and expanded for stronger operator binding
energies. (B) The effective Hill coefficient generally decreases with
increasing repressor copy number, indicating a flatter normalized response. The
maximum possible Hill coefficient is roughly 1.75 for all repressor-operator
binding energies. The [Python code (`ch6_figS18-S19.py`)](https://github.com/gchure/phd/blob/master/src/chapter_06/code/ch6_figS18-S19.py) used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch6_figS19){#fig:ec50_ehill short-caption="[EC$_{50}$] and
effective Hill coefficient depend strongly on repressor copy number and operator
binding energy."}

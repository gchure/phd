## Results

### Experimental Design

We test our model by predicting the induction profiles for an array of
strains that could be made using previously characterized repressor copy
numbers and DNA binding energies. Our approach contrasts with previous
studies that have parameterized induction curves of simple repression motifs,
as these have relied on expression systems where proteins are expressed from
plasmids, resulting in highly variable and unconstrained copy numbers
[@murphy2007; @murphy2010; @daber2009a; @daber2011a; @sochor2014]. Instead,
our approach relies on a foundation of previous work as depicted in Fig.
@fig:induction*types (C). This includes work from our laboratory that used
*E. coli* constructs based on components of the *lac* system to demonstrate
how the Lac repressor (LacI) copy number $R$ and operator binding energy
\$\Delta\varepsilon*{RA}\$ affect gene expression in the absence of inducer
[@garcia2011]. @rydenfelt2014b extended the theory used in that work to the
case of multiple promoters competing for a given transcription factor, which
was validated experimentally by @brewster2014, who modified this system to
consider expression from multiple-copy plasmids as well as the presence of
competing repressor binding sites.

The present study extends this body of work by introducing three additional
biophysical parameters – $\Delta\varepsilon_{AI}$, $K_A$, and $K_I$ – which
capture the allosteric nature of the transcription factor and complement the
results shown by @garcia2011 and @brewster2014. Although the current work
focuses on systems with a single site of repression, in the Materials and
Methods, we utilize data from @brewster2014, in which multiple sites of
repression are explored, to characterize the allosteric free energy
difference $\Delta\varepsilon_{AI}$ between the repressor’s active and
inactive states. This additional data set is critical because multiple
degenerate sets of parameters can characterize an induction curve equally
well, with the $\Delta\varepsilon_{AI}$ parameter compensated by the inducer
dissociation constants $K_A$ and $K_I$ (see Chapter 7, Supplementary
Information). After fixing $\Delta\varepsilon_{AI}$ as described in the
Materials and Methods, we can use data from single-site simple repression
systems to determine the values of $K_A$ and $K_I$.

We determine the values of $K_A$ and $K_I$ by fitting to a single
induction profile using Bayesian inferential methods [@sivia2006]. We then use
Eq. @eq:foldchange_full to
predict gene expression for any concentration of inducer, repressor copy
number, and DNA binding energy and compare these predictions against
experimental measurements. To obtain induction profiles for a set of
strains with varying repressor copy numbers, we used modified *lacI*
ribosomal binding sites from @garcia2011 to generate strains with mean repressor
copy number per cell of $R = 22 \pm 4$, $60 \pm 20$, $124 \pm 30$,
$260 \pm 40$, $1220 \pm 160$, and $1740 \pm 340$, where the error
denotes standard deviation of at least three replicates as measured by @garcia2011.
We note that $R$ refers to the number of repressor dimers in the cell,
which is twice the number of repressor tetramers reported by @garcia2011; since
both heads of the repressor are assumed to always be either specifically
or non-specifically bound to the genome, the two repressor dimers in
each LacI tetramer can be considered independently. Gene expression was
measured using a Yellow Fluorescent Protein (YFP) gene, driven by a
*lacUV5* promoter. Each of the six repressor copy number variants were
paired with the native O1, O2, or O3 *lac* operator [@oehler1994] placed at the YFP
transcription start site, thereby generating eighteen unique strains.
The repressor-operator binding energies (O1 \$\Delta\varepsilon*{RA} =
-15.3 \pm 0.2~k*BT$, O2 $\Delta\varepsilon*{RA} = -13.9~k*BT \pm 0.2$, and O3
$\Delta\varepsilon*{RA} = -9.7 \pm 0.1~k*BT$) were previously inferred
by measuring the fold-change of the *lac* system at different repressor
copy numbers, where the error arises from model fitting [@garcia2011]. Additionally,
we were able to obtain the value $\Delta \varepsilon*{AI} = 4.5\, k*BT$ by
fitting to previous data as discussed in the Materials and Methods. We
measure fold-change over a range of known IPTG concentrations $c$,
using $n=2$ inducer binding sites per LacI dimer and approximating the
number of non-specific binding sites as the length in base-pairs of the
*E. coli* genome, $N_{NS} = 4.6 \times 10^6$.

Our experimental pipeline for determining fold-change using flow
cytometry is shown in Fig. @fig:flowchart. Briefly, cells were grown to exponential phase,
in which gene expression reaches steady state [@scott2010], under concentrations of
the inducer IPTG ranging between 0 and $5000\, \mu$M. We measure YFP
fluorescence using flow cytometry and automatically gate the data to
include only single-cell measurements (see Materials and Methods). To validate
the use of flow cytometry, we also measured the fold-change of a subset
of strains using the established method of single-cell microscopy (see
Supplemental Information, Chapter 7). We found that
the fold-change measurements obtained from microscopy were
indistinguishable from that of flow-cytometry and yielded values for the
inducer binding constants $K_A$ and $K_I$ that were within error.

![An experimental pipeline for high-throughput fold-change
measurements. Cells are grown to exponential steady state and
their fluorescence is measured using flow cytometry. Automatic gating
methods using forward- and side-scattering are used to ensure that all
measurements come from single cells (see Methods). Mean expression is
then quantified at different IPTG concentrations (top, blue histograms)
and for a strain without repressor (bottom, green histograms), which
shows no response to IPTG as expected. Fold-change is computed by
dividing the mean fluorescence in the presence of repressor by the mean
fluorescence in the absence of repressor.](ch2_fig3){#fig:flowchart
short-caption="An experimental pipeline for high-throughput fold-change
measurements."}


## Determination of the *in vivo* MWC Parameters

The three parameters that we tune experimentally are shown in , leaving the
three allosteric parameters ($\Delta \varepsilon_{AI}$, $K_A$, and $K_I$) to
be determined by fitting. We used previous LacI fold-change data
[@brewster2014] to infer that $\Delta\varepsilon_{AI} = 4.5\, k_BT$ (see
Materials \& Methods). Rather than fitting $K_A$ and $K_I$ to our entire data
set of eighteen unique constructs, we performed Bayesian parameter estimation
on data from a single strain with $R=260$ and an O2 operator
[$\Delta\varepsilon_{RA}=-13.9\,k_BT$ [@garcia2011]] shown in Fig.
[@fig:induction_predictions](D, orange points). Using Markov Chain Monte
Carlo, we determine the most likely parameter values to be
$K_A=139^{+29}_{-22} \times 10^{-6} \, \text{M}$ and
$K_I=0.53^{+0.04}_{-0.04} \times 10^{-6}\, \text{M}$, which are the modes of
their respective distributions, where the superscripts and subscripts
represent the upper and lower bounds of the $95^\text{th}$ percentile of the
parameter value distributions [see Fig. [@fig:induction_predictions] (B)].
Unfortunately, we are not able to make a meaningful value-for-value
comparison of our parameters to those of earlier studies [@daber2009a;
@daber2011a] because of uncertainties in both gene copy number and
transcription factor copy numbers in these studies, as illustrated in
Supplemental Chapter 7. We then predicted the fold-change for the remaining
seventeen strains with no further fitting [see Fig.
[@fig:induction_predictions] (C - E)] together with the specific phenotypic
properties described in and discussed in detail below [see Fig.
[@fig:induciton_predictions] (F - J)]. The shaded regions in Fig.
[@fig:induction_predictions] (C - E) denote the 95% credible regions. Factors
determining the width of the credible regions are explored in Supplemental
Chapter 7.

We stress that the entire suite of predictions in is based upon the induction
profile of a single strain. Our ability to make such a broad range of
predictions stems from the fact that our parameters of interest such as the
repressor copy number and DNA binding energy - appear as distinct physical
parameters within our model. While the single data set in Fig.
@fig:induction_predictions could also be fit using a Hill function, such an
analysis would be unable to predict any of the other curves in the figure
(see Materials \& Methods). Phenomenological expressions such as the Hill
function can describe data, but lack predictive power and are thus unable to
build our intuition, help us design *de novo* input-output functions, or
guide future experiments [@kuhlman2007; @murphy2007].

![**Predicting induction profiles for different biological control
parameters.** We can quantitatively tune $R$ via ribosomal binding site (RBS)
modifications, $\Delta\varepsilon_{RA}$ by mutating the operator sequence,
and $c$ by adding different amounts of IPTG to the growth medium. Previous
experiments have characterized the $R$, $N_{NS}$, $\Delta\varepsilon_{RA}$,
and $\Delta\varepsilon_{AI}$ parameters (see ), leaving only the unknown
dissociation constants $K_A$ and $K_I$ between the inducer and the repressor
in the active and inactive states, respectively. These two parameters can be
inferred using Bayesian parameter estimation from a single induction curve.
Predicted IPTG titration curves for different repressor copy numbers and
operator strengths. Titration data for the O2 strain [white circles in Panel
] with $R=260$, $\Delta\varepsilon_{RA} = -13.9\, k_BT$, $n=2$, and
$\Delta\varepsilon_{AI}=4.5~k_BT$ can be used to determine the thermodynamic
parameters $K_A=139^{+29}_{-22} \times 10^{-6} , \text{M}$ and
$K_I=0.53^{+0.04}_{-0.04} \times 10^{-6}, \text{M}$ (orange line). The
remaining solid lines predict the fold-change for all other combinations of
repressor copy numbers (shown in the legend) and repressor-DNA binding
energies corresponding to the O1 operator ($-15.3\, k_B T$), O2 operator
($-13.9\, k_B T$), and O3 operator ($-9.7\, k_B T$). Error bars of
experimental data show the standard error of the mean (eight or more
replicates) when this error is not smaller than the diameter of the data
point. The shaded regions denote the 95% credible region, although the
credible region is obscured when it is thinner than the curve itself. To
display the measured fold-change in the absence of inducer, we alter the
scaling of the $x$-axis between $0$ and $10^{-7}$ M to linear rather than
logarithmic, as indicated by a dashed line. Additionally, our model allows us
to investigate key phenotypic properties of the induction profiles (see ).
Specifically, we show predictions for the leakiness, saturation, dynamic
range, $[EC_{50}]$, and effective Hill coefficient of the induction
profiles.](ch2_fig4.png){@fig:induction_predictions
short-caption="Predicting induction profiles for different biological control
parameters."}


## Comparison of Experimental Measurements with Theoretical Predictions

We tested the predictions shown in Fig. @fig:induction_predictions by
measuring fold-change induction profiles in strains with a broad range of
repressor copy numbers and repressor binding energies as characterized in
@garcia2011. With a few notable exceptions, the results shown in Fig.
@fig:induction_experiments demonstrate
agreement between theory and experiment. We note that there was an apparently
systematic shift in the O3 $\Delta\varepsilon_{RA} = -9.7\ k_BT$ strains [Fig.
@fig:induction_experiments (C)]  and all of the $R=1220$ and $R =1740$ strains. This may be partially due to
imprecise previous determinations of their $\Delta\varepsilon_{RA}$ and $R$
values. By performing a global fit where we infer all parameters including
the repressor copy number $R$ and the binding energy
$\Delta\varepsilon_{RA}$, we found better agreement for these strains,
although a discrepancy in the steepness of the response for all O3 strains
remains (see Materials \& Methods). We considered a number of hypotheses to explain
these discrepancies such as including other states (e.g. non-negligible
binding of the inactive repressor), relaxing the weak promoter approximation,
and accounting for variations in gene and repressor copy number throughout
the cell cycle, but none explained the observed discrepancies. As an
additional test of our model, we considered strains using the synthetic Oid
operator which exhibits an especially strong binding energy of
$\Delta\varepsilon_{RA}=-17\,k_B T$ [@garcia2011]. The global fit agrees well with the Oid
microscopy data, though it asserts a stronger Oid binding energy of
$\Delta\varepsilon_{RA}=-17.7\,k_B T$ (see Supplemental Chapter 7).

To ensure that the agreement between our predictions and data is not an
accident of the strain we used to perform our fitting, we also inferred
$K_A$ and $K_I$ from each of the other strains. As discussed in 
the Materials \& Methods and Fig. @fig:induction_predictions, the inferred values of $K_A$ and $K_I$ depend
minimally upon which strain is chosen, indicating that these parameter
values are highly robust. We also performed a global fit using the data
from all eighteen strains in which we fitted for the inducer
dissociation constants $K_A$ and $K_I$, the repressor copy number
$R$, and the repressor DNA binding energy $\Delta\varepsilon_{RA}$
(see Materials \& Methods). The resulting parameter values were nearly identical
to those fitted from any single strain. For the remainder of the text we
continue using parameters fitted from the strain with $R=260$
repressors and an O2 operator.

![**Comparison of predictions against measured and inferred data.** Flow
cytometry measurements of fold-change over a range of IPTG
concentrations for O1, O2, and O3 strains at varying repressor copy
numbers, overlaid on the predicted responses. Error bars for the
experimental data show the standard error of the mean (eight or more
replicates). As discussed in Fig. @fig:induction_predictions, all of the predicted induction curves
were generated prior to measurement by inferring the MWC parameters
using a single data set (O2 $R=260$, shown by white circles in Panel
B). The predictions may therefore depend upon which strain is used to
infer the parameters. The inferred parameter values of the dissociation
constants $K_A$ and $K_I$ using any of the eighteen strains instead
of the O2 $R=260$ strain. Nearly identical parameter values are
inferred from each strain, demonstrating that the same set of induction
profiles would have been predicted regardless of which strain was
chosen. The points show the mode, and the error bars denote the $95%$
credible region of the parameter value distribution. Error bars not
visible are smaller than the size of the
marker.](ch2_fig5){@fig:induction_experiments short-caption="Comparison of predictions against measured and inferred data."}

## Predicting the Phenotypic Traits of the Induction Response

A subset of the properties shown in Fig. @fig:inducible_types (i.e. the leakiness, saturation,
dynamic range, $[EC_{50}]$, and effective Hill coefficient) are of
significant interest to synthetic biology. For example, synthetic
biology is often focused on generating large responses (i.e. a large
dynamic range) or finding a strong binding partner (i.e. a small
$[EC_{50}]$) [@brophy2014; @shis2014]. While these properties are all individually
informative, when taken together they capture the essential features of
the induction response. We reiterate that a Hill function approach
cannot predict these features *a priori* and furthermore requires
fitting each curve individually. The MWC model, on the other hand,
enables us to quantify how each trait depends upon a single set of
physical parameters as shown by Fig. @fig:induction_predictions (F-J).

We define these five phenotypic traits using expressions derived from the
model presented in Eq. @eq:foldchange_full. These results build upon
extensive work by @martins2011, who computed many such properties for
ligand-receptor binding within the MWC model. We begin by analyzing the
leakiness, which is the minimum fold-change observed in the absence of
ligand, given by

$$
\begin{aligned}
\text{leakiness} &= \text{fold-change}(c=0) \nonumber\\
&= \left(
1+\frac{1}{1+e^{-\beta \Delta \varepsilon_{AI} }}\frac{R}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}} \right)^{-1},
\end{aligned}
$${#eq:leakiness}

and the saturation, which is the maximum fold change observed in the
presence of saturating ligand,

$$
\begin{aligned}
\text{saturation} &= \text{fold-change}(c \to \infty) \nonumber\\
&= \left(
1+\frac{1}{1+e^{-\beta \Delta \varepsilon_{AI} } \left(\frac{K_A}{K_I}\right)^n
}\frac{R}{N_{NS}}e^{-\beta \Delta\varepsilon_{RA}} \right)^{-1}.
\end{aligned}
$${#eq:saturation}

Systems that minimize leakiness repress strongly in the absence of
effector while systems that maximize saturation have high expression in
the presence of effector. Together, these two properties determine the
dynamic range of a system’s response, which is given by the difference

$$
\text{dynamic range} = \text{saturation} - \text{leakiness}
$${#eq:dynamic_range}

These
three properties are shown in Fig. @fig:induction_predictions (F-H). We discuss these properties in greater
detail in the Materials \& Methods section. Fig. @fig:properties_experiment shows that the measurements of these three
properties, derived from the fold-change data in the absence of IPTG and
the presence of saturating IPTG, closely match the predictions for all
three operators.

![**Predictions and experimental measurements of key properties of
induction profiles.** Data for the leakiness, saturation, and dynamic
range are obtained from fold-change measurements in in the absence of
IPTG and at saturating concentrations of IPTG. The three
repressor-operator binding energies in the legend correspond to the O1
operator ($-15.3~k_B T$), O2 operator ($-13.9~k_B T$), and O3
operator ($-9.7~k_B T$). Both the $[EC_{50}]$ and effective Hill
coefficient are inferred by individually fitting each operator-repressor
pairing in - separately to in order to smoothly interpolate between the
data points. Error bars for - represent the standard error of the mean
for eight or more replicates; error bars for - represent the 95%
credible region for the parameter found by propagating the credible
region of our estimates of $K_A$ and $K_I$ into
.](ch2_fig6){#eq:properties_experiment short-caption="Predictions and experimental measurements of key properties of
induction profiles."}

<span id="fig_properties_data" label="fig_properties_data">\[fig_properties_data\]</span>

Two additional properties of induction profiles are the $[EC_{50}]$
and effective Hill coefficient, which determine the range of inducer
concentration in which the system’s output goes from its minimum to
maximum value. The $[EC_{50}]$ denotes the inducer concentration
required to generate a system response halfway between its minimum and
maximum value, \[\label{ec50}
\foldchange(c = [EC_{50}]) = \frac{\text{leakiness} + \text{saturation}}{2}.\]
The effective Hill coefficient $h$, which quantifies the steepness of
the curve at the $[EC_{50}]$ , is given by \[\label{effectiveHill}
h = \left( 2 \frac{d}{d \log c} \left[ \log \left( \frac{ \foldchange(c) - \text{leakiness}}{\text{dynamic range}} \right) \right] \right)_{c = [EC_{50}]}.\]

- shows how the $[EC_{50}]$ and effective Hill coefficient depend on
  the repressor copy number. In , Section “”, we discuss the analytic
  forms of these two properties as well as their dependence on the
  repressor-DNA binding energy.

\- shows the estimated values of the $[EC_{50}]$ and the effective
Hill coefficient overlaid on the theoretical predictions. Both
properties were obtained by fitting to each individual titration curve
and computing the $[EC_{50}]$ and effective Hill coefficient using and
, respectively. We find that the predictions made with the single strain
fit closely match those made for each of the strains with O1 and O2
operators, but the predictions for the O3 operator are markedly off. In
the , Section “”, we show that the large, asymmetric error bars for the
O3 $R=22$ strain arise from its nearly flat response, where the lack
of dynamic range makes it impossible to determine the value of the
inducer dissociation constants $K_A$ and $K_I$, as can be seen in
the uncertainty of both the $[EC_{50}]$ and effective Hill
coefficient. Discrepancies between theory and data for O3 are improved,
but not fully resolved, by performing a global fit or fitting the MWC
model individually to each curve (see , Sections “” and “”). It remains
an open question how to account for discrepancies in O3, in particular
regarding the significant mismatch between the predicted and fitted
effective Hill coefficients.

## Data Collapse of Induction Profiles

Our primary interest heretofore was to determine the system response at
a specific inducer concentration, repressor copy number, and
repressor-DNA binding energy. However, the cell does not necessarily
“care about” the precise number of repressors in the system or the
binding energy of an individual operator. The relevant quantity for
cellular function is the fold-change enacted by the regulatory system.
This raises the question: given a specific value of the fold-change,
what combination of parameters will give rise to this desired response?
In other words, what trade-offs between the parameters of the system
will give rise to the same mean cellular output? These are key questions
both for understanding how the system is governed and for engineering
specific responses in a synthetic biology context. To address these
questions, we follow the data collapse strategy used in a number of
previous studies , and rewrite as a Fermi function,
\[\label{eq*free_energy_definition} \foldchange= \frac{1}{1+e^{-
F(c)}},\] where $F(c)$ is the free energy of the repressor binding to
the operator of interest relative to the unbound operator state in
$k_B T$ units , which is given by
\[\label{eq_free_energy_MWC_parameters}
F(c) = \frac{\Delta\varepsilon*{RA}}{k*BT} - \log \frac{\left(1+\frac{c}{K_A}\right)^n}{\left(1+\frac{c}{K_A}\right)^n+e^{-\beta \Delta\varepsilon*{AI} }\left(1+\frac{c}{K*I}\right)^n} - \log \frac{R}{N*{NS}}.\]
The first term in $F(c)$ denotes the repressor-operator binding
energy, the second the contribution from the inducer concentration, and
the last the effect of the repressor copy number. We note that
elsewhere, this free energy has been dubbed the Bohr parameter since
such families of curves are analogous to the shifts in hemoglobin
binding curves at different pHs known as the Bohr effect .

Instead of analyzing each induction curve individually, the free energy
provides a natural means to simultaneously characterize the diversity in
our eighteen induction profiles. demonstrates how the various induction
curves from - all collapse onto a single master curve, where points from
every induction profile that yield the same fold-change are mapped onto
the same free energy. shows this data collapse for the 216 data points
in -, demonstrating the close match between the theoretical predictions
and experimental measurements across all eighteen strains.

There are many different combinations of parameter values that can
result in the same free energy as defined in . For example, suppose a
system originally has a fold-change of 0.2 at a specific inducer
concentration, and then operator mutations increase the
$\Delta\varepsilon_{RA}$ binding energy . While this serves to
initially increase both the free energy and the fold-change, a
subsequent increase in the repressor copy number could bring the cell
back to the original fold-change level. Such trade-offs hint that there
need not be a single set of parameters that evoke a specific cellular
response, but rather that the cell explores a large but degenerate space
of parameters with multiple, equally valid paths.

![**Fold-change data from a broad collection of different strains
collapse onto a single master curve.** Any combination of parameters can
be mapped to a single physiological response (i.e. fold-change) via the
free energy, which encompasses the parametric details of the model.
Experimental data from collapse onto a single master curve as a function
of the free energy . The free energy for each strain was calculated from
using $n=2$, $\Delta\varepsilon_{AI}=4.5~k_BT$,
$K_A=139 \times 10^{-6} , \text{M}$,
$K_I=0.53 \times 10^{-6}, \text{M}$, and the strain-specific $R$
and $\Delta\varepsilon_{RA}$. All data points represent the mean, and
error bars are the standard error of the mean for eight or more
replicates.](main_figs/fig7.pdf)

<span id="fig_datacollapse" label="fig_datacollapse">\[fig_datacollapse\]</span>

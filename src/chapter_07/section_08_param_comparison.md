## Comparing Parameter Values To The Literature

In this section, we compare and contrast the biophysical parameter
values we use to characterize the wild-type Lac repressor with the rich
literature that consists of *in vitro* and *in vivo* experiments. This
section has an accompanying interactive figure available on the [paper
website](http://rpgroup.caltech.edu/mwc_mutants) which allows the reader
to examine different combinations of parameter values and their
agreement or disagreement with data taken from 
@garcia2011; @brewster2014 and @razo-mejia2018.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;While the mutations used in this work and those in @daber2011a
are the same, we report significantly different values for the inducer
binding, DNA binding parameters, and the relative energy difference
between active and inactive states of the mutant repressors. The
apparent disagreement of parameter values between the present work and
those presented in @daber2011a in part stem from different
treatments of the values for the wild-type Lac repressor. Since its
isolation by Gilbert and Müller-Hill in the 1960's @gilbert1966, the
Lac repressor has been the subject of intense biochemical and structural
study. Many measurements of the inducer and DNA binding kinetics of the
repressor *in vitro* [such as @ogorman1980] and their values
have informed the fitting of other parameters from measurements *in
vivo* (such as @daber2011a and  @daber2009a). All of these
measurements, however, do not *directly* measure the DNA- or
inducer-binding kinetics nor the equilibrium constant between the active
and inactive states of the repressor. To properly estimate the
parameters, one must either have direct measurement of a subset of the
parameters or make assumptions regarding the states of the system.
Examples of the estimated allosteric parameter values of the wild-type
LacI repressor from our previous work [@razo-mejia2018], that of @daber2011a, and *in vitro* measurements from @ogorman1980 are given in Table 7.3
The theoretical predictions for the
fold-change in gene expression, along with the values reported in
@daber2009a, can be seen using the interactive figure on the [paper
website](https://rpgroup.caltech.edu/mwc_mutants), where the reader can
also enter their own parameter values and independently assess the
agreement or lack thereof with the data.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It is notable that differences between the various references shown in
Table 7.3 can be drastic, in some cases differing by
almost an order of magnitude. Of particular note is the disagreement in
the energy difference between the active and inactive states of the
repressor, $\Delta\varepsilon_{AI}$. @daber2011a determines a
negative value of $\Delta\varepsilon_{AI}$ meaning that the inactive
state of the repressor is energetically favorable to the active state.
In stark contrast is the value reported in our previous work of
$+4.5\, k_BT$, implying that the active state is significantly more
stable than the inactive state. The now seminal *in vitro* measurements
reported in @ogorman1980 suggest that the two states are
nearly energetically equivalent.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The wide range of these reported values demonstrate that such
thermodynamic models are highly degenerate, meaning that many
combinations of parameter values can result in nearly equally good
descriptions of the data. To illustrate this point, we estimated $K_A$,
$K_I$, and DNA binding energy $\Delta\varepsilon_{RA}$ for each operator
to the data reported in 
@razo-mejia2018; @garcia2011; @brewster2014, using the three values of
$\Delta\varepsilon_{AI}$ shown in Table 7.4. The resulting fits can be seen in @Fig:lit_degeneracy.
Despite the drastically different values of $\Delta\varepsilon_{AI}$ it
is possible to generate adequate fits by modulating the other
parameters. The parameter values of these fits, reported in Table 7.3
further illustrate this point as they differ significantly from one another.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The theoretically predicted fold-change in the presence of multiple
promoters, shown in @Fig:lit_degeneracy (E) is perhaps the most informative test of
the parameter values. The theoretical advancements made in @weinert2014
provide a means to mathematically grasp the
intricacies of plasmid-borne expression through calculation of the
chemical potential. This formalism transforms the intimidating
combinatorics of $R$ repressors and $N$ specific binding sites (i.e.
plasmid reporter genes) into a two-state system where one can compute an
*effective* repressor copy number regulating a single promoter. Using
this formalism, the input-output function for the fold-change in gene
expression can be written as 
$$
\text{fold-change} = {1\over 1 +
    \lambda_{R}(c)e^{-\beta\Delta\varepsilon_{RA}}},
$${#eq:fc_fug}
where the effective repressor copy number (termed
the fugacity) is denoted as $\lambda_{R}(c)$. In Chapter 6 of this thesis, we show that the
inflection point of 
@Eq:fc_fug, whose curves shown in the bottom right-hand plot
of @Fig:lit_degeneracy, is located where the number of specific
binding sites for the repressor is approximately equal to the number of
repressors in the active state. Using this key feature, one can infer
$\Delta\varepsilon_{AI}$ given prior knowledge of
$\Delta\varepsilon_{RA}$ and the total number of repressors per cell. As
shown in @Fig:lit_degeneracy, using values of $\Delta\varepsilon_{AI}$
from @razo-mejia2018; @ogorman1980 approximately agree with the
measurements whereas the predicted curves using $\Delta\varepsilon_{AI}$
from @daber2011a overestimates the fold-change,
even though these values accurately describe the simple-repression data
shown in the other panels of @Fig:lit_degeneracy.

![**Degenerate fits of data from 
@razo-mejia2018; @brewster2014; @garcia2011 using different values for
active/inactive state energy difference $\Delta\varepsilon_{AI}$**. In
all plots, the solid, dashed, and dotted lines correspond to the
best-fit curves conditioned on the data using parameter values for
$\Delta\varepsilon_{AI}$ of $4.5\, k_BT$, $-1.75\,
k_BT$, and $0.35\, k_BT$, respectively. Induction profiles from @razo-mejia2018
for operators O1 (A), O2 (B), and O3 (C) are
shown as points and errors which correspond to the mean and standard
error of at least 10 biological replicates. (D) Leakiness measurements
of the simple repression motif with one unique regulated reporter gene.
Data shown are from @garcia2011; @brewster2014. (E) The
transcription factor titration effect. For gene expression measurements
on plasmids, the fold-change as a function of repressor copy number
exhibits strong nonlinearities. We used this effect as a way to
independently infer the parameter $\Delta\varepsilon_{AI}$ and it can be
seen that this breaks the degeneracy between different
parameters. The [Python code                                                
(`ch7_figS20.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS20.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).  ](ch7_figS20){#fig:lit_degeneracy short-caption="Degenerate fits of data using parameter values from the literature."} 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Without some direct *in vivo* measurements of
these parameters, one must make assumptions about the system to make any
quantitative progress. We chose to use the parameter values defined in our
laboratory as the repressor fugacity provides us with an independent, albeit
indirect, measurement of $\Delta\varepsilon_{AI}$ which other works such as
@daber2011a and @ogorman1980 do not. Both of these works determine all
of parameter values simultaneously by fitting to a single set of
measurements. While these measurements accurately describe their data, their
parameter values are less successful in accounting for data from @brewster2014; @garcia2011; @razo-mejia2018. In the context
of this work, we emphasize that we make many of the same qualitative
conclusions as in @daber2011a with respect to the effects of the
mutations using our particular set of parameter values.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;While we use different values for
$\Delta\varepsilon_{AI}$, the qualitative results between this work and that
of @daber2011a are in agreement. For example, both works determine that
mutations in the DNA binding domain alter only the DNA binding affinity
whereas the mutations in the inducer binding pocket affect only the
allosteric parameters. Because the biological variables such as repressor and
reporter gene copy number are tightly controlled in our system, we are able
to more precisely measure features of the induction profiles such as the
leakiness in gene expression. This ability allows us to detect changes in the
active/inactive equilibrium which were masked in @daber2011a
by the experimental design.While the precise parameter
values may be different between publications, the exploration of free energy
differences resulting from mutations are parameter-value independent and any
parameter disagreements do not change our theoretical model. Fig. 2 of the
main text is presented with no knowledge of parameter values -- it simply
shows the mathematics of the model. While the value of $\Delta F$ will
ultimately depend on the parameter values, the formalism of this work remains
agnostic to the parameter values and can be a useful tool for classifying
mutations and couple the sequence-level variation to the systems-level
response.

| **Parameter** | **Value** | **Reference**|
|:--:|:--:|:--:|
| $\Delta\varepsilon_{AI}$ | $\geq 4.5\, k_BT$ | @razo-mejia2018 |
| |$-1.7\, k_BT$ | @daber2011a |
| |$0.35\, k_BT$ | @ogorman1980 |
| $K_A$ |$139^{+22}_{-20}\, \mu$M  | @razo-mejia2018 |
| | $16\,\mu$M  | @daber2011a |
| |$133\,\mu$M  | @ogorman1980 |
| $K_I$ | $0.53^{+0.01}_{-0.01}\, \mu$M | @razo-mejia2018 |
| | $2\,\mu$M | @daber2011a |
| | $4\,\mu$M | @ogorman1980 |
: Thermodynamic parameter values of wild-type LacI from the
  literature.


|**Parameter** | **Value** | $\Delta\varepsilon_{AI}$ **Reference Value** |
|:--:|:--:|:--:| 
| $\Delta\varepsilon_{RA}$ (O1 operator) |  $-15.1^{+0.1}_{-0.1}\, k_BT$ |$4.5\, k_BT$ [@razo-mejia2018]|
|                                         |  $-17.1_{-0.1}^{+0.1}\, k_BT$ | $-1.75\, k_BT$ [@daber2011a]|
|                                         |  $-15.7_{-0.1}^{+0.1}\, k_BT$ | $0.35\, k_BT$ [@ogorman1980]|
| $\Delta\varepsilon_{RA}$ (O1 operator)  | $-15.1^{+0.1}_{-0.1}\, k_BT$ |$4.5\, k_BT$ [@razo-mejia2018]| 
|                                         |  $-17.1_{-0.1}^{+0.1}\, k_BT$| $-1.75\, k_BT$ [@daber2011a]|
|                                         |  $-15.7_{-0.1}^{+0.1}\, k_BT$| $0.35\, k_BT$ [@oforman1980]|
|  $\Delta\varepsilon_{RA}$ (O2 operator) |  $-13.4_{-0.1}^{+0.1}\, k_BT$   |   $4.5 k_BT$ [@razo-mejia2018]| 
|                                         |  $-15.4_{-0.1}^{+0.1}\, k_BT$   |   $-1.75\, k_BT$ [@daber2011a]|
|                                         |  $-14.0_{-0.1}^{+0.1}\, k_BT$   |   $0.35\, k_BT$ [@ogorman1980]|
|  $\Delta\varepsilon_{RA}$ (O3 operator) |  $-9.21^{+0.06}_{-0.06}\, k_BT$ |   $4.5 k_BT$ [@razo-mejia2018]|
|                                         |  $-11.29^{+0.06}_{-0.06}\, k_BT$|   $-1.75\, k_BT$ [@daber2011a]|
|                                         |  $-9.85^{+0.06}_{-0.05}\, k_BT$ |   $0.35\, k_BT$ [@ogorman1980]|
|  $K_A$                                  |  $225^{+10}_{-10}\, \mu$M       | $4.5\, k_BT$ [@razo-mejia2018]|
|                                         |  $290_{-20}^{+20}\, \mu$M       |   $-1.75\, k_BT$ [@daber2011a]|
|                                         |  $270_{-20}^{+20}\, \mu$M       |   $0.35\, k_BT$ [@ogorman1980]|
|  $K_I$                                  |  $0.81_{-0.05}^{+0.05}\, \mu$M  | $4.5\, k_BT$ [@razo-mejia2018]|
|                                         |  $8.2_{-0.5}^{+0.5}\, \mu$M     |   $-1.75\, k_BT$ [@daber2011a]|
|                                         |  $5.5_{-0.3}^{+0.5}\, \mu$M     |   $0.35\, k_BT$ [@ogorman1980]|

: Estimated parameters from global fits of data from literature.



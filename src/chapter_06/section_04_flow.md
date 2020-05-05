## Flow Cytometry

In this section, we provide information regarding the equipment used to
make experimental measurements of the fold-change in gene expression in
the interests of transparency and reproducibility. We also provide a
summary of our unsupervised method of gating the flow cytometry
measurements for consistency between experimental runs.

### Equipment
Due to past experience using the Miltenyi Biotec MACSQuant flow
cytometer during the Physiology summer course at the Marine Biological
Laboratory, we used the same flow cytometer for the formal measurements
in this work graciously provided by the Pamela Björkman lab at Caltech.
All measurements were made using an excitation wavelength of
$488\,\text{nm}$ with an emission filter set of 525/$50\,\text{nm}$.
This excitation wavelength provides approximately 40% of the maximum YFP
absorbance , and this was found to be sufficient for the purposes of
these experiments. A useful feature of modern flow cytometry is the
high-sensitivity signal detection through the use of photomultiplier
tubes (PMT) whose response can be tuned by adjusting the voltage. Thus,
the voltage for the forward-scatter (FSC), side-scatter (SSC), and gene
expression measurements were tuned manually to maximize the dynamic
range between autofluorescence signal and maximal expression without
losing the details of the population distribution. Once these voltages
were determined, they were used for all subsequent measurements.
Extremely low signal producing particles were discarded before data
storage by setting a basal voltage threshold, thus removing the majority
of spurious events. The various instrument settings for data collection
are given in Table 6.1.

|**Laser** | **Channel** | **Sensor Voltage**|
|:--|:--|:--:|
| 488 nm | Forward-Scatter (FSC) | 423 V |
| 488 nm | Side-Scatter (SSC) | 537 V |
| 488 nm | Intensity (B1 Filter, 525/50 nm) | 790 V|
| 488 nm | Trigger (debris threshold) | 24.5V|
Table: Instrument settings for data collectuion using the Miltenyi Biotec MACSQuant
flow cytometer.


### Unsupervised Gating
Flow cytometry data will frequently include a number of spurious events
or other undesirable data points such as cell doublets and debris. The
process of restricting the collected data set to those data determined
to be “real” is commonly referred to as gating. These gates are
typically drawn manually  and restrict the data set to those points
which display a high degree of linear correlation between their
forward-scatter (FSC) and side-scatter (SSC). The development of
unbiased and unsupervised methods of drawing these gates is an active
area of research [@aghaeepour2013].

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For this study, we used an automatic unsupervised gating
procedure to filter the flow cytometry data based on the front and
side-scattering values returned by the MACSQuant flow cytometer. We assume
that the region with highest density of points in these two channels
corresponds to single-cell measurements. Everything extending outside of this
region was discarded in order to exclude sources of error such as cell
clustering, particulates, or other spurious events.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In order to define the gated region we fit a two-dimensional Gaussian
function to the $\log_{10}$ forward-scattering (FSC) and the
$\log_{10}$ side-scattering (SSC) data. We then kept a fraction
$\alpha \in [0, 1]$ of the data by defining an elliptical region given
by 
$$
\left(\boldsymbol{x} - \boldsymbol{\mu} \right)^T \boldsymbol{\Sigma}^{-1}
\left(\boldsymbol{x} - \boldsymbol{\mu} \right) \leq \chi^2_\alpha(p),
$${#eq:gauss_gate}

where $\boldsymbol{x}$ is the $2 \times 1$ vector containing the
$\log(\text{FSC})$ and $\log(\text{SSC})$, $\boldsymbol{\mu}$ is
the $2 \times 1$ vector representing the mean values of $\log(\text{FSC})$ and
$\log(\text{SSC})$ as obtained from fitting a two-dimensional Gaussian
to the data, and $\boldsymbol{\Sigma}$ is the $2\times 2$ covariance
matrix also obtained from the Gaussian fit. $\chi^2_\alpha(p)$ is the
quantile function for probability $p$ of the chi-squared distribution
with two degrees of freedom. @Fig:gate_contours shows an example of different gating
contours that would arise from different values of $\alpha$ in @Eq:gauss_gate. In
this work, we chose $\alpha = 0.4$ which we deemed was a sufficient
constraint to minimize the noise in the data.  The specific code where this
gating is implemented can be found in [GitHub repository](https://github.com/RPGroup-PBoC/mwc_induction/blob/master/code/analysis/unsupervised_gating.ipynb).

![**Representative unsupervised gating contours of flow-cytometry data.** Points
indicate individual flow cytometry measurements of forward scatter and side
scatter. Colored contours indicate arbitrary gating contours ranging from 100\%
($\alpha = 1$) to 5\% ($\alpha=0.05$). All measurements shown in Chapters 2 and
3 in this work were made by computing the mean fluorescence from the
40$^\text{th}$ percentile ($\alpha = 0.4)$. The [Python code
(`ch6_figS8.py`)](https://github.com/gchure/phd/blob/master/src/chapter_06/code/ch6_figS8.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch6_figS8){#fig:gate_contours
short-caption="Representative unsupervised gating contours of flow-cytometry data."}

### Comparison of Flow Cytometry with Other Methods

Previous work from the Phillips' lab experimentally determined fold-change
for similar simple repression constructs using a variety of different
measurement methods [@garcia2011b]. Garcia and Phillips used the same
background strains as the ones used in this work, but gene expression was
measured with Miller assays based on colorimetric enzymatic reactions with
the LacZ protein [@garcia2011]. The experiments in @brewster2014 (as well as
in Chapter 4 of this dissertation) used a LacI dimer with the tetramerization
region replaced with an mCherry tag, where the fold-change was measured as
the ratio of the gene expression rate rather than a single snapshot of the
gene output.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;@Fig:new_old_comparison shows the comparison of these methods along
with the flow cytometry method used in this work. The consistency of these
three readouts validates the quantitative use of flow cytometry and
unsupervised gating to determine the fold-change in gene expression. However,
one important caveat revealed by this figure is that the sensitivity of flow
cytometer measurements is not sufficient to accurately determine the
fold-change for the high repressor copy number strains in O1 without
induction. Instead, a method with a large dynamic range such as the Miller
assay is needed to accurately resolve the fold-change at such low levels of
expression.

![**Comparison of experimental methods to determine the fold-change.** The
fold-change in gene expression ofr equivalent simple-repression constructs has
been determined using three independent methods: flow cytometry (Chapter 2),
colorimetric Miller assays (@garcia2011), and video microscopy (@brewster2014).
All three methods give consistent results, although flow cytometry meeasurements
lose accuracy for fold-change less than 0.01. Note that the repressor-DNA
binding energies $\Delta\varepsilon_{RA}$ used for the theoretical predictions
were determined in @garcia2011. The [Python code (`ch6_figS9.py`)](https://github.com/gchure/phd/blob/master/src/chapter_06/code/ch6_figS9.py) used to generate this figure can be found on the thesis [GitHub repository](https://github.com/gchure/phd). ](ch6_figS9){#fig:new_old_comparison
short-caption="Comparison of experimental methods to determine the fold-change
in gene expression."}

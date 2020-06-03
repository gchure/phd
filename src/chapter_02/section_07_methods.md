## Materials \& Methods

### Bacterial Strains and DNA Constructs

All strains used in these experiments were derived from *E. coli* K12
MG1655 with the *lac* operon removed, adapted from those created and
described in @garcia2011. Briefly, the operator variants and YFP reporter gene were
cloned into a pZS25 background which contains a *lacUV5* promoter that
drives expression as is shown schematically in  @Fig:states_weights. These constructs
carried a kanamycin resistance gene and were integrated into the *galK*
locus of the chromosome using $\lambda$ Red recombineering [@sharan2009]. The
*lacI* gene was constitutively expressed via a
P$_\text{LtetO-1}$ promoter [@lutz1997], with ribosomal binding site
mutations made to vary the LacI copy number as described in @salis2009a using
site-directed mutagenesis (Quickchange II; Stratagene), with further
details in @garcia2011. These *lacI* constructs carried a
chloramphenicol resistance gene and were integrated into the *ybcN*
locus of the chromosome. Final strain construction was achieved by
performing repeated P1 transduction [@thomason2007]  of the different operator and
*lacI* constructs to generate each combination used in this work.
Integration was confirmed by PCR amplification of the replaced
chromosomal region and by sequencing. Primers and final strain genotypes
are listed in supplemental Chapter 6.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It is important to note that the rest of the *lac* operon (*lacZYA*) was
never expressed. The LacY protein is a transmembrane protein which
actively transports lactose as well as IPTG into the cell. As LacY was
never produced in our strains, we assume that the extracellular and
intracellular IPTG concentration was approximately equal due to
diffusion across the membrane into the cell as is suggested by previous
work [@fernandez-castane2012].

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To make this theory applicable to transcription factors with any number
of DNA binding domains, we used a different definition for repressor
copy number than has been used previously. We define the LacI copy
number as the average number of repressor dimers per cell, whereas in @garcia2011,
the copy number is defined as the average number of repressor tetramers
in each cell. To motivate this decision, we consider the fact that the
LacI repressor molecule exists as a tetramer in *E. coli*  [@lewis1996] in which a
single DNA binding domain is formed from dimerization of LacI proteins,
so that wild-type LacI might be described as dimer of dimers. Since each
dimer is allosterically independent (i.e. either dimer can be
allosterically active or inactive, independent of the configuration of
the other dimer) [@daber2009a], a single LacI tetramer can be treated as two
functional repressors. Therefore, we have simply multiplied the number
of repressors reported in @garcia2011 by a factor of two. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A subset of strains in these experiments were measured using fluorescence
microscopy for validation of the flow cytometry data and results. To aid in
the high-fidelity segmentation of individual cells, the strains were modified
to constitutively express an mCherry fluorophore. This reporter was cloned
into a pZS4\*1 backbone [@lutz1997] in which mCherry is driven by the
*lacUV5* promoter. All microscopy and flow cytometry experiments were
performed using these strains.

### Growth Conditions for Flow Cytometry Measurements

All measurements were performed with *E. coli* cells grown to
mid-exponential phase in standard M9 minimal media (M9 5X Salts,
Sigma-Aldrich M6030; $2\,\text{mM}$ magnesium sulfate, Mallinckrodt
Chemicals 6066-04; $100\,\mu\text{M}$ calcium chloride, Fisher
Chemicals C79-500) supplemented with 0.5% (w/v) glucose. Briefly,
$500\,\mu\text{L}$ cultures of *E. coli* were inoculated into Lysogeny
Broth (LB Miller Powder, BD Medical) from a 50% glycerol frozen stock
(-80$^\circ$C) and were grown overnight in a $2\,\text{mL}$
96-deep-well plate sealed with a breathable nylon cover (Lab Pak - Nitex
Nylon, Sefar America Inc. Cat. No. 241205) with rapid agitation for
proper aeration. After approximately $12$ to $15$ hours, the
cultures had reached saturation and were diluted 1000-fold into a second
$2\,\text{mL}$ 96-deep-well plate where each well contained
$500\,\mu\text{L}$ of M9 minimal media supplemented with 0.5% w/v
glucose (anhydrous D-Glucose, Macron Chemicals) and the appropriate
concentration of IPTG (Isopropyl $\beta$-D-1 thiogalactopyranoside
Dioxane Free, Research Products International). These were sealed with a
breathable cover and were allowed to grow for approximately eight hours.
Cells were then diluted ten-fold into a round-bottom 96-well plate
(Corning Cat. No. 3365) containing $90\,\mu\text{L}$ of M9 minimal
media supplemented with 0.5% w/v glucose along with the corresponding
IPTG concentrations. For each IPTG concentration, a stock of 100-fold
concentrated IPTG in double distilled water was prepared and partitioned
into $100\,\mu\text{L}$ aliquots. The same parent stock was used for
all experiments described in this work.

### Flow Cytometry

All fold-change measurements were collected
on a Miltenyi Biotec MACSquant Analyzer 10 Flow Cytometer graciously
provided by the Pamela Björkman lab at Caltech. Detailed information
regarding the voltage settings of the photo-multiplier detectors can be
found in the supplemental Chapter 6. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Prior to each
day’s experiments, the analyzer was calibrated using MACSQuant
Calibration Beads (Cat. No. 130-093-607) such that day-to-day
experiments would be comparable. All YFP fluorescence measurements were
collected via $488\,\text{nm}$ laser excitation coupled with a
525/$50\,\text{nm}$ emission filter. Unless otherwise specified, all
measurements were taken over the course of two to three hours using
automated sampling from a 96-well plate kept at approximately
$4^\circ \, \hbox{-} \,
10^\circ$C on a MACS Chill 96 Rack (Cat. No. 130-094-459). Cells were
diluted to a final concentration of approximately $4\times 10^{4}$
cells per $\mu\text{L}$ which corresponded to a flow rate of
2,000-6,000 measurements per second, and acquisition for each well was
halted after 100,000 events were detected. Once completed, the data were
extracted and immediately processed using the following methods.

### Unsupervised Gating of Flow Cytometry Data

Flow cytometry data will frequently include a number of spurious events
or other undesirable data points such as cell doublets and debris. The
process of restricting the collected data set to those data determined
to be “real” is commonly referred to as gating. These gates are
typically drawn manually  and restrict the data set to those points
which display a high degree of linear correlation between their
forward-scatter (FSC) and side-scatter (SSC). The development of
unbiased and unsupervised methods of drawing these gates is an active
area of research [@lo2008; @aghaeepour2013]. For our purposes, we assume that the fluorescence
level of the population should be log-normally distributed about some
mean value. With this assumption in place, we developed a method that
allows us to restrict the data used to compute the mean fluorescence
intensity of the population to the smallest two-dimensional region of
the $\log(\mathrm{FSC})$ vs. $\log(\mathrm{SSC})$ space in which 40%
of the data is found. This was performed by fitting a bivariate Gaussian
distribution and restricting the data used for calculation to those that
reside within the 40th percentile. This procedure is described in more
detail in the supplemental Chapter 6.

### Experimental Determination of Fold-Change
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For each strain and IPTG concentration, the fold-change in gene
expression was calculated by taking the ratio of the population mean YFP
expression in the presence of LacI repressor to that of the population
mean in the absence of LacI repressor. However, the measured
fluorescence intensity of each cell also includes the autofluorescence
contributed by the weak excitation of the myriad protein and small
molecules within the cell. To correct for this background, we computed
the fold change as
$$
\text{fold-change} = \frac{\langle I_{R > 0} \rangle - \langle
I_\text{auto}\rangle}{\langle I_{R = 0} \rangle - \langle I_\text{auto}\rangle},
$${#eq:induction_image_def}
where $\langle I_{R > 0}\rangle$ is the average cell YFP intensity in
the presence of repressor, $\langle I_{R = 0}\rangle$ is the average
cell YFP intensity in the absence of repressor, and
$\langle I_\text{auto} \rangle$ is the average cell autofluorescence
intensity, as measured from cells that lack the *lac*-YFP construct.

### Bayesian Parameter Estimation
In this work, we determine the most likely parameter values for the
inducer dissociation constants $K_A$ and $K_I$ of the active and
inactive state, respectively, using Bayesian methods. We compute the
probability distribution of the value of each parameter given the data
$D$, which by Bayes’ theorem is given by 
$$P(K_A, K_I \,\vert\, D) = \frac{P(D \,\vert\, K_A, K_I)P(K_A, K_I)}{P(D)},
$${#eq:induction_bayes}
where $D$ is all the data composed of independent variables (repressor
copy number $R$, repressor-DNA binding energy
$\Delta\varepsilon_{RA}$, and inducer concentration $c$) and one
dependent variable (experimental fold-change). $P(D
\mid K_A, K_I)$ is the likelihood of having observed the data given the
parameter values for the dissociation constants, $P(K_A, K_I)$
contains all the prior information on these parameters, and $P(D)$
serves as a normalization constant, which we can ignore in our parameter
estimation. @Eq:fold_change_full assumes a deterministic relationship between the parameters
and the data, so in order to construct a probabilistic relationship as
required by @Eq:induction_bayes, we assume that the experimental fold-change for the
$i^\text{th}$ datum given the parameters is of the form
$$
\text{fold-change}_\text{exp}^{(i)} = \left( 1 + \frac{\left(1 +
\frac{c^{(i)}}{K_A}\right)^2}{\left( 1 + \frac{c^{(i)}}{K_A}\right)^2 +
e^{-\beta \Delta \varepsilon_{AI}} \left(1 + \frac{c^{(i)}}{K_I} \right)^2} \frac{R^{(i)}}{N_{NS}} e^{-\beta
\Delta \varepsilon_{RA}^{(i)}}\right)^{-1} + \epsilon^{(i)},
$${#eq:induction_fold_change_experimental}
where $\epsilon^{(i)}$ represents the
departure from the deterministic theoretical prediction for the
$i^\text{th}$ data point. If we assume that these $\epsilon^{(i)}$
errors are normally distributed with mean zero and standard deviation
$\sigma$, the likelihood of the data given the parameters is of the
form 
$$
\begin{aligned}
P(D \vert K_A, K_I, \sigma) &=
\frac{1}{(2\pi\sigma^2)^{\frac{n}{2}}}\times\\
&\prod\limits_{i=1}^n \exp 
\left[-\frac{(\text{fold-change}^{(i)}_\text{exp} - \text{fc}^\text{(theo)}(K_A, K_I, R^{(i)},
    \Delta\varepsilon_{RA}^{(i)}, c^{(i)}))^2}{2\sigma^2}\right],
\end{aligned}
$${#eq:likelihood}
where $\text{fold-change}^{(i)}_{\text{exp}}$ is the experimental fold-change
and $\text{fc}^\text{(theo)}(\,\cdots)$ is the theoretical prediction. The product
$\prod\limits_{i=1}^n$ captures the assumption that the $n$ data points are
independent. Note that the likelihood and prior terms now include the extra
unknown parameter $\sigma$. In applying  @Eq:likelihood, a choice of $K_A$
and $K_I$ that provides better agreement between theoretical fold-change
predictions and experimental measurements will result in a more probable
likelihood.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Both mathematically and numerically, it is convenient to define
$\tilde{k}_A = -\log \frac{K_A}{1\,\mu\text{M}}$ and
$\tilde{k}_I = -\log \frac{K_I}{1\,\mu\text{M}}$ and fit for these
parameters on a log scale. Dissociation constants are scale invariant,
so that a change from $10\,\mu\text{M}$ to $1\,\mu\text{M}$ leads to
an equivalent increase in affinity as a change from $1\,\mu\text{M}$
to $0.1\,\mu\text{M}$. With these definitions, we assume for the prior
$P(\tilde{k}_A, \tilde{k}_I, \sigma)$ that all three parameters are
independent. In addition, we assume a uniform distribution for
$\tilde{k}_A$ and $\tilde{k}_I$ and a Jeffreys prior  for the scale
parameter $\sigma$. This yields the complete prior
$$P(\tilde{k}_A, \tilde{k}_I, \sigma) \equiv \frac{1}{(\tilde{k}_A^{\max} -
\tilde{k}_A^{\min})} \frac{1}{(\tilde{k}_I^{\max} -
\tilde{k}_I^{\min})}\frac{1}{\sigma}.
$${#eq:induction_prior} 

These priors are maximally uninformative meaning that they imply no prior knowledge of the
parameter values. We defined the $\tilde{k}_A$ and $\tilde{k}_A$
ranges uniform on the range of $-7$ to $7$, although we note that
this particular choice does not affect the outcome provided the chosen
range is sufficiently wide.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Putting all these terms together we can now sample from $P(\tilde{k}_A,
\tilde{k}_I, \sigma \mid D)$ using Markov chain Monte Carlo  to compute the most
likely parameter as well as the error bars (given by the 95\% credible region)
for $K_A$ and $K_I$.

### Data Curation

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;All of the data used in this work as well as all relevant code can be
found at the [website](http://rpgroup-pboc.github.io/mwc_induction) associated
with the publication mentioned at the beginning of this chapter. Data were
collected, stored, and preserved using the Git version control software
in combination with off-site storage and hosting website GitHub. Code
used to generate all figures and complete all processing step as and
analyses are available on the GitHub repository. Many analysis files are
stored as instructive Jupyter Notebooks. The scientific community is
invited to fork our repositories and open constructive issues on the
[GitHub repository](https://www.github.com/rpgroup-pboc/mwc_induction).

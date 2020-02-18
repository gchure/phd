## Materials \& Methods

### Bacterial Strains and DNA Constructs

All wild-type strains from which the mutants were derived were generated in
previous work from the Phillips group [@garcia2011; @razo-mejia2018].
Briefly, mutations were first introduced into the *lacI* gene of our
pZS3\*1-lacI plasmid [@garcia2011] using a combination of overhang PCR
Gibson assembly as well as QuickChange mutagenesis (Agligent Technologies).
The oligonucleotide sequences used to generate each mutant as well as the
method are provided in the supplemental Chapter 7.

For mutants generated through overhang PCR and Gibson assembly,
oligonucleotide primers were purchased containing an overhang with the
desired mutation and used to amplify the entire plasmid. Using the homology
of the primer overhang, Gibson assembly was performed to circularize the DNA
prior to electroporation into MG1655 *E. coli* cells. Integration of LacI
mutants was performed with $\lambda$ Red recombineering as described in
@sharan2009 and @garcia2011.

The mutants studied in this work were chosen from data reported in
@daber2011a. In selecting mutations, we looked for mutants which suggested
moderate to strong deviations from the behavior of the wild-type repressor. We
note that the variant of LacI used in this work has an additional three amino
acids (Met-Val-Asn) added to the N-terminus than the canonical LacI sequence
reported in @farabaugh1978. To remain consistent with the field, we have
identified the mutations with respect to their positions in the canonical
sequence and those in @daber2011a. However, their positions in the raw data files correspond to
that of our LacI variant and is noted in the `README` files associated
with the data. 

### Flow Cytometry
All fold-change measurements were performed on a MACSQuant flow cytometer as
described in @razo-mejia2018. Briefly, saturated
overnight cultures 500 $\mu$L in volume were grown in deep-well 96 well plates 
covered with  a breathable nylon cover (Lab Pak - Nitex
Nylon, Sefar America, Cat. No. 241205). After approximately 12 to 15 hr, the
cultures reached saturation and were diluted 1000-fold into a second 2 mL
96-deep-well plate where each well contained 500 $\mu$L of M9 minimal media
supplemented with 0.5\% w/v glucose (anhydrous D-Glucose, Macron Chemicals)
and the appropriate concentration of IPTG (Isopropyl
$\beta$-D-1-thiogalactopyranoside, Dioxane Free, Research Products International).
These were sealed with a breathable cover and were allowed to grow for
approximately 8 hours until the OD$_\text{600nm} \approx 0.3$. 
Cells were then diluted ten-fold into a round-bottom 96-well plate 
(Corning Cat. No. 3365) containing 90 $\mu$L of M9 minimal media supplemented 
with 0.5\% w/v glucose along with the corresponding IPTG concentrations.

The flow cytometer was calibrated prior to use with MACSQuant Calibration
Beads (Cat. No. 130-093-607). During measurement, the cultures were held at
approximately 4$^\circ$ C by placing the 96-well plate on a MACSQuant ice block.
All fluorescence measurements were made using a 488 nm excitation wavelength
with a 525/50 nm emission filter. The photomultiplier tube voltage settings for
the instrument are the same as those used in @razo-mejia2018 and are listed in
supplemental Chapter 6.

The data was processed using an automatic unsupervised gating procedure based
on the front and side-scattering values, where we fit a two-dimensional
Gaussian function to the $\log_{10}$ forward-scattering (FSC) and the
$\log_{10}$ side-scattering (SSC) data. Here we assume that the region with
highest density of points in these two channels corresponds to single-cell
measurements and consider data points that fall within 40\% of the highest
density region of the two-dimensional Gaussian function. We direct the reader
to Reference [@razo-mejia2018] for further detail and comparison of flow
cytometry with single-cell microscopy.

### Bayesian Parameter Estimation
We used a Bayesian definition of probability in the statistical analysis of
all mutants in this work. In the SI text, we derive
in detail the statistical models used for the various parameters as well as
multiple diagnostic tests. Here, we give a generic description of our
approach. To be succinct in notation, we consider a generic parameter
$\theta$ which represents $\Delta\varepsilon_{RA}$, $K_A$, $K_I$, and/or
$\Delta\varepsilon_{AI}$ depending on the specific LacI mutant.

As prescribed by Bayes' theorem, we are interested in the posterior probability
distribution

$$
g(\theta\,\vert\, y) \propto {f(y\,\vert\,\theta)g(\theta)},
$${#eq:bayes_generic}

where we use $g$ and $f$ to represent probability densities over parameters and
data, respectively, and $y$ to represent a set of fold-change measurements. The
likelihood of observing our dataset $y$ given a value of $\theta$ is captured by
$f(y\,\vert\,\theta)$. All prior information we have about the possible values
of $\theta$ are described by $g(\theta)$. 

In all inferential models used in this work, we assumed that all experimental measurements at a
given inducer concentration were normally distributed about a mean value $\mu$
dictated by Eq. @eq:foldchange with a variance $\sigma^2$, 

$$
f(y\,\vert\, \theta) = {1 \over (2\pi\sigma^2)^{N/2}}\prod\limits_i^N \exp\left[-{(y_i - \mu(\theta))^2 \over 2\sigma^2}\right], 
$${#eq:generic_likelihood}

where $N$ is the number of measurements in the data set $y$.

This choice of likelihood is justified as each individual measurement at a given
inducer concentration is a biological replicate and independent of all other
experiments. By using a Gaussian likelihood, we introduce another parameter
$\sigma$. As $\sigma$ must be positive and greater than zero, we define as a prior
distribution  a half-normal distribution with a standard deviation $\phi$, 

$$
g(\sigma) = { {1 \over \phi}\sqrt{2 \over \pi} }\exp\left[-{x \over 2\phi^2}\right]\,;\, x \geq 0,
$${#eq:sigma_prior}

where $x$ is a given range of values for $\sigma$. A standard deviation of
$\phi=0.1$ was chosen given our knowledge of the scale of our measurement
error from other experiments. As the absolute measurement of fold-change is
restricted between $0$ and $1.0$, and given our knowledge of the sensitivity
of the experiment, it is reasonable to assume that the error will be closer
to $0$ than to $1.0$. Further justification of this choice of prior through
simulation based methods are given in the supplemental Chapter 7. The prior distribution for
$\theta$ is dependent on the parameter and its associated physical and
physiological restrictions. Detailed discussion of our chosen prior
distributions for each model can also be found in the supplemental Chapter 7.

All statistical modeling and parameter inference was performed using Markov
chain Monte Carlo (MCMC). Specifically, Hamiltonian Monte Carlo sampling was
used as is implemented in the Stan probabilistic programming language [@carpenter2017]. All
statistical models saved as `.stan` models and can be accessed at the
[GitHub repository](https://www.github.com/rpgroup-pboc/mwc_mutants)
associated with this work (DOI: 10.5281/zenodo.2721798) or can be downloaded directly from the
[paper website](https://www.rpgroup.caltech.edu/mwc_mutants). 

### Inference of Free Energy From Fold-Change Data

While the fold-change in gene expression is restricted to be between 0 and 1, experimental
noise can generate fold-change measurements beyond these bounds. To determine
the free energy for a given set of fold-change measurements (for one unique
strain at a single inducer concentration), we modeled the observed fold-change
measurements as being drawn from a normal distribution with a mean $\mu$
and standard deviation $\sigma$. Using Bayes' theorem, we can write
the posterior distribution as 
$$
g(\mu, \sigma\,\vert y) \propto
g(\mu)g(\sigma){1\over(2\pi\sigma^2)^{N/2}}\prod\limits_i^N
\exp\left[{-(y_i - \mu)^2 \over 2\sigma^2}\right]
$${#eq:post_fc_mu}

where $y$ is a collection of fold-change measurements. The prior distribution
for $\mu$ was chosen to be uniform between 0 and 1 while the prior on $\sigma$
was chosen to be half normal, as written in Eq. @{eq:sigma_prior. The
posterior distribution was sampled independently for each set of fold-change
measurements using MCMC. The `.stan` model for this inference is available on the
[paper website](http://www.rpgroup.caltech.edu/mwc_mutants). 

For each MCMC sample of $\mu$, the free energy was calculated as 
$$
F = -\log\left(\mu^{-1} - 1\right)
$${#eq:empirical_F}

which is simply the rearrangement of Eq. @eq:collapse. Using simulated
data, we determined that when $\mu < \sigma$ or $(1 - \mu) < \sigma$, the
mean fold-change in gene expression was over or underestimated for the lower
and upper limit, respectively. This means that there are maximum and minimum
levels of fold-change that can be detected using flow cytometry which are set
by the distribution of fold-change measurements resulting from various
sources of day-to-day variation. This results in a systematic error in the
calculation of the free energy, making proper inference beyond these limits
difficult. This bounds the range in which we can confidently infer this
quantity with flow cytometry. We hypothesize that more sensitive methods,
such as single cell microscopy, colorimetric assays, or direct counting of
mRNA transcripts via Fluorescence *In Situ* Hybridization (FISH) would
improve the measurement of $\Delta F$. We further discuss details of this
limitation in the supplemental Chapter 7.

### Data and Code Availability

All data was collected, stored, and preserved using the Git version control
software. Code for data processing, analysis, and figure generation is
available on the GitHub repository
([https://www.github.com/rpgroup-pboc/mwc_mutants}{https://www.github.com/rpgroup-pboc/mwc\_mutants](https://www.github.com/rpgroup-pboc/mwc_mutants}{https://www.github.com/rpgroup-pboc/mwc\_mutants))
or can be accessed via the
[paper website](http://www.rpgroup.caltech.edu/mwc_mutants). Raw flow
cytometry data is stored on the CaltechDATA data repository and can be
accessed via DOI 10.22002/D1.1241.


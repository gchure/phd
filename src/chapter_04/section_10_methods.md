## Materials & Methods


### Bacterial Strains and Growth Media


Three genotypes were primarily used in this work, all in the genetic
background of *Escherichia coli* MG1655-K12 and all derived from those
used in [@brewster2014. The three genotypes are as follows. For
each experiment, an autofluorescence control was used which contained no
fluorescent reporters [except for a CFP volume marker used for
segmentation in @brewster2014] which had the *lacI* and *lacZYA*
genes deleted from the chromosome. The constitutive expression strain
($\Delta lacI; \Delta lacZYA$) included a YFP reporter gene integrated
into the *galK* locus of the chromosome along with a kanamycin
resistance cassette. The experimental strains in which LacI expression
was controlled contained a *lacI-mCherry* fluorescent fusion integrated
into the *ybcN* locus of the chromosome along with a chloramphenicol
resistance cassette. This cassette was later deleted from the chromosome
using FLP/FRT recombination [@schlake1994; @zhu1995]. The strain was
then transformed with plasmid [pZS3-pN25-tetR following notation
described in @lutz1997] constitutively expressing the TetR
repressor along with a chloramphenicol resistance cassette. All
bacterial strains and plasmids used in this work are reported in the SI
Text.

### Bacterial Growth Curves

Bacterial growth curves were measured in a multi-well plate reader
(BioTek Cytation5) generously provided by the David Van Valen lab at
Caltech. Cells constitutively expressing YFP were grown overnight to
saturation in LB broth (BD Medical) at 37$^\circ$ C with aeration. Once
saturated, cells were diluted 1000 fold into 50 mL of the desired growth
medium and were allowed to grow at the appropriate experimental
temperature with aeration for several hours until an
OD$_{600nm} \approx 0.3$ was reached. Cells were then diluted 1:10 into
the desired growth media at the desired temperature. After being
thoroughly mixed, 500 $\mu$L aliquots were transferred to a black-walled
96-well plate (Brooks Automation Incorporated, Cat No. MGB096-1-2-LG-L),
leaving two rows and two columns of wells on each side of the plate
filled with sterile growth medium to serve as blanks and buffer against
temperature variation. The plate was then transferred to the pre-warmed
plate reader. OD$_{600nm}$ measurements were made every five minutes for
12 to 24 hours until cultures had saturated. In between measurements,
the plate incubated at the appropriate temperature with a linear shaking
mode. We found that double-orbital shaking modes led to the formation of
cell aggregates which gave inconsistent measurements.

### Estimation of Bacterial Growth Rate

Non-parametric estimation of the maximum growth rate was performed using
the FitDeriv Python package as described in @swain2016. Using
this approach, the bacterial growth curve is modeled as a Gaussian
process in which the measured growth at a given time point is modeled as
a Gaussian distribution whose mean is dependent on the mean of the
neighboring time points. This allows for smooth interpolation between
adjacent measurements and calculation of second derivatives without an
underlying parametric model. The reported growth rates are the maximum
value inferred from the exponential phase of the experimental growth
curve.

### Growth Conditions

Parent strains (autofluorescence control, $\Delta lacI$ constitutive
control, and the ATC-inducible LacI-mCherry strain) were grown in LB
Miller broth (B.D. Medical, Cat. No. 244620 ) at 37$^\circ$ C with
vigorous aeration until saturated. Cells were then diluted between 1000
and 5000 fold into 3 mL of M9 minimal medium (B.D. Medical, Cat. No.
248510). The bacterial strain expressing the tetracycline-inducible
LacI-mCherry was diluted into six 3 mL cultures with differing
concentrations of ATC (Chemodex, Cat. No. CDX-A0197-T78) ranging from
0.1 to 10 ng / mL to induce expression of the transcription factor.
These concentrations were reached by dilution from 1 $\mu$g / mL stock
in 50\% ethanol. All cultures were shielded from ambient light using
either aluminum foil or via an enclosure and were grown at the
appropriate experimental temperature with aeration until an OD$_{600nm}$
of approximately $0.25 - 0.35$. Due to pipetting errors, cultures
reached OD$_{600nm} \approx 0.3$ at slightly different points in time.
To ensure that strains could be directly compared, all strains were back
diluted by several fold until equivalent. When all samples reached the
appropriate OD$_{600nm}$, the cells were harvested for imaging.

### Imaging Sample Preparation

Prior to the preparation of cell cultures for imaging, a 2\% (w/v)
agarose substrate (UltraPure, Thermo Scientific) was prepared and
allowed to reach room temperature. For experiments conducted at
42$^\circ C$, 4\% (w/v) agarose substrates were prepared. Briefly, the
agarose was mixed with the appropriate growth medium in a 50 mL conical
polystyrene tube and then microwaved until molten. A 300 to 500 $\mu$L
aliquot was then sandwiched between two glass coverslips to ensure a
flat imaging surface. Once solidified, the agarose pads were cut into
squares approximately $0.5$ cm per side.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To determine the calibration factor between fluorescence and protein
copy number, production the fluorophore in question must be halted such
that all differences in intensity between two daughter cells result from
binomial partitioning of the fluorophores and not from continuing
expression. This was achieved by removing the anhydrous tetracycline
inducer from the growth medium through several washing steps as outlined
in @brewster2014 ]. Aliquots of 100 $\mu$L from each ATC-induced
culture were combined and pelleted at 13000$\times g$ for 2 minutes. The
supernatant (containing ATC) was aspirated and replaced with 1 mL of M9
growth medium without ATC. The pellet was resuspended and pelleted at
13000$\times g$. This wash step was repeated three times to ensure
residual ATC had been removed and LacI-mCherry production was halted.
After the final wash, the cell pellet was resuspended in 1 mL of M9
medium and diluted ten-fold. A 1 $\mu$L aliquot of the diluted mixture
was then spotted onto an agarose substrate containing the appropriate
carbon source.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The remaining bacterial cultures
(autofluorescence control, constitutive expression control, and the
ATC-induced samples) were diluted ten-fold into sterile M9 medium.This
dilution was thoroughly mixed and 1 $\mu$L aliquots were spotted onto agarose
substrates lacking the carbon source.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Once the spotted cells had dried onto the
agarose substrates (about 5 to 10 minutes after deposition), the agarose pads
were inverted and pressed onto a glass bottom dish (Electron Microscopy
Sciences, Cat. No. 70674-52) and sealed with parafilm. Strips of double stick
tape were added to the edge of the dish to help immobilize the sample on the
microscope stage and minimize drift.

### Microscopy

All imaging was performed on a Nikon Ti-Eclipse inverted microscope
outfitted with a SOLA LED fluorescence illumination system. All images
were acquired on a Andor Zyla 5.5 sCMOS camera (Oxford Instruments
Group). The microscope body and stage was enclosed in a plexiglass
incubation chamber (Haison, approximately 1$^\circ$ C regulation
control) connected to an external heater. Temperature of the stage was
measured via a thermometer which controlled heating of the system.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;All static images (i.e. images from which
fold-change and repressor counts were calculated) were measured in an
identical manner. Ten to fifteen fields of view containing on average 25
cells were imaged using phase contrast and fluorescence excitation.
Fluorescence exposures were each 5 seconds while the phase contrast exposure
time was between 75 ms and 150 ms. This procedure was repeated for each
unique strain and ATC induction concentration.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To compute the calibration factor for that day
of imaging, time-lapse images were taken of a separate agarose pad covered in
cells containing various levels of LacI-mCherry. Fifteen to twenty positions
were marked, choosing fields of view containing 20 to 50 cells. Cells were
allowed to grow for a period of 90 to 120 minutes (depending on the
medium-dependent growth rate) with phase contrast images taken every 5 to 10
minutes. At the final time-point, both phase contrast and fluorescence images
were acquired using the same settings for the snapshots. Once the experiment
was completed, images were exported to `.tif` format and transferred to cold
storage and a computational cluster for analysis.

### Lineage Tracking

Cells were segmented and lineages reconstructed using the
deep-learning-based bacterial segmentation software SuperSegger `v1.1`
[@stylianidou2016a; @cass2017] operated through MATLAB (Mathworks,
`v2017b`). We found that the pre-trained network constants `100XEcM9`
(packaged with the SuperSegger software) worked well for all growth
conditions tested in this work. The generated files (`clist.mat`)
associated with each sample and position were parsed using bespoke
Python scripts to reconstruct lineages and apply appropriate filtering
steps before calculating the fluorescence calibration factor. We direct
the reader to the SI text for details of our data validation procedure
to ensure proper lineage tracking.

### Calculation of the Calibration Factor

To estimate the calibration factor $\alpha$, we used a Bayesian
definition of probability to define a posterior distribution of $\alpha$
conditioned on intensity measurements of sibling cells after division.
We direct the reader to the SI text for a detailed discussion of the
inferential procedures and estimation of systematic error. We give a
brief description of the inference below.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We are interested in determining the
fluorescence of a single LacI-mCherry repressor dimer given a set of
intensity measurements of sibling cells, $[I_1, I_2]$. The intensity of a
given cell $I$ is related to the number of LacI-mCherry dimers it is
expressing by a multiplicative factor $\alpha$ which can be enumerated
mathematically as 
$$
I = \alpha N,
$${#eq:growth_ian}
where $N$ is the total number of LacI-mCherry dimers. We can define the
posterior probability distribution of $\alpha$ conditioned on the intensity
measurements using
Bayes' theorem as
$$
g(\alpha\,\vert\, [I_1, I_2]) = \frac{f([I_1, I_2]\,\vert\,
\alpha)g(\alpha)}{f([I_1, I_2])}.
$${#eq:growth_bayes_generic}
where we have used $g$ and $f$ as probability
densities over parameters and data, respectively. The denominator of
this expression (the evidence) is equivalent to the first term of the
numerator (the likelihood) marginalized over $\alpha$. In this work,
this term serves as normalization factor and can be neglected.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Assuming that no more LacI-mCherry dimers are produced during cell
division, the number of repressors that each sibling cell receives after
division of the parent cell is binomially distributed with a probability
$p$. We can make the approximation that partitioning of the repressors
is even such that $p =
1/2$. The validity of this approximation is discussed in detail in the
SI text. Using @Eq:growth_ian and the change-of-variables formula, we can define
the likelihood $g([I_1, I_2]\,\vert\, \alpha)$ as
$$
g([I_1, I_2]\,\vert\,\alpha) =
\frac{1}{\alpha^k}\prod_i^k\frac{\Gamma\left(\frac{ {I_1}_i + {I_2}_i}{\alpha} +
1\right)}{\Gamma\left(\frac{ {I_1}_i}{\alpha} +
1\right)\Gamma\left(\frac{ {I_2}_i}{\alpha} + 1\right)}2^{-\frac{ {I_1}_i +
{I_2}_i}{\alpha}},
$${#eq:growth_binom_like}
where $k$ is the total number of sibling pairs measured.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With a likelihood defined, we must now define a functional form for
$g(\alpha)$ which describes all prior information known about the
calibration factor knowing nothing about the actual measurements.
Knowing that we design the experiments such that only $\approx 2/3$ of
the dynamic range of the camera is used and $\alpha$ cannot be less than
or equal to zero, we can define a half-normal distribution with a
standard deviation of $\sigma$ as
$$
g(\alpha) = \sqrt{2 \over \pi\sigma^2}\exp\left[{-\alpha^2 \over
2\sigma^2}\right]\,;\, \forall \alpha > 0.
$${#eq:alpha_prior}
where the standard deviation is large, for
example, $\sigma = 500$ a.u. / fluorophore. We evaluated the posterior
distribution using Markov chain Monte Carlo (MCMC) as is implemented in
the Stan probabilistic programming language [@carpenter2017]. The
`.stan` file associated with this model along with the Python code used
to execute it can be accessed on the [paper
website](https://rpgroup.caltech.edu/mwc_growth).

### Counting Repressors 
Given an estimation for $\alpha$ for each experiment, we calculate the
total number of repressors per cell from a
$$
R = \frac{I_{mCherry}}{\alpha}.
$${#eq:nia}
However, as discussed in detail in Chapter 8, a
systematic error in the repressor count is introduced due to division in
the asynchronous culture between the cessation of LacI-mCherry
production and the actual imaging. The entire sample preparation
procedure is $\approx 30 - 60$ min, during which time some cells
complete a cell division, thereby diluting the total repressor count. To
ensure that the measured number of repressors corresponded to the
measured YFP intensity, we restricted our dataset for all experiments to
cells that had a pole-to-pole length $\ell \geq 3.5\,\mu$m, indicating
that they had likely not undergone a division during the sample
preparation.

### Code and Data Availability {#code-and-data-availability .unnumbered}

All code used in this work is available on the [paper
website](https://github.com/rpgroup-pboc/mwc_growth) and [associated
GitHub repository](https://www.github.com/rpgroup-pboc/mwc_growth)(DOI:
0.5281/zenodo.3560369). This work also used the open-source software
tools
[SuperSegger](http://mtshasta.phys.washington.edu/website/SuperSegger.php)
`v.1.1`[@stylianidou2016a; @cass2017] for lineage tracking and
[FitDeriv](http://swainlab.bio.ed.ac.uk/software/fitderiv/) `v.1.03`
[@swain2016] for the nonparametric estimation of growth rates. Raw image
files are large (1.8 Tb) and are therefore available upon request. The
clist.mat files used to calculate fold-change and to assign sibling
cells can be accessed via the associated [GitHub
repository](https://github.com/rpgroup-pboc/mwc_growth) via ([DOI:
0.5281/zenodo.3560369](http://doi.org/10.5281/zenodo.3560369)) or
through Caltech DATA under the DOI: 10.22002/D1.1315.

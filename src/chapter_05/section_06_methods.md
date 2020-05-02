## Materials & Methods

### Bacterial strains and growth conditions

The bacterial strains are described in Table 9.1. The parent strain for the
mutants used in this study was MJF641 [@edwards2012], a strain which had all seven
mechanosensitive channels deleted. The MscL-sfGFP coding region from MLG910
[@bialecka-fornal2012] was integrated into MJF641 by P1 transduction, creating the strain
D6LG-Tn10. Selection pressure for MscL integration was created by
incorporating an osmotic shock into the transduction protocol, which favored
the survival of MscL-expressing stains relative to MJF641 by ~100-fold.
Screening for integration candidates was based on fluorescence expression of
plated colonies. Successful integration was verified by sequencing. Attempts
to transduce RBS-modified MscL-sfGFP coding regions became increasingly
inefficient as the targeted expression level of MscL was reduced. This was
due to the decreasing fluorescence levels and survival rates of the
integration candidates. Consequently, Shine-Dalgarno sequence modifications were made by
inserting DNA oligos with lambda Red-mediated homologous recombination, i.e.,
recombineering [@sharan2009]. The oligos had a designed mutation (Fig. @fig:mscl_boxplot)
flanked by ~25 base pairs that matched the targeted MscL region (Table S2). A
two-step recombineering process of selection followed by counter selection
using a *tetA-sacB* gene fusion cassette [@li2013] was chosen because of its
capabilities to integrate with efficiencies comparable to P1 transduction and
not leave antibiotic resistance markers or scar sequences in the final
strain. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; To prepare the strain D6LG-Tn10 for this scheme, the Tn10 transposon
containing the *tetA* gene needed to be removed to avoid interference with the
*tetA-sacB* cassette. Tn10 was removed from the middle of the *ycjM* gene with
the primer Tn10delR (Table S2) by recombineering, creating the strain D6LG
(SD0). Counter selection against the *tetA* gene was promoted by using agar
media with fusaric acid [@bochner1980; @li2013]. The *tetA-sacB* cassette was
PCR amplified out of the strain XTL298 using primers MscLSPSac and MscLSPSacR
(Table S2). The cassette was integrated in place of the spacer region in front
of the MscL start codon of D6LG (SD0) by recombineering, creating the
intermediate strain D6LTetSac. Positive selection for cassette integration
was provided by agar media with tetracycline. Finally, the RBS modifying
oligos were integrated into place by replacing the *tetA-sacB* cassette by
recombineering. Counter selection against both *tetA* and *sacB* was ensured by
using agar media with fusaric acid and sucrose [@li2013], creating the Shine-Dalgarno
mutant strains used in this work.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Strain cultures were grown in 5 mL of LB-Lennox media with antibiotic
(apramycin) overnight at 37$^\circ$C. The next day, 50 $\mu$L of overnight culture was
inoculated into 5 mL of LB-Lenox with antibiotic and the culture was grown to
OD$_\text{600nm} \sim .25$. Subsequently, 500 $\mu$L of that culture was inoculated into 5 mL of
LB-Lennox supplemented with 500mM of NaCl and the culture was regrown to OD$_\text{600nm}$
~0.25. A 1 mL aliquot was taken and used to load the flow cell.

### Flow cell

All experiments were conducted in a home-made
flow cell as is shown in @Fig:flow_cell(A). This flow cell has two inlets which allow
media of different osmolarity to be exchanged over the course of the
experiment. The imaging region is approximately 10 mm wide and 100 $\mu$m in
depth. All imaging took place within 1 – 2 cm of the outlet to avoid imaging
cells within a non-uniform gradient of osmolarity. The interior of the flow
cell was functionalized with a 1:400 dilution of polyethylenimine prior to
addition of cells with the excess washed away with water. A dilute cell
suspension in LB Lennox with 500 mM NaCl was loaded into one inlet while the
other was connected to a vial of LB medium with no NaCl. This hypotonic
medium was clamped during the loading of the cells.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Once the cells had adhered to the polyethylenimine
coated surface, the excess cells were washed away with the 500 mM NaCl growth
medium followed by a small (\~20 $\mu$L) air bubble. This air bubble forced
the cells to lay flat against the imaging surface, improving the time-lapse
imaging. Over the observation period, cells not exposed to an osmotic shock
were able to grow for 4 – 6 divisions, showing that the flow cell does not
directly impede cell growth.

### Imaging conditions

All imaging was performed in a flow cell held at
30$^\circ$C on a Nikon Ti-Eclipse microscope outfitted with a Perfect Focus system
enclosed in a Haison environmental chamber (approximately 1$^\circ$C regulation
efficiency). The microscope was equipped with a 488 nm laser excitation
source (CrystaLaser) and a 520/35 laser optimized filter set (Semrock). The
images were collected on an Andor iXon EM+ 897 EMCCD camera and all microscope
and acquisition operations were controlled via the open source $\mu$Manager
microscope control software [@edelstein2014]. Once cells were securely mounted onto the
surface of the glass coverslip, between 15 and 20 positions containing 5 to
10 cells were marked and the coordinates recorded. At each position, a phase
contrast and GFP fluorescence image was acquired for segmentation and
subsequent measurement of channel copy number. To perform the osmotic shock,
LB media containing no NaCl was pulled into the flow cell through a syringe
pump. To monitor the media exchange, both the high salt and no salt LB media
were supplemented with a low-affinity version of the calcium-sensitive dye
Rhod-2 (250 nM; TEF Labs) which fluoresces when bound to Ca^2+^. The no salt
medium was also supplemented with 1$\mu$M CaCl$_2$ to make the media mildly
fluorescent and the exchange rate was calculated by measuring the
fluorescence increase across an illuminated section of one of the positions.
These images were collected in real time for the duration of the shock. The
difference in measured fluorescence between the pre-shock images and those at
the end of the shock set the scale of a 500 mM NaCl down shock. The rate was
calculated by fitting a line to the middle region of this trace. Further
details regarding this procedure can be found in @bialecka-fornal2015.

### Image Processing

Images were processed using a combination of
automated and manual methods. First, expression of MscL was measured via
segmenting individual cells or small clusters of cells in phase contrast and
computing the mean pixel value of the fluorescence image for each segmented
object. The fluorescence images were passed through several filtering
operations which reduced high-frequency noise as well as corrected for uneven
illumination of the excitation wavelength.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Survival or death classification was performed
manually using the CellProfiler plugin for ImageJ software (NIH). A survivor
was defined as a cell which was able to undergo at least two division events
after the osmotic down shock. Cell death was recognized by stark changes in
cell morphology including loss of phase contrast through ejection of
cytoplasmic material, structural decomposition of the cell wall and membrane,
and the inability to divide. To confirm that these morphological cues
corresponded with cell death, we probed cell viability on a subset of our
strains after osmotic shock through staining with propidium iodide, a DNA
intercalating dye commonly used to identifying dead cells (LIVE/DEAD BacLight
Bacterial Cell Viability Assay, Thermo Fisher). We found that our
classification based on morphology agreed with that based off of staining
within 1%. More information regarding these experiments can be found in 
Chapter 9. Cells which detached from the surface during the post-shock
growth phase or those which became indistinguishable from other cells due to
clustering were not counted as survival or death and were removed from the
dataset completely. A region of the cell was manually marked with 1.0
(survival) or 0.0 (death) by clicking on the image. The `xy` coordinates of
the click as well as the assigned value were saved as an `.xml` file for that
position.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The connection between the segmented cells and
their corresponding manual markers was automated. As the manual markings were
made on the first phase contrast image after the osmotic shock, small shifts
in the positions of the cell made one-to-one mapping with the segmentation
mask non-trivial. The linkages between segmented cell and manual marker were
made by computing all pairwise distances between the manual marker and the
segmented cell centroid, taking the shortest distance as the true pairing.
The linkages were then inspected manually and incorrect mappings were
corrected as necessary.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;All relevant statistics about the segmented
objects as well as the sample identity, date of acquisition, osmotic shock
rate, and camera exposure time were saved as `.csv` files for each individual
experiment. A more in-depth description of the segmentation procedure as well
as the relevant code can be accessed as a Jupyter Notebook at
(`http://rpgroup.caltech.edu/mscl_survival`).

### Calculation of effective channel copy number
To compute the MscL channel copy number, we relied
on measuring the fluorescence level of a bacterial strain in which the mean
MscL channel copy number was known via fluorescence microscopy
[@bialecka-fornal2012]. *E. coli* strain MLG910, which expresses the MscL-sfGFP fusion
protein from the wild-type SD sequence, was grown under identical conditions
to those described in Bialecka-Fornal et al. 2015 in LB Miller medium (BD Medical Sciences) to an OD~600nm~ of ~0.3. The cells
were then diluted ten fold and immobilized on a rigid 2% agarose substrate
and placed onto a glass bottom petri dish and imaged in the same conditions
as described previously.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Images were taken of six biological replicates
of MLG910 and were processed identically to those in the osmotic shock
experiments. A calibration factor between the average cell fluorescence level
and mean MscL copy number was then computed. We assumed that all measured
fluorescence (after filtering and background subtraction) was derived from
the MscL-sfGFP fusion, $$ \langle I_\text{tot}\rangle = \alpha \langle N
\rangle, $${#eq:ian} in which $\alpha$ is the calibration factor and $\langle
N \rangle$ is the mean cellular MscL-sfGFP copy number as reported in
@bialecka-fornal2012. To correct for errors in segmentation, the intensity
was computed as an areal density $\langle I_A \rangle$ and was multiplied by
the average cell area $\langle A \rangle$ of the population. The calibration
factor was therefore computed as $$ \alpha = {\langle I_A \rangle \langle A
\rangle \over \langle N \rangle}. $${#eq:calibration_factor}
We used Bayesian inferential methods to compute this
calibration factor taking measurement error and replicate-to-replicate
variation into account. The resulting average cell area and calibration
factor was used to convert the measured cell intensities from the osmotic
shock experiments to cell copy number. The details of this inference are
described in depth in the supplemental Chapter 9.

### Logistic regression

We used Bayesian inferential methods to find the
most probable values of the coefficients $\beta_0$ and $\beta_1$ and the
appropriate credible regions and is described in detail in the supplemental information (*Logistic Regression*).
Briefly, we used Markov chain Monte Carlo (MCMC) to sample from the log
posterior distribution and took the most probable value as the mean of the
samples for each parameter. The MCMC was performed using the Stan
probabilistic programming language [@carpenter2017] and all models can be
found on the GitHub repository (`http://github.com/rpgroup-pboc/mscl_survival`).

### Calculation of survival probability error

The vertical error bars for the points shown in @Fig:survival represent
our uncertainty in the survival probability given our measurement of $n$
survivors out of a total $N$ single-cell measurements. The probability
distribution of the survival probability $p_s$ given these measurements can
be written using Bayes' theorem as
$$
g(p_s\,\vert\, n, N) = {f(n\,\vert\,p_s, N)g(p_s) \over f(n\,\vert\, N)},
$${#eq:probability_bayes}

where $g$ and $f$ represent probability density functions over parameters and
data, respectively. The likelihood $f(n\,\vert p_s, N)$ represents the
probability of measuring $n$ survival events, given a total of $N$
measurements each with a probability of survival $p_s$. This matches the
story for the Binomial distribution and can be written as
$$
f(n\,\vert\,p_s, N) = {N! \over n!(N - n)!}p_s^n(1 - p_s)^{N - n}.
$${#eq:binomial}

To maintain maximal ignorance we can assume that any value for $p_s$ is
valid, such that is in the range [0, 1]. This prior knowledge, represented by
$g(p_s)$, can be written as
$$
g(p_s) = \begin{cases}1 & 0\leq p_s\leq 1 \\
0 & \text{otherwise} \end{cases}.
$${#eq:uniform_prob}
We can also assume maximal ignorance for the total number of survival events
we could measure given $N$ observations, $f(n\, \vert\, N)$. Assuming all
observations are equally likely, this can be written as
$$
f(n\,\vert\, N) = {1 \over N + 1}
$${#eq:evidence}

where the addition of one comes from the possibility of observing zero
survival events. Combining @Eq:binomial, @Eq:uniform_prob, and @Eq:evidence, the
posterior distribution $g(p_s\,\vert\, n, N)$ is
$$
g(p_s\,\vert\, n, N) = {(N+1)! \over n!(N - n)!}p_s^{n}(1 - p_s)^{N - n}.
$${#eq:probability_posterior}

The most probable value of $p_s$, where the posterior probability
distribution given by @Eq:probability_posterior is maximized, can be
found by computing the point at which derivative of the log posterior with
respect to $p_s$ goes to zero,
$$
{d\log g(p_s\,\vert\,n, N) \over d p_s} = {n \over p_s} - {N - n  \over 1 - p_s} = 0.
$${#eq:deriv_ps}
Solving @Eq:deriv_ps for $p_s$ gives the most likely value for the
probability,
$$
p_s^* = {n \over N}.
$${#eq:most_prob_prob}

So long as $N >> np_s^*$, @Eq:probability_posterior can be approximated as
a Gaussian distribution with a mean $p_s^*$ and a variance $\sigma_{p_s}^2$.
By definition, the variance of a Gaussian distribution is computed as the
negative reciprocal of the second derivative of the log posterior evaluated
at $p_s = p_s^*$,
$$
\sigma_{p_s}^2 = - \left({d^2 \log g(p_s\,\vert\, n, N) \over
dp_s^2}\Bigg\vert_{p_s=p_s^*}\right)^{-1}.
$${#eq:variance_def}

Evaluating @Eq:variance_def yields
$$
\sigma_{p_s}^2 = {n(N-n)\over N^3}.
$${#eq:prob_variance}

Given @Eq:most_prob_prob and @Eq:prob_variance, the most-likely
survival probability and estimate of the uncertainty can be expressed as
$$
p_s = p_s^* \pm \sigma_{p_s}.
$${#eq:}

### Data and software availability

All raw image data is freely available and is stored on the CaltechDATA
Research Data Repository. The raw Markov chain Monte Carlo
samples are stored as `.csv` files on CaltechDATA. All processed
experimental data, Python, and Stan code used in this work are freely
available through the  paper [GitHub repository](http://github.com/rpgroup-pboc/mscl_survival) accessible
through DOI: 10.5281/zenodo.1252524. The scientific community is invited to
fork our repository and open constructive issues.


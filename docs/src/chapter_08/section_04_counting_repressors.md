## Counting Repressors


In this section, we expand upon the theoretical and experimental
implementation of the fluorescence calibration method derived in 
@rosenfeld2005. We cover several experimental data validation steps as
well as details regarding the parameter inference. Finally, we comment
on the presence of a systematic error in the repressor counts due to
continued asynchronous division between sample preparation and imaging.

### Theoretical Background of the Binomial Partitioning Method {#sec:cal_factor}

A key component of this work is the direct measurement of the repressor copy
number in each growth condition using fluorescence microscopy. To translate
between absolute fluorescence and protein copy number, we must be able to
estimate the average brightness of a single fluorophore or, in other words,
determine a calibration factor $\alpha$ that permits translation from copy
number to intensity or vice versa. Several methods have been used over the
past decade to estimate this factor, such as measuring single-molecule
photobleaching steps [@garcia2011c; @bialecka-fornal2012], measurement of *in
vivo* photobleaching rates [@nayak2011; @kim2016], and through measuring the
partitioning of fluorescent molecules between sibling cells after cell
division [@rosenfeld2005; @rosenfeld2006; @brewster2014]. In this work, we
used the latter method to estimate the brightness of a single LacI-mCherry
dimer. Here, we derive a simple expression which allows the determination of
$\alpha$ from measurements of the fluorescence intensities of a collection of
sibling cells.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In the absence of measurement error, the fluorescence intensity of a
given cell is proportional to the total number of fluorescent proteins
$N_\text{prot}$ by some factor $\alpha$,
$$
I_\text{cell} = \alpha N_\text{prot}.
$${#eq:growth_si_ian}
Assuming that no fluorophores are produced or
degraded over the course of the cell cycle, the fluorescent proteins
will be partitioned into the two siblings cells such that the intensity
of each sibling can be computed as
$$
I_1 = \alpha N_1\, ;\, I_2 = \alpha (N_\text{ tot} - N_1),
$${#eq:ian_division} 
where $N_\text{ tot}$ is the total
number of proteins in the parent cell and $N_1$ and $N_2$ correspond to
the number of proteins in sibling cell 1 and 2, respectively. This
explicitly states that fluorescence is conserved upon a division,
$$
I_{tot} = I_1 + I_2.
$${#eq:fluo_cons}
As the observed intensity is directly
proportional to the number of proteins per cell, measuring the variance
in intensity between sibling cells provides some information as to how
many proteins were there to begin with. We can compute these
fluctuations as the squared intensity difference between the two
siblings as
$$
\langle \left(I_1 - I_2\right)^2 \rangle = \langle \left(2I_1 -
    I_\text{ tot}\right)^2\rangle.
$${#eq:step1}
We can relate @Eq:step1 in terms of the number of proteins using @Eq:growth_si_ian
as
$$
\langle \left(I_1 - I_2\right)^2 \rangle = 4\alpha^2 \langle N_1^2
    \rangle - 4\alpha^2\langle N_1 \rangle N_\text{ tot} +
    \alpha^2N_\text{ tot}^2,
$${#eq:step2}
where the squared fluctuations are now cast in
terms of the first and second moment of the probability distribution for
$N_1$.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Without any active partitioning of the proteins into the sibling cells,
one can model the probability distribution $g(N_1)$ of finding $N_1$
proteins in sibling cell 1 as a binomial distribution,
$$
g(N_1\vert\, N_\text{ tot}, p) = \frac{N_\text{
    tot}!}{N_1!(N_\text{ tot} - N_1)!}p^{N_1}(1 - p)^{N_\text{
    tot} - N_1},
$${#eq:growth_si_binomial}
where $p$ is the probability of a protein
being partitioned into one sibling over the other. With a probability
distribution for $N_1$ in hand, we can begin to simplify 
@Eq:step2. Recall that the mean and variance of a binomial distribution are
$$
\langle N_{1} \rangle = N_\text{ tot} p,
$${#eq:binomial_mean}
and
$$
\langle N_1^2 \rangle - \langle N_1\rangle^2 = N_\text{ tot}p(1-p),
$${#eq:binom_variance}
respectively. With knowledge of the mean and variance, we can solve for the
second moment as
$$
\langle N_1^2 \rangle = N_\text{ tot}p(1 - p) + N_\text{
tot}^2p^2 = N_\text{ tot}p\left(1 - p + N_\text{
tot}p\right).
$${#eq:binomial_mom2}
By plugging @Eq:binomial_mean and @Eq:binomial_mom2 into our expression for the fluctuations
(@Eq:step2), we arrive at
$$
\langle \left( I_1 - I_2\right)^2 \rangle =
    4\alpha^2\left[(N_\text{ tot}p[1 - p - N_\text{ tot}p]) -
    N_\text{ tot}^2\left(p + \frac{1}{4}\right)\right],
$${#eq:p_agnostic_result}
which is now defined in terms of the total number of proteins present in the parent cell. Assuming that the
proteins are equally partitioned $p=1/2$, @Eq:p_agnostic_result reduces to
$$
\langle \left(I_1 - I_2\right)^2\rangle = \alpha^2N_\text{ tot} =
    \alpha I_\text{ tot}.
$${#eq:fluct_Itot}
Invoking our assertion that fluorescence is
conserved (@Eq:fluo_cons), $I_{tot}$ is equivalent to the sum total
fluorescence of the siblings,
$$
\langle \left(I_1 - I_2\right)^2 \rangle = \alpha (I_1 + I_2).
$${#eq:dilution}
Thus, given snapshots of cell intensities and
information of their lineage, one can compute how may arbitrary
fluorescence units correspond to a single fluorescent protein.

### Cell Husbandry and Time-Lapse Microscopy 

The fluorescence calibration method  was first described and implemented by
@rosenfeld2005 followed by a more in-depth approach on the statistical
inference of the calibration factor in @rosenfeld2006. In both of
these works, the partitioning of a fluorescent protein was tracked
across many generations from a single parent cell, permitting inference
of a calibration factor from a single lineage. @brewster2014
applied this method in a slightly different manner by
quantifying the fluorescence across a large number of *single* division
events. Thus, rather than examining the partitioning of fluorescence
down many branches of a single family tree, it was estimated from an
array of single division events where the fluorescence intensity of the
parent cell was variable. In the present work, we take a similar
approach to that of @brewster2014 and examine the partitioning of
fluorescence among a large number of independent cell divisions.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A typical experimental work-flow is shown in
@Fig:growth_si_microscopy_workflow. For each experiment, the strains were
grown in varying concentrations of ATC to tune the expression of the
repressor. Once the cells had reached exponential phase growth (OD$_{600nm}
\approx 0.3$), the cells were harvested and prepared for imaging. This
involved two separate sample handling procedures, one for preparing samples
for lineage tracking and estimation of the calibration factor and another for
taking snapshots of cells from each ATC induction condition for the
calculation of fold-change.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To prepare cells for the calibration factor
measurement, a 100 $\mu$L aliquot of each ATC induction condition was
combined and mixed in a 1.5 mL centrifuge tube. This cell mixture was then
pelleted at 13000$\times g$ for 1 -- 2 minutes. The supernatant containing
ATC was then aspirated and the pellet was resuspended in an equal volume of
sterile growth media without ATC. This washing procedure was repeated three
times to ensure that any residual ATC had been removed from the culture and
that expression of LacI-mCherry had ceased. Once washed and resuspended, the
cells were diluted ten fold into sterile M9 medium and then imaged on a rigid
agarose substrate. Depending on the precise growth condition, a variety of
positions were imaged for 1.5 to 4 hours with a phase contrast image acquired
every 5 to 15 minutes to facilitate lineage tracking. On the final image of
the experiment, an mCherry fluorescence image was acquired of every position.
The experiments were then transferred to a computing cluster and the images
were computationally analyzed, as described in the next section.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;During the washing steps, the remaining ATC
induced samples were prepared for snapshot imaging to determine the repressor
copy number and fold-change in gene expression for each ATC induction
condition. Without mixing the induction conditions together, each ATC induced
sample was diluted 1:10 into sterile M9 minimal medium and vigorously mixed.
Once mixed, a small aliquot of the samples were deposited onto rigid agarose
substrates for later imaging. While this step of the experiment was
relatively simple, the total preparation procedure typically lasted between
30 and 60 minutes. As is discussed later, the continued growth of the
asynchronously growing culture upon dilution into the sterile medium results
in a systematic error in the calculation of the repressor copy number.

![**An experimental workflow for time-lapse imaging.** (A) the series of steps
followed in a given experiment. Cells are grown in various concentrations of ATC
(shaded red cultures) to an OD$_{600nm} \approx 0.3$. Equal aliquots of each
ATC-induced culture are mixed into a single eppendorf tube and pelleted via
centrifugation. The supernatant containing ATC is aspirated and replaced with an
equal volume of sterile, ATC-free growth medium. This washing procedure is
repeated three times to ensure that residual ATC is removed from the culture and
expression of the LacI-mCherry fusion is ceased. After a final resuspension in
sterile ATC-free medium, the cell mixture is diluted 10 fold to reduce cell
density.  small aliquot of this mixture is then mounted and imaged at 100x
magnification until at least one cell division has occurred. (B) Two
representative microcolonies from a time-lapse growth experiment. The time point
is provided above each image. After at least one division has occurred, a final
mCherry fluorescence image is acquired and
quantified.](ch8_figS3){#fig:growth_si_microscopy_workflow short-caption="An experimental
workflow for determination of a fluorescence calibration factor."}

### Lineage Tracking and Fluorescence Quantification

Segmentation and lineage tracking of both the fluorescence snap shots
and time-lapse growth images were performed using the SuperSegger
`v1.03` [@stylianidou2016] software using MATLAB R2017B (MathWorks,
Inc). The result of this segmentation is a list of matrices for each
unique imaged position with identifying data for each segmented cell
such as an assigned ID number, the ID of the sibling cell, the ID of the
parent cell, and various statistics. These files were then analyzed
using Python 3.7. All scripts and software used to perform this analysis
can be found on the [associated paper website](https://www.rpgroup.caltech.edu/mwc_growth) and [GitHub
repository](https://github.com/rpgroup-pboc/mwc_growth).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Using the ID numbers assigned to each cell in a
given position, we matched all sibling pairs present in the last frame of the
growth movie when the final mCherry fluorescence image was acquired. These
cells were then filtered to exclude segmentation artifacts (such as
exceptionally large or small cells) as well as any cells which the
SuperSegger software identified as having an error in segmentation. Given the
large number of cells tracked in a given experiment, we could not manually
correct these segmentation artifacts, even though it is possible using the
software. To err on the side of caution, we did not consider these edge cases
in our analysis.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With sibling cells identified, we performed a
series of validation checks on the data to ensure that both the experiment
and analysis behaved as expected. Three validation checks are illustrated in
@Fig:sanity_check. To make sure that the computational pairing of sibling
cells was correct, we examined the intensity distributions of each sibling
pair. If siblings were being paired solely on their lineage history and not
by other features (such as size, fluorescence, etc), one would expect the
fluorescence distributions between the two sibling cells to be identical.
@Fig:sanity_check(A) shows the nearly identical intensity distributions of
all siblings from a single experiment, indicating that the pairing of
siblings is independent of their fluorescence.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Furthermore, we examined the partitioning of the fluorescence between
the siblings. In @Eq:binomial_mean, we defined the mean number of proteins
inherited by one sibling. This can be easily translated into the
language of intensity as
$$
\langle I_1 \rangle = \alpha \left(I_1 + I_2\right)p,
$${#eq:binom_mean_int}
where we assume that fluorescence is conserved during a division event such that
$I_\text{ tot} = I_1 + I_2$. As described previously, 
we make the simplifying assumption that partitioning of the fluorophores
between the two sibling cells is fair, meaning that $p = 0.5$. We can
see if this approximation is valid by computing the fractional intensity
of each sibling cell as
$$
p = \frac{\langle I_1 \rangle}{\langle I_1 + I_2 \rangle}.
$${#eq:frac_intensity}
@Fig:sanity_check(B) shows the distribution of the
fractional intensity for each sibling pair. The distribution is
approximately symmetric about 0.5, indicating that siblings are
correctly paired and that the partitioning of fluorescence is
approximately equal between siblings. Furthermore, we see no correlation
between the cell volume immediately after division and the observed
fractional intensity. This suggests that the probability of partitioning
to one sibling or the other is not dependent on the cell size. An
assumption (backed by experimental measurements
[@garcia2011; @phillips2019]) in our thermodynamic model is that all
repressors in a given cell are bound to the chromosome, either
specifically or nonspecifically. As the chromosome is duplicated and
partitioned into the two siblings without fail, our assumption of
repressor adsorption implies that partitioning should be independent of
the size of the respective sibling. The collection of these validation
statistics give us confidence that both the experimentation and the
analysis are properly implemented and not introducing bias into our
estimation of the calibration factor.

![**Experimental sanity checks and inference of a fluorescence calibration
factor.** (A) Intensity distributions of sibling cells after division.
Arbitrarily labeled "sibling 1" and "sibling 2" distributions are shown in
purple and orange, respectively. Similarity of the distributions illustrates
lack of intensity bias on sibling pair assignment. (B) the fractional intensity
of sibling 1 upon division. For each sibling pair, the fractional intensity is
computed as the intensity of sibling 1 $I_2$ divided by the summed intensities
of both siblings $I_1 + I_2$. (C) Partitioning intensity as a function of cell
volume. The fractional intensity of every sibling cell is plotted against its
estimated newborn volume. The [Python code                                                
(`ch8_figS4.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS4.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch8_figS4){#fig:sanity_check
short-caption="Experimental sanity checks and inference of a fluorescence
calibration factor."}

### Statistical Inference of the Fluorescence Calibration Factor

As is outlined in the Materials \& Methods section of the Chapter 4, we
took a Bayesian approach towards our inference of the calibration factor
given fluorescence measurements of sibling cells. Here, we expand in
detail on this statistical model and its implementation.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To estimate the calibration factor $\alpha$ from a set of lineage
measurements, we assume that fluctuations in intensity resulting from
measurement noise is negligible compared to that resulting from binomial
partitioning of the repressors upon cell division. In the absence of
measurement noise, the intensity of a given cell $I_{cell}$ can be
directly related to the total number of fluorophores $N$ through a
scaling factor $\alpha$, such that 
$$
I_{cell} = \alpha N.
$${#eq:no_err}
Assuming no fluorophores are produced or degraded over the course of a
division cycle, the fluorescence of the parent cell before division is equal
to the sum of the intensities of the sibling cells,
$$
I_{parent} = I_1 + I_2 = \alpha\left(N_1 + N_2\right),
$${#eq:conservation_law}
where subscripts $1$ and $2$ correspond to arbitrary labels of the two sibling
cells. We are ultimately interested in knowing the probability of a given value
of $\alpha$ which, using Bayes' theorem, can be written as
$$
g(\alpha\, \vert\, I_1, I_2) = \frac{f(I_1, I_2\, \vert\, \alpha) g(\alpha)}
    {f(I_1, I_2)},
$${#eq:app_bayes_generic}
where we have used $g$ and $f$ to denote probability densities over parameters and data, respectively. The
first quantity in the numerator $f(I_1, I_2\,\vert\, \alpha)$ describes
the likelihood of observing the data $I_1, I_2$ given a value for the
calibration factor $\alpha$. The term $g(\alpha)$ captures all prior
knowledge we have about what the calibration factor could be, remaining
ignorant of the collected data. The denominator $f(I_1, I_2)$ is the
likelihood of observing our data $I_1, I_2$ irrespective of the
calibration factor and is a loose measure of how well our statistical
model describes the data. As it is difficult to assign a functional form
to this term and serves as a multiplicative constant, it can be
neglected for the purposes of this work.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Knowing that the two observed sibling cell
intensities are related, conditional probability allows us to rewrite the
likelihood as
$$
f(I_1\,I_2\,\vert\,\alpha) = f(I_1\,\vert\,I_2, \alpha)f(I_2\,\vert\,\alpha),
$${#eq:arb_unlabeled_eq}
where $f(I_2\,\vert\, \alpha)$ describes the likelihood of observing
$I_2$ given a value of $\alpha$. As $I_2$ can take any value with equal
probability, this term can be treated as a constant. Through change of
variables and noting that $I_2 = \alpha N_2$, we can cast the likelihood
$f(I_1\,\vert\,I_2, \alpha)$ in terms of the number of proteins $N_1$
and $N_2$ as
$$
f(I_1\,\vert\, I_2, \alpha) = f(N_1\,\vert\,N_2, \alpha)\bigg\vert {d N_1
\over dI_1}\bigg\vert= \frac{1}{\alpha}f(N_1\,\vert\,N_2,\alpha).
$${#eq:change_of_var}
Given that the proteins are binomially
distributed between the sibling cells with a probability $p$ and that
the intensity is proportional to the number of fluorophores, this
likelihood becomes
$$
f(I_1\,\vert\,I_2, \alpha, p) = \frac{1}{\alpha}{\left(\frac{I_1 +
I_2}{
\alpha}\right)! \over \left({I_1 \over \alpha}\right)!\left({I_2 \over
\alpha}\right)!}p^{{I_1 \over \alpha}}(1 - p)^{I_2 \over \alpha}
$${#eq:binom_like}
However, The quantity $I / \alpha$ is not
exact, making calculation of its factorial undefined. These factorials
can therefore be approximated by a gamma function as
$n! = \Gamma(n + 1)$.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Assuming that partitioning of a protein between
the two sibling cells is a fair process ($p = 1/2$), @Eq:binom_like can be generalized to a set of lineage
measurements $[I_1, I_2]$ as
$$
f([I_1]\,\vert\,[I_2], \alpha) = {1 \over \alpha^k}
\prod\limits_{i}^k {\Gamma\left({I_1 + I_2 \over \alpha} + 1\right) \over
\Gamma\left({I_1 \over \alpha} + 1\right)\Gamma\left({I_2 \over \alpha} + 1
\right)}2^{-{I_1 + I_2 \over \alpha}}, 
$${#eq:full_likelihood}
where $k$ is the number of division events observed.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With a likelihood in place, we can now assign a
functional form to the prior distribution for the calibration factor
$g(\alpha)$. Though ignorant of data, the experimental design is such that
imaging of a typical highly-expressing cell will occupy 2/3 of the dynamic
range of the camera. We can assume it's more likely that the calibration
factor will be closer to 0 a.u. than the bit depth of the camera (4095 a.u.)
or larger. We also know that it is physically impossible for the fluorophore
to be less than 0 a.u., providing a hard lower-bound on its value. We can
therefore impose a weakly informative prior distribution as a half normal
distribution,
$$
g(\alpha) = \sqrt{2 \over \pi\sigma^2}\exp\left[{-\alpha^2 \over
   2\sigma^2}\right]\,;\, \forall \alpha > 0. 
$${#eq:growth_alpha_prior}
where the standard deviation is large, for example, $\sigma = 500$ a.u.
/ fluorophore. We evaluated the posterior distribution using Markov
chain Monte Carlo (MCMC) as is implemented in the Stan probabilistic
programming language [@carpenter2017]. The `.stan` file associated with
this model along with the Python code used to execute it can be accessed
on the [paper website](https://www.rpgroup.caltech.edu/mwc_growth) and
   [GitHub repository](https://github.com/rpgroup-pboc/mwc_growth).
@Fig:cal_factor_posts  shows the posterior probability
distributions of the calibration factor estimated for several biological
replicates of the glucose growth condition at 37$^\circ$ C. For each
posterior distribution, the mean and standard deviation was used as the
calibration factor and uncertainty for the corresponding data set.

![**Representative posterior distributions for the calibration factor.** (A)
Empirical cumulative distribution functions and (B) histograms of the posterior
distribution for the calibration factor $\alpha$ for several biological
replicates. Different colors indicate different biological replicates of
experiments performed in glucose supplemented medium at 37$^\circ$ C.
The [Python code                                                
(`ch8_figS5.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS5.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch8_figS5){#fig:cal_factor_posts short-caption="Representative posterior distributions
for inference of the fluorescence calibration factor."}

### Correcting for Systematic Experimental Error

While determination of the calibration factor relies on time-resolved
measurement of fluorescence partitioning, we computed the repressor copy
number and fold-change in gene expression from still snapshots of each
ATC induction condition and two control samples, as is illustrated in
@Fig:growth_si_microscopy_workflow
Given these snapshots, individual cells were segmented again using the SuperSegger software in MATLAB
R2017b. The result of this analysis is an array of single-cell
measurements of the YFP and mCherry fluorescence intensities. With these
values and a calibration factor estimated by @Eq:full_likelihood 
and @Eq:growth_alpha_prior, we can compute the estimated number of
repressors per cell in every condition.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;However, as outlined previously
and in @Fig:growth_si_microscopy_workflow, the sample preparation steps for
these experiments involve several steps which require careful manual
labor. This results in an approximately 30 to 60 minute delay from when
production of the LacI-mCherry construct is halted by removal of ATC to
actual imaging on the microscope. During this time, the diluted cultures
are asynchronously growing, meaning that the cells of the culture are at
different steps in the cell cycle. Thus, at any point in time, a subset
of cells will be on the precipice of undergoing division, partitioning
the cytosolic milieu between the two progeny. As LacI-mCherry is no
longer being produced, the cells that divided during the dwell time from
cell harvesting to imaging will have reduced the number of repressors by
a factor of 2 on average. This principle of continued cell division is
shown @Fig:correction_origin.


![**Continued division results in a systematic error in repressor counts.**
During the sample preparation steps, the asynchronous culture continues to
divide though the length of the sample preparation step is less than a cell
doubling time. Because production of LacI-mCherry is halted, cells that complete
a division cycle during this time will partition their repressors between the
progeny. The total number of repressors in parent cells just before division is,
on average, twice that of the newborn cells.](ch8_figS6){#fig:correction_origin
short-caption="Origin of a systematic error in repressor counting due to
continued growth."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;How does this partitioning affect our
calculation and interpretation of the fold-change in gene expression? Like
the LacI-mCherry fusion, the YFP reporter proteins are also partitioned
between the progeny after a division event such that, on average, the total
YFP signal of the newborn cells is one-half that of the parent cell. As the
maturation time of the mCherry and YFP variants used in this work are
relatively long in *E. coli* [@balleza2018; @nagai2002], we can make the
assumption any newly-expressed YFP molecules after cells have divided are not
yet visible in our experiments. Thus, the fold-change in gene expression of
the average parent cell can be calculated given knowledge of the average
expression of the progeny.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The fold-change in gene expression is a relative measurement to a
control which is constitutively expressing YFP. As the latter control
sample is also asynchronously dividing, the measured YFP intensity of
the newborn constitutively expressing cells is on average 1/2 that of
the parent cell. Therefore, the fold-change in gene expression can be
calculated as 
$$
\langle \text{ fold-change}_{parent}\rangle =
\frac{2 \times \langle I^\text{(YFP)}_\text{ newborn}(R > 0) \rangle}{2
\times \langle I^\text{(YFP)}_\text{ newborn}(R = 0)\rangle} =
\frac{\langle I^\text{(YFP)}_\text{ newborn}(R > 0) \rangle}{\langle
I^\text{(YFP)}_\text{ newborn}(R = 0)\rangle}.
$${#eq:fc_correction}
Thus, when calculating the fold-change in gene expression one does not
need to correct for any cell division that occurs between sample
harvesting and imaging as it is a relative measurement. However, the
determination of the repressor copy number is a *direct* measurement and
requires a consideration of unknown division events.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To address this source of systematic
experimental error, we examined the cell length distributions of all
segmented cells from the snapshots as well as the distribution of newborn
cell lengths from the time-lapse measurements. @Fig:length_distributions
shows that a significant portion of the cells from the snapshots (colored distributions) overlap with the
distribution of newborn cell lengths (grey distributions). For each
condition, we partitioned the cells from the snapshots into three bins
based on their lengths -- "small" cells had a cell length less than
2.5 $\mu$m, "medium" cells had lengths between 2.5 $\mu$m and 3.5
$\mu$m, and "large" cells being longer than 3.5 $\mu$m. These
thresholds were chosen manually by examining the newborn cell-size
distributions and are shown as red vertical lines in @Fig:length_distributions. 
Under this partitioning, we consider all "small" cells to have divided between cessation of
LacI-mCherry production and imaging, "medium"-length cells to have a
mixture of long newborn cells (from the tail of the newborn cell length
distribution) and cells that haven't divided, and cells in the \"long\"
group to be composed entirely of cells which did not undergo a division
over the course of sample preparation.

![**Cell length distributions of fluorescence snapshots and newborn cells.** Top
row shows the distribution of cell lengths (pole-to-pole, colored distributions)
off cells imaged for calculation of the repressor copy number and fold-change in
gene expression. The distribution of newborn cell lengths for that given
condition is shown in grey. The vertical red lines correspond to the cell length
threshold of 2.5$\mu$m and 3.5$\mu$m, from left to right, respectively. Cels to
the left of the first vertical line were identified as "small", cells in between
the two vertical lines to be "medium" sized, and "long" cells tot he right of
the second vertical line. Cells below the 2.4$\mu$m threshold were treated as
cells who divided after production of lacI-mCherry had been halted. Bottom row
shows the same data as the top row but as the empirical cumulative
distribution. The [Python code                                                
(`ch8_figS7.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS7.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch8_figS7){#fig:length_distributions short-caption="Cell length
distributions of fluorescence snapshots and newborn cells."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Given this coarse delineation of cell age by
length, we examined how correction factors could be applied to correct for
the the undesired systematic error due to dilution of repressors. We took the
data collected in this work and compared the results to the fold-change in
gene expression reported in the literature for the same regulatory
architecture. Without correcting for undesired cell division, the observed
fold-change in gene expression falls below the prediction and does not
overlap with data from the literature (@Fig:correction_factor(A), light
purple). Using the uncorrected measurements, we estimated the DNA binding
energy to be $\Delta\varepsilon_R \approx -15\, k_BT$ which does not agree
with the value for the O2 operator reported in @garcia2011 or with the
inferred DNA binding energies from the other data sources (@Fig:correction_factor (B)).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;These results emphasize the need to correct for undesired dilution of
repressors through cell division during the sample preparation period,
and we now consider several different manners of applying this
correction. We first consider the extreme case where all cells of the
culture underwent an undesired division after LacI-mCherry production
was halted. This means that the average repressor copy number measured
from all cells is off by a factor of 2. The result of applying a factor
of 2 correction to all measurements can be seen in @Fig:correction_factor
as dark red points. Upon applying
this correction, we find that the observed fold-change in gene
expression agrees with the prediction and data from other sources in the
literature. The estimated DNA binding energy $\Delta\varepsilon_R$ from
these data is also in agreement with other data sources
(@Fig:correction_factor(B), dark red). 
This result suggests
that over the course of sample preparation, a non-negligible fraction of
the diluted culture undergoes a division event before being imaged.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We now begin to relax assumptions as to what fraction of the measured
cells underwent a division event before imaging. As described above, in
drawing distinctions between "small", "medium", and "large" cells,
we assume the latter represent cells which *did not* undergo a division
between the harvesting and imaging of the samples. Thus, the repressor
counts of these cells should require no correction. The white-faced
points in @Fig:correction_factor(A) shows the fold-change in gene
expression of *only* the large cell fraction, which falls within error
of the theoretical prediction. Furthermore, the inferred DNA binding
energy falls within error of that inferred from data of @garcia2011
and that inferred from data assuming all
cells underwent a division event [@Fig:correction_factor(B)], though it does not fall within
error of the binding energy reported in @garcia2011.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The most realistic approach that can be taken
to avoid using only the "large" bin of cells is to to assume that all cells
with a length $\ell < 2.5\, \mu$m have undergone a division, requiring a
two-fold correction to their average repressor copy number. The result of
this approach can be seen in @Fig:correction_factor as purple points. The
inferred DNA binding energy from this correction approach falls within error
of that inferred from only the large cells [white-faced points in
@Fig:correction_factor(B)] as well as overlapping with the estimate treating
all cells as having undergone a division (red).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;There are several experimental techniques that
could be implemented to avoid needing to apply a correction factor as
described here. In @brewster2014, the fold-change in gene expression was
measured by tracking the production rate of a YFP reporter before and after a
single cell division event coupled with direct measurement of the repressor
copy number using the same binomial partitioning method. This implementation
required an extensive degree of manual curation of segmentation as well as
correcting for photobleaching of the reporter, which in itself is a
non-trivial correction [@garcia2011c]. The experimental approach presented
here sacrifices a direct measure of the repressor copy number for each cell
via the binomial partitioning method, but permits the much higher throughput
needed to assay the variety of environmental conditions. Ultimately, the
inferred DNA binding energy for all of the scenarios described above agree
within $1\, k_BT$, a value smaller than the natural variation in the DNA
binding energies of the three native *lacO*, which is $\approx 6\, k_BT$. For
the purposes of this work, we erred on the side of caution and only used the
cells deemed "large" for the measurements reported in Chapter 4 and the
remaining sections of this chapter. 

![**Influence of a correction factor on fold-change and the DNA binding energy.**
(A) Fold-change in gene expression measurements from [@brewster2014;
@garcia2011; @razo-mejia2018] along with data rom this work. Data from this work
shown are with no correction (purple), correcting for small cells only (dark
purple), using only the large cell fraction (black faced point), and treating
all cells as newborn cells (red). Where visible, errors correspond to the
standard error of 5 to 10 biological replicates. (B) estimated DNA binding
energy from each data set. Points are the median of the posterior distribution
over the DNA binding energy. Horizontal lines indicate the width of the
95\%
credible region of the posterior distribution. Grey lines in (A) and (B)
correspond to the theoretical prediction and estimated binding energy from
@garcia2011, respectively. The [Python code                                                
(`ch8_figS8.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS8.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch8_figS8){#fig:correction_factor
short-caption="Influence of a correction factor on fold-change and the DNA
binding energy."}



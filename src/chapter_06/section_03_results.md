## Results

### Quantifying the single-cell MscL copy number ###

&nbsp;&nbsp;&nbsp;&nbsp;The principal goal of this work is to examine the
contribution of a single mechanosensitive channel species to cell survival
under a hypo-osmotic shock. While this procedure could be performed for any
species of channel, we chose MscL as it is the most well characterized and
one of the most abundant species in *E. coli*. To probe the contribution of
MscL alone, we integrated an *mscL* gene encoding an MscL super-folder GFP
(sfGFP) fusion into a strain in which all seven known mechanosensitive
channel genes were deleted from the chromosome [@edwards2012]. Chromosomal
integration imposes strict control on the gene copy number compared to
plasmid borne expression systems, which is important to minimize variation in
channel expression across the population and provide conditions more
representative of native cell physiology. Abrogation of activity, mislocalization, or cytotoxicity are all inherent risks associated with
creating chimeric reporter constructs. In Supplement A, we carefully dissect
the functionality of this protein through electrophysiology (Fig. S1),
measure the rate of fluorophore maturation (Fig. S2), and quantify
potential aggregates (Figs. S3 and S4). To the best of our knowledge, the
MscL-sfGFP fusion protein functions identically to the wild-type, allowing us
to confidently draw conclusions about the physiological role this channel
plays in wild-type cells.

&nbsp;&nbsp;&nbsp;&nbsp;To modulate the number of MscL channels per cell, we
developed a series of mutants which were designed to decrease the expression
relative to wild-type. These changes involved direct alterations of the
Shine-Dalgarno sequence as well as the inclusion of AT hairpins of varying
length directly upstream of the start codon which influences the translation rate and hence the number of MscL proteins produced ([@Fig:boxplot]A). The six
Shine-Dalgarno sequences used in this work were chosen using the RBS binding
site strength calculator from the Salis Laboratory at the Pennsylvania State
University [@espahborujeni2014; @salis2009]. While the designed
Shine-Dalgarno sequence mutations decreased the expression relative to
wild-type as intended, the distribution of expression is remarkably wide
spanning an order of magnitude.

&nbsp;&nbsp;&nbsp;&nbsp;To measure the number of MscL channels per cell, we
determined a fluorescence calibration factor to translate arbitrary
fluorescence units per cell to protein copy number. While there have been
numerous techniques developed over the past decade to directly measure this
calibration factor, such as quantifying single-molecule photobleaching
constants or measuring the binomial partitioning of fluorescent proteins upon
cell division [@bialecka-fornal2012; @elowitz2002], we used *a priori*
knowledge of the mean MscL-sfGFP expression level of a particular *E. coli*
strain to estimate the average fluorescence of a single channel. In
Bialecka-Fornal et al. 2012 [@bialecka-fornal2012], the authors used
single-molecule photobleaching and quantitative Western blotting to probe the
expression of MscL-sfGFP under a wide range of growth conditions. To compute
a calibration factor, we used the strain MLG910 (*E. coli* K12 MG1655
$\phi$(mscL-sfGFP)) as a "standard candle", highlighted in yellow in
[@Fig:boxplot]B. This standard candle strain was grown and imaged in
identical conditions in which the MscL count was determined through fluorescence microscopy. The calibration
factor was computed by dividing the mean total cell fluorescence by the known
MscL copy number, resulting in a measure of arbitrary fluorescence units per
MscL channel. Details regarding this calculation and appropriate propagation
of error as well as its sensitivity to varying growth media can be found in the Materials & Methods as well as Supplement B (Fig. S5 - S8).

&nbsp;&nbsp;&nbsp;&nbsp;While it is seemingly straightforward to use this calibration
factor to determine the total number of channels per cell for wild-type or
highly expressing strains, the calculation for the lowest expressing strains
is complicated by distorted cell morphology. We observed that as the channel
copy number decreases, cellular morphology becomes increasingly aberrant with
filamentous, bulging, and branched cells becoming more abundant (Fig. S7A).
This morphological defect has been observed when altering the abundance of
several species of mechanosensitive channels, suggesting that they play an
important role in general architectural stability [@bialecka-fornal2012;
@bialecka-fornal2015]. As these aberrant morphologies can vary widely in size
and shape, calculating the number of channels per cell becomes a more nuanced
endeavor. For example, taking the total MscL copy number for these cells
could skew the final calculation of survival probability as a large but
severely distorted cell would be interpreted as having more channels than a
smaller, wild-type shaped cell (Fig. S7B). To correct for this pathology, we
computed the average expression level per unit area for each cell and
multiplied this by the average cellular area of our standard candle strain
which is morphologically indistinguishable from wild-type *E. coli*, allowing
for the calculation of an effective channel copy number. The effect of this
correction can be seen in Fig. S7C and D, which illustrate that there is no
other correlation between cell area and channel expression.

&nbsp;&nbsp;&nbsp;&nbsp;Our calculation of the effective channel copy number
for our suite of Shine-Dalgarno mutants is shown in [@Fig:boxplot]B. The
expression of these strains cover nearly three orders of magnitude with the
extremes ranging from approximately four channels per cell to nearly one
thousand. While the means of each strain are somewhat distinct, the
distributions show a large degree of overlap, making one strain nearly
indistinguishable from another. This variance is a quantity that is lost in
the context of bulk scale experiments but can be accounted for via
single-cell methods.

![Control of MscL expression and calculation of channel copy number. (A)
Schematic view of the expression modifications performed in this work. The
beginning portion of the native *mscL* sequence is shown with the
Shine-Dalgarno sequence, spacer region, and start codon shaded in red, green,
and blue, respectively. The Shine-Dalgarno sequence was modified through the
Salis lab Ribosomal Binding Strength calculator [@espahborujeni2014;
@salis2009]. The wild-type sequence (MLG910) is shown in black with mutations
for the other four Shine-Dalgarno mutants highlighted in red. Expression was
further modified by the insertion of repetitive `AT` bases into the spacer
region, generating hairpins of varying length which acted as a thermodynamic
barrier for translation initiation. (B) Variability in effective channel copy
number is computed using the standard candle. The boxes represent the
interquartile region of the distribution, the center line displays the
median, and the whiskers represent 1.5 times the maximum and minimum of the
interquartile region. Individual measurements are denoted as black points.
The strain used for calibration of channel copy number (MLG910) is
highlighted in yellow. ](../figs/fig2.pdf){#fig:boxplot}


### Performing a single-cell hypo-osmotic challenge assay ###

 &nbsp;&nbsp;&nbsp;&nbsp;To measure the channel copy number of a single cell
 and query its survival after a hypo-osmotic shock, we used a custom-made
 flow cell in which osmotic shock and growth can be monitored in real time
 using video microscopy ([@Fig:flow_cell]A). The design and characterization
 of this device has been described in depth previously and is briefly
 described in the Materials & Methods [@bialecka-fornal2015]. Using this
 device, cells were exposed to a large hypo-osmotic shock by switching
 between LB Lennox medium supplemented with 500 mM NaCl and LB Lennox media alone. All six Shine-Dalgarno modifications shown in [@Fig:boxplot]B
 (excluding MLG910) were subjected to a hypo-osmotic shock at controlled
 rates while under observation. After the application of the osmotic shock,
 the cells were imaged every sixty seconds for four to six hours. Each cell
 was monitored over the outgrowth period and was manually scored as either a
 survivor, fatality, or inconclusive observation. The criteria used for scoring death were the same as those previously described in Bialecka-Fornal et al. 2015 [@bialecka-fornal2015]. Survivors were defined as cells that underwent multiple divisions  post-shock. To qualify as survivors, cells must undergo at least two divisions, although more typically,  four to eight divisions are observed without any signs of slowing down. Imaging is stopped when the survivors cells begin to go out of focus or overlap each other.  Survivors do not show any sign of ceasing division. More information regarding this  classification can be found in the Materials and Methods as well as Supplement C (Fig. S9 - S10 and Table S1 - S2). The brief experimental protocol can be seen in  [@Fig:flow_cell]B.

![Experimental approach to measuring survival probability. (A) Layout of
a home-made flow cell for subjecting cells to osmotic shock. Cells are
attached to a polyethylenimine functionalized surface of a glass coverslip
within the flow chamber by loading a dilute cell suspension through one of
the inlets. (B) The typical experimental procedure. Cells are loaded into a
flow chamber as shown in (A) and mounted to the glass coverslip surface.
Cells are subjected to a hypo-osmotic shock by flowing hypotonic medium into
the flow cell. After shock, the cells are monitored for several hours and
surviving cells are identified.](../figs/fig3.pdf){#fig:flow_cell}

&nbsp;&nbsp;&nbsp;&nbsp;Due to the extensive overlap in expression between
the different Shine-Dalgarno mutants (see [@Fig:boxplot]B), computing the
survival probability by treating each mutant as an individual bin obfuscates
the relationship between channel abundance and survival. To more thoroughly
examine this relationship, all measurements were pooled together with each
cell being treated as an individual experiment. The hypo-osmotic shock
applied in these experiments was varied across a range of 0.02 Hz (complete
exchange in 50 s) to 2.2 Hz (complete exchange in 0.45 s). Rather than
pooling this wide range of shock rates into a single data set, we chose to
separate the data into “slow shock” ( &lt; 1.0 Hz) and “fast shock” ($\geq
1.0$ Hz) classes. Other groupings of shock rate were explored and are
discussed in Supplement D (Fig. S11 and S12). The
cumulative distributions of channel copy number separated by survival are
shown in [@Fig:survival_dists]. In these experiments, survival was never
observed for a cell containing less than approximately 100 channels per cell,
indicated by the red stripe in [@Fig:survival_dists]. This suggests that
there is a minimum number of channels needed for survival on the order of 100
per cell. We also observe a slight shift in the surviving fraction of the
cells towards higher effective copy number, which matches our intuition that
including more mechanosensitive channels increases the survival probability.

![Distributions of survival and death as a function of effective channel
number. (A) Empirical cumulative distributions of channel copy number
separated by survival (green) or death (purple) after a slow ($< 1.0$ Hz)
osmotic shock. (B) The empirical cumulative distribution for a fast ($\geq
1.0$ Hz) osmotic shock. Shaded green and purple regions represent the 95%
credible region of the effective channel number calculation for each cell.
Shaded red stripe signifies the range of channels in which no survival was
observed.](../figs/fig4.pdf){#fig:survival_dists}

### Prediction of survival probability as a function of channel copy number

&nbsp;&nbsp;&nbsp;&nbsp;There are several ways by which the survival
probability can be calculated. The most obvious approach would be to group
each individual Shine-Dalgarno mutant as a single bin and compute the average
MscL copy number and the survival probability. Binning by strain is the most
frequently used approach for such measurements and has provided valuable
insight into the qualitative relationship of survival on other physiological
factors [@bialecka-fornal2015; @vandenberg2016]. However the copy number
distribution for each Shine-Dalgarno mutant ([@Fig:boxplot]B) is remarkably wide and
overlaps with the other strains. We argue that this coarse-grained binning
negates the benefits of performing single-cell measurements as two strains
with different means but overlapping quartiles would be treated as distinctly
different distributions. 

&nbsp;&nbsp;&nbsp;&nbsp;Another approach would be to pool all data together,
irrespective of the Shine-Dalgarno mutation, and bin by a defined range of channels.
Depending on the width of the bin, this could allow for finer resolution of
the quantitative trend, but the choice of the bin width is arbitrary with the
*a priori* knowledge that is available. Drawing a narrow bin width can easily
restrict the number of observed events to small numbers where the statistical
precision of the survival probability is lost. On the other hand, drawing
wide bins increases the precision of the estimate, but becomes further
removed from a true single-cell measurement and represents a population mean,
even though it may be a smaller population than binning by the Shine-Dalgarno sequence
alone. In both of these approaches, it is difficult to extrapolate the
quantitative trend outside of the experimentally observed region of channel
copy number. Here, we present a method to estimate the probability of
survival for any channel copy number, even those that lie outside of the
experimentally queried range.

&nbsp;&nbsp;&nbsp;&nbsp;To quantify the survival probability while
maintaining single-cell resolution, we chose to use a logistic regression
model which does not require grouping data into arbitrary bins and treats
each cell measurement as an independent experiment. Logistic regression is
an inferential method to model the probability of a Boolean or categorical
event (such as survival or death) given one or several predictor variables
and is commonly used in medical statistics to compute survival rates and
dose response curves [@anderson2003; @mishra2016]. The primary assumption of
logistic regression is that the log-odds probability of survival $p_{s}$ is
linearly dependent on the predictor variable, in our case the log
channels per cell $N_{c}$ with a dimensionless intercept $\beta_0$ and slope
$\beta_1$,
$$
\log{p_s \over 1 - p_s} = \beta_0 + \beta_1 \log N_c.
$${#eq:linear_channel_logit}
Under this assumption of linearity, $\beta_0$ is the log-odds probability of
survival with no MscL channels. The slope $\beta_1$ represents the change in
the log-odds probability of survival conveyed by a single channel. As the
calculated number of channels in this work spans nearly three orders of
magnitude, it is better to perform this regression on $\log N_c$ as
regressing on $N_c$ directly would give undue weight for lower channel copy
numbers due to the sparse sampling of high-copy number cells. The functional
form shown in [@Eq:linear_channel_logit] can be derived directly from Bayes’
theorem and is shown in Supplement E.
If one knows the values of $\beta_0$ and $\beta_1$, the survival probability
can be expressed as
$$
p_s = \frac{1}{1 + N_c^{-\beta_1}e^{-\beta_0}}.
$${#eq:prob}
In this analysis, we used Bayesian inferential methods to determine the
most likely values of the coefficients and is described in detail in the
Supplement E (Fig. S13 and S14).

&nbsp;&nbsp;&nbsp;&nbsp; The results of the logistic regression are shown in
[@Fig:survival]. We see a slight rightward shift the survival probability
curve under fast shock relative to the slow shock case, reaffirming the
conclusion that survival is also dependent on the rate of osmotic shock
[@bialecka-fornal2015]. This rate dependence has been observed for cells
expressing MscL alongside other species of mechanosensitive channels, but not
for MscL alone. This suggests that MscL responds differently to different
rates of shock, highlighting the need for further study of rate dependence
and the coordination between different species of mechanosensitive channels.
[@Fig:survival] also shows that several hundred channels are required to
provide appreciable protection from osmotic shock. For a survival probability of 80%, a cell must have
approximately 500 to 700 channels per cell for a fast and slow shock,
respectively. The results from the logistic regression are showed as continuous colored curves. The individual
cell measurements separated by survival and death are shown at the top and
bottom of each plot, respectively, and are included to provide a sense of
sampling density. 

&nbsp;&nbsp;&nbsp;&nbsp;Over the explored range of MscL copy number, we observed a
maximum of 80\% survival for any binning method. The remaining 20\% survival
may be attained when the other species of mechanosensitive channels are
expressed alongside MscL. However, it is possible that the flow cell method
performed in this work lowers the maximal survival fraction as the cells are
exposed to several, albeit minor, mechanical stresses such as loading into
the flow cell and chemical adherence to the glass surface. To ensure that the
results from logistic regression accurately describe the data, we can compare
the survival probabilities to those using the binning methods described
earlier (red and black points, [@Fig:survival]). Nearly all binned data fall
within error of the prediction (see Materials & Methods for definition of
error bar on probability), suggesting that this approach accurately reflects
the survival probability and gives license to extrapolate the estimation of
survival probability to regions of outside of our experimentally explored
copy number regime.

&nbsp;&nbsp;&nbsp;&nbsp;Thus far, we’ve dictated that for a given rate of
osmotic shock (i.e. "fast" or "slow"), the survival probability is dependent
only on the number of channels. In Fig. S13, we show the result of including
other predictor variables, such as area and shock rate alone. In such cases,
including other predictors resulted in pathological curves showing that
channel copy number is the most informative out of the available predictor
variables.

![**Probability of survival as a function of MscL copy number.** (A)
Estimated survival probability for survival under slow shock as a function of
channel copy number. (B) The estimated survival probability of survival under
a fast shock as a function of channel copy number. Solid curves correspond to
the most probable survival probability from a one-dimensional logistic
regression. Shaded regions represent the 95\% credible regions. Points at the
top and bottom of plots represent individual cell measurements which survived
and perished, respectively. The red and black points correspond to the survival
probability estimated via binning by Shine-Dalgarno sequence and binning by groups of 50
channels per cell, respectively. Horizontal error bars represent the standard
error of the mean from at least 25 measurements. Vertical error bars
represent the certainty of the probability estimate given $n$ survival events
from $N$ total observations.](../figs/fig5.pdf){#fig:survival}

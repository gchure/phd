### Experimental Setup

Seminal studies from the burgeoning years of bacterial physiology have
demonstrated a strong dependence of the total cellular protein content on the
growth rate [@schaechter1958], a relationship which has been rigorously
quantified in recent years using mass spectrometry [@schmidt2016] and
ribosomal profiling [@li2014a]. Their combined results illustrate that
modulation of the growth rate, either through controlling the available
carbon source or the temperature of the growth medium, significantly alters
the physiological state of the cell, triggering the reallocation of resources
to prioritize expression of ribosome-associated genes. @Eq:foldchange_simple
has no explicit dependence on the available carbon source but does depend on
the temperature through the energetic parameters $\Delta\varepsilon_{R}$ and
$\Delta\varepsilon_{AI}$ which are defined relative to the thermal energy,
$k_BT$. With this parametric knowledge, we are able to draw quantitative
predictions of the fold-change in gene expression in these physiologically
distinct states.

We modulated growth of *Escherichia coli* by varying either the quality of
the available carbon source (differing ATP yield per C atom) or the
temperature of the growth medium [@Fig:growth_control (A)]. All experiments
were performed in a defined M9 minimal medium supplemented with one of three
carbon sources – glucose, glycerol, or acetate – at concentrations such that
the total number of carbon atoms available to the culture remained the same.
These carbon sources have been shown to drastically alter growth rate and
gene expression profiles , indicating changes in the proteomic composition
and distinct physiological states. These carbon sources yield an approximate
four-fold modulation of the growth rate with doubling times ranging from
$\approx$ 220 minutes to $\approx$ 65 minutes in an acetate or glucose
supplemented medium, respectively [@Fig:growth_control (B) and (C)]. While
the growth temperature was varied over 10$^\circ$ C, both 32$^\circ$ and
42$^\circ$ C result in approximately the same doubling time of $\approx$ 90
min, which is 1.5 times slower than the optimal temperature of 37$^\circ$ C
[@Fig:growth_control (B) and (C)\].

![**Control of physiological state via growth rate through environmental
factors.** (A) Bacterial growth can be controlled by varying the available
carbon source (top panel) or temperature (bottom panel). (B) Bulk bacterial
growth curves under all conditions illustrated in (A). The y-axis is the
optical density measurements at 600 nm relative to the initial value.
Interval between points is $\approx$ 6 min. Points and errors represent the
mean and standard deviation of three to eight biological replicates. (C)
Inferred maximum growth rate of each biological replicate for each condition.
Points represent the doubling time computed from the maximum growth rate.
Error bars correspond to the standard deviation of the inferred growth rate.
Where not visible, error bars are smaller than the
marker.](ch4_fig1){#fig:growth_control short-caption="Controlling
physiological state of *E. coli* via growth rate modulation by environmental
factors."}

The growth rate dependence of the proteome composition suggests that
changing physiological conditions could change the total repressor copy
number of the cell. As shown in our previous work , it can be difficult
to differentiate between a change in repressor copy number $R$ and the
allosteric energy difference $\Delta\varepsilon_{AI}$ as there are
many combinations of parameter values that yield the same fold-change.
To combat this degeneracy, we used a genetically engineered strain of
*E. coli* in which the expression of the repressor copy number and its
regulated gene product (YFP) can be simultaneously measured. This
strain, used previously to interrogate the transcription factor
titration effect [@brewster2014], is diagrammed in 
@Fig:dilution_circuit (A). A dimeric form of the LacI repressor
N-terminally tagged with an mCherry fluorophore is itself regulated
through the action of the TetR repressor whose level of activity can be
modulated through the addition of the allosteric effector anhydrous
tetracycline (ATC). This dual repression genetic circuit allows for the
expression of the LacI repressor to be tuned over several orders of
magnitude. This is demonstrated in Fig.
@Fig:dilution_circuit (B) where a titration of ATC in the
growth medium results in a steady increase in the expression of the
LacI-mCherry gene product (red lines and points) which in turn represses
expression of the YFP reporter (yellow lines and points).

While the mCherry fluorescence is proportional to the repressor copy number,
it is not a direct measurement as the fluorescence of a single LacI-mCherry
dimer is unknown *a priori*. Using video microscopy, we measure the
partitioning statistics of the fluorescence intensity into two sibling cells
after division [ @Fig:dilution_circuit (C)]. This method, described in detail
in the Materials \& Methods as well as in @rosenfeld2002, @rosenfeld2005, and
@brewster2014 reveals a linear relationship between the variance in intensity
between two sibling cells and the intensity of the parent cell, the slope of
which is equal to the brightness of a single LacI repressor. Since this
measurement is performed simultaneously with measurement of the expression of
the YFP reporter, this calibration factor was determined for each unique
experimental replicate. We direct the reader to the supplemental Chapter 8
for a more thorough discussion of this inference.


![**Control and quantification of repressor copy number. (A) The dual
repression expression system.** The inducible repressor TetR (purple blob) is
expressed from a low-copy-number plasmid in the cell and represses expression
of the LacI-mCherry repressor by binding to its cognate operator (tetO). In
the presence of anhydrous tetracycline (ATC, green sphere), the inactive
state of TetR becomes energetically favorable, permitting expression of the
LacI-mCherry construct (red). This in turn binds to the lacO operator
sequence repressing the expression of the reporter Yellow Fluorescent Protein
(YFP, yellow lightbulb). (B) An ATC titration curve showing anticorrelated
YFP (yellow) and mCherry (red) intensities. Reported values are scaled to
the maximum mean fluorescence for each channel. Points and errors correspond
to the mean and standard error of eight biological replicates. (C)
Determination of a fluorescence calibration factor. After cessation of
LacI-mCherry expression, cells are allowed to divide, partitioning the
fluorescently tagged LacI repressors into the two daughter cells (left
panel). The total intensity of the parent cell is equivalent to the summed
intensities of the daughters. The squared fluctuations in intensity of the
two sibling cells is linearly related to the parent cell with a slope $\alpha$,
which is the fluorescence signal measured per partitioned repressor (right
panel). Black points represent single divisions and red points are the means
of 50 division events. Line corresponds to linear fit to the black points
with a slope of $\alpha = 740 \pm 40$ a.u. per LacI.](ch4_fig2){#fig:dilution_circuit
short-caption="Control and quantification of repressor copy number through the
binomial partitioning method."}

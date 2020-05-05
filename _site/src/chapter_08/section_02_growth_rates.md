## Non-parametric Inference of Growth Rates

In this section, we discuss the measurement of the bacterial growth
curves as well as our strategy for estimating the growth rates for all
experimental conditions in this work.

### Experimental growth curves

As is described in the Materials and Methods section of Chapter 4,
we measured the growth of *E. coli* strains constitutively expressing
YFP using a BioTek Cytation 5 96-well plate reader generously provided
by the lab of Prof. David Van Valen at Caltech. Briefly, cells were
grown overnight in a nutrient rich LB medium to saturation and were
subsequently diluted 1000 fold into the appropriate growth medium. Once
these diluted cultures reached an OD$_{600nm}$ of $\approx 0.3$, the
cells were again diluted 100 fold into fresh medium preheated to the
appropriate temperature . Aliquots of 300 $\mu$L of this dilution were
then transferred to a 96-well plate leaving two-rows on all sides of the
plate filled with sterile media to serve as blank measurements. Once
prepared, the plate was transferred to the plate reader and OD$_{600nm}$
measurements were measured every $\approx$ 5 - 10 min for 12 to 18
hours. Between measurements, the plates were agitated with a linear
shaking mode to avoid sedimentation of the culture. A series of
technical replicates of the growth curve in glycerol supplemented medium
at 37$^\circ$ C is shown in @Fig:gp_growth_curve (A).

### Inference of maximum growth rate

The phenomenon of collective bacterial growth has been the subject of
intense research for the better part of a century
[@schaechter1958; @jun2018] yielding many parametric descriptions of the
bulk growth rates each with varying degrees of detail
[@jun2018; @allen2018]. For this scope of this work, we are not
particularly interested in estimating numerous parameters that describe
the phenomenology of the growth curves. Rather, we are interested in
knowing a single quantity -- the maximum growth rate -- and the degree
to which it is tuned across the different experimental conditions. To
avoid forcing the bacterial growth curves into a parametric form, we
treated the observed growth curves using Gaussian process modeling using
as implemented in the FitDeriv software described in Ref. [@swain2016].
This method permits an estimation of the most-likely OD$_{600nm}$ value
at each point in time given knowledge of the adjacent measurements. The
weighting given to all points in the series of measurements is defined
by the covariance kernel function and we direct the reader to Ref.
[@swain2016] for a more detailed discussion on this kernel choice and
overall implementation of Gaussian process modeling of time-series data.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As this approach estimates the probability of a given OD$_{600nm}$ at
each time point, we can compute the mean and standard deviation. @Fig:gp_growth_curve (B) show the raw measurements and the
mean estimated value in green and light, respectively. With the high
temporal resolution of the OD$_{600nm}$ measurements, modeling the
entire growth curve becomes a computationally intensive task.
Furthermore, as we are interested in only the maximum growth rate, there
is no benefit in analyzing the latter portion of the experiment where
growth slows and the population reaches saturation. We therefore
manually restricted each analyzed growth curve to a region capturing
early and mid exponential phase growth, illustrated by the shaded region
in @Fig:gp_growth_curve (A).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With a smooth description of the OD$_{600nm}$
measurements as a function of time with an appropriate measure of
uncertainty, we can easily compute the time derivative which is equivalent to
the growth rate. A representative time derivative is shown in
@Fig:gp_growth_curve (C). Here, the dark green curve is the mean value of
the time derivative and the shaded region is $\pm$ one standard deviation.
The maximum value of this inferred derivative is the reported maximum growth
rate of that experimental condition.

![**Non-parametric characterization of bacterial growth curves and estimation of
the growth rate.** (A) Representative biological replicate growth curve of a
$\Delta$*lacIZYA E. coli* strain in glycerol supplemented M9 minimal medium at
37$^\circ$ C. Points represent individual optical density measurements across
eight technical replicates. The shaded region illustrates the time window from
which the maximum growth rate was inferred. (B) Green points are those from the
shaded region in (A).  Line is the estimated value of the optical density
resulting from the gaussian process modeling. (C) The value of the first
derivative of the optical density as a function of time estimated via gaussian
process modeling. The first derivative is taken as the growth rate at each time
point. The solid line is the estimated mean value and the shaded region
represents the standard deviation of the posterior
The [Python code                                                
(`ch8_figS1.py`)](https://github.com/gchure/phd/blob/master/src/chapter_08/code/ch8_figS1.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).
distribution.](ch8_figS1){#fig:gp_growth_curve short-caption="Non-parametric
estimation of maximum growth rate from bacterial growth curves."}

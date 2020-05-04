
## Classification Of Shock Rate

&nbsp;&nbsp;&nbsp;&nbsp;Its been previously shown that the rate of
hypo-osmotic shock dictates the survival probability [@bialecka-fornal2015].
To investigate how a single channel contributes to survival, we queried
survival at several shock rates with varying MscL copy number. In the main
text of this work, we separated our experiments into arbitrary bins of "fast"
($\geq$ 1.0 Hz) and "slow" ($<$ 1.0 Hz) shock rates. In this section, we
discuss our rationale for coarse graining our data into these two groupings.

&nbsp;&nbsp;&nbsp;&nbsp;As is discussed in Chapter 5, we used a bin-free method of
estimating the survival probability given the MscL channel copy number as a
predictor variable. While this method requires no binning of the data, it
requires a data set that sufficiently covers the physiological range of
channel copy number to accurately allow prediction of survivability.
@Fig:indiv_shock_group shows the results of the logistic regression
treating each shock rate as an individual data set. The most striking feature
of the plots shown in @Fig:indiv_shock_group is the inconsistent behavior
of the predicted survivability from shock rate to shock rate. The appearance
of bottle necks in the credible regions for some shock rates (0.2Hz, 0.5Hz,
2.00Hz, and 2.20 Hz) appear due to a high density of measurements within a
narrow range of the channel copy number at the narrowest point in the bottle
neck. While this results in a seemingly accurate prediction of the survival
probability at that point, the lack of data in other copy number regimes
severely limits our extrapolation outside of the copy number range of that
data set. Other shock rates (0.018 Hz, 0.04 Hz, and 1.00 Hz) demonstrate
completely pathological survival probability curves due to either complete
survival or complete death of the population.

&nbsp;&nbsp;&nbsp;&nbsp;Ideally, we would like to have a wide range of MscL
channel copy numbers at each shock rate shown in @Fig:indiv_shock_group.
However, the low-throughput nature of these single-cell measurements
prohibits completion of this within a reasonable time frame. It is also
unlikely that thoroughly dissecting the shock rate dependence will change the
overall finding from our work that several hundred MscL channels are needed
to convey survival under hypo-osmotic stress.

![**Binning by individual shock rates.** Survival probability estimates from
logistic regression (red lines) and the computed survival probability for all
SD mutants subjected to that shock rate (blue points). Black points at top
and bottom of each plot correspond to single cell measurements of survival
(top) and death (bottom). Red shaded regions signify the 95\% credible region
of the logistic regression. Horizontal error bars of blue points are the
standard error of the mean channel copy number. Vertical error bars of blue
points correspond to the uncertainty in survival probability by observing $n$
survival events from $N$ single-cell
measurements. The [Python code                                                
(`ch9_figS11.py`)](https://github.com/gchure/phd/blob/master/src/chapter_09/code/ch9_figS11.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch9_figS11){#fig:indiv_shock_group short-caption="Binning by
individual shock rates."}


&nbsp;&nbsp;&nbsp;&nbsp; Given the data shown in @Fig:indiv_shock_group, we
can try to combine the data sets into several bins. @Fig:three_shock_group
shows the data presented in @Fig:indiv_shock_group separated into "slow"
($<$ 0.5 Hz, A), "intermediate" (0.5 - 1.0 Hz, B), and "fast" ($>$ 1.0 Hz, C)
shock groups. Using these groupings, the full range of MscL channel copy
numbers are covered for each case, with the intermediate shock rate sparsely
sampling copy numbers greater than 200 channels per cell. In all three of
these cases, the same qualitative story is told -- several hundred channels
per cell are necessary for an appreciable level of survival when subjected to
an osmotic shock. This argument is strengthened when examining the predicted
survival probability by considering all shock rates as a single group, shown
in @Fig:three_shock_group (D). This treatment tells nearly the same
quantitative and qualitative story as the three rate grouping shown in this
section and the two rate grouping presented in the main text. While there
does appear to be a dependence on the shock rate for survival when only MscL
is expressed, the effect is relatively weak with overlapping credible regions
for the logistic regression across the all curves. To account for the sparse
sampling of high copy numbers observed in the intermediate shock group, we
split this set and partitioned the measurements into either the "slow" ($<$
1.0 Hz) or "fast" ($\geq$ 1.0 Hz) groups presented in the main text of this
work.

![**Coarse graining shock rates into different groups.** Estimated survival
probability curve for slow (A), intermediate (B), and fast (C) shock rates.
(D) Estimated survival probability curve from pooling all data together,
ignoring varying shock rates. Red shaded regions correspond to the 95\%
credible region of the survival probability estimated via logistic
regression. Black points at top and bottom of each plot represent single-cell
measurements of cells which survived and died, respectively. Black points and
error bars represent survival probability calculations from bins of 50
channels per cell. Blue points represent the survival probability for a given
Shine-Dalgarno mutant. Horizontal error bars are the standard error of the
mean with at least 25 measurements and vertical error bars signifies the
uncertainty in the survival probability from observing $n$ survival events
out of $N$ total measurements.](ch9_figS12){#fig:three_shock_group
short-caption="Coarse graining shock rates into different groups."}

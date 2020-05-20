## Comparison of Survival Probability with van den Berg et al. (2016)

&nbsp;&nbsp;&nbsp;&nbsp;In @vandenberg2016, the authors report a
100% survival rate at approximately 100 channels per cell. While the number
of mechanosensitive channels per cell was quantified at the level of single
cells, the survival probability was measured in bulk using ensemble plating
assays. The results of these experiments considering the contribution of MscL
to survival is shown in Figure 5 of their work, although without displayed
uncertainty in the survival probability. Figure S6B of their work shows the
approximate error in survival probability through ensemble plating assays for
three different strains (@Fig:poolman_comparison (A)), which is approximately
30%. Using this approximate error and the data shown in their Figure 5B, we
have reproduced this plot with error bars in both measured dimensions
(@Fig:poolman_comparison (B)). This plot shows that even when the mean
survival probability is 100\%, the variation in the measured survival
probability is large, extending as low as $\approx$ 70\%. This variation is likely born
from a multitude of experimental steps including time of outgrowth, variation
in shock rate, plating efficiency, and counting errors. As our experimental
approach directly measures the survival/death of individual cells, we remove
many sources of error that would arise from an ensemble approach, albeit at
lower throughput. While it is possible that the discrepancy between @vandenberg2016
and the work presented in Chapter 5 could arise from other unknown
factors, we believe that single-cell experiments introduce the fewest sources
of error.

![**MscL abundance vs survival data reported in @vandenberg2016 with
included error.** (A) Reported survival probabilities of a strain lacking all
mechanosensitive channels (“no plasmid”), plasmid borne MscL-mEos3.2, and
plasmid borne MscS-mEos3.2. Approximate reported errors for MscL-mEos3.2
survival probability is 30\%. (B) The measurement of survival probability as a
function of MscL channel copy number was obtained from Figure 5B in @vandenberg2016. Errors in channel copy number represent the standard
deviation of several biological replicates (present in original figure) while
the error in survival probability is taken as  $\approx$
30\%. The [Python code (`ch9_figS12.py`)](https://github.com/gchure/phd/blob/master/src/chapter_09/code/ch9_figS12.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch9_figS12){#fig:poolman_comparison
short-caption="MscL abundance versus survival data reported in van den Berg et
al. 2016 with included error."}

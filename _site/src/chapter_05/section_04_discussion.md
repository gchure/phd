
## Discussion

&nbsp;&nbsp;&nbsp;&nbsp;One of the most challenging endeavors in the
biological sciences is linking the microscopic details of cellular components
to the macro-scale physiology of the organism. This formidable task has been
met repeatedly in the recent history of biology, especially in the era of DNA
sequencing and single molecule biochemistry. For example, the scientific
community has been able to connect sickle-cell anemia to a single amino acid
substitution in Hemoglobin which promotes precipitation under a change in
O$_2$ partial pressure [@feeling-taylor2004; @finch1973; @perutz1950]. Others
have assembled a physical model that quantitatively describes chemosensation
in bacteria [@berg1977] in which the arbiter of sensory adaptation is the
repeated methylation of chemoreceptors [@colin2017; @krembel2015a;
@krembel2015; @sourjik2002b]. In the past ~50 years alone, numerous biological
and physical models of the many facets of the central dogma have been
assembled that give us a sense of the interplay between the genome and
physiology. For example, the combination of biochemical experimentation and
biophysical models have given us a picture of how gene dosage affects furrow
positioning in *Drosophila* [@liu2013], how recombination of V(D)J gene
segments generates an extraordinarily diverse antibody repertoire
[@lovely2015; @schatz2004; @schatz2011], and how telomere shortening through
DNA replication is intrinsically tied to cell senescence [@herbig2004;
@victorelli2017], to name just a few of many such examples.

&nbsp;&nbsp;&nbsp;&nbsp; By no means are we ``finished” with any of these
topics. Rather, it's quite the opposite in the sense that having a handle on
the biophysical knobs that tune the behavior opens the door to a litany of
new scientific questions. In the case of mechanosenstaion and osmoregulation,
we have only recently been able to determine some of the basic facts that
allow us to approach this fascinating biological phenomenon biophysically.
The dependence of survival on mechanosensitive channel abundance is a key
quantity that is missing from our collection of critical facts. To our
knowledge, this work represents the first attempt to quantitatively control
the abundance of a single species of mechanosensitive channel and examine the
physiological consequences in terms of survival probability at single-cell
resolution. Our results reveal two notable quantities. First, out of the
several hundred single-cell measurements, we never observed a cell which had
less than approximately 100 channels per cell and survived an osmotic shock,
irrespective of the shock rate. The second is that between 500 and 700
channels per cell are needed to provide $\geq 80\%$ survival, depending on
the shock rate.

&nbsp;&nbsp;&nbsp;&nbsp;Only recently has the relationship between the MscL
copy number and the probability of survival been approached experimentally.
In @vandenberg2016, the authors examined the contribution of MscL
to survival in a genetic background where all other known mechanosensitive
channels had been deleted from the chromosome and plasmid-borne expression of
an MscL-mEos3.2 fusion was tuned through an IPTG inducible promoter
[@vandenberg2016]. In this work, they measured the single-cell channel
abundance through super-resolution microscopy and queried survival through
bulk assays. They report a nearly linear relationship between survival and
copy number, with approximately 100 channels per cell conveying 100%
survival. This number is significantly smaller than our observation of
approximately 100 channels as the *minimum* number needed to convey an  y
observable degree of survival.


The disagreement between the numbers reported in this work and in van den
Berg et al. may partially arise from subtle differences in the experimental
approach. The primary practical difference is the magnitude of the osmotic
shock. van den Berg et al. applied an approximately 600 mOsm downshock in
bulk whereas we applied a 1 Osm downshock, which would lead to lower survival
[@levina1999]. In their work, the uncertainty in both the MscL channel count
and survival probability is roughly 30% (Fig. S14). Given this uncertainty,
it is reasonable to interpret that the number of channels needed for complete
protection from osmotic downshock is between 100 and 250 per cell. The
uncertainty in determining the number of channels per cell is consistent with
the observed width of the channel number distribution of the Shine-Dalgarno
sequence mutants used in this work [Fig. @fig:mscl_boxplot(B)]. A unique
property of the single-cell measurements performed in this work allow is the
direct observation of survival or death of individual cells. We find that
morphological classification and classification through a propidium iodide
staining agree within 1% (Supplement C). The bulk plating assays, as are used
in van den Berg et al., rely on colony formation and outgrowth to determine
survival probability. As is reported in their supplemental information, the
precision in this measurement is around 30% (Fig. S14). Accounting for this
uncertainty brings both measurements within a few fold where we still
consistently observe lower survival for a given channel number. This
remaining disagreement may be accounted for by systematic uncertainty in both
experimental methods.

For example, variation in the length of outgrowth, variable shock rate, and
counting statistics could bias towards higher observed survival rates in
ensemble plating assays. During the outgrowth phase, the control sample not
exposed to an osmotic shock is allowed to grow for approximately 30 minutes
in a high-salt medium before plating. The shocked cells, however, are allowed
to grow in a low-salt medium. We have found that the difference between the
growth rates in these two conditions can be appreciable (approximately 35
minutes versus 20 minutes, respectively) as can be seen in Fig S2. Cells that
survived an osmotic shock may have a growth advantage relative to the control
sample if the shock-induced lag phase is less than the outgrowth, leading to
higher observed survival rates [@levina1999]. This is one possible
explanation for the survival rates which are reported in excess of 100%.
Cells that survived an osmotic shock may have a growth advantage relative to
the normalization sample if the shock-induced lag phase is less than the
outgrowth, leading to higher observed survival rates, even surpassing 100%.
We have performed these assays ourselves and have observed survival rates
above of 100% (ranging from 110% to 125%) with an approximate 30% error (see
Fig. S3 in Bialecka-Fornal et al. 2012 [@bialecka-fornal2012]) which we
concluded to arise from differences in growth rate. We also note that
survival rates greater than 100% are observed in van den Berg et al. (Fig.
S14). For strains that have survival rates between 80% and 100% the
uncertainty is typically large, making it difficult to make precise
statements regarding when full survival is achieved.
 
It has been shown that there is a strong inverse relationship between the
rate of osmotic shock and survival probability [@bialecka-fornal2015]. Any
experiment in which the shock was applied more slowly or quickly than another
would bias toward higher or lower survivability, respectively. The shocks
applied in bulk assays are often performed manually which can be highly
variable. We note that in our experiments, we frequently observe cells which
do not separate and form chains of two or more cells (Fig. S9 and S10). In
plating assays, it is assumed that colonies arise from a single founding cell
however a colony formed by a cluster of living and dead cells would be
interpreted as a single surviving cell, effectively masking the death of the
others in the colony forming unit. This too could bias the measurement toward
higher survival rates. Single-cell shock experiments can also have systematic
errors which can bias the results towards lower survival rates. Such errors
are associated with handling of the cells such as shear damage from loading
into the flow cell, adhering the cells to the coverslip, and any chemical
perturbations introduced by the dye used to measure the shock rate.

Despite these experimental differences, the results of this work and van den
Berg et al., are in agreement that MscL must be present at the level of 100
or more channels per cell in wild-type cells to convey appreciable survival.
As both of these works were performed in a strain in which the only
mechanosensitive channel was MscL, it remains unknown how the presence of the
other channel species would alter the number of MscL needed for complete
survival. In our experiments, we observed a maximum survival probability of
approximately 80\% even with close to 1000 MscL channels per cell. It is
possible that the combined effort of the six other mechanosensitive channels
would make up for some if not all of the remaining 20\%. To explore the
contribution of another channel to survival, van den Berg et al. also queried
the contribution of MscS, another mechanosensitive channel, to survival in
the absence of any other species of mechansensitive channel. It was found
that over the explored range of MscS channel copy numbers, the maximum
survival rate was approximately 50\%, suggesting that different
mechanosensitive channels have an upper limit to how much protection they can
confer. Both van den Berg et al. and our work show that there is still much
to be learned with respect to the interplay between the various species of
mechanosensitive channel as well as their regulation.

Recent work has shown that both magnitude and the rate of osmotic down shock
are important factors in determining cell survival [@bialecka-fornal2015]. In
this work, we show that this finding holds true for a single species of
mechanosensitive channel, even at high levels of expression. One might
naïvely expect that this rate-dependent effect would disappear once a certain
threshold of channels had been met. Our experiments, however, show that even
at nearly 1000 channels per cell the predicted survival curves for a slow ($<
1.0$ Hz) and fast ($\geq 1.0$ Hz) are shifted relative to each other with the
fast shock predicting lower rates of survival. This suggests either we have
not reached this threshold in our experiments or there is more to understand
about the relationship between abundance, channel species, and the shock
rate.

Some experimental and theoretical treatments suggest that only a few copies
of MscL or MscS should be necessary for 100% protection given our knowledge
of the conductance and the maximal water flux through the channel in its open
state [@louhivuori2010; @booth2014]. However, recent proteomic studies have
revealed average MscL copy numbers to be in the range of several hundred per
cell, depending on the condition, as can be seen in Table 1 [@li2014a;
@schmidt2016; @soufi2015]. Studies focusing solely on MscL have shown similar
counts through quantitative Western blotting and fluorescence microscopy
[@bialecka-fornal2012]. Electrophysiology studies have told another story
with copy number estimates ranging between 4 and 100 channels per cell
[@blount1999; @stokes2003a; @booth2005]. These measurements, however, measure
the active number of channels. The factors regulating channel activity in
these experiments could be due to perturbations during the sample preparation
or reflect some unknown mechanism of regulation, such as the presence or
absence of interacting cofactors [@schumann2010]. The work described here, on
the other hand, measures the *maximum* number of channels that could be
active and may be able to explain why the channel abundance is higher than
estimated by theoretical means. There remains much more to be learned about
the regulation of activity in these systems. As the *in vivo* measurement of
protein copy number becomes accessible through novel single-cell and
single-molecule methods, we will continue to collect more facts about this
fascinating system and hopefully connect the molecular details of
mechanosensation with perhaps the most important physiological response --
life or death.

| Reported channels per cell | Method  | Reference |
|:---|:---|:---:|
| 480 $\pm$ 103 | Western blotting   | @bialecka-fornal2012 |
|  560\* |  Ribosomal profiling  | @li2014a|
| 331\* | Mass spectrometry  | @schmidt2016|
| 583\* | Mass spectrometry  | @soufi2015|
| 4 - 5 | Electrophysiology  | @stokes2003a|
| 10 - 100 | Electrophysiology  | @booth2005|
|10 - 15 | Electrophysiology | @blount1999|
Table: Measured cellular copy numbers of MscL. Asterisk (\*) Indicates
inferred MscL channel copy number from the total number of detected MscL
peptides.

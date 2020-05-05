
## Introduction

Changes in the extracellular osmolarity can be a fatal event for the bacterial
cell. Upon a hypo-osmotic shock, water rushes into the cell across
the membrane, leaving the cell with no choice but to equalize the pressure.
This equalization occurs either through damage to the cell membrane
(resulting in death) or through the regulated flux of water molecules through
transmembrane protein channels (Fig 1A). Such proteinaceous pressure release
valves have been found across all domains of life, with the first bacterial
channel being described in 1987 [@martinac1987]. Over the past thirty years,
several more channels have been discovered, described, and (in many cases)
biophysically characterized. *E. coli*, for example, has seven of these
channels (one MscL and six MscS homologs) which have varied conductance,
gating mechanisms, and expression levels. While they have been the subject of
much experimental and theoretical dissection, much remains a mystery with
regard to the roles their abundance and interaction with other cellular
processes play in the greater context of physiology [@bavi2016;
@bialecka-fornal2012; @bialecka-fornal2015; @edwards2012; @naismith2012;
@ursell2008; @vandenberg2016].

&nbsp;&nbsp; &nbsp; &nbsp;Of the seven channels in *E. coli*, the
mechanosensitive channel of large conductance (MscL) is one of the most
abundant and the best characterized. This channel has a large conductance (3
nS) and mediates the flux of water molecules across the membrane via a ~3 nm
wide pore in the open state [@cruickshank1997; @haswell2011]. Molecular
dynamics simulations indicate that a single open MscL channel permits the
flux of $4 \times 10^9$ water molecules per second, which is an order of
magnitude larger than a single aquaporin channel (BNID 100479)
[@louhivuori2010; @milo2010]. This suggests that having only a few channels
per cell could be sufficient to relieve even large changes in membrane
tension. Electrophysiological experiments have suggested a small number 
of channels per cell [@booth2005; @hase1997], however, more recent approaches
using quantitative Western blotting, fluorescence microscopy, and proteomics
have measured several hundred MscL per cell [@bialecka-fornal2012;
@schmidt2016; @soufi2015]. To further complicate matters, the expression
profile of MscL appears to depend on growth phase, available carbon source,
and other environmental challenges [@bialecka-fornal2012, @schmidt2016;
@soufi2015; @stokes2003a]. While there are likely more than just a few
channels per cell, why cells seem to need so many and the biological
rationale behind their condition-dependent expression both remain a mystery.

&nbsp;&nbsp;&nbsp; &nbsp; While their biochemical and biophysical
characteristics have received much attention, their connection to cell survival is
understudied. Drawing such a direct connection between channel copy number
and survival requires quantitative *in vivo* experiments. To our knowledge,
the work presented in van den Berg et al. 2016 [@vandenberg2016] is the first
attempt to simultaneously measure channel abundance and survivability for a
single species of mechanosensitive channel. While the measurement of channel
copy number was performed at the level of single cells using super-resolution
microscopy, survivability after a hypo-osmotic shock was assessed in bulk
plating assays which rely on serial dilutions of a shocked culture followed
by counting the number of resulting colonies after incubation. Such bulk
assays have long been the standard for querying cell viability after an
osmotic challenge. While they have been highly informative, they reflect only
the mean survival rate of the population, obfuscating the variability in
survival of the population. The stochastic nature of gene expression results
in a noisy distribution of MscL channels rather than a single value, meaning
those found in the long tails of the distribution have quite different
survival rates than the mean but are lost in the final calculation of
survival probability.

&nbsp; &nbsp; &nbsp; &nbsp;In this work, we present an experimental system to
quantitatively probe the interplay between MscL copy number and survival at
single-cell resolution, as is seen in @Fig:overview(B). We generated an *E.
coli* strain in which all seven mechanosensitive channels had been deleted
from the chromosome followed by a chromosomal integration of a single gene
encoding an MscL-super-folder GFP (sfGFP) fusion protein. To explore copy
number regimes beyond those of the wild-type expression level, we modified
the Shine-Dalgarno sequence of this integrated construct, allowing us to
cover nearly three decades of MscL copy number. To probe survivability, we
exposed cells to a large hypo-osmotic shock at controlled rates in a flow
cell under a microscope, allowing the observation of the single-cell channel
copy number and the resulting survivability of single cells. With this large
set of single cell measurements, we approach the calculation of survival
probability in a manner that is free of binning bias which allows the
reasonable extrapolation of survival probability to copy numbers outside of
the observed range. In addition, we show that several hundred channels are
needed to convey high rates of survival and observe a minimum number of
channels needed to permit any degree of survival.

![**Role of mechanosensitive channels during hypo-osmotic shock.** (A) A
hypo-osmotic shock results in a large difference in the osmotic strength
between the intracellular and extracellular spaces. As a result, water rushes
into the cell to equalize this gradient increasing the turgor pressure and
tension in the cell membrane. If no mechanosensitive channels are present and
membrane tension is high (left panel), the membrane ruptures releasing
intracellular content into the environment resulting in cell death. If
mechanosensitive channels are present (right panel) and membrane tension is
beyond the gating tension, the mechanosensitive channel MscL opens, releasing
water and small intracellular molecules into the environment thus relieving
pressure and membrane tension. (B) The experimental approach undertaken in
this work. The number of mechanosensitive channels tagged with a fluorescent
reporter is tuned through modification of the Shine-Dalgarno sequence of the
*mscL* gene. The cells are then subjected to a hypo-osmotic shock and the
number of surviving cells are counted, allowing the calculation of a survival
probability.](ch5_fig1){#fig:overview short-caption="Role of mechanosensitive
channels during hypo-osmotic shock."}

## Experimental validation of MscL-sfGFP
Despite revolutionizing modern cell biology, tagging proteins with
fluorophores can lead to myriad deleterious effects such as mislocalization,
abrogation of functionality, or even cytotoxicity. In this section, we
examine the stability and functionality of the MscL-sfGFP construct used in
this work.

### Comparing functionality of wild-type and fluorescently tagged MscL
&nbsp;&nbsp;&nbsp;&nbsp; To quantitatively compare the functionality between
the wild-type MscL and MscL-sfGFP, patch-clamp electrophysiology experiments
were conducted on each channel. Patch-clamp recordings were performed 
on membrane patches derived from giant protoplasts which were prepared as
previously described [@blount1999]. In brief, cells were grown in
Luria-Bertani (LB) medium with 0.06 mg/ml cephalexin for 2.5 hours. The
elongated cells were then collected by centrifugation and digested by 0.2
mg/ml lysozyme to form giant protoplasts.

&nbsp;&nbsp;&nbsp;&nbsp; Excised, inside-out patches were analyzed at a
membrane potential of -20 mV with pipette and bath solutions containing 200
mM KCl, 90 mM MgCl$_2$, 10 mM CaCl$_2$, and 5 mM HEPES buffer at pH 7. All
data were acquired at a sampling rate of 50 kHz with 5 kHz filtration using
an AxoPatch 200B amplifier and pClamp software (Molecular Devices). The
pressure threshold for activation a single MscS channel (blue stripe in
[@Fig:ephys]) was compared to that of single MscL channels (orange strip in
[@Fig:ephys]). The pressure threshold for activation of the MscL channels was
referenced against the activation threshold of MscS to determine the pressure
ratio (PL:PS) for gating as previously described [@blount1996]. Recordings of
the transmembrane current were made of three individual patches with an
average PL:PS ratio of 1.56 for MscL-sfGFP. This ratio quantitatively agrees
with the PL:PS ratio of 1.54 measured in a strain (MJF429 from the Booth
laboratory) which expresses the wild-type MscL protein from the chromosome.
The average transient current change from MscL openings (@Fig:ephys shaded
orange region) is 75 pA, corresponding to a single channel conductance of
3.7 nS, comparable to the reported values of wild-type MscL. The agreement
between these two strains indicates that there is negligible difference in
functionality between MscL and MscL-sfGFP, allowing us to make physiological
conclusions of the wild-type channel from our experiments.

![**Characteristic MscL-sfGFP conductance obtained through patch-clamp
electrophysiology**. Top panel presents a characteristic measurement of
channel current obtained through a patch-clamp electrophysiology measurement
of bacterial protoplasts. The bottom panel shows the applied pressure through
the micropipette to facilitate opening of the mechanosensitive channels. The
blue shaded region indicates opening of the mechanosensitive channel of small
conductance (MscS). The shaded orange region represents opening of single
MscL channels. These regions were used to compute the PL:PS ratio. The [Python
code
(`ch9_figS1.py`)](https://github.com/gchure/phd/blob/master/src/chapter_09/code/ch9_figS1.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch9_figS1){#fig:ephys short-caption="Characteristic MscL-sfGFP conductance
obtained through patch-clamp electrophysiology"}

### Maturation time of MscL-sfGFP

&nbsp;&nbsp;&nbsp;&nbsp; Reliable quantification of the channel copy number
is paramount to this work. As such, it is important to verify that the
detected fluorescence per cell accurately represents the total cellular MscL
copy number. We have made the assumption that the total fluorescence per
represents all MscL-sfGFP channels present. However, it is possible that
there are more channels present per cell but are not detected as the
fluorophores have not properly matured. This potential error becomes more
significant with longer maturation times of the fluorophore as the mean
expression level changes with the growth phase of the culture. With a
maturation time much longer than the typical cell division time, it is
possible that the measured channel copy number represents only a fraction of
the total number inherited over generations.

&nbsp;&nbsp;&nbsp;&nbsp; In our earlier work, we quantified the MscL-sfGFP
channel copy number using fluorescence microscopy as well as with
quantitative Western blotting. We found that these two methods agreed within
20% of the mean value, often with the counts resulting from microscopy being
slightly larger than those measured through Western blotting
[@bialecka-fornal2012]. This strongly suggests that a negligible amount of
channels are not observed due to inactive fluorophores.

&nbsp;&nbsp;&nbsp;&nbsp;Despite these suggestive data, we directly measured
the maturation time of the superfolder GFP protein. We constructed a
chromosomal integration of sfGFP expressed from a promoter under regulation
from plasmid-borne TetR (*E. coli* MG1655 K12
*$\Delta$lacIZYA ybcN::sfGFP*) . These cells were allowed to grow in LB
supplemented with 500 mM NaCl held at 37$^\circ$C to an OD$_{600\text{nm}}$ of approximately
0.3. At this time, transcription and translation of the sfGFP gene was
induced by addition of 10 ng/mL of anhydrous tetracycline. This expression
was allowed to occur for three minutes before the addition of 100 $\mu$g/mL of
kanamycin, ceasing proper protein synthesis. Three minutes of expression was
chosen to provide enough time for transcription and translation. The sfGFP
variant used in this work is 1155 base pairs. We can assume that the rate for
transcription is 42 nucleotides per second (BNID 108488)[@milo2010], meaning
approximately 28 seconds are needed to transcribe the gene. The translation
rate is on the order of 10 amino acids per second, (12 - 42 amino acids / s,
BNID 100059)[@milo2010]. This means that 39 seconds are needed to complete
translation. In total, approximately one minute is needed to complete
expression of the genes. These numbers are not known for LB supplemented with
500 mM NaCl but may be reduced. For this reason, we extended the length of
induction to three minutes before translation was ceased.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The excess anhydrous tetracycline was removed from the culture through
centrifugation and washing with one volume of LB supplemented with 500 mM
NaCl and 100 $\mu$g/mL kanamycin at 37$^\circ$C. The maturation of sfGFP was then
monitored through flow cytometry by measuring the mean expression of 100,000
cells every 60 to 90 seconds. The result of these measurements are shown in
@Fig:maturation_time.

&nbsp;&nbsp;&nbsp;&nbsp;We observe complete maturation of the protein within
20 minutes after translation of the sfGFP gene was ceased. While the growth
rate in LB + 500mM NaCl varies depending on the expression of MscL-sfGFP, we
typically observe doubling times between 30 and 40 minutes, as indicated by a
orange stripe in @Fig:maturation_time (A). To examine the ``best case”
scenario for cell growth in this medium, we measured the growth rate of the
same *E. coli* strain used to measure the fluorophore maturation time
([@Fig:maturation_time] B). We observed a doubling time of 35 $\pm$ 1 min,
which falls in the middle of the yellow stripe shown in
[@Fig:maturation_time] A. These data, coupled with our previous
quantification of MscL copy number using independent methods, suggests that
the fluorescence measurements made in this work reflect the total amount of
MscL protein expressed.

![**Measurement of sfGFP maturation as a function of time through flow
cytometry.** (A) Measurement of sfGFP fluorescence intensity as a function of
time after cessation of protein translation. Points and connected lines
indicate means of gated flow cytometry intensity distributions. Yellow stripe
indicates the range of doubling times observed for the various RBS mutant
strains described in this work (B) Growth curve of *E. coli* MG1655 cells in
LB + 500mM NaCl. Red points indicate individual absorbance measurements. Line
of best fit is shown in black with the uncertainty shown in shaded gray. The
measured doubling time was 35 $\pm$ 1 min. The [Python
code (`ch9_figS2.py`)](https://github.com/gchure/phd/blob/master/src/chapter_09/code/ch9_figS2.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch9_figS2){#fig:maturation_time
short-caption="Measurement of sfGFP maturation as a function of time through
flow cytometry."}


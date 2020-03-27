## Single-Cell Microscopy

In this section, we detail the procedures and results from single-cell
microscopy verification of our flow cytometry measurements. Our previous
measurements of fold-change in gene expression have been measured using
bulk-scale Miller assays  or through single-cell microscopy . In this
work, flow cytometry was an attractive method due to the ability to
screen through many different strains at different concentrations of
inducer in a short amount of time. To verify our results from flow
cytometry, we examined two bacterial strains with different
repressor-DNA binding energies ($\Delta\varepsilon_{RA}$) of
$-13.9~k_BT$ and $-15.3~k_BT$ with $R = 260$ repressors per cell
using fluorescence microscopy and estimated the values of the parameters
$K_A$ and $K_I$ for direct comparison between the two methods. For a
detailed explanation of the Python code implementation of the processing
steps described below, please see this paper’s [GitHub
repository](https://rpgroup-pboc.github.io/mwc_induction/code/notebooks/unsupervised_gating.html).


### Strains and Growth Conditions

Cells were grown in an identical manner to those used for measurement
via flow cytometry (see Materials \& Methods of Chapter 2). Briefly, cells
were grown overnight (between 10 and 13 hours) to saturation in rich media
broth (LB) with $100\,\mu\text{g} \cdot \text{mL}^{-1}$ spectinomycin in a
deep-well 96 well plate at $37^\circ \text{C}$. These cultures were then
diluted 1000-fold into $500\,\mu\text{L}$ of M9 minimal medium supplemented
with 0.5% glucose and the appropriate concentration of the inducer IPTG.
Strains were allowed to grow at $37^\circ \text{C}$ with vigorous aeration
for approximately 8 hours. Prior to mounting for microscopy, the cultures
were diluted 10-fold into M9 glucose minimal medium in the absence of IPTG.
Each construct was measured using the same range of inducer concentration
values as was performed in the flow cytometry measurements (between
$100\,\text{nM}$ and $5\,\text{mM}$ IPTG). Each condition was measured in
triplicate in microscopy whereas approximately ten measurements were made
using flow cytometry.

### Imaging Procedure

During the last hour of cell growth, an agarose mounting substrate was
prepared containing the appropriate concentration of the IPTG inducer.
This mounting substrate was composed of M9 minimal medium supplemented
with 0.5% glucose and 2% agarose (Life Technologies UltraPure Agarose,
Cat. No. 16500100). This solution was heated in a microwave until molten
followed by addition of the IPTG to the appropriate final concentration.
This solution was then thoroughly mixed and a $500\,\mu\text{L}$
aliquot was sandwiched between two glass coverslips and was allowed to
solidify.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Once solid, the agarose substrates were cut
into approximately $10\,\text{mm}\times 10\,\text{mm}$ squares. An aliquot of
one to two microliters of the diluted cell suspension was then added to each
pad. For each concentration of inducer, a sample of the autofluorescence
control, the $\Delta lacI$ constitutive expression control, and the
experimental strain was prepared yielding a total of thirty-six agarose
mounts per experiment. These samples were then mounted onto two glass-bottom
dishes (Ted Pella Wilco Dish, Cat. No. 14027-20) and sealed with parafilm.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;All imaging was performed on a Nikon Ti-Eclipse
inverted fluorescent microscope outfitted with a custom-built laser
illumination system and operated by the open-source MicroManager control
software . The YFP fluorescence was imaged using a CrystaLaser
$514\,\text{nm}$ excitation laser coupled with a laser-optimized (Semrock
Cat. No. LF514-C-000) emission filter.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For each sample, between fifteen and twenty
positions were imaged allowing for measurement of several hundred cells. At each
position, a
phase contrast image, an mCherry image, and a YFP image were collected
in that order with exposures on a time scale of ten to twenty
milliseconds. For each channel, the same exposure time was used across
all samples in a given experiment. All images were collected and stored
in `ome.tiff` format. All microscopy images are available on the
CaltechDATA online repository under DOI: 10.22002/D1.229.

### Image Processing

#### Correcting Uneven Illumination

The excitation laser has a two-dimensional gaussian profile. To minimize
non-uniform illumination of a single field of view, the excitation beam
was expanded to illuminate an area larger than that of the camera
sensor. While this allowed for an entire field of view to be
illuminated, there was still approximately a 10% difference in
illumination across both dimensions. This nonuniformity was corrected
for in post-processing by capturing twenty images of a homogeneously
fluorescent plastic slide (Autofluorescent Plastic Slides, Chroma Cat.
No. 920001) and averaging to generate a map of illumination intensity at
any pixel $I_\text{YFP}$. To correct for shot noise in the camera
(Andor iXon+ 897 EMCCD), twenty images were captured in the absence of
illumination using the exposure time used for the experimental data.
Averaging over these images produced a map of background noise at any
pixel $I_\text{dark}$. To perform the correction, each fluorescent
image in the experimental acquisition was renormalized with respect to
these average maps as 
$$
I_\text{flat} = \frac{I - I_\text{dark}}{I_\text{YFP} - I_\text{dark}}\langle
I_\text{YFP} - I_\text{dark} \rangle,
$${#eq:flatfield_corr}
where $I_\text{flat}$ is the renormalized image and $I$ is the
original fluorescence image. 

### Cell Segmentation

Each bacterial strain constitutively expressed an mCherry fluorophore
from a low copy-number plasmid. This served as a volume marker of cell
mass allowing us to segment individual cells through edge detection in
fluorescence. We used the Marr-Hildreth edge detector  which identifies
edges by taking the second derivative of a lightly Gaussian blurred
image. Edges are identified as those regions which cross from highly
negative to highly positive values or vice-versa within a specified
neighborhood. Bacterial cells were defined as regions within an intact
and closed identified edge. All segmented objects were then labeled and
passed through a series of filtering steps.

To ensure that primarily single cells were segmented, we imposed area
and eccentricity bounds. We assumed that single cells projected into two
dimensions are roughly $2\,\mu\text{m}$ long and $1\,\mu\text{m}$
wide, so that cells are likely to have an area between
$0.5\,\mu\text{m}^2$ and $6\,\mu\text{m}$. To determine the
eccentricity bounds, we assumed that the a single cell can be
approximated by an ellipse with semi-major ($a$) and semi-minor
($b$) axis lengths of $0.5\,\mu\text{m}$ and $0.25\,\mu\text{m}$,
respectively. The eccentricity of this hypothetical cell can be computed
as 
$$
\text{eccentricity} = \sqrt{1 - \left(\frac{b}{a}\right)^2},
$${eq:eccentricity}
yielding a value of approximately 0.8. Any objects with an eccentricity
below this value were not considered to be single cells. After imposing
both an area and eccentricity filter, the remaining objects were
considered cells of interest and the mean fluorescence intensity of
each cell was extracted.

### Calculation of Fold-Change

Cells exhibited background fluorescence even in the absence of an
expressed fluorophore. We corrected for this autofluorescence
contribution to the fold-change calculation by subtracting the mean YFP
fluorescence of cells expressing only the mCherry volume marker from
each experimental measurement. The fold-change in gene expression was
therefore calculated as
$$
\text{fold-change} = \frac{\langle I_{R > 0} \rangle - \langle I_\text{auto}
\rangle}{\langle I_{R = 0} \rangle - \langle I_\text{auto} \rangle},
$${#eq:microscopy_fc}
where $\langle I_{R > 0}\rangle$ is the mean fluorescence intensity of cells
expressing LacI repressors, $\langle I_\text{auto}\rangle$ is the mean
intensity of cells expressing only the mCherry volume marker, and $\langle
I_{R = 0}\rangle$ is the mean fluorescence intensity of cells in the absence
of LacI. These fold-change values were very similar to those obtained through
flow cytometry and were well described using the thermodynamic parameters
used in the main text. With these experimentally measured fold-change values,
the best-fit parameter values of the model were inferred and compared to
those obtained from flow cytometry.

## Parameter Estimation and Comparison

![**Comparison of measured fold-change between flow cytometry and
single-cell microscopy.** Experimentally measured fold-change values
obtained through single-cell microscopy and flow cytometry are shown as
white filled and solid colored circles, respectively. Solid and dashed
lines indicate the predicted behavior using the most likely parameter
values of $K_A$ and $K_I$ inferred from flow cytometry data and
microscopy data, respectively. The red and blue plotting elements
correspond to the different operators O1 and O2 with binding energies
$\Delta\varepsilon_{RA}$ of $-13.9~k_BT$ and $-15.3~k_BT$,
respectively . The marginalized posterior distributions for $K_A$ and
$K_I$ are shown in the top and bottom panel, respectively. The
posterior distribution determined using the microscopy data is wider
than that computed using the flow cytometry data due to a smaller fig
collection of data sets (three for microscopy and ten for flow
cytometry).](ch6_figS10){#fig:flow_microscopy_comparison
short-caption="Comparison of measured fold-change between flow cytometry and
single-cell microscopy."}



To confirm quantitative consistency between flow cytometry and
microscopy, the parameter values of $K_A$ and $K_I$ were also
estimated from three biological replicates of IPTG titration curves
obtained by microscopy for strains with $R=260$ and operators O1 and
O2. shows the data from these measurements (orange circles) and the ten
biological replicates from our flow cytometry measurements (blue
circles), along with the fold-change predictions from each inference. In
comparison with the values obtained by flow cytometry, each parameter
estimate overlapped with the 95\% credible region of our flow cytometry
estimates, as shown in @Fig:flow_microscopy_comparison. Specifically, these
values were $K_A=142^{+40}_{-34}\,\mu\text{M}$ and
$K_I=0.6^{+0.1}_{-0.1}\,\mu\text{M}$ from microscopy and
$K_A = 149^{+14}_{-12}\,\mu\text{M}$ and $K_I =
0.57^{+0.03}_{-0.02}\,\mu\text{M}$ from the flow cytometry data. We
note that the credible regions from the microscopy data shown in are
much broader than those from flow cytometry due to the fewer number of
replicates performed.
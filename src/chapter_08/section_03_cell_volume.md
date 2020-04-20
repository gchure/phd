## Approximating Cell Volume


In Fig. 4.2 of Chapter 4, we make reference to the volume of the cells
grown in various conditions. Here, we illustrate how we approximated
this estimate using measurements of the individual cell segmentation
masks.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Estimation of bacterial cell volume and
its dependence on the total growth rate has been the target of numerous
quantitative studies using a variety of methods including microscopy
[@pilizota2012; @pilizota2014; @taheri-araghi2015a; @schmidt2016;
@schaechter1958] and microfluidics [@kubitschek1986a] revealing fascinating
phenomenology of growth at the single cell level . Despite the high precision
and extensive calibration of these methods, it is not uncommon to have
different methods yield different estimates, indicating that it is not a
trivial measurement to make. In the present work, we sought to estimate the
cell volume and compare it to the well-established empirical results of the
field of bacterial physiology to ensure that our experimental protocol does
not alter the physiology beyond expectations. As the bulk of this work is
performed using single-cell microscopy, we chose to infer the approximate the
cell volume from the segmentation masks produced by the SuperSegger MATLAB
software [@stylianidou2016a] which reported the cell length and width in
units of pixels which can be converted to meaningful units given knowledge of
the camera interpixel distance.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We approximated each segmented cell as a cylinder of length $a$ and
radius $r$ capped on each end by hemispheres with a radius $r$. With
these measurements in hand, the total cell volume was computed as
$$
V_\text{cell} =\pi r^2 \left(a + \frac{4}{3}r\right).
$${#eq:cell_volume}
The output of the SuperSegger segmentation
process is an individual matrix for each cell with a variety of
fluorescence statistics and information regarding the cell shape. Of the
latter category, the software reports in pixels the total length $\ell$
and width $w$ of the cell segmentation mask. Given these measurements,
we computed the radius $r$ of the spherocylinder as 
$$
r = \frac{w}{2}
$${#eq:radius}
and the cylinder length $a$ as 

$$
a = \ell - w.
$${#eq:length}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;@Fig:cell_shape shows the validity of modeling the
segmentation masks as a spherocylinder in two dimensions. Here, the thin
colored lines are the contours of a collection of segmentation masks and
the black solid lines are spherocylinders using the average length and
width of the segmentation masks, as calculated by @Eq:radius and @Eq:length.
It appears that this simple approximation is reasonable for the purposes of this work.

![**Growth-rate dependence of cell shape and spherocylinder approximation.**
Contours of segmentation masks for a single experiment of each condition are
shown as thin colored lines for (A) carbon quality variation and (B) temperature
variation.](ch8_figS2){#fig:cell_shape short-caption="Growth-rate dependence of
cell shape and approximation as a spherocylinder."}

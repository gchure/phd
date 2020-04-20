## Results

### Thermodynamic model
This chapter builds off the theoretical details presented in Chapters 2-3 of
this thesis. Here, we once again consider the simple repression motif in which
expression of a target gene is regulated by the action of a single allosteric
repressor. The key measurable quantity of the following work is the fold-change
in expression of the regulated gene. Thermodynamic models described previously
[@garcia2011; @razo-mejia2018; @phillips2019] and in Chapters 2-3 of this work
result in a succinct input-output function to quantitatively describe the
fold-change in gene expression and is of the form


$$
\text{fold-change} = \left(1 + p_{act}(c)\frac{R}{N_{NS}}e^{-\Delta\varepsilon_{RA}
/ k_BT}\right)^{-1},
$${#eq:foldchange_simple} 


where $R$ is the total number of allosteric repressors per cell, $N_{NS}$ is
the number of nonspecific binding sites for the repressor,
$\Delta\varepsilon_{RA}$ is the repressor-DNA binding energy, and $k_BT$ is the
thermal energy of the system. The prefactor $p_{act}(c)$ defines the
probability of the repressor being in the active state at a given concentration
of inducer $c$. In the absence of inducer, $p_{act}(c = 0)$ can be written as

$$
p_{act}(c = 0) = \left(1 + e^{-\Delta\varepsilon_{AI} / k_BT}\right)^{-1},
$${#eq:growth_pact}

where $\Delta\varepsilon_{AI}$ is the energy difference between the active and
inactive states. Conditioned on only a handful of experimentally accessible
parameters, this model has been verified using the well-characterized LacI
repressor of *Escherichia coli* where parameters such as the repressor copy
number and DNA binding affinity [@garcia2011], copy number of the regulated
promoter [@brewster2014], and the concentration of an extracellular inducer
[@razo-mejia2018, Chapter 2] can be tuned over orders of magnitude. Chapter
3 and the associated publication [@chure2019] illustrated that this model
permits the mapping of mutations within the repressor protein directly to
biophysical parameters in a manner that permits accurate prediction of double
mutant phenotypes . All of these applications, however, have been performed in
a single physiological state where cells are grown in a glucose-supplemented
minimal medium held at 37$^\circ$ C with aeration. In this work, we challenge
this model by changing the environmental conditions away from this
gold-standard condition, perturbing the physiological state of the cell.


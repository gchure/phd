## Discussion

Since the early work by Monod, Wyman, and Changeux [@monod1963; @monod1965],
an array of biological phenomena has been tied to the existence of
macromolecules that switch between inactive and active states. Examples can
be found in a wide variety of cellular processes, including ligand-gated ion
channels [@auerbach2012a], enzymatic reactions [@velyvis2007; @einav2016],
chemotaxis [@keymer2006], quorum sensing [@swem2008], G-protein
coupled receptors [@canals2012], physiologically important proteins [@milo2007;
@levantino2012], and beyond. One of the most ubiquitous examples of allostery is
in the context of gene expression, where an array of molecular players bind to
transcription factors to influence their ability to regulate gene activity
[@huang2011; @li2014a]. A number of studies have focused on developing a
quantitative understanding of allosteric regulatory systems. @martins2011 and
@marzen2013 analytically derived fundamental properties of the MWC model,
including the leakiness and dynamic range described in this work, noting the
inherent trade-offs in these properties when tuning the model’s parameters. Work
in the Church and Voigt labs, among others, has expanded on the availability of
allosteric circuits for synthetic biology [@lutz1997; @moon2012; @rogers2015;
@rohlhill2017]. Somewhat recently, 
@daber2009a theoretically explored the induction of simple repression
within the MWC model and experimentally measured how mutations alter the
induction profiles of transcription factors @daber2011a. @vilar2013 analyzed a
variety of interactions in inducible *lac*-based systems including the
effects of oligomerization and DNA folding on transcription factor induction
. Other work has attempted to use the *lac* system to reconcile *in vitro*
and *in vivo* measurements [@tungtur2011].

Although this body of work has done much to improve our understanding of
allosteric transcription factors, there have been few attempts to
explicitly connect quantitative models to experiments. Here, we generate
a predictive model of allosteric transcriptional regulation and then
test the model against a thorough set of experiments using
well-characterized regulatory components. Specifically, we used the MWC
model to build upon a well-established thermodynamic model of
transcriptional regulation [@bintu2005a; @garcia2011], allowing us to compose the model from a
minimal set of biologically meaningful parameters. This model combines
both theoretical and experimental insights; for example, rather than
considering gene expression directly we analyze the fold-change in
expression, where the weak promoter approximation (see ) circumvents
uncertainty in the RNAP copy number. The resulting model depended upon
experimentally accessible parameters, namely, the repressor copy number,
the repressor-DNA binding energy, and the concentration of inducer. We
tested these predictions on a range of strains whose repressor copy
number spanned two orders of magnitude and whose DNA binding affinity
spanned 6 $k_BT$. We argue that one would not be able to generate such
a wide array of predictions by using a Hill function, which abstracts
away the biophysical meaning of the parameters into phenomenological
parameters [@forsen1995]. Furthermore, our model reveals systematic
relationships between behaviors that previously were only determined empirically.

One such property is the dynamic range, which is of considerable interest when designing or
characterizing a genetic circuit, is revealed to have an interesting
property: although changing the value of $\Delta \varepsilon_{RA}$
causes the dynamic range curves to shift to the right or left, each
curve has the same shape and in particular the same maximum value. This
means that strains with strong or weak binding energies can attain the
same dynamic range when the value of $R$ is tuned to compensate for
the binding energy. This feature is not immediately apparent from the
IPTG induction curves, which show very low dynamic ranges for several of
the O1 and O3 strains. Without the benefit of models that can predict
such phenotypic traits, efforts to engineer genetic circuits with
allosteric transcription factors must rely on trial and error to achieve
specific responses [@rogers2015, @rohlhill2017]. Other calculable properties,
such as leakiness, saturation, [EC$_{50}$] and the effective Hill coefficient,
agree well with experimental measurement. One exception is the titration profile
of the weakest operator, O3. While performing a global fit for all model
parameters marginally improves the prediction of all properties for O3 (see
supplemental Chapter 6), a noticeable difference remains when inferring the
effective Hill coefficient or the [EC$_{50}$]. WE further tried including
additional states (such as allowing the inactive repressor to bind to the
operator), relaxing the weak promoter approximation, accounting for changes in
gene and repressor copy number throughout the cell cycle [@jones2014], and
refitting the original binding energies from @garcia2011, but such
generalizations were unable to account for the O3 data. It remains an open
question as to how the discrepancy between the theory and measurements for O3
can be reconciled.

Despite the diversity observed in the induction profiles of each of our
strains, our data are unified by their reliance on fundamental
biophysical parameters. In particular, we have shown that our model for
fold-change can be rewritten in terms of the free energy, which
encompasses all of the physical parameters of the system. This has
proven to be an illuminating technique in a number of studies of
allosteric proteins . Although it is experimentally straightforward to
observe system responses to changes in effector concentration $c$,
framing the input-output function in terms of $c$ can give the
misleading impression that changes in system parameters lead to
fundamentally altered system responses. Alternatively, if one can find
the “natural variable" that enables the output to collapse onto a single
curve, it becomes clear that the system’s output is not governed by
individual system parameters, but rather the contributions of multiple
parameters that define the natural variable. When our fold-change data
are plotted against the respective free energies for each construct,
they collapse cleanly onto a single curve (see ). This enables us to
analyze how parameters can compensate each other. For example, rather
than viewing strong repression as a consequence of low IPTG
concentration $c$ or high repressor copy number $R$, we can now
observe that strong repression is achieved when the free energy
$F(c) \leq -5 k_BT$, a condition which can be reached in a number of
ways.

While our experiments validated the theoretical predictions in the case
of simple repression, we expect the framework presented here to apply
much more generally to different biological instances of allosteric
regulation. For example, we can use this model to study more complex
systems such as when transcription factors interact with multiple
operators . We can further explore different regulatory configurations
such as corepression, activation, and coactivation, each of which are
found in *E. coli* (see Appendix supplemental Chapter 6). This work can also
serve as a springboard to characterize not just the mean but the full
gene expression distribution and thus quantify the impact of noise on
the system . Another extension of this approach would be to
theoretically predict and experimentally verify whether the
repressor-inducer dissociation constants $K_A$ and $K_I$ or the
energy difference $\Delta \varepsilon_{AI}$ between the allosteric
states can be tuned by making single amino acid substitutions in the
transcription factor . Finally, we expect that the kind of rigorous
quantitative description of the allosteric phenomenon provided here will
make it possible to construct biophysical models of fitness for
allosteric proteins similar to those already invoked to explore the
fitness effects of transcription factor binding site strengths and
protein stability .

To conclude, we find that our application of the MWC model provides an
accurate, predictive framework for understanding simple repression by
allosteric transcription factors. To reach this conclusion, we analyzed
the model in the context of a well-characterized system, in which each
parameter had a clear biophysical meaning. As many of these parameters
had been measured or inferred in previous studies, this gave us a
minimal model with only two free parameters which we inferred from a
single data set. We then accurately predicted the behavior of seventeen
other data sets in which repressor copy number and repressor-DNA binding
energy were systematically varied. In addition, our model allowed us to
understand how key properties such as the leakiness, saturation, dynamic
range, $[EC_{50}]$, and effective Hill coefficient depended upon the
small set of parameters governing this system. Finally, we show that by
framing inducible simple repression in terms of free energy, the data
from all of our experimental strains collapse cleanly onto a single
curve, illustrating the many ways in which a particular output can be
targeted. In total, these results show that a thermodynamic formulation
of the MWC model supersedes phenomenological fitting functions for
understanding transcriptional regulation by allosteric proteins.

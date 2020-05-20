### Inducer Binding Domain Mutants

Much as in the case of the DNA binding mutants, we cannot safely assume
*a priori* that a given mutation in the inducer binding domain affects
only the inducer binding constants $K_A$ and $K_I$. While it is easy
to associate the inducer binding constants with the inducer binding
domain, the critical parameter in our allosteric model
$\Delta\varepsilon_{AI}$ is harder to restrict to a single spatial
region of the protein. As $K_A$, $K_I$, and
$\Delta\varepsilon_{AI}$ are all parameters dictating the allosteric
response, we consider two hypotheses in which inducer binding mutations
alter either all three parameters or only $K_A$ and $K_I$.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We made four point mutations within the inducer binding domain of LacI
(F161T, Q291V, Q291R, and Q291K) that have been shown previously to
alter binding to multiple allosteric effectors [@daber2009a]. In contrast to the DNA
binding domain mutants, we paired the inducer binding domain mutations
with the three native LacI operator sequences (which have various
affinities for the repressor) and a single ribosomal binding site
sequence. This ribosomal binding site sequence, as reported in @garcia2011, 
expresses the wild-type LacI repressor to an average copy number of
approximately $260$ per cell. As the free energy differences resulting
from point mutations in the DNA binding domain can be described solely
by changes to $\Delta\varepsilon_{RA}$, we continue under the
assumption that the inducer binding domain mutations do not
significantly alter the repressor copy number.

![**Induction profiles and free-energy differences of inducer binding domain
mutants.** Open points represent the strain to which the parameters were
fit — namely, the O2 operator sequence. Each column corresponds to the mutant
highlighted at the top of the figure. All strains have $R = 260$ per cell. (A)
The fold change in gene expression as a function of the inducer concentration
for three operator sequences of varying strength. Dashed lines correspond to the
curve of best fit resulting from fitting $K_A$ and $K_I$ alone. Shaded curves
correspond to the 95\% credible region of the induction profile determined
from fitting $K_A$, $K_I$, and $\Delta\varepsilon_{AI}$. Points correspond to the mean measurement of
6–12 biological replicates. Error bars are the SEM. (B) Points in A collapsed
as a function of the free energy calculated from redetermining $K_A$, $K_I$, and
$\Delta\varepsilon_{AI}$. (C) Change in free energy resulting from each mutation as a function of
the inducer concentration. Points correspond to the median of the posterior
distribution for the free energy. Error bars represent the upper and lower
bounds of the 95\% credible region. Shaded curves are the predictions. IPTG
concentration is shown on a symmetric log scaling axis with the linear region
spanning from 0 to $10^{-2}$ $\mu$M and log scaling
elsewhere. The [Python code
(`ch3_fig4.py`)](https://github.com/gchure/phd/blob/master/src/chapter_03/code/ch3_fig4.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch3_fig4){#fig:IND_muts short-caption="Induction profiles and
free-energy differences of inducer binding domain mutants."}


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The induction profiles for these four mutants
are shown in @Fig:IND_muts (A). Of the mutations chosen, Q291R and Q291K
appear to have the most significant impact, with Q291R abolishing the
characteristic sigmoidal titration curve entirely. It is notable that both
Q291R and Q291K have elevated expression in the absence of inducer compared
to the other two mutants paired with the same operator sequence. Panel (A) in
@Fig:deltaF_theory illustrates that if only $K_A$ and $K_I$ were being
affected by the mutations, the fold-change should be identical for all
mutants in the absence of inducer. This discrepancy in the observed leakiness
immediately suggests that more than $K_A$ and $K_I$ are affected for Q291K
and Q291R.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Using a single induction profile for each mutant
(shown in @Fig:IND_muts as white-faced circles), we inferred
the parameter combinations for both hypotheses and drew predictions for
the induction profiles with other operator sequences. We find that the
simplest hypothesis (in which only $K_A$ and $K_I$ are altered) does
not permit accurate prediction of most induction profiles. These curves,
shown as dotted lines in @Fig:IND_muts (A),
fail spectacularly in the case of Q291R and Q291K, and undershoot the
observed profiles for F161T and Q291V, especially when paired with the
weak operator sequence O3. The change in the leakiness for Q291R and
Q291K is particularly evident as the expression at $c = 0$ should be
identical to the wild-type repressor under this hypothesis. Altering
only $K_A$ and $K_I$ is not sufficient to accurately predict the
induction profiles for F161T and Q291V, but not to the same degree as
Q291K and Q291R. The disagreement is most evident for the weakest
operator O3 [green lines in 
@Fig:IND_muts (A)], though we have discussed
previously that the induction profiles for weak operators are difficult
to accurately describe and can result in comparable disagreement for the
wild-type repressor [@razo-mejia2018].

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Including $\Delta\varepsilon_{AI}$ as a perturbed parameter in
addition to $K_A$ and $K_I$ improves the predicted profiles for all
four mutants. By fitting these three parameters to a single strain, we
are able to accurately predict the induction profiles of other operators
as seen by the shaded lines in 
@Fig:IND_muts (A). With these modified parameters,
all experimental measurements collapse as a function of their free
energy as prescribed by @Eq:mut_collapse (@Fig:IND_muts (B)). All four mutations
significantly diminish the binding affinity of both states of the
repressor to the inducer, as seen by the estimated parameter values
reported in Table 3.1. As evident in
the data alone, Q291R abrogates inducibility outright
($K_A \approx K_I$). For Q291K, the active state of the repressor can
no longer bind inducer whereas the inactive state binds with weak
affinity. The remaining two mutants, Q291V and F161T, both show
diminished binding affinity of the inducer to both the active and
inactive states of the repressor relative to the wild-type.

| Mutant | $K_A$ | $K_I$ | $\Delta\varepsilon_{AI}$ [$k_BT$] |  Reference |
| :----- | :----: | :----: | :-----------------------------: | ---------: |
| WT     |  $139^{+29}_{-22}\,\mu$M  | $0.53^{+0.04}_{-0.04}\,\mu$M |   4.5  |  @razo-mejia2018  |
| F161T  |  $165^{+90}_{-65}\,\mu$M  | $3^{+6}_{-3}\,\mu$M          |  $1^{+5}_{-2}$ | This study |
| Q291V  | $650^{+450}_{-250}\,\mu$M | $8^{+8}_{-8}\,\mu$M          |  $3^{+6}_{-3}$ | This study |
| Q291K  | $> 1$ mM                  | $310^{+70}_{-60}\,\mu$M      |  $-3.11^{+0.07}_{-0.07}$ | This study |
| Q291R  | $9_{-9}^{+20}\,\mu$M      | $8^{+20}_{-8}\,\mu$M         |  $-2.35^{+0.01}_{-0.09}$ | This study |
:  Inferred values of $K_A$, $K_I$, and $\Delta\varepsilon_{AI}$ for inducer binding mutants 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Given the collection of fold-change
measurements, we computed the $\Delta F$ relative to the wild-type strain
with the same operator and repressor copy number. This leaves differences in
$p_{act}(c)$ as the sole contributor to the free energy difference, assuming
our hypothesis that $K_A$, $K_I$, and $\Delta\varepsilon_{AI}$ are the only
perturbed parameters is correct. The change in free energy can be seen in
@Fig:IND_muts (C). For all mutants, the free energy difference inferred from
the observed fold-change measurements falls within error of the predictions
generated under the hypothesis that $K_A$, $K_I$, and
$\Delta\varepsilon_{AI}$ are all affected by the mutation (shaded curves in
@Fig:IND_muts (C)). The profile of the free energy change exhibits some of the
rich phenomenology illustrated in @Fig:deltaF_theory (A) and (B). Q291K,
F161T, and Q291V exhibit a non-monotonic dependence on the inducer
concentration, a feature that can only appear when $K_A$ and $K_I$ are
altered. The non-zero $\Delta F$ at $c=0$ for Q291R and Q291K coupled with an
inducer concentration dependence is a telling sign that
$\Delta\varepsilon_{AI}$ must be significantly modified. This shift in
$\Delta F$ is positive in all cases, indicating that $\Delta\varepsilon_{AI}$
must have decreased, and that the inactive state has become more
energetically favorable for these mutants than for the wild-type protein.
Indeed the estimates for $\Delta\varepsilon_{AI}$ (Table 3.1) reveal both
mutations Q291R and Q291K make the inactive state more favorable than the
active state. Thus, for these two mutations, only $\approx 10\%$ of the
repressors are active in the absence of inducer, whereas the basal active
fraction is $\approx 99\%$ for the wild-type repressor [@razo-mejia2018].

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We note that the parameter values reported here disagree with those
reported in @daber2011. This disagreement stems from different assumptions
regarding the residual activity of the repressor in the absence of
inducer and the parametric degeneracy of the MWC model without a
concrete independent measure of $\Delta\varepsilon_{AI}$. A detailed
discussion of the difference in parameter values between our previous
work [@garcia2011; @razo-mejia2018], that of @daber2011, and those of other
seminal works  can be found in the supplemental Chapter 7. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Taken together, these parametric changes
diminish the response of the regulatory architecture as a whole to changing
inducer concentrations. They furthermore reveal that the parameters which
govern the allosteric response are interdependent and no single parameter is
insulated from the others. However, as *only* the allosteric parameters are
changed, one can say that the allosteric parameters as a whole are insulated
from the other components which define the regulatory response, such as
repressor copy number and DNA binding affinity.

## Additional Characterization of Inducer Binding Domain Mutants 

To predict the induction profiles of the inducer binding mutants, we
used only the induction profile of each mutant paired with the native O2
*lac* operator to infer the parameters. Here, we examine the influence
the choice of fit strain has on the predictions of the induction
profiles and $\Delta F$ for each mutant.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In Chapter 3, we dismissed the hypothesis that only $K_A$ and $K_I$
were changing due to the mutation and based the fit to a single
induction profile. In @Fig:KaKi_only_pairwise_comparison, the fits and predictions
for each mutant paired with each operator sequence queried. Here, the
rows correspond to the operator sequence of the fit strain while the
columns correspond to the operator sequence of the predicted strain. The
diagonals show the fit induction profiles and the
corresponding data. Regardless of the choice of fit strain, the
predicted induction profiles of the repressor paired with the O3
operator are poor, with the leakiness in each case being significantly
underestimated. We also see that fitting the allosteric parameters to O3 results in poor
predictions with incredibly wide credible regions for the other two
operators. In @razo-mejia2018 and in Chapter 6, we also found
that fitting $K_A$ and $K_I$ to the induction profile of O3 generally
resulted in poor predictions of the other strains with comparably wide
credible regions.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;When $\Delta\varepsilon_{AI}$ is included as a parameter, however, the
predictive power is improved for all three operators, as can be seen in
@Fig:KaKi_epAI_pairwise_comparison. While the credible
regions are still wide when fit to the O3 operator, they are much
narrower than under the first hypothesis. We emphasize that we are able
to accurately predict the leakiness of nearly every strain by
redetermining $\Delta\varepsilon_{AI}$ whereas the leakiness was not
predicted when only $K_A$ and $K_I$ were considered. Thus, we conclude
that all three allosteric parameters $K_A$, $K_I$, and
$\Delta\varepsilon_{AI}$ are modified for these four inducer binding
domain mutations. The values of the inferred parameters are reported in
Table 7.3.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We also examined the effect the choice of
fit strain has on the predicted $\Delta F$, shown in
@Fig:IND_deltaF_comparison. We find that the predictions agree with the data
regardless of the choice of fit strain. One exception is the prediction of
the Q291K $\Delta F$ when the parameters fit to the O3 induction profile are
used. As the induction profile for Q291K paired with O3 is effectively flat
at a fold-change of 1, it is difficult to properly estimate the parameters of
our sigmoidal function. We note that all measurements of $\Delta F$ for Q291K are
described by using the parameters fit to either O1 or O3 induction
profiles, suggesting that the choice of fit strain makes little difference.

![**Pairwise comparison of fit strain versus predictions assuming only
$K_A$ and $K_I$ are influenced by the mutation.** Rows correspond to the
operator sequence of the strain used for the parameter inference.
Columns correspond to the operator sequence of the predicted strain.
Colors identify the mutation. Diagonal positions show
the induction fit strain and
profiles. The [Python code                                                
(`ch7_figS17-S18.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS17-S18.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS17){#fig:KaKi_only_pairwise_comparison
short-caption="Pairwise comparison of fit strain versus predictions assuming
only $K_A$ and $K_I$ are influenced by mutations in the inducer binding domain."}


![**Pairwise comparison of fit strain versus predictions assuming all
allosteric parameters are affected by the mutation.** Rows correspond to
the operator of the strain used to fit the parameters. Columns
correspond to the operator of the strains whose induction profile is
predicted. Mutants are identified by color. Diagonals (gray background)
show the induction profiles of the strain to which the parameters were
fit. The [Python code                                                
(`ch7_figS17-S18.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS17-S18.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).  ](ch7_figS18){#fig:KaKi_epAI_pairwise_comparison short-caption="Pairwise
comparison of fit strain versus predictions assuming all allosteric parameters
are affected by the mutation in the inducer binding domain."}

![**Comparison of choice of fit strain on predicted $\Delta F$ profiles.**
Rows correspond to the operator of the strain to which the parameters
were fit. Columns correspond to mutations. Points are colored by their
operator sequence. The data corresponding to the operator of the fit
strain are shown as white-faced
points. The [Python code                                                
(`ch7_figS19.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS19.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch7_figS19){#fig:IND_deltaF_comparison short-caption="Comparison of
choice of fit strain on predicted $\Delta F$ profiles for inducer binding domain
mutants."}

|**Mutant** | **Operator** |  $\mathbf{K_A}$ **\[$\mu$M\]** |  $\mathbf{K_I}$ **\[$\mu$M\]** | $\Delta\varepsilon_{AI}$ **\[$\mathbf{k_BT}$\]**|  
|:--:|:--:|:--:|:--:|:--:|
|F164T| O1 | $290_{-56}^{+60}$ | $1^{+4}_{-0.98}$ |$4^{+5}_{-3}$|
| | O2| $165^{+90}_{-65}$ | $3^{+6}_{-3}$ | $1_{-2}^{+5}$|
| | O3 | $110_{-105}^{+700}$ | $7_{-4}^{+5}$ | $-0.9_{-0.3}^{+0.4}$|
| Q291K | O1 | $> 1000$ | $410^{+150}_{-100}$ | $-3.2^{+0.1}_{-0.1}$|
| | O2 | $> 1000$ | $310_{-60}^{+70}$ | $-3.11^{+0.07}_{-0.07}$|
| | O3 | $10_{-10}^{+200}$ | $1^{+9}_{-1}$ | $-7_{-5}^{+3}$ |
| Q291R |  O1 | $3^{+27}_{-3}$ | $2^{+20}_{-2}$ | $-1.9_{-0.3}^{+0.4}$ |
| | O2 | $9_{-9}^{+20}$ | $8_{-8}^{+20}$ | $-2.32_{-0.09}^{+0.01}$|
| | O3 | $6_{-6}^{+24}$ | $9_{-9}^{+30}$ | $-2.6_{-0.5}^{+0.4}$|
| Q291V | O1 | $> 1000$ | $3_{-3}^{+13}$ | $6_{-4}^{+4}$|
| | O2 | $650_{-250}^{+450}$| $8_{-8}^{+8}$ | $3_{-3}^{+6}$|
| | O3 | $100_{-90}^{+400}$ | $22_{-18}^{+33}$ | $0.1_{-0.6}^{+0.8}$|                          
  : Inferred values of $K_A$, $K_I$, and $\Delta\varepsilon_{AI}$ for
  inducer binding domain mutants. Values reported are the mean of the
  posterior distribution with the upper and lower bounds of the 95%
  credible region.

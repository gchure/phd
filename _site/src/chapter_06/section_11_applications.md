## Applications to Other Regulatory Architectures

In this section, we discuss how the theoretical framework presented in this
work is sufficiently general to include a variety of regulatory architectures
outside of simple repression by LacI. We begin by noting that the exact same
formula for fold-change given in can also describe corepression. We then
demonstrate how our model can be generalized to include other architectures,
such as a coactivator binding to an activator to promote gene expression. In
each case, we briefly describe the system and describe its corresponding
theoretical description. For further details, we invite the interested reader
to read @bintu2005 and @marzen2013.

### Corepression

Consider a regulatory architecture where binding of a transcriptional
repressor occludes the binding of RNAP to the DNA. A corepressor
molecule binds to the repressor and shifts its allosteric equilibrium
towards the active state in which it binds more tightly to the DNA,
thereby decreasing gene expression (in contrast, an inducer shifts the
allosteric equilibrium towards the inactive state where the repressor
binds more weakly to the DNA). As in the main text, we can enumerate the
states and statistical weights of the promoter and the allosteric states
of the repressor. We note that these states and weights exactly match those in
Fig. 2.2  and yield the same fold-change equation, 
$$
\text{fold-change} \approx \left(1 + {\left(1 + {c \over K_A}\right)^n \over
\left(1 + {c \over K_A}\right)^n + e^{\beta\Delta\varepsilon_{AI} }\left(1 + {c
\over K_I}\right)^n} {R \over
N_{NS} }e^{-\beta\Delta\varepsilon_{RA} }\right)^{-1},
$${#eq:corepression_fc}
where $c$ now represents the concentration of the corepressor
molecule. Mathematically, the difference between these two architectures
can be seen in the relative sizes of the dissociation constants $K_A$
and $K_I$ between the inducer and repressor in the active and inactive
states, respectively. The corepressor is defined by $K_A < K_I$, since
the corepressor favors binding to the repressor’s active state; an
inducer must satisfy $K_I < K_A$, as was found in Chapter 2. Much as was
performed in Chapter 2, we can make some predictions about the how the response
of a corepressor. In @Fig:applications (A), we show how varying the repressor
copy number $R$ and the repressor-DNA binding energy $\Delta\varepsilon_{RA}$
influence the response. We draw the reader’s attention to the decrease in
fold-change as the concentration of effector is increased.

### Activation

We now turn to the case of activation. While this architecture was not
studied in this work, we wish to demonstrate how the framework presented
here can be extended to include transcription factors other than
repressors. To that end, we consider a transcriptional activator which
binds to DNA and aids in the binding of RNAP through energetic
interaction term $\varepsilon_{AP}$. Note that in this architecture,
binding of the activator does not occlude binding of the polymerase.
Binding of a coactivator molecule shifts its allosteric equilibrium
towards the active state ($K_A < K_I$), where the activator is more
likely to be bound to the DNA and promote expression. Enumerating all of
the states and statistical weights of this architecture and making the
approximation that the promoter is weak generates a fold-change equation
of the form 
$$
\text{fold-change} = {1 + \frac{\left(1 + \frac{c}{K_A}\right)^n}{\left( 1 + \frac{c}{K_A}\right)^n + e^{\beta\Delta\varepsilon_{AI} }\left(1 + \frac{c}{K_I}\right)^n}\frac{A}{
N_{NS} }e^{-\beta\Delta\varepsilon_{AA} }e^{-\beta\varepsilon_{AP} } \over 1 +
{\left(1 + {c \over K_A}\right)^n \over \left(1 + {c \over K_A}\right)^n +
e^{\beta\Delta\varepsilon_{AI} }\left(1 + {c \over K_I}\right)^n}{A \over
N_{NS} }e^{-\beta\Delta\varepsilon_{AA} } },
$${#eq:foldchange_activation}
where $A$ is the total number of activators per cell, $c$ is the
concentration of a coactivator molecule, $\Delta\varepsilon_{AA}$ is
the binding energy of the activator to the DNA in the active allosteric
state, and $\varepsilon_{AP}$ is the interaction energy between the
activator and the RNAP. Unlike in the cases of induction and
corepression, the fold-change formula for activation includes terms from
when the RNAP is bound by itself on the DNA as well as when both RNAP
and the activator are simultaneously bound to the DNA. explores
predictions of the fold-change in gene expression by manipulating the
activator copy number, DNA binding energy, and the polymerase-activator
interaction energy. Note that with this activation scheme, the
fold-change must necessarily be greater than one. An interesting feature
of these predictions is the observation that even small changes in the
interaction energy ($< 0.5\, k_BT$) can result in dramatic increase in
fold-change.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As in the case of induction, the is straightforward to generalize. For
example, the relative values of $K_I$ and $K_A$ can be switched such
that $K_I < K_A$ in which the secondary molecule drives the activator
to assume the inactive state represents induction of an activator. While
these cases might be viewed as separate biological phenomena,
mathematically they can all be described by the same underlying
formalism.

![**Representative fold-change predictions for allosteric
corepression and activation.** (A) Contrary to the case of
induction described in the main text, addition of a corepressor
decreases fold-change in gene expression. The left and right panels
demonstrate how varying the values of the repressor copy number $R$
and repressor-DNA binding energy $\Delta\varepsilon_{RA}$,
respectively, change the predicted response profiles. (B) In the case of
inducible activation, binding of an effector molecule to an activator
transcription factor increases the fold-change in gene expression. Note
that for activation, the fold-change is greater than 1. The left and
center panels show how changing the activator copy number $A$ and
activator-DNA binding energy $\Delta\varepsilon_{AA}$ alter response,
respectively. The right panel shows how varying the polymerase-activator
interaction energy $\varepsilon_{AP}$ alters the fold-change.
Relatively small perturbations to this energetic parameter drastically
changes the level of activation and plays a major role in dictating the
dynamic range of the system.](ch6_figS20){#fig:applications
short-caption="Representative fold-change predictions for allosteric
corepression and activation."}

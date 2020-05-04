## Non-Monotonic Behavior of $\Delta F$ Under Changing $K_A$ and $K_I$

In Chapter 3, we illustrated that perturbations only to the allosteric
parameters $K_A$ and $K_I$ relative to the wild-type values can result in a
non-monotonic dependence of $\Delta F$ on the inducer concentration $c$. In
this section, we prove that when the ratio of $K_A$ to $K_I$ is the same
between the mutant and wild-type proteins, the function must be monotonic.
This section is paired with an interactive figure available on the [paper
website](https://www.rpgroup.caltech.edu/mwc_mutants) which illustrates how
scaling $K_A$ and $K_I$ relative to the wild-type value results in
non-monotonic behavior.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We define a monotonic function as a continuous
function whose derivative does not change sign across the domain upon which
it is defined. To show that $\Delta F$ is non-monotonic when $K_A$ and $K_I$
are perturbed, we can compute the derivative of $\Delta F$ with respect to
the inducer concentration $c$ and evaluate the sign of the derivative at the
limits of inducer concentration. If the sign of the derivative is different
at the limits of $c = 0$ and $c \gg 0$, we can see that the function is
non-monotonic. However, if the sign is the same in both limits, we can not
say conclusively if it is non-monotonic and must consider other diagnostics.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The free energy difference between a mutant and
wild-type repressor when all parameters other than $K_A$ and $K_I$ are
unperturbed can be written as
$$
\beta\Delta F(c) = -\log\left({\left[1 + e^{-\beta\Delta\varepsilon_{AI} }\left({1
+ {c \over K_I^\text{(mut)} } \over 1 + {c \over
K_A^\text{(mut)} }}\right)^2\right]^{-1} \over \left[1 +
e^{-\beta\Delta\varepsilon_{AI} }\left({1 + {c \over K_I^\text{(wt)} } \over 1
+ {c \over K_A^\text{(wt)} }}\right)^2\right]^{-1} }\right),
$${#eq:delF_pact}
in which $\Delta\varepsilon_{AI}$ is the energy difference between the active
and inactive states of the repressor, $c$ is the inducer concentration, and
$\beta = 1 / k_BT$ where $k_B$ is the Boltzmann constant and $T$ is the
temperature. The derivative with respect to $c$, which we determined using
Mathematica's (Wolfram Research, version 11.2) symbolic computing ability, is
given as
$$
\begin{gathered}
{\partial \beta\Delta F(c) \over \partial c} =
2e^{-\beta\Delta\varepsilon_{AI} }\left(\frac{ {K_A^\text{(mut)} }^2 \left(K_A^\text{(mut)}
- K_I^\text{(mut)}\right)\left(c + K_I^\text{(mut)}\right)}{\left(c +
K_A^\text{(mut)}\right) \left[\left(c +
K_A^\text{(mut)}\right)^2{K_I^\text{(mut)} }^2 +
e^{-\beta\Delta\varepsilon_{AI} }{K_A^\text{(mut)} }^2\left(c +
K_I^\text{(mut)}\right)^2\right]} \right. \\ - \left.
\frac{ {K_A^\text{(wt)} }^2\left(K_A^\text{(wt)} - K_I^\text{(wt)}\right)\left(c +
K_I^\text{(wt)}\right)}{\left(c + K_A^\text{(wt)}\right) \left[\left(c +
K_A^\text{(wt)}\right)^2{K_I^\text{(wt)} }^2 +
e^{-\beta\Delta\varepsilon_{AI} }{K_A^\text{(wt)} }^2\left(c +
K_I^\text{(wt)}\right)^2\right]} \right).
\end{gathered}
$${#eq:partial_explicit} 
This unwieldy expression can be simplified by defining the values of
$K_A^\text{(mut)} = \theta K_A^\text{(wt)}$ and
$K_I^\text{(mut)} = \theta K_I^\text{(wt)}$ as relative changes to the
wild-type values where $\theta$ is a scaling parameter. While we can
permit $K_A^\text{(mut)}$ and $K_I^\text{(mut)}$ to vary by different
degrees, we will consider the case in which they are equally perturbed
such that the ratio of $K_A$ to $K_I$ is the same between the mutant and
wild-type versions of the repressor. While the equations become more
cumbersome when one permits the dissociation constants to vary by
different amounts (i.e. $\theta_{K_A}, \theta_{K_I}$), one arrives at
the same conclusion. This definition allows us to rewrite 
@Eq:partial_explicit in the form of 
$$
\begin{gathered}
{\partial \beta\Delta F(c) \over \partial c} = 2{K_A^\text{(wt)} }e^{-\beta\Delta\varepsilon_{AI} }\left({\theta^3 \left({K_A^\text{(wt)} }- {K_I^\text{(wt)} }\right)\left(c + \theta {K_I^\text{(wt)} }\right) \over \left(c + \theta {K_A^\text{(wt)} }\right)\left[\theta^2 {K_I^\text{(wt)} }^2\left(c + \theta {K_A^\text{(wt)} }\right)^2 + e^{-\beta\Delta\varepsilon_{AI} }\theta^2{K_A^\text{(wt)} }^2\left(c + \theta {K_I^\text{(wt)} }\right)^2\right]} \right. \\ - \left. {\left({K_A^\text{(wt)} } - {K_I^\text{(wt)} }\right)\left(c + {K_I^\text{(wt)} }\right) 
\over \left(c + {K_A^\text{(wt)} }\right) \left[\left(c +
{K_A^\text{(wt)} }\right)^2{K_I^\text{(wt)} }^2 +
e^{-\beta\Delta\varepsilon_{AI} }{K_A^\text{(wt)} }^2\left(c +
{K_I^\text{(wt)} }\right)^2\right]} \right).
\end{gathered}
$${#eq:partial_theta} 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With this derivative in hand, we can examine the limits of inducer
concentration. As discussed in the main text, the free energy difference
between the mutant and wild-type repressors when $c = 0$ should be equal to
$0$. However, the derivative at $c =0$ will be different between the
wild-type and the mutant. In this limit, @Eq:partial_theta simplifies to
$$
{\partial \beta\Delta F(c) \over \partial c}\bigg\vert_{c = 0} =
{2e^{-\beta\Delta\varepsilon_{AI} }\left({K_A^\text{(wt)} } - {K_I^\text{(wt)} }\right) \over {K_A^\text{(wt)} }{K_I^\text{(wt)} }\left(1
+ e^{-\beta\Delta\varepsilon_{AI} }\right)}\left({1 \over \theta} - 1\right).
$${#eq:partial_c0} 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;When $\theta < 1$, meaning that the affinity of
the active and inactive states of the repressor to the inducer is increased
relative to wild-type, the derivative is positive. Thus, the repressor bound
state of the promoter becomes less energetically favorable than the repressor
bound state. Similarly, if $\theta > 1$, binding of the inducer to the mutant
repressor is weaker than the wild-type repressor, making $\partial\beta\Delta
F(c) / \partial c < 0$, meaning the repressor bound state becomes more
energetically favorable than the repressor unbound state of the promoter.

With an intuition for the sign of the derivative when $c = 0$, we can
compute the derivative at another extreme where $c \gg 0$. Here,
@Eq:partial_theta reduces to 
$$
{\partial \beta\Delta F(c) \over \partial c}\bigg\vert_{c \gg 0} \approx
{2e^{-\beta\Delta\varepsilon_{AI} }{K_A^\text{(wt)} }^2({K_A^\text{(wt)} } -
{K_I^\text{(wt)} }) \over c^2\left({K_I^\text{(wt)} }^2 +
e^{-\beta\Delta\varepsilon_{AI} }{K_A^\text{(wt)} }^2\right)} \left(\theta -
1\right).
$${#eq:partial_cinf} 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;When $\theta > 1$, @Eq:partial_cinf is positive. This is the opposite sign of
the derivative when $c = 0$ when $\theta > 1$. When $\theta < 1$, @Eq:partial_cinf
becomes negative whereas @Eq:partial_c0 is positive. As the derivative of $\Delta F$
with respect to $c$ changes signs across the defined range of inducer
concentrations, we can say the function is non-monotonic.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;@Fig:nonmono shows the non-monotonic behavior
of $\Delta F$ when $K_A$ and $K_I$ change by the same factor $\theta$
(maintaining the wild-type ratio, @Fig:nonmono (A)) and when $K_A$ and $K_I$
change by different factors (@Fig:nonmono (B)). In both cases, non-monotonic
behavior is observed with the peak difference in the free energy covering
several $k_BT$. We have hosted an interactive figure similar to @Fig:nonmono
on the [paper website](https://rpgroup.caltech.edu/mwc_mutants) where the
reader can modify how $K_A$ and $K_I$ are affected by a mutation and examine
how the active probability, free energy difference, and $\partial \beta
\Delta F / \partial c$ are tuned.

![**Non-monotonic behavior of $\Delta F$ with changes in $K_A$ and $K_I$.**
Middle column shows the allosteric contribution of free energy $F$ plotted as
a function of the inducer concentration. Right column shows the free energy
difference $\Delta F$ as a function of inducer concentration, revealing
non-monotonicity. (A) Behavior of $F$ and $\Delta F$ when the values of $K_A$
and $K_I$ change relative to wild-type, but maintain the same ratio. $\theta$
is the scaling factor for both inducer dissociation constants. (B) Behavior
of $F$ and $\Delta F$ when the values of $K_A$ and $K_I$ change relative to
the wild-type, but by different factors. In both panels, the wild-type
parameter values were taken to be $K_A = 200\,\mu$M, $K_I = 1\,\mu$M and
$\Delta\varepsilon_{AI} = 4.5\,k_BT$. An interactive version of this figure
is available on the [paper website.](http://rpgroup.caltech.edu/mwc_mutants) The [Python code                                                
(`ch7_figS1.py`)](https://github.com/gchure/phd/blob/master/src/chapter_07/code/ch7_figS1.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/gchure/phd).](ch7_figS1){#fig:nonmono
short-caption="Non-monotonic behavior of $\Delta F$ with changes in $K_A$ and
$K_I$."}

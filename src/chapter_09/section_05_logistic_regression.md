
## Logistic Regression

&nbsp;&nbsp;&nbsp;&nbsp;In this work, we were interested in computing the
survival probability under a large hypo-osmotic shock as a function of MscL
channel number. As the channel copy number distributions for each
Shine-Dalgarno sequence mutant were broad and overlapping, we chose to
calculate the survival probability through logistic regression -- a method
that requires no binning of the data providing the least biased estimate of
survival probability. Logistic regression is a technique that has been used
in medical statistics since the late 1950's to describe diverse phenomena
such as dose response curves, criminal recidivism, and survival probabilities
for patients after treatment [@anderson2003; @mishra2016; @stahler2013]. It
has also found much use in machine learning to tune a binary or categorical
response given a continuous input [@cheng2009; @dreiseitl2002].

&nbsp;&nbsp;&nbsp;&nbsp;In this section, we derive a statistical model for
estimating the most-likely values for the coefficients $\beta_0$ and
$\beta_1$ and use Bayes' theorem to provide an interpretation for the
statistical meaning.

### Bayesian parameter estimation of $\beta_0$ and $\beta_1$ 

&nbsp;&nbsp;&nbsp;&nbsp;The central challenge of this work is to estimate the
probability of survival $p_s$ given only a measure of the total number of
MscL channels in that cell. In other words, for a given measurement of $N_c$
channels, we want to know likelihood that a cell would survive an osmotic
shock. Using Bayes' theorem, we can write a statistical model for the
survival probability as
$$
g(p_s\,\vert\, N_c) = {f(N_c\,\vert\, p_s)g(p_s) \over f(N_c)},
$${#eq:bayes_surv_prob}

where $g$ and $f$ represent probability density functions over parameters and
data, respectively. The posterior probability distribution $g(p_s\,\vert\,
N_c)$ describes the probability of $p_s$ given a specific number of channels
$N_c$. This distribution is dependent on the likelihood of observing $N_c$
channels assuming a value of $p_s$ multiplied by all prior knowledge we have
about knowing nothing about the data, $g(s)$. The denominator $f(N_c)$ in
@Eq:bayes_surv_prob captures all knowledge we have about the available
values of $N_c$, knowing nothing about the true survival probability. As this
term acts as a normalization constant, we will neglect it in the following
calculations for convenience.
&nbsp;&nbsp;&nbsp;&nbsp;To begin, we must come up with a statistical model
that describes the experimental measurable in our experiment -- survival or
death. As this is a binary response, we can consider each measurement as a
Bernoulli trial with a probability of success matching our probability of
survival $p_s$, $$ f(s\, \vert\, p_s) = p_s^s(1 - p_s)^{1-s},
$${#eq:bernoulli} where $s$ is the binary response of $1$ or $0$ for survival
and death, respectively. As is stated in the introduction to this section, we
decided to use a logistic function to describe the survival probability. We
assume that the log-odds of survival is linear with respect to the effective
channel copy number $N_c$ as $$ \log{p_s \over 1 - p_s} = \beta_0 + \beta_1
N_c, $${#eq:logit} where $\beta_0$ and $\beta_1$ are coefficients which
describe the survival probability in the absence of channels and the increase
in log-odds of survival conveyed by a single channel. The rationale behind
this interpretation is presented in the following section, *A Bayesian
interpretation of $\beta_0$ and $\beta_1$*. Using this assumption, we can
solve for the survival probability $p_s$ as, $$ p_s = {1 \over 1 +
e^{-\beta_0 -\beta_1 N_c}}. $${#eq:logistic_prob}

With a functional form for the survival probability, the likelihood stated in
@Eq:bayes_surv_prob can be restated as
$$
f(N_c, s\,\vert\,\beta_0,\beta_1) = \left({1 \over 1 + e^{-\beta_0 - \beta_1 N_c}}\right)^s\left(1 - {1 \over 1 + e^{-\beta_0 - \beta_1 N_c}}\right)^{1 - s}.
$${#eq:bernoulli_likelihood}
As we have now introduced two parameters, $\beta_0$, and $\beta_1$, we must
provide some description of our prior knowledge regarding their values. As is
typically the case, we know nothing about the values for $\beta_0$ and
$\beta_1$. These parameters are allowed to take any value, so long as it is a
real number. Since all values are allowable, we can assume a flat
distribution where any value has an equally likely probability. This value of
this constant probability is not necessary for our calculation and is
ignored. For a set of $k$ single-cell measurements, we can write the
posterior probability distribution stated in @Eq:bayes_surv_prob as
$$
g(\beta_0, \beta_1\,\vert\, N_c, s) = \prod\limits_{i=1}^n\left({1 \over 1 + e^{-\beta_0 - \beta_1 N_c^{(i)}}}\right)^{s^{(i)}}\left(1 - {1 \over 1 + e^{-\beta_0 - \beta_1 N_c^{(i)}}}\right)^{1 - s^{(i)}}
$${#eq:known_nc_post}

&nbsp;&nbsp;&nbsp;&nbsp;Implicitly stated in @Eq:known_nc_post is absolute
knowledge of the channel copy number $N_c$. However, as is described in
*Standard Candle Calibration*, we must convert from a measured areal sfGFP
intensity $I_A$ to a effective channel copy number,
$$
N_c = {I_A \tilde{\langle A \rangle} \over \tilde{\alpha}},
$${#eq:standard_candle}
where $\tilde{\langle A \rangle}$ is the average cell area of the standard
candle strain and $\tilde{\alpha}$ is the most-likely value for the
calibration factor between arbitrary units and protein copy number. In
*Standard Candle Calibration*, we detailed a process for generating an
estimate for the most-likely value of $\tilde{\langle A \rangle}$ and
$\tilde{\alpha}$. Given these estimates, we can include an informative prior
for each value. From the Markov chain Monte Carlo samples shown in
[@Fig:posterior_samples], the posterior distribution for each parameter is
approximately Gaussian. By approximating them as Gaussian distributions, we
can assign an informative prior for each as
$$
g(\alpha\,\vert\,\tilde{\alpha}, \tilde{\sigma}_\alpha) \propto {1 \over \tilde{\sigma}_\alpha^k}\prod\limits_{i=1}^k\exp\left[-{(\alpha_i - \tilde{\alpha})^2 \over  2\tilde{\sigma}_\alpha^2}\right]
$${#eq:alpha_prior}
for the calibration factor for each cell and 
$$
g(\langle A \rangle\,\vert\,\tilde{\langle A \rangle},\tilde{\sigma}_{\langle A \rangle}) = {1 \over \tilde{\sigma}_{\langle A \rangle}^k}\prod\limits_{i=1}^k\exp\left[-{(\langle A \rangle_i - \tilde{\langle A \rangle})^2 \over 2\tilde{\sigma}_{\langle A \rangle}^2}\right],
$${#eq:area_prior}
where $\tilde{\sigma}_\alpha$ and $\tilde{\sigma}_{\langle A \rangle}$
represent the variance from approximating each posterior as a Gaussian. The
proportionality for each prior arises from the neglecting of normalization
constants for notational convenience.

&nbsp;&nbsp;&nbsp;&nbsp;Given [@Eq:bernoulli_likelihood] through
[@Eq:area_prior], the complete posterior distribution for estimating the most
likely values of $\beta_0$ and $\beta_1$ can be written as
$$
\begin{aligned}
g(\beta_0, &\beta_1\,\vert\,[I_A, s],\tilde{\langle A \rangle}, \tilde{\sigma}_{\langle A \rangle}, \tilde{\alpha}, \tilde{\sigma}_\alpha) \propto{1 \over (\tilde{\sigma}_\alpha\tilde{\sigma}_{\langle A \rangle})^k}\prod\limits_{i=1}^k\left(1 + \exp\left[-\beta_0 - \beta_1 {{I_A}_i \langle A \rangle_i \over \alpha_i}\right]\right)^{-s_i}\,\times\,\\
&\left(1 - \left(1 + \exp\left[-\beta_0 - \beta_1 {{I_A}_i\langle A \rangle_i \over \alpha_i}\right]\right)^{-1}\right)^{1 - s_i}
\exp\left[-{(\langle A \rangle_i - \tilde{\langle A \rangle})^2 \over 2\tilde{\sigma}_{\langle A \rangle}} - {(\alpha_i - \tilde{\alpha})^2\over 2\tilde{\sigma}_\alpha^2}\right]
\end{aligned}.
$${#eq:logistic_posterior}

As this posterior distribution is not solvable analytically, we used Markov
chain Monte Carlo to draw samples out of this distribution, using the log of
the effective channel number as described in the main text. The posterior
distributions for $\beta_0$ and $\beta_1$ for both slow and fast shock rate
data can be seen in [@Fig:logistic_posterior_distributions]

![**Posterior distributions for logistic regression coefficients evaluated for
fast and slow shock rates.** (A) Kernel density estimates of the posterior
distribution for $\beta_0$ for fast (blue) and slow (purple) shock rates. (B)
Kernel density estimates of posterior distribution for
$\beta_1$.](../figs/figS4.pdf ){#fig:logistic_posterior_distributions
short-caption="Posterior distributions for logistic regression coefficients
evaluated for fast and slow shock rates."}

### A Bayesian interpretation of $\beta_0$ and $\beta_1$ 

The assumption of a linear relationship between the log-odds of survival and
the predictor variable $N_c$ appears to be arbitrary and is presented without
justification. However, this relationship is directly connected to the manner
in which Bayes' theorem updates the posterior probability distribution upon
the observation of new data. In following section, we will demonstrate this
connection using the relationship between survival and channel copy number.
However, this description is general and can be applied to any logistic
regression model so long as the response variable is binary. This connection
was shown briefly by Allen Downey in 2014 and has been expanded upon in this
work [@downey2014].

&nbsp;&nbsp;&nbsp;&nbsp;The probability of observing a survival event $s$
given a measurement of $N_c$ channels can be stated using Bayes' theorem as
$$
g(s\,\vert\, N_c) = {f(N_c\,\vert\, s)g(s) \over f(N_c)}.
$${#eq:survival_bayes}
where $g$ and $f$ represent probability density functions over parameters and data respectively. The posterior distribution $g(s\,\vert\, N_c)$ is the quantity of interest and implicitly related to the probability of survival. The likelihood $g(N_c\,\vert\, s)$ tells us the probability of observing $N_c$ channels in this cell given that it survives. The quantity $g(s)$ captures all *a priori* knowledge we have regarding the probability of this cell surviving and the denominator $f(N_c)$ tells us the converse -- the probability of observing $N_c$ cells irrespective of the survival outcome.

&nbsp;&nbsp;&nbsp;&nbsp;Proper calculation of [@Eq:survival_bayes] requires
that we have knowledge of $f(N_c)$, which is difficult to estimate. While we
are able to give appropriate bounds on this term, such as a requirement of
positivity and some knowledge of the maximum membrane packing density, it is
not so obvious to determine the distribution between these bounds. Given this
difficulty, it's easier to compute the odds of survival
$\mathcal{O}(s\,\vert\, N_c)$, the probability of survival $s$ relative to
death $d$,
$$
\mathcal{O}(s\,\vert\, N_c) = {g(s\,\vert\,N_c) \over g(d\,\vert\, N_c)} = {f(N_c\,\vert\, s)g(s) \over f(N_c\,\vert\,d)g(d)},
$${#eq:odds_definition}
where $f(N_c)$ is cancelled. The only stipulation on the possible value of
the odds is that it must be a positive value. As we would like to equally
weigh odds less than one as those of several hundred or thousand, it is more
convenient to compute the log-odds, given as
$$
\log \mathcal{O}(s\,\vert\,N_c)= \log {g(s) \over g(d)} + \log {f(N_c \,\vert\, s )\over f(N_c\,\vert\, d)}.
$${#eq:log_odds}
Computing the log-transform reveals two interesting quantities. The first
term is the ratio of the priors and tells us the *a priori* knowledge of the
odds of survival irrespective of the number of channels. As we have no reason
to think that survival is more likely than death, this ratio goes to unity.
The second term is the log likelihood ratio and tells us how likely we are to
observe a given channel copy number $N_c$ given the cell survives relative to
when it dies.

&nbsp;&nbsp;&nbsp;&nbsp;For each channel copy number, we can evaluate
@Eq:log_odds to measure the log-odds of survival. If we start with zero
channels per cell, we can write the log-odds of survival as
$$
\log \mathcal{O}(s\,\vert\,N_c=0) = \log {g(s) \over g(d)} + \log {f(N_c=0\,\vert\, s) \over f(N_c=0\,\vert\, d)}.
$${#eq:nc_0}
For a channel copy number of one, the odds of survival is
$$
\log \mathcal{O}(s\,\vert\,N_c=1) = \log{g(s) \over g(d)} + \log{f(N_c=1\,\vert\, s) \over f(N_c=1\,\vert\, d)}.
$${#eq:nc_1}
In both @Eq:nc_0 and @Eq:nc_1, the log of our *a priori* knowledge of
survival versus death remains. The only factor that is changing is log
likelihood ratio. We can be more general in our language and say that the
log-odds of survival is increased by the difference in the log-odds conveyed
by addition of a single channel. We can rewrite the log likelihood ratio in a
more general form as
$$
\log {f(N_c\, \vert\, s) \over f(N_c\,\vert\, d)} = \log{f(N_c = 0\,\vert\,s) \over f(N_c=0\,\vert\, d)} + N_c \left[\log{f(N_c=1 \,\vert\,s) \over f(N_c=1\,\vert\, d)} - \log{f(N_c=0\,\vert\, s) \over f(N_c=0\,\vert\, d)}\right],
$${#eq:generalized_LLR}
where we are now only considering the case in which $N_c \in [0, 1]$.
The bracketed term in @Eq:generalized_LLR is the log of the odds of survival given a single channel relative to the odds of survival given no channels. Mathematically, this odds-ratio can be expressed as
$$
\log\mathcal{OR}_{N_c}(s) = \log{{f(N_c=1\,\vert\,s)g(s)\over f(N_c=1\,\vert\,d)g(d)}\over {f(N_c=0\,\vert\,s)g(s)\over f(N_c=0\,\vert\,d)g(d)}} = \log{f(N_c=1\,\vert\,s) \over f(N_c=1\,\vert\,d)} - \log{f(N_c=0\,\vert\,s)\over f(N_c=0\,\vert\,d)} .
$${#eq:odds_ratio}
@Eq:odds_ratio is mathematically equivalent to the bracketed term shown in @Eq:generalized_LLR. 

&nbsp;&nbsp;&nbsp;&nbsp;We can now begin to staple these pieces together to
arrive at an expression for the log odds of survival. Combining
@Eq:generalized_LLR with @Eq:log_odds yields
$$
\log \mathcal{O}(s\,\vert\,N_c) = \log{g(s) \over g(d)} + \log {f(N_c=0\,\vert\, s) \over f(N_c=0\,\vert\, d)} + N_c\left[{f(N_c=1\,\vert\,s)\over f(N_c=1\,\vert\,d)} - \log{f(N_c=0\,\vert\,s) \over f(N_c=0\,\vert\,d)}\right].
$${#eq:combined_v1}
Using our knowledge that the bracketed term is the log odds-ratio and the
first two times represents the log-odds of survival with $N_c=0$, we conclude
with
$$
\log\mathcal{O}(s\,\vert\,N_c) = \log \mathcal{O}(s\,\vert\,N_c=0) + N_c \log \mathcal{OR}_{N_c}(s).
$${#eq:bayesian_logit}
This result can be directly compared to Eq. 1 presented in the main text,
$$
\log {p_s \over 1 - p_s} = \beta_0 + \beta_1 N_c,
$${#eq:}
which allows for an interpretation of the seemingly arbitrary coefficients
$\beta_0$ and $\beta_1$. The intercept term, $\beta_0$, captures the log-odds
of survival with no MscL channels. The slope, $\beta_1$, describes the log
odds-ratio of survival which a single channel relative to the odds of
survival with no channels at all. While we have examined this considering
only two possible channel copy numbers ($1$ and $0$), the relationship
between them is linear. We can therefore generalize this for any MscL copy
number as the increase in the log-odds of survival is constant for addition
of a single channel.

### Other properties as predictor variables ###

&nbsp;&nbsp;&nbsp;&nbsp; The previous two sections discuss in detail the
logic and practice behind the application of logistic regression to cell
survival data using only the effective channel copy number as the predictor
of survival. However, there are a variety of properties that could rightly be
used as predictor variables, such as cell area and shock rate. As is
stipulated in our standard candle calibration, there should be no correlation
between survival and cell area. [@Fig:alternative_predictor_variables]A and B
show the logistic regression performed on the cell area. We see for both slow
and fast shock groups,there is little change in survival probability with
changing cell area and the wide credible regions allow for both positive and
negative correlation between survival and area. The appearance of a bottle
neck in the notably wide credible regions is a result of a large fraction of
the measurements being tightly distributed about a mean value.
[@Fig:alternative_predictor_variables]C shows the predicted survival
probability as a function of the the shock rate. There is a slight decrease
in survivability as a function of increasing shock rate, however the width of
the credible region allows for slightly positive or slightly negative
correlation. While we have presented logistic regression in this section as a
one-dimensional method, [@Eq:logit] can be generalized to $n$ predictor
variables $x$ as
$$
\log {p_s \over 1 - p_s} = \beta_0 + \sum\limits_{i}^n \beta_ix_i.
$${#eq:generalized_logit}
Using this generalization, we can use both shock rate and the effective
channel copy number as predictor variables. The resulting two-dimensional
surface of survival probability is shown in
@Fig:alternative_predictor_variables (D). As is suggested by
@Fig:alternative_predictor_variables (C), the magnitude of change in
survivability as the shock rate is increased is smaller than that along
increasing channel copy number, supporting our conclusion that for MscL
alone, the copy number is the most important variable in determining
survival.


![**Survival probability estimation using alternative predictor variables.** (A)
Estimated survival probability as a function of cell area for the slow shock
group. (B) Estimated survival probability as a function of cell area for the
fast shock group. (C) Estimated survival probability as a function shock
rate. Black points at top and bottom of plots represent single-cell
measurements of cells who survived and perished, respectively. Shaded regions
in (A) - (C) represent the 95\% credible region. (D) Surface of estimated
survival probability using both shock rate and effective channel number as
predictor variables. Black points at left and right of plot represent
single-cell measurements of cells which survived and died, respectively,
sorted by shock rate. Points at top and bottom of plot represent survival and
death sorted by their effective channel copy number. Labeled contours
correspond to the survival
probability.](../figs/figSX_alternative_predictors.png){#fig:alternative_predictor_variables
short-caption="Survival probability estimation using alternative predictor variables."}
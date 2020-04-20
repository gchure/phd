## Induction of Simple Repression with Multiple Promoters or Competitor Sites
We made the choice to perform all of our experiments using strains in which
a single copy of our simple repression construct had been integrated into the
chromosome. This stands in contrast to the methods used by a number of other
studies [@oehler1994; @setty2003; @oehler2006; @daber2009a; @daber2011a;
@vilar2013; @shis2014; @sochor2014], in which reporter constructs are placed on
plasmid, meaning that the number of constructs in the cell is not precisely
known. It is also common to express repressor on plasmid to boost its copy
number, which results in an uncertain value for repressor copy number. Here we
show that our treatment of the MWC model has broad predictive power beyond the
single-promoter scenario we explore experimentally, and indeed can account for
systems in which multiple promoters compete for the repressor of interest.
Additionally, we demonstrate the importance of having precise control over
these parameters, as they can have a significant effect on the induction
profile.

### Chemical Potential Formulation to Calculate Fold-Change

In this section, we discuss a simple repression construct which we
generalize in two ways from the scenario discussed in the text. First,
we will allow the repressor to bind to $N_S$ identical specific
promoters whose fold-change we are interested in measuring, with each
promoter containing a single repressor binding site ($N_S = 1$ in Chapter 2). Second, we consider $N_C$ identical competitor sites which
do not regulate the promoter of interest, but whose binding energies are
substantially stronger than non-specific binding ($N_C = 0$ in Chapter 2). As in Chapter 2, we assume that the rest of the genome
contains $N_{NS}$ non-specific binding sites for the repressor. Using the formalism described in the previous section, we can write the
fold-change in the grand canonical ensemble as 
$$
\text{fold-change} = \frac{1}{1 + \lambda_r e^{-\beta \Delta \varepsilon_{RA}}},
$${#eq:fc_fugacity}
where $\lambda_r$ is the fugacity of the repressor and $\Delta
\varepsilon_{RA}$ represents the energy difference between the
repressor’s binding affinity to the specific operator of interest
relative to the repressor’s non-specific binding affinity to the rest of
the genome.

We now expand our definition of the total number of repressors in the
system, $R_{\text{tot}}$, so that it is given by
$$
R_\text{tot} = R_S + R_{NS} + R_C,
$${#eq:fugacity_rtot}
where $R_S$, $R_{NS}$, and $R_C$ represent the number of repressors bound to the specific promoter, a non-specific binding site, or to a competitor binding site, respectively. The value of $R_S$ is given by 
$$
R_S = N_S \frac{\lambda_r e^{-\beta \Delta \varepsilon_{RA}}}{1 + \lambda_r e^{-\beta \Delta \varepsilon_{RA}}},
$${#eq:fugacity_rs}
where $N_S$ is the number of specific binding sites in the cell. The
value of $R_{NS}$ is similarly given by 
$$
R_{NS} = N_{NS} \frac{\lambda_r}{1 + \lambda_r},
$${#eq:fugacity_rns}
where $N_{NS}$ is the number of non-specific sites in the cell (recall that we use
$N_{NS} = 4.6 \times 10^6$ for *E. coli*), and $R_C$ is given by
$$
R_C = N_C \frac{\lambda_r e^{-\beta \Delta \varepsilon_C}}{1 + \lambda_r e^{-\beta \Delta \varepsilon_C}},
$${#eq:fugacity_rnc}
where $N_C$ is the number of competitor sites in the cell and $\Delta
\varepsilon_C$ is the binding energy of the repressor to the competitor
site relative to its non-specific binding energy to the rest of the
genome.

To account for the induction of the repressor, we replace the total
number of repressors $R_{\text{tot}}$ in @Eq:fugacity_rtot by the number of active repressors in the cell, $p_\text{act}(c) R_{\text{tot}}$. Here, $p_\text{act}$ denotes the probability that the repressor is in the active state (Eq. 2.4),
$$
p_\text{act}(c)=\frac{\left(1+\frac{c}{K_A}\right)^n}{\left(1+\frac{c}{K_A}\right)^n+e^{-\beta \Delta \varepsilon_{AI} }\left(1+\frac{c}{K_I}\right)^n}.
$${#eq:fugacity_pact}
Substituting @Eq:fugacity_rs - @Eq:fugacity_rc into the modified @Eq:fugacity_rtot yields
$$
p_\text{active}(c) R_{\text{tot}} = N_S \frac{\lambda_r e^{-\beta \Delta \varepsilon_{RA}}}{1 + \lambda_r e^{-\beta \Delta \varepsilon_{RA}}} + N_{NS} \frac{\lambda_r}{1 + \lambda_r} + N_C \frac{\lambda_r e^{-\beta \Delta \varepsilon_C}}{1 + \lambda_r e^{-\beta \Delta \varepsilon_C}}.
$${#eq:fugacity_reff}
For systems where the number of binding sites $N_S$, $N_{NS}$, and
$N_C$ are known, together with the binding affinities
$\Delta \varepsilon_{RA}$ and $\Delta \varepsilon_C$, we can solve
numerically for $\lambda_r$ and then substitute it into to obtain a
fold-change at any concentration of inducer $c$. In the following
sections, we will theoretically explore the induction curves given by
@Eq:fugacity_reff for a number of different combinations of simple repression binding sites, thereby predicting how the system would behave if additional
specific or competitor binding sites were introduced.

## Variable Repressor Copy Number ($\boldsymbol{R}$) with Multiple Specific Binding Sites ($\boldsymbol{N_S > 1}$)

In Chapter 2, we consider the induction profiles of strains with
varying $R$ but a single, specific binding site $N_S = 1$.
Here we predict the induction profiles for similar strains in which
$R$ is varied, but $N_S > 1$, as shown in @Fig:fugacity_r. The top row shows
induction profiles in which $N_S = 10$ and the bottom row shows
profiles in which $N_S = 100$, assuming three different choices for
the specific operator binding sites given by the O1, O2, and O3
operators. These values of $N_S$ were chosen to mimic the common
scenario in which a promoter construct is placed on either a low or high
copy number plasmid. A few features stand out in these profiles. First,
as the magnitude of $N_S$ surpasses the number of repressors $R$,
the leakiness begins to increase significantly, since there are no
longer enough repressors to regulate all copies of the promoter of
interest. Second, in the cases where $\Delta \varepsilon_{RA} = -15.3\,k_B T$ for the O1 operator or $\Delta \varepsilon_{RA} = -13.9\,k_B T$ for the O2 operator, the profiles where $N_S = 100$ are notably sharper than the profiles where $N_S = 10$, and it is possible to achieve dynamic ranges approaching 1. Finally, it is interesting to note that the profiles for the O3 operator where
$\Delta \varepsilon_{RA} = -9.7~k_B T$ are nearly indifferent to the
value of $N_S$.

![**Induction with variable $\boldsymbol{R}$ and multiple specific
binding sites.** Induction profiles are shown for strains with variable
R and $\Delta \varepsilon_{RA} = -15.3$, $-13.9$, or $-9.7~k_B T$. The number
of specific sites, $N_S$, is held constant at 10 as $R$ and $\Delta \varepsilon_{RA}$ are varied. $N_S$ is held constant at 100 as $R$ and $\Delta \varepsilon_{RA}$ are varied. These situations mimic the common scenario in which a promoter construct is placed on either a low or high copy number plasmid.](ch6_figS3){#fig:fugacity_r short-caption="Induction with variable $\boldsymbol{R}$ and multiple specific binding sites."}


## Variable Number of Specific Binding Sites $\boldsymbol{N_S}$ with Fixed Repressor Copy Number ($\boldsymbol{R}$)

The second set of scenarios we consider is the case in which the repressor copy
number $R=260$ is held constant while the number of specific promoters $N_S$ is
varied (see @Fig:fugacity_Ns). Again, we see that leakiness is increased
significantly when $N_S > R$, though all profiles for $\Delta \varepsilon_{RA}
= -9.7\, k_BT$ exhibit high leakiness, making the effect less dramatic for this
operator. Additionally, we find again that adjusting the number of specific
sites can produce induction profiles with maximal dynamic ranges. In
particular, the O1 and O2 profiles with $\Delta \varepsilon_{RA} = -15.3$ and
$-13.9\, k_BT$, respectively, have dynamic ranges approaching 1 for $N_S = 50$
and 100.

![**Induction with variable specific sites and fixed
$\boldsymbol{R}$.** Induction profiles are shown for strains with
$R=260$ and $\Delta \varepsilon_{RA} = -15.3~k_BT$, $\Delta \varepsilon_{RA} = -13.9\,k_BT$, or $\Delta \varepsilon_{RA} = -9.7~k_BT$. The number
of specific sites $N_S$ is varied from 1 to 500.](ch6_figS4){#fig:fugacity_Ns short-caption="Induction with variable specific sites and fixed
$\boldsymbol{R}$."}


## Competitor Binding Sites

An intriguing scenario is presented by the possibility of competitor
sites elsewhere in the genome. This serves as a model for situations in
which a promoter of interest is regulated by a transcription factor that
has multiple targets. This is highly relevant, as the majority of
transcription factors in *E. coli* have at least two known binding
sites, with approximately 50 transcription factors having more than ten
known binding sites [@rydenfelt2014; @schmidt2016]. If the number of competitor sites and their average binding energy is known, however, they can be accounted for in the model. Here, we predict the induction profiles for strains in which
$R=260$ and $N_S=1$, but there is a variable number of competitor
sites $N_C$ with a strong binding energy $\Delta \varepsilon_C = -17.0~k_BT$. In the presence of such a strong competitor, when $N_C > R$ the leakiness is greatly increased, as many repressors are siphoned into the pool of competitor sites. This is most dramatic for the case where $\Delta \varepsilon_{RA} = -9.7~k_B T$, in which it appears that no repression occurs at all when $N_C = 500$.
Interestingly, when $N_C < R$ the effects of the competitor are not
especially notable.

![**Induction with variable competitor sites, a single specific site,
and fixed $\boldsymbol{R}$.** Induction profiles are shown for strains
with $R=260$, $N_s=1$, and $\Delta \varepsilon_{RA} = -15.3~k_B T$
for the O1 operator, $\Delta \varepsilon_{RA} = -13.9~k_B T$ for the
O2 operator, or $\Delta \varepsilon_{RA} = -9.7~k_B T$ for the O3
operator. The number of specific sites, $N_C$, is varied from 1 to
500. This mimics the common scenario in which a transcription factor has
multiple binding sites in the
genome.<span label="fig:Nc"></span>](ch6_figS5){#fig:Nc short-caption="Induction with variable competitor sites, a single specific site, and fixed $\boldsymbol{R}$."}

## Properties of the Induction Response

As discussed in the main body of the paper, our treatment of the MWC
model allows us to predict key properties of induction responses. Here,
we consider the leakiness, saturation, and dynamic range (diagrammed in Fig. 2.1) by numerically solving @Eq:fugacity_reff in the absence of inducer, $c=0$, and in the presence of saturating inducer $c \to \infty$. Using @Eq:fugacity_pact, the former case is given by
$$
R_\text{tot} \frac{1}{1 + e^{-\beta \Delta \varepsilon_{AI}}} = N_S \frac{\lambda_r e^{-\beta \Delta \varepsilon_{RA}}}{1 + \lambda_r e^{-\beta \Delta \varepsilon_{RA}}} + N_{NS} \frac{\lambda_r}{1 + \lambda_r} + N_C \frac{\lambda_r e^{-\beta \Delta \varepsilon_C}}{1 + \lambda_r e^{-\beta \Delta \varepsilon_C}},
$${#eq:fugacity_leakiness}
whereupon substituting in the value of $\lambda_r$ into @Eq:fugacity_foldchange will yield the leakiness. Similarly, the limit of saturating inducer is found by
determining $\lambda_r$ from the form
$$
R_{\text{tot}} \frac{1}{1 + e^{-\beta \Delta \varepsilon_{AI}} \left(\frac{K_A}{K_I} \right)^2} = N_S \frac{\lambda_r e^{-\beta \Delta \varepsilon_{RA}}}{1 + \lambda_r e^{-\beta \Delta \varepsilon_{RA}}} + N_{NS} \frac{\lambda_r}{1 + \lambda_r} + N_C \frac{\lambda_r e^{-\beta \Delta \varepsilon_C}}{1 + \lambda_r e^{-\beta \Delta \varepsilon_C}}.
$${#eq:fugacity_saturation}

In @Fig:fugacity_properties_Ns, we show how the leakiness, saturation, and dynamic range vary with $R$ and $\Delta \varepsilon_{RA}$ in systems with $N_S
=10$ or $N_S = 100$. An inflection point occurs where $N_S = R$,
with leakiness and dynamic range behaving differently when $R < N_S$
than when $R > N_S$. This transition is more dramatic for
$N_S = 100$ than for $N_S = 10$. Interestingly, the saturation
values consistently approach 1, indicating that full induction is easier
to achieve when multiple specific sites are present. Moreover, dynamic
range values for O1 and O2 strains with
$\Delta \varepsilon_{RA} = -15.3$ and $-13.9~k_B T$ approach 1 when
$R > N_S$, although when $N_S = 10$ there is a slight downward dip
owing to saturation values of less than 1 at high repressor copy
numbers.

![**Phenotypic properties of induction with multiple specific binding
sites.** The leakiness, saturation, and dynamic range are shown for
systems with number of specific binding sites $N_S = 10$ or
$N_S = 100$ . The vertical white line indicates the point at which
$N_S = R$.](ch6_figS6){#fig:fugacity_properties_Ns short-caption="Phenotypic properties of induction with multiple specific binding sites."}

In @Fig:fugacity_properties_Nc, we similarly show how the leakiness, saturation, and dynamic range
vary with $R$ and $\Delta \varepsilon_{RA}$ in systems with
$N_S =1$ and multiple competitor sites $N_C = 10$ or $N_C = 100$.
Each of the competitor sites has a binding energy of
$\Delta \varepsilon_C = -17.0~k_BT$. The phenotypic profiles are very similar to those for multiple specific sites shown in , with sharper transitions at
$R = N_C$ due to the greater binding strength of the competitor site.
This indicates that introducing competitors has much the same effect on
the induction phenotypes as introducing additional specific sites, as in
either case the influence of the repressors is dampened when there are
insufficient repressors to interact with all of the specific binding
sites.

![**Phenotypic properties of induction with a single specific site and
multiple competitor sites.** The leakiness, saturation, and dynamic
range are shown for systems with a single specific binding site
$N_S = 1$ and a number of competitor sites $N_C = 10$ or
$N_C = 100$. All competitor sites have a binding energy of
$\Delta \varepsilon_C = -17.0~k_BT$. The vertical white line
indicates the point at which $N_C = R$.](ch6_figS7){#fig:fugacity_properties_Nc short-caption="Phenotypic properties of induction with a single specific site and multiple competitor sites."}


This section of the appendix gives a quantitative analysis of the
nuances imposed on induction response in the case of systems involving
multiple gene copies as are found in the vast majority of studies on
induction. In these cases, the intrinsic parameters of the MWC model get
entangled with the parameters describing gene copy number.

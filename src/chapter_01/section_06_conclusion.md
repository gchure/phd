## On Molecular Biophysics and Evolutionary Dynamics

This thesis as a whole presents an attempt to understand how adaptive
processes operate in biological systems at a mechanistic level beyond
qualitative description. The thermodynamic model derived and explored in
Chapter 2 presents a concrete theoretical framework through which we can
understand how mutations and environmental perturbations influence the output
of a simple genetic circuit with quantitative precision. While the work here
specifically explores the *mean* level of gene expression of a population,
I have had the privilege to be involved in several projects which explore the
complete distribution of gene expression of various regulatory motifs using
non-equilibrium models [@razo-mejia2020; @laxhuber2020]. Both
equilibrium and non-equilibrium approaches, while differing in their
fundamental assumptions of the system, can be used to understand how the
regulation of gene expression occurs *in vivo* and should be viewed as
complimentary rather than adversarial approaches.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A combination of these types of approaches will be necessary to attack what I
believe is the next great frontier of biological physics -- predicting
evolution. While this thesis is  primarily focused 
on a single type of regulatory architecture regulating a single
promoter via a single species of transcription factor, it is worth remembering
that systems-level phenotypes are often complex and result from the concerted
action of an array of biological processes. As was mentioned in our discussion
on physiological adaptation, it has been known for nearly a century that the
bacterial growth rate is directly correlated to the total protein content of
the cell, with recent works illustrating rich phenomenology in the structure of
the bacterial proteome as a whole [@li2014a; @schmidt2016; @hui2015;
@scott2010a; @klumpp2014]. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In collaboration again with Nathan M.
Belliveau, we have begun to explore how the composition of the bacterial
proteome is structured at the single-protein level. @Fig:proteomics (A) shows
data compiled from several independent proteomic data sets (using either
quantitative mass spectrometry [@schmidt2016; @peebo2015; @valgepea2013a] or
ribosomal profiling [@li2014a]) where the abundance of different molecular
constituents of the bacterial proteome is plotted as a function of the
growth rate. These components, broken down by their functional designation
according to their Cluster of Orthologous Groups (COG) annotation
[@galperin2015a], reveal varied dependencies on the growth rate. Of note are
the COG classes "cellular processes and signaling," "metabolism," and
"information storage and processing" which all appear to have a 
correlation between the cellular growth rate and the total mass of that
proteome sector. However, when plotted as the total
*mass fraction* of the proteome instead of the total mass, a striking result
is observed. @Fig:proteomics (B) reveals a very strong, negative correlation
between the mass fraction of the proteome dedicated to information storage
and processing (including ribosomal and transcriptional machinery) and the
proteome fraction dedicated to metabolism. This direct competition for
resources between the proteins involved in translation (ribosomes, elongation
factors, etc.) and metabolic networks has been shown previously [@klumpp2008;
@scott2010a; @hui2015] and suggests a strong evolutionary constraint on how
resources can be optimally partitioned.

![**Allocation of cellular resources induces compositional structure in the *E.
coli* proteome.** (A) The total proteome mass of the five major annotated COG
categories is shown as a function of the experimental growth rate. Different
marker shapes represent different data sets. (B) The fraction by mass of the
proteome dedicated to metabolic machinery plotted as a function of the total
proteome mass dedicated to the processes of the central dogma. Different shapes
correspond to the different data sets shown in (A). Color indicates increasing
growth rate from yellow to dark blue. Data shown in this figure come from
@peebo2015 (inverted triangles),  @li2014a (triangles), @schmidt2016 (circles),
and @valgepea2013a (diamonds). The [Python code
(`ch1_fig13.py`)](https://github.com/gchure/phd/blob/master/src/chapter_01/code/ch1_fig13.py)
used to generate this figure can be accessed via the thesis [GitHub
repository](https://github.com/gchure/phd). ](ch1_fig13){#fig:proteomics
short-caption="Allocation of cellular resources induces compositional structure
in the *E. coli* proteome."}

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; As of this writing, our understanding of the
cellular resource allocation visible in @Fig:proteomics remains largely
phenomenological [@scott2014]. This is in part due to the tremendously
high-dimensional nature of systems-level organization. Our understanding of
systems with such huge numbers of degrees-of-freedom have classically benefited enormously
from the application of statistical mechanics as this thesis shows in the
context of transcriptional regulation. The quantitative framework derived and
carefully dissected in this thesis, I believe, lays the groundwork to understand
how phenomena such as that shown in @Fig:proteomics (B) arise, and perhaps more
importantly, evolve. A recent work from Michael Lässig, Ville
Mustonen, and Aleksandra Walczak entitled *Predicting evolution* [@lassig2017]
describes what the future of evolutionary theory may look like given these types
of models. Recent technological advancements in sequencing, microscopy, and
computation coupled with theoretical advancements in the biophysics of gene
regulation present an opportunity for a rich theoretical dialogue between
molecular biophysics and evolutionary dynamics coupled with experimental
dissection (@Fig:molbiophys_evoldynam).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The first 20 years of the 21$^\text{st}$ century
brought a paradigm shift in our
understanding of noise in biological networks, illustrating how
cross-disciplinary approaches to scientific discovery can solve and
even create new fields of biological inquiry [@ciechonska2016;@eldar2010b]. I hope that
some of the material described in the coming chapters can help contribute to
a systems-biology approach to evolution.

![**The coming interplay between molecular biophysics and evolutionary
dynamics.** (A) Recent progress in our understanding of the structure and function of
biological networks has resulted in many examples where high-dimensional
biological phenomena can be boiled down to effective phenomena. Future work will
draw from our understanding of these networks to place them in an evolutionary
perspective (B) where the connection between perturbations at the level of
nodes in biological networks can be drawn to fitness and evolutionary
trajectories can be predicted.](ch1_fig14){#fig:molbiophys_evoldynam
short-caption="The coming interplay between molecular biophysics and
evolutionary dynamics."}















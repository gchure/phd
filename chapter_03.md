---
# Page settings
layout: default
keywords:
comments: true
image: switchboard.png

# Hero section
title: Chapter 3
subtitle: >
    Unknown Knowns, Known Unknowns, and Unforseen Consequences: Using Free Energy Shifts To Predict Mutant Phenotypes
    
# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >

        Here, I present a biophysical model of allosteric transcriptional regulation
        that directly links the location of a mutation within a repressor to the
        biophysical parameters that describe its behavior. Here, we explore the
        phenotypic space of a repressor with mutations in either the inducer
        binding or DNA binding domains. Using the LacI repressor in <i>Escherichia</i>
        <i>coli</i>, we make sharp, falsifiable predictions and use this framework to
        generate a null hypothesis for how double mutants behave, given knowledge of
        the single mutants. Linking mutations to the parameters which govern the
        system allows for quantitative predictions of how the free energy of the
        system changes as a result, permitting coarse graining of high-dimensional
        data into a single-parameter description of the mutational consequences.
        
# Page navigation
page_nav:
    prev:
        content: Chapter II 
        url: chapter_02
    next:
        content: Chapter IV
        url: chapter_04
prefix: chapter_03
contents:
    - section_01_header
    - section_02_abstract
    - section_03_introduction
    - section_04_model
    - section_05_DNA_mutants
    - section_06_IND_mutants
    - section_07_DBL_mutants
    - section_08_discussion
    - section_09_methods
---

**Published as ...**
<hr/>
{% if page.contents %}
{% for val in page.contents %}
{% if jekyll.environment == production %}
{% include_relative {{site.doks.baseurl}}src/{{page.prefix}}/{{val}}.md %}
{% else %}
{% include_relative src/{{page.prefix}}/{{val}}.md %}
{% endif %}
{% endfor %}
{% endif %}

## References
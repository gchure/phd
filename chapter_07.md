---
# Page settings
layout: default
keywords:
comments: true
image: waves.jpg

# Hero section
title: Chapter VII 
subtitle: > 
    Supplemental Information For Chapter III: Predictive Shifts in Free Energy Couple Mutations to Their Phenotypic Consequences

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false

page_nav:
    prev:
        content: Chapter VI
        url: phd/chapter_06
    next:
        content: Chapter VIII
        url: phd/chapter_08
prefix: chapter_07
contents:
    - section_01_header
    - section_02_nonmonotonicity
    - section_03_DNA_epRA_inference
    - section_04_delta_F_inference
    - section_05_DNA_mutants
    - section_06_IND_inference
    - section_07_IND_mutants
    - section_08_param_comparison
    - section_09_global_fit
    - section_10_generality
    - section_11_strains
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
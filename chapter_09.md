---
# Page settings
layout: default
keywords:
comments: true
image: waves.jpg

# Hero section
title: Chapter IX
subtitle: > 
    Supplemental Information For Chapter V: How Bacteria Survive Osmotic Shocks

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false

page_nav:
    prev:
        content: Chapter VII
        url: chapter_07
    next:
        content: Chapter IX
        url: chapter_09
prefix: chapter_09
contents:
    - section_01_header
    - section_02_electrophysiology
    - section_03_cal_factor_derivation
    - section_04_survival_classification
    - section_05_logistic_regression
    - section_06_shock_classification
    - section_07_vandenberg_comparison
    - section_08_strains
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
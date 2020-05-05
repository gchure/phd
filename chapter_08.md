---
# Page settings
layout: default
keywords:
comments: true
image: prismatic.jpg

# Hero section
title: Chapter VIII
subtitle: > 
    Supplemental Information For Chapter IV: The Physiological Adaptability of a
    Simple Genetic Circuit 

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
prefix: chapter_08
contents:
    - section_01_header
    - section_02_growth_rates
    - section_03_cell_volume
    - section_04_counting_repressors
    - section_05_carbon_inference
    - section_06_entropy
    - section_07_strains_recipes
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
---
# Page settings
layout: default
keywords:
comments: true
image: prismatic.jpg

# Hero section
title: Chapter IV
subtitle: > 
   On The Physiological Adaptability of a Simple Genetic Circuit 
# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
        Cells adapt to changing environmental conditions by repressing or activating
        gene expression from enormous fractions of their genome, drastically changing
        the molecular composition of the cell. This requires the concerted adaptation
        of transcription factors to the environmental signals, leading to binding or
        releasing of their cognate sequences. Here, we dissect a well characterized
        genetic circuit in a number of physiological states, make predictions of the
        response, and measure how the copy number of a regulator and its gene target
        are affected. We find the parameters defining the regulators behavior are
        remarkably robust to changes in the nutrient availability, but are 
        susceptible to temperature changes. We quantitatively explore these two effects and
        discuss how they challenge current models of transcriptional regulation.
# Page navigation
page_nav:
    prev:
        content: Chapter III
        url: chapter_03
    next:
        content: Chapter 2
        url: chapter_05
prefix: chapter_04
contents:
    - section_01_header
    - section_02_abstract
    - section_03_introduction
    - section_04_model
    - section_05_experimental_design
    - section_06_scaling
    - section_08_temperature
    - section_09_discussion
    - section_10_methods
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
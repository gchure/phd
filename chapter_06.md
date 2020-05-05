---
# Page settings
layout: default
keywords:
comments: true
image: waves.jpg

# Hero section
title: Chapter VI 
subtitle: > 
    Supplemental Information For Chapter I: Signal Processing Via Allosteric Transcription Factors 

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false

page_nav:
    prev:
        content: Chapter V
        url: phd/chapter_05
    next:
        content: Chapter VII
        url: phd/chapter_07
prefix: chapter_06
contents:
    - section_01_header
    - section_02_epAI
    - section_03_fugacity
    - section_04_flow
    - section_05_microscopy
    - section_06_sensitivity
    - section_07_global_fit
    - section_08_Oid
    - section_09_comparison
    - section_10_properties
    - section_11_applications
    - section_12_tables
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
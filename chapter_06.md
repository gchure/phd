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
        content: Published Content
        url: '{{site.doks.baseurl}}/published_content'
    next:
        content: Chapter 2
        url: '{{site.doks.baseurl}}/chapter_02'
prefix: chapter_06
contents:
    - section_01_header
    - section_02_epAI
    - section_03_fugacity
    - section_04_flow
    - section_05_microscopy
    - section_06_sensitivity
    - section_07_global_fit
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
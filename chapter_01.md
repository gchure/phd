---
# Page settings
layout: default
keywords:
comments: true
image: migration.jpg
# Hero section
title:  Chapter I
subtitle: >  
    The Phenomenon of Adaptation Across The Biological Scales 

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
       All forms of life are united in their obedience to the forces of evolution.
       The ability to adapt -- either at the level of molecules or at the level of
       populations -- is central to evolution. This chapter outlines the focus of
       my dissertation and places the results in a historical context.
# Page navigation
page_nav:
    prev:
        content: Abstract
        url: abstract
    next:
        content: Chapter 2
        url: chapter_02
prefix: chapter_01
contents:
    - section_01_introduction
    - section_02_molecular_adaptation
    - section_03_evolutionary_adaptation
    - section_04_physiological_adaptation
    - section_05_survival
    - section_06_conclusion
---

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
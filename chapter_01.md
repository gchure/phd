---
# Page settings
layout: default
keywords:
comments: true
image: aspens.jpg

# Hero section
title:  Chapter I
subtitle: >  
    Janus-Faced Molecules and Adaptation Across Biological Scales     

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
       All forms of life are united in their obedience to the forces of evolution.
       The ability to adapt -- either at the level of molecules or at the level of
       populations -- is central to evolution. Blah Blah 
# Page navigation
page_nav:
    prev:
        content: Published Content
        url: '{{site.doks.baseurl}}/published_content'
    next:
        content: Chapter 2
        url: '{{site.doks.baseurl}}/chapter_2'
prefix: chapter_01
contents:
    - section_01_introduction
    - section_02_monod
    - section_03_molecular_adaptation
    - section_04_data_collapse
    - section_05_evolutionary_adaptation
    - section_06_physiological_adaptation
    - section_07_survival
    - section_08_conclusion
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


Blank placeholder
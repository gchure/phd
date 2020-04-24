---
# Page settings
layout: default
keywords:
comments: true
image: waves.jpg

# Hero section
title: Chapter V
subtitle: > 
    'Water, Water Everywhere, Nor Any Drop to Drink': How Bacteria Adapt To
    Changes in Osmolarity

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
        Mechanosensitive (MS) channels are transmembrane protein complexes which open
        and close in response to changes in membrane tension as a result of osmotic
        shock. Despite extensive biophysical characterization, the contribution of
        these channels to cell survival remains largely unknown. In this work, we use
        quantitative video microscopy to measure the abundance of a single species of
        MS channel in single cells followed by their survival after a large osmotic
        shock. We observe total death of the population with less than ~100 channels
        per cell and determine that approximately 500 - 700 channels are needed for
        80% survival. The number of channels we find to confer nearly full survival
        is consistent with the counts of the number of channels in wild type cells in
        several earlier studies. These results prompt further studies to dissect the
        contribution of other channel species to survival.
# Page navigation
page_nav:
    prev:
        content: Published Content
        url: '{{site.doks.baseurl}}/published_content'
    next:
        content: Chapter 2
        url: '{{site.doks.baseurl}}/chapter_2'
prefix: chapter_05
contents:
    - section_01_header
    - section_02_abstract
    - section_03_introduction
    - section_04_results
    - section_05_discussion
    - section_06_methods
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
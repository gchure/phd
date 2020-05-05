---
# Page settings
layout: default
keywords:
comments: true
image: switchboard.png

# Hero section
title: Chapter 2
subtitle: >
    Through The Intramolecular Grapevine: Cellular Decision Making Via Allosteric
    Transcription factors

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
        Allosteric regulation is found across all domains of life, yet we still lack
        simple, predictive theories that directly link the experimentally tunable
        parameters of a system to its input-output response. This chapter presents
        a general theory of allosteric transcriptional regulation using the
        Monod-Wyman-Changeux model. We rigorously test this model using the
        ubiquitous simple repression motif in bacteria by first predicting the
        behavior of strains that span a large range of repressor copy numbers and DNA
        binding strengths and then constructing and measuring their response.

# Page navigation
page_nav:
    prev:
        content: Chapter I
        url: phd/chapter_01
    next:
        content: Chapter III
        url: phd/chapter_03
prefix: chapter_02
contents:
    - section_01_header
    - section_02_abstract
    - section_03_introduction
    - section_04_model   
    - section_05_results
    - section_06_discussion
    - section_07_methods
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
# Visual Compositional Data Analytics for Spatial Transcriptomics
David Hägele, Yuxuan Tang, Daniel Weiskopf
***
This repository contains the code for the visual analytics system prototype that was proposed as a solution for the [2024 Bio+MedVis Challenge](http://biovis.net/2024/biovisChallenges_vis/).
The software is a [Bokeh server](https://docs.bokeh.org/en/latest/docs/user_guide/server.html) application written in Python that provides a web front-end.

## Challenge Description
Spatial transcriptomics technology can detect cell types at different locations on cellular tissue. 
Due to limited resolution, a mixture of cell types is detected per *spot*.
However, the proportions of the individual cell types for an individual spot can be determined, resulting in a composition such as `[80% type1, 10% type2, 4% type3, ...]`.
This data can be visualized using pie chart glyphs superimposed onto a histological image, so analysts can relate the cell type proportions to locations on the tissue and observe areas of similar patterns.
Such a visualization has certain limitations and therefore the challenge asks for a redesign.

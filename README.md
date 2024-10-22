# Visual Compositional Data Analytics for Spatial Transcriptomics
David Hägele, Yuxuan Tang, Daniel Weiskopf
***
This repository contains the code for the visual analytics system prototype that was proposed as a solution for the [2024 Bio+MedVis Challenge](http://biovis.net/2024/biovisChallenges_vis/).
The software is a [Bokeh server](https://docs.bokeh.org/en/latest/docs/user_guide/server.html) application written in Python that provides a web front-end.

| Challenge Description | Visualization to Redesign |
| --- | --- |
|Spatial transcriptomics technology can detect cell types at different locations on cellular tissue. Due to limited resolution, a mixture of cell types is detected per *spot*. However, the proportions of the individual cell types for an individual spot can be determined, resulting in a composition such as `[80% type1, 10% type2, 4% type3, ...]`. This data can be visualized using pie chart glyphs superimposed onto a histological image, so analysts can relate the cell type proportions to locations on the tissue and observe areas of similar patterns. Such a visualization has certain limitations and therefore the challenge asks for a redesign. | ![scatter-pies-celltypes](https://github.com/user-attachments/assets/4619661d-b2fe-4eb3-b05c-ba4525f1e956) |


## Proposed Redesign
We propose a visual analytics system to explore the cell type compositions and relate them to the histological image of the tissue.
![image](https://github.com/user-attachments/assets/fe4ad2fb-e042-4732-af0d-20e5e03a61a7)
There are 3 views that support brushing and linking, i.e., selections made in one view are reflected in the other views.
- Histological image view - shows tissue and locations of spots (toggleable) which can be highlighted when a selection of spots is issued.
- Stacked bar chart of cell type mixtures - shows the cell type proportions of selected spots. 
- Dimensionality reduction (similarity) of cell type mixtures - using PCA of the proportions in Aitchison gemotry shows similar mixtures being grouped into blobs.
  - additional k-means clustering for color coding.



### Setup Instructions
To set up the project you need an up to date Python 3 installation.
Then you can use the bash scripts to set up and run the server application.
```bash
# set up a virtual python environment
./setup_venv.sh
# install the dependencies
./setup_dependencies.sh
# start the server
./start_server.sh
```

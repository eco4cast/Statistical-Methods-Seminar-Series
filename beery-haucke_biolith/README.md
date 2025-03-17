# Occupancy Modelling in Python

*Presented by [Sara Beery](https://beerys.github.io/) and [Timm Haucke](https://timm.haucke.xyz/)*

Incorporating machine learning into statistical ecological models is crucial in order to effectively and efficiently leverage large datasets that are infeasible to label manually. However, there is a disconnect between the tools statistical ecologists use and those used for machine learning. To help bridge disconnect, we propose a new toolbox for Bayesian ecological modeling called [Biolith](https://github.com/timmh/biolith). Biolith is written entirely in Python and supports a range of occupancy models. Using Biolith, we present an exemplary workflow to automatically classify species in a camera trap dataset, calibrate a confidence threshold, and perform occupancy modeling on the resulting detection / non-detection data.

## Running the Code

We provide our code as a Jupyter notebook that you can [download](beery-haucke_biolith/OccupancyModellingInPython.ipynb) and run on your own computer, provided you have Python and [JupyterLab](https://jupyter.org/) setup. Alternatively, you can run the provided notebook without downloading anything on Google Colab, which already provides the required environment.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1828fk-7DEsDL9reK5oYSOrsYA68cim-W)

# Pathway Enrichment Analysis

This repository contains the Python code for conducting pathway enrichment analysis using custom pathway annotations. The analysis is focused on identifying regulated pathways based on gene set enrichment analysis (GSEA) and correlation matrices.

## Getting Started

### Prerequisites

The code is written in Python and requires the following libraries:
- pandas
- gseapy
- matplotlib
- seaborn
- numpy

Jupyter Notebooks are also provided for interactive analysis:

- `Overall_Pathway_Enrichment_mean_033122.ipynb`: For pathway enrichment analysis using WormPaths.
- `Overall_pathway_enrichment_ClusterSets_062822_3_2.ipynb`: For enrichment analysis using smaller cluster sets (deepSplit=2, minclustersize=3).
- `Overall_pathway_enrichment_ClusterSets_062822_6_3.ipynb`: For enrichment analysis using larger cluster sets (deepSplit=3, minclustersize=6).

To open a notebook, run:

```bash
jupyter notebook jupyter_notebooks/<notebook_name>.ipynb

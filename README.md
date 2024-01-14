# **Warning- This repository is under development as we are in the process of cleaning up and organising codes. There might be missing files.**
# Systems-level Transcriptional Regulation of Metabolism in *C.elegans*

## Table of Contents
1. [Introduction](#Introduction)
2. [System Components](#System-Components)
3. [Features](#Features)
4. [Requirements](#Requirements)
5. [Installation](#Installation)
6. [Authors](#Authors)
7. [Contributing and Contact](#Contributing-and-Contact)


## Introduction

Welcome to our project! We've constructed a computational pipeline aimed at elucidating the transcriptional regulation of metabolism in *C.elegans* at a systems level. Primarily written in Python, this pipeline also utilizes components of MATLAB and shell script. It has been rigorously developed for a development dataset, tissue dataset, and a compendium of 177 expression datasets. Importantly, our pipeline can be applied to any expression dataset, whether it's RNA-Seq or microarray.

## System Components

Our pipeline comprises the following key modules:

### 1) Determining the Extent of Transcriptional Regulation of Metabolism
* During Development - Utilizing Variation Score (VS)
* Across Tissues - Utilizing Coefficient of Variation (CV)
* Across a Compendium of Expression Datasets - Utilizing CV 

### 2) Identifying the Prevalence of Transcriptional Regulation at the Pathway Level
* Supervised Approach
* [Unsupervised Approach](https://github.com/WalhoutLab/WormClust/blob/master/2_b_UnsupervisedApproach/Unsupervised_Analysis.md)

### 3) Uncovering Activation/Repression Conditions of Metabolic Sub-Pathways

### 4) WormClust: A Gene Query Web Application
* Allows gene-by-gene queries of all *C.elegans* genes to associate them with metabolic (sub)-pathways.
* For all iCEL genes in the metabolic network model, it generates a clustered heatmap of the query gene with other closely associated metabolic network genes based on coflux and coexpression.
* For all non-iCEL genes, it identifies the pathway enrichment of closely associated metabolic network model genes.

## Features

* **Python-Centric:** Our pipeline is primarily written in Python, ensuring high readability and maintainability.
* **Broad Applicability:** The methodologies implemented in this pipeline can be extrapolated to any organism for which large gene expression profile compendia and high-quality metabolic network models are available, including humans.
* **Interactive Web Tool - WormClust:** This feature enables users to assess the association of a specific gene with the metabolic network based on similarities in gene expression. It is consistent with the methodology proposed by [Nanda et al., 2023](https://www.embopress.org/doi/full/10.15252/msb.202211443).

## Requirements

Please refer to the attached `requirements.txt` file for software dependencies.

## Installation

The repository can be cloned and the dependencies installed as follows:

1. Clone the repository: `git clone https://github.com/WalhoutLab/WormClust.git`
2. Navigate to the cloned directory: `cd WormClust`
3. Install the required dependencies: `pip install -r requirements.txt`
4. Alternatively, Create a Conda environment from the requirements.txt file: ```shell conda create --name myenv --file requirements.txt


## Authors

* Shivani Nanda - [GitHub](https://github.com/shivani710)
* Safak Yilmaz - Email: LutfuSSafak.Yilmaz@umassmed.edu

## Contributing and Contact

For bug reports, contributions, or further queries, please reach out to us. Email communication is preferred. Contact Shivani Nanda (shivani.nanda@umassmed.edu) or Safak Yilmaz (LutfuSSafak.Yilmaz@umassmed.edu). 

Project Link: [https://github.com/WalhoutLab/WormClust](https://github.com/WalhoutLab/WormClust)

WormClust Web Server Link: [http://wormflux.umassmed.edu/WormClust/wormclust.php](http://wormflux.umassmed.edu/WormClust/wormclust.php)


  
  

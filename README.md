# <b>Systems-level transcriptional regulation of metabolism in *C.elegans*
  ## **Introduction**

We developed a computational pipeline to unravel transcriptional regulation of metabolism in *C.elegans* at a systems-level. It is written mainly in Python, with parts written in MATLAB and shell script. It has been developed for a development dataset, tissue dataset and a compendium of 177 expression datasets. However, this pipeline is applicable to any expression dataset- whether RNA-Seq or microarray. It consists of following major parts:
  
  ### 1) Finding extent of transcriptional regulation of metabolism
    a) During Development- Using Variation Score (VS)
    b) Across tissues- Using Coefficient of Variation (CV)
    c) Across compendium of expression datasets - Using CV 
    
  ### 2) Finding prevalence of transcriptional regulation at pathway level
    a) Supervised approach
    b) Unsupervised approach
    
  ### 3) Finding activation/ repression conditions of metabolic sub-pathways
  ### 4) WormClust: a web application that enables gene-by-gene query of all *C.elegans* genes to associate them to metabolic (sub)-pathways. 
    a)For all iCEL genes that are part of metabolic network model, obtain a clustered heatmap of the query gene with other closely associated metabolic network genes based on coflux and coexpression
    b) For all non-iCEL genes, find the pathway enrichment to closely associated metabolic network model genes
  
  
  ## **Features**
### 1) Written in python
### 2) Provides WormClust- an online tool to evaluate the association of a given gene with the metabolic network based on similarities in gene expression, according to Nanda et al., 2022 (in review). 

  
  
  ## **Requirements**
  numpy (>=1.12.1)
  scipy (>=0.19.0)
  pandas (>=0.20.1)
  matplotlib (>=2.0.2)
  scikit-learn (>=0.18.1)
  seaborn
  
  ## **Contact or Contributions**
  Please contact us for reporting bugs or contributing purposes. Email contacting is encouraged. Please send to Shivani Nanda or Safak Yilmaz.
  Project link: https://github.com/WalhoutLab/WormClust
  

  
  
  

# <b>Systems-level transcriptional regulation of metabolism in *C.elegans*

 ## ** Warning: This repository is still under development. There might be missing files as we are in the process of organising and cleaning up codes.
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
    a) For all iCEL genes that are part of metabolic network model, obtain a clustered heatmap of the query gene with other closely associated metabolic network genes based on coflux and coexpression
    b) For all non-iCEL genes, find the pathway enrichment to closely associated metabolic network model genes
  
  
  ## **Features**
### 1) Written in python
### 2) Approches can be applied to any other organism, for which large gene expression profile compendia and high-quality metabolic network models are available, including humans.
### 3) Provides WormClust- an online tool to evaluate the association of a given gene with the metabolic network based on similarities in gene expression, according to Nanda et al., 2022 (in review). 

  
  ## **Requirements**
  Please find attached the requirements.txt file.
  
  
  ## **Authors**
  Shivani Nanda- https://github.com/shivani710
  
  
  ## **Contact or Contributions**
  Please contact us for reporting bugs or contributing purposes. Email contacting is encouraged. Please send to Shivani Nanda(shivani.nanda@umassmed.edu) or Safak Yilmaz(LutfuSSafak.Yilmaz@umassmed.edu).<br/>Project link: https://github.com/WalhoutLab/WormClust 
  <br> WormClust web server link- http://wormflux.umassmed.edu/WormClust/wormclust.php
  
  

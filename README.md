# <b>Systems-level transcriptional regulation of metabolism in *C.elegans*
  ## **Introduction**

We developed a computational pipeline to unravel transcriptional regulation of metabolism in *C.elegans* at a systems-level. It is written mainly in Python, with parts written in MATLAB and shell script. It has been developed for a development dataset, tissue dataset and a compendium of 177 expression datasets. However, this pipeline is applicable to any expression dataset- whether RNA-Seq or microarray. It consists of following major parts:
  
  ### 1) Finding extent of transcriptional regulation of metabolism
    a) During Development- Using Variation Score (VS)
    b) Across tissues- Using Coefficient of Variation (CV)
    c) Across compendium of expression datasets - Using CV 

  ## **Features**
  
  
  ## **Requirements**
  numpy (>=1.12.1)
  scipy (>=0.19.0)
  pandas (>=0.20.1)
  matplotlib (>=2.0.2)
  scikit-learn (>=0.18.1)
  seaborn

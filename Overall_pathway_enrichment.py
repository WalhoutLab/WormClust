#!/usr/bin/env python
# coding: utf-8

# ## Importing modules

# In[59]:


import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
from gseapy.parser import Biomart
import os
import numpy as np
import seaborn as sns
from gseapy.plot import gseaplot


# ## Setting base directory

# In[60]:


Base_dir='/data/nandas/Combined_coexp/Pathway_enrichment'
os.chdir(Base_dir)


# ## Reading required files: GeneSets(gmt), PathwayToGenes and Gene Correlations 

# In[61]:


pathway_filename = 'Pathway2Gene_2.gmt';
metabolic_corr_df=pd.read_csv("normalized_metaboliccorr_df.csv",index_col=0,header='infer')
Pathway_df=pd.read_csv(pathway_filename,index_col=0,sep='\t')


# ### Setting default coregulated state of pathway

# In[62]:


Pathway_df['IsRegulated'] = False


# In[63]:


Pathway_df['IsRegulated']


# ## PreRank Gene set enrichment analyses for custom pathway annotations

# In[64]:


def PreRank(genes, outdir):
#     print("Genes: {}".format(genes));
#     print("Length of genes:{}".format(len(genes)))
    intersection_list = list(set(metabolic_corr_df.index).intersection(set(genes)))
    missing_genes=list(set(genes).difference(set(intersection_list)))
#     print("IntersectionList: {}".format(intersection_list));
#     print("Length of intersection list:{}".format(len(intersection_list)))
    print('Missing genes:{}\n{}'.format(len(missing_genes),missing_genes))
    Combined=metabolic_corr_df[intersection_list];
    Mean=Combined.mean(axis=1)
    rnk=Mean.sort_values(ascending=False)
#     print("Rank: {}".format(rnk))    
    pre_res = gp.prerank(rnk=rnk, gene_sets=pathway_filename, processes=4,min_size=1, outdir=outdir, format='png', weighted_score_type=1,verbose=True)
    return pre_res

def _is_regulated_pathway_(pre_res, pathway):
    pathway_pre_res = pre_res.res2d.loc[pathway];
    is_regulated_pathway = pathway_pre_res.es >= 0.70 and pathway_pre_res.fdr <= 0.05
#     is_regulated_pathway =  pathway_pre_res.fdr <= 0.05
    return is_regulated_pathway;

def PlotEnrichment(pre_res,pathway, outdir):
    Sorted_values=pre_res.res2d.sort_values(ascending=False,by=['es'])[0:15]
    fig = plt.figure(figsize=(12,12))
    df = pd.DataFrame({'Enrichment Score': Sorted_values.es,
                   'p-value': Sorted_values.pval,'FDR':Sorted_values.fdr}, index=Sorted_values.index)
    ax = df.plot.barh(rot=0)
    plt.legend(loc='best', bbox_to_anchor=(1, 1))
    plt.savefig("{}/{}_plot.png".format(outdir, pathway))
    plt.show()
def PlotGSEA(pre_res, pathway, outdir):
    terms = pre_res.res2d.sort_values(by=['es'],ascending=False).index
    fig=gseaplot(rank_metric=pre_res.ranking, term=pathway, **pre_res.results[pathway],
                 ofname='{}/{}_gsea.png'.format(outdir,pathway))
    
    
    


# In[65]:


Pathway_df_withoutIsRegulated = Pathway_df.drop(['IsRegulated'], axis=1);


# In[66]:


New_df = pd.DataFrame([])
for pathway in Pathway_df.index:
    print(pathway)
#     pathway = 'ALA_ASP_AND_GLU_METABOLISM';
    genes = list(Pathway_df_withoutIsRegulated.loc[pathway].dropna());
    pre_res = PreRank(genes, pathway);
    Pathway_df.at[pathway, 'IsRegulated'] = _is_regulated_pathway_(pre_res, pathway);
    print("{} is regulated:{}".format(pathway,_is_regulated_pathway_(pre_res, pathway)))
    PlotEnrichment(pre_res, pathway, outdir=pathway)
    PlotGSEA(pre_res, pathway,pathway)
    gsea_result_df=pre_res.res2d.loc[pathway];
    New_df=New_df.append(gsea_result_df)
#     print(New_df)
Pathway_df.to_csv("Pathway_Regulation_status.csv")
New_df.to_csv("Final_pathway_gsea.csv")




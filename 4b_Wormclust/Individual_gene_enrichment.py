#!/usr/bin/env python
# coding: utf-8

# ## Importing modules

# In[10]:


import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
from gseapy.parser import Biomart
import os
import numpy as np
import seaborn as sns
from gseapy.plot import gseaplot
from scipy import stats as st
import fnmatch


# ## Setting base directory

# In[11]:


Base_dir='/data/nandas/Combined_coexp/Pathway_enrichment/ORgenes_ignore_negative'
os.chdir(Base_dir)


# ## Reading required files: GeneSets(gmt), PathwayToGenes and Gene Correlations 

# In[12]:


pathway_filename = '../Pathway2Gene_2.gmt';
metabolic_corr_df=pd.read_csv("../normalized_metaboliccorr_df.csv",index_col=0,header='infer')
Pathway_df=pd.read_csv(pathway_filename,index_col=0,sep='\t')


# ### Setting default coregulated state of pathway

# In[13]:


Genes=metabolic_corr_df
Genes ['IsRegulated']=False
Genes=Genes['IsRegulated']


# In[14]:


Genes=pd.DataFrame(Genes)


# In[15]:


# Pathway_df['IsRegulated'] = False


# In[16]:


# Pathway_df['IsRegulated']


# ## PreRank Gene set enrichment analyses for custom pathway annotations

# In[17]:


def PreRank(gene, outdir):
    Combined=metabolic_corr_df[gene];
    Combined.dropna(inplace=True)
#     print(Mean)
    rnk=Combined.sort_values(ascending=False)
    #print("Rank: {}".format(rnk))    
    pre_res = gp.prerank(rnk=rnk, gene_sets=pathway_filename, processes=4,min_size=2, outdir=outdir, format='png', weighted_score_type=1,verbose=True)
    return pre_res

def _is_regulated_pathway_(pre_res, gene):
#     print('Hello There: {}'.format(pre_res));
#     print('Shivani Here: {}'.format(pre_res.res2d))
    pathway_pre_res = pre_res.res2d.sort_values(by=['fdr'],ascending=True);
#     is_regulated_pathway = pathway_pre_res.es >= 0.70 and pathway_pre_res.fdr <= 0.05
    is_regulated_pathway=False;
    for index in pathway_pre_res.index:
        val = pathway_pre_res.loc[index];
#         print(val)
        if(val.fdr <= 0.05 and val.nes > 0 and val.es >0):
            is_regulated_pathway = True;
#             print("Changed IsRegulatedValue to True");
#             print("FDR: {}---NES: {}---: ES:{}".format(val.fdr, val.nes, val.es))
            break;
    print("OutsideLoop: {}".format(is_regulated_pathway));
    return is_regulated_pathway;

def PlotEnrichment(pre_res,gene, outdir):
    Sorted_values=pre_res.res2d.sort_values(ascending=False,by=['nes'])[0:20]
    Sorted_values.to_csv("Pathways_enriched_{}.csv".format(gene))
    fig = plt.figure()
    df = pd.DataFrame({'Enrichment Score': Sorted_values.es,
                   'p-value': Sorted_values.pval,'FDR':Sorted_values.fdr}, index=Sorted_values.index)
    df2=pd.DataFrame(Sorted_values.nes,index=Sorted_values.index)
    ax = df.plot.barh(rot=0)
    ax.invert_yaxis()
    plt.legend(loc='best', bbox_to_anchor=(1, 1))
    plt.savefig("{}/{}_plot.png".format(outdir, gene))
    plt.show()
    plt.close()
    ax2=df2.plot.barh(rot=0)
    ax2.invert_yaxis()
    plt.legend(loc='best', bbox_to_anchor=(1, 1))
    plt.xlabel('Normalized enrichment score')
    plt.savefig("{}/{}_nes.png".format(outdir, gene))
    plt.show()
    plt.close()
    
def PlotGSEA(pre_res, gene, outdir):
    terms = pre_res.res2d.sort_values(by=['nes'],ascending=False).index
    fig=gseaplot(rank_metric=pre_res.ranking, term=terms[0], **pre_res.results[terms[0]])
    plt.savefig('{}/{}_gsea.png'.format(outdir,gene))
    
def ConvertPairsToMatrix(bayesian_metabol_df):
    a = np.unique(bayesian_metabol_df['Gene1'])
    b = np.unique(bayesian_metabol_df['Gene2'])
    c = np.union1d(a,b);
    data = np.zeros((len(c), len(c)));
    output_df = pd.DataFrame(data, index=c, columns=c)
    for values in bayesian_metabol_df.values: 
        output_df[values[0]][values[1]] = values[2];
        output_df[values[1]][values[0]]=values[2];
    np.fill_diagonal(output_df.values,1)
    return output_df
    
    


# In[18]:


New_df = pd.DataFrame([])
for genes in metabolic_corr_df.index:
#     pathway = 'PROPIONATE_SHUNT';
    print(genes)
#     pathway = 'ALA_ASP_AND_GLU_METABOLISM';
#     genes = list(Pathway_df_withoutIsRegulated.loc[pathway].dropna());
    file_exist = False;
    for file in os.listdir('./'):
        if fnmatch.fnmatch(file, "Pathways_enriched_{}.csv".format(genes)):
            print("File: {} found, skipping!!!".format(file))
            file_exist = True;
#     pre_res = PreRank(genes, genes);
    if(not file_exist):
        pre_res = PreRank(genes, genes);
        if(pre_res is None):
            continue;
        else:
            Genes.at[genes, 'IsRegulated'] = _is_regulated_pathway_(pre_res, genes);
            print("{} is regulated:{}".format(genes,_is_regulated_pathway_(pre_res, genes)))
            PlotEnrichment(pre_res, genes, outdir=genes)
            PlotGSEA(pre_res=pre_res,gene=genes,outdir=genes)
Genes.to_csv("Genes_Regulation_status.csv")
New_df.to_csv("Final_pathway_gsea.csv")




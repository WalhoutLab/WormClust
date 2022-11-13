#!/usr/bin/env python
# coding: utf-8

# ## Importing modules

# In[1]:


import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
from gseapy.parser import Biomart
import os
import numpy as np
import seaborn as sns
from gseapy.plot import gseaplot
from scipy import stats as st
import random


# ## Setting base directory

# In[2]:


Base_dir='/data/nandas/Combined_coexp/Pathway_enrichment/Gmean/test'
os.chdir(Base_dir)


# ## Reading required files: GeneSets(gmt), PathwayToGenes and Gene Correlations 

# In[3]:


pathway_filename = '../../Pathway2Gene_2.gmt';
metabolic_corr_df=pd.read_csv("../../normalized_metaboliccorr_df.csv",index_col=0,header='infer')
Pathway_df=pd.read_csv(pathway_filename,index_col=0,sep='\t')


# ### Setting default coregulated state of pathway

# In[4]:


Pathway_df['IsRegulated'] = False


# In[ ]:


Pathway_df['IsRegulated']


# ## PreRank Gene set enrichment analyses for custom pathway annotations

# In[ ]:


def shuffle(df):
    a = metabolic_corr_df.index
    # print(a)
    # print(len(a))
    indices = random.sample(range(len(a)), len(a));
    metabolic_corr_df['NewCol']=indices
    metabolic_corr_df['NewIndices']=metabolic_corr_df.iloc[indices].index
    shuffled_df=metabolic_corr_df.set_index('NewIndices')
    shuffled_df.drop(columns=['NewCol'],inplace=True)
    # shuffled_df = metabolic_corr_df.iloc[indices,indices]
    # return shuffled_df
    shuffled_df=shuffled_df.transpose()
    shuffled_df['NewCol']=indices
    shuffled_df['NewIndices']=metabolic_corr_df.iloc[indices].index
    shuffled_df.set_index('NewIndices',inplace=True)
    shuffled_df.drop(columns=['NewCol'],inplace=True)
    return shuffled_df

def PreRank(genes, outdir):
    
    print("Length of genes:{}".format(len(genes)))
    intersection_list = list(set(metabolic_corr_df.index).intersection(set(genes)))
    missing_genes=list(set(genes).difference(set(intersection_list)))
#     print("IntersectionList: {}".format(intersection_list));
#     print("Length of intersection list:{}".format(len(intersection_list)))
    print('Missing genes:{}\n{}'.format(len(missing_genes),missing_genes))
    if(len(missing_genes) > 0.6* len(genes)):
        return;
    Combined=metabolic_corr_df[intersection_list];
#     print(Combined)
    Combined=shuffle(Combined)
#     print(Combined)
    Mean=np.exp(np.log(Combined.prod(axis=1))/Combined.notna().sum(1))
    Mean.dropna(inplace=True)
#     print(Mean)
    rnk=Mean.sample(frac=1)
    rnk=rnk.sort_values(ascending=False)
#     print("Rank: {}".format(rnk))    
    pre_res = gp.prerank(rnk=rnk, gene_sets=pathway_filename, processes=4,min_size=2, outdir=outdir, format='png', weighted_score_type=1,verbose=True)
    return pre_res

def _is_regulated_pathway_(pre_res, pathway):
#     print('Hello There: {}'.format(pre_res));
#     print('Shivani Here: {}'.format(pre_res.res2d))
    if(pathway not in pre_res.res2d.index):
        return "NaN"
    pathway_pre_res = pre_res.res2d.loc[pathway];
#     is_regulated_pathway = pathway_pre_res.es >= 0.70 and pathway_pre_res.fdr <= 0.05
    is_regulated_pathway =  pathway_pre_res.fdr <= 0.05
    return is_regulated_pathway;

def PlotEnrichment(pre_res,pathway, outdir):
    Sorted_values=pre_res.res2d.sort_values(ascending=False,by=['nes'])[0:20]
    fig = plt.figure()
    df = pd.DataFrame({'Enrichment Score': Sorted_values.es,
                   'p-value': Sorted_values.pval,'FDR':Sorted_values.fdr}, index=Sorted_values.index)
    ax = df.plot.barh(rot=0)
    plt.legend(loc='best', bbox_to_anchor=(1, 1))
    plt.savefig("{}/{}_plot.png".format(outdir, pathway))
    plt.show()
def PlotGSEA(pre_res, pathway, outdir):
    terms = pre_res.res2d.sort_values(by=['es'],ascending=False).index
    fig=gseaplot(rank_metric=pre_res.ranking, term=pathway, **pre_res.results[pathway])
    plt.savefig('{}/{}_gsea.png'.format(outdir,pathway))
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
    
    


# In[ ]:


Pathway_df_withoutIsRegulated = Pathway_df.drop(['IsRegulated'], axis=1);


# In[ ]:


New_df = pd.DataFrame([])
for pathway in Pathway_df.index:
#     pathway = 'PROPIONATE_SHUNT';
    print(pathway)
#     pathway = 'ALA_ASP_AND_GLU_METABOLISM';
    genes = list(Pathway_df_withoutIsRegulated.loc[pathway].dropna());
    for i in range(0,100):
        pre_res = PreRank(genes, pathway);
        if(pre_res is None):
            continue; 
        Number=i
        Pathway_df.at[pathway, i] = _is_regulated_pathway_(pre_res, pathway);
        print("{} is regulated:{}".format(pathway,_is_regulated_pathway_(pre_res, pathway)))
#     PlotEnrichment(pre_res, pathway, outdir=pathway)
#     if(pathway in pre_res.res2d.index):
#         PlotGSEA(pre_res, pathway,pathway)
#         gsea_result_df=pre_res.res2d.loc[pathway];
#         New_df=New_df.append(gsea_result_df)
Pathway_df.to_csv("Pathway_Regulation_status.csv")
# New_df.to_csv("Final_pathway_gsea.csv")




# In[ ]:


Pathway_df


# In[ ]:


metabolic_corr_df.loc['aap-1']['ace-1']


# In[ ]:


shuffled_df.loc['aap-1']['ace-1']


# In[ ]:


status=pd.read_csv("Pathway_Regulation_status.csv",index_col=0)


# In[ ]:


status=status['IsRegulated']


# In[ ]:


status.replace('False','Non-coregulated')


# In[ ]:


status.IsRegulated.groupby(status.IsRegulated).count()


# In[ ]:


my_labels='Non-Regulated','Regulated'
sums = status.IsRegulated.groupby(status.IsRegulated).count()
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct


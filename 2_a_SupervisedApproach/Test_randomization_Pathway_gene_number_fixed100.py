#!/usr/bin/env python
# coding: utf-8

# ## Importing modules

# In[84]:


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

# In[85]:


Base_dir='/data/nandas/Combined_coexp/Pathway_enrichment/Gmean/Pathway_gene_randomization2'
os.chdir(Base_dir)


# ## Reading required files: GeneSets(gmt), PathwayToGenes and Gene Correlations 

# In[86]:


pathway_filename = '../../Pathway2Gene_2.gmt';
metabolic_corr_df=pd.read_csv("../../normalized_metaboliccorr_df.csv",index_col=0,header='infer')
Pathway_df=pd.read_csv(pathway_filename,index_col=0,sep='\t')
Genes=pd.read_excel("/data/nandas/Combined_coexp_TFplusMetabolic/Pathway_centric/PATHWAYS AND CATEGORIES APRIL 17 2020.xlsx",
                    sheet_name="Gene2Pathway",index_col=0,header='infer')


# In[87]:


master_gene_list = [];
for gene in Genes.index:
    pathways = Genes.loc[gene]['Categories/Pathways'].split(';')
    master_gene_list += len(pathways) * [gene];


# In[88]:


print(master_gene_list)


# In[89]:


Pathway_df


# In[90]:


number_of_genes = len(master_gene_list)
new_pathway_df = Pathway_df;
random_sample = random.sample(range(number_of_genes), number_of_genes);
for pathway in new_pathway_df.index:
    pathway_list = new_pathway_df.loc[pathway];
    pathway_list_dropNA = pathway_list.dropna();
    sample_subset = random_sample[:len(pathway_list_dropNA)];
    random_sample = random_sample[len(pathway_list_dropNA):]
    sampled_gene_list = [master_gene_list[x] for x in sample_subset]
    pathway_list.iloc[1:len(sampled_gene_list)+1] = sampled_gene_list;
    new_pathway_df.loc[pathway] = pathway_list
   

    


# ### Setting default coregulated state of pathway

# In[91]:


Pathway_df['IsRegulated'] = False


# In[92]:


Pathway_df['IsRegulated']


# ## PreRank Gene set enrichment analyses for custom pathway annotations

# In[93]:


def generate_random_pathway_df_V1(dataframe, master_gene_list):
    number_of_genes = len(master_gene_list)
    random_sample = random.sample(range(number_of_genes), number_of_genes);
    new_pathway_df = dataframe;
    for pathway in new_pathway_df.index:
        pathway_list = new_pathway_df.loc[pathway];
        pathway_list_dropNA = pathway_list.dropna();
        sample_subset = random_sample[:len(pathway_list_dropNA)];
        
        random_sample = random_sample[len(pathway_list_dropNA):]
        sampled_gene_list = [master_gene_list[x] for x in sample_subset]
    
        pathway_list.iloc[1:len(sampled_gene_list)+1] = sampled_gene_list;
        new_pathway_df.loc[pathway] = pathway_list
        pathway=new_pathway_df
        new_pathway_df.to_csv("Pathway_gene_{}.gmt".format(i),sep='\t')
        pathway_filename="Pathway_gene_{}.gmt".format(i)
    return new_pathway_df,pathway_filename


def generate_random_pathway_df(dataframe):
    new_pathway_df = dataframe;
    for pathway in new_pathway_df.index:
        pathway_list = new_pathway_df.loc[pathway];

        pathway_list_dropNA = pathway_list.dropna();
        random_sample = random.sample(range(number_of_genes), number_of_genes);
        sampled_gene_list = list(Genes.iloc[random_sample[:len(pathway_list_dropNA)]].index)
    #     print(sampled_gene_list)
        pathway_list.iloc[1:len(sampled_gene_list)+1] = sampled_gene_list;
        new_pathway_df.loc[pathway] = pathway_list
        pathway=new_pathway_df
        new_pathway_df.to_csv("Pathway_gene_{}.gmt".format(i),sep='\t')
        pathway_filename="Pathway_gene_{}.gmt".format(i)
    return new_pathway_df,pathway_filename

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
    #Combined=shuffle(Combined)
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
    fig = plt.figure(figsize=(11,12))
    df = pd.DataFrame({'Enrichment Score': Sorted_values.es,
                   'p-value': Sorted_values.pval,'FDR':Sorted_values.fdr}, index=Sorted_values.index)
    ax = df.plot.barh(rot=0)
    plt.legend(loc='best', bbox_to_anchor=(1, 1))
    plt.savefig("{}/{}_plot.png".format(outdir, pathway))
    plt.close()
    
def PlotGSEA(pre_res, pathway, outdir):
    terms = pre_res.res2d.sort_values(by=['es'],ascending=False).index
    fig=gseaplot(rank_metric=pre_res.ranking, term=pathway, **pre_res.results[pathway])
    plt.savefig("'{}/{}_gsea.png'.format(outdir,pathway)")
    
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
    
    


# In[94]:


Pathway_df_withoutIsRegulated = Pathway_df.drop(['IsRegulated'], axis=1);


# In[95]:


New_df = pd.DataFrame([])
Status_df=Pathway_df['IsRegulated']
Status_df=pd.DataFrame(Status_df)
for i in range(0,100):
    print(i)
    New_df,pathway_filename=generate_random_pathway_df_V1(Pathway_df, master_gene_list)
    #print(New_df)
    print(pathway_filename)
    for pathway in New_df.index:
#     pathway = 'PROPIONATE_SHUNT';
        print(pathway)
    #     pathway = 'ALA_ASP_AND_GLU_METABOLISM';
        genes = list(Pathway_df_withoutIsRegulated.loc[pathway].dropna());
        pre_res = PreRank(genes, pathway);
        if(pre_res is None):
            continue; 
        Status_df.at[pathway, i] = _is_regulated_pathway_(pre_res, pathway);
        print("{} is regulated:{}".format(pathway,_is_regulated_pathway_(pre_res, pathway)))
#     PlotEnrichment(pre_res, pathway, outdir=pathway)
#     if(pathway in pre_res.res2d.index):
#         PlotGSEA(pre_res, pathway,pathway)
#         gsea_result_df=pre_res.res2d.loc[pathway];
#         New_df=New_df.append(gsea_result_df)
#     break;
Status_df.to_csv("Random10_Pathway_Regulation_status.csv")
# New_df.to_csv("Final_pathway_gsea.csv")




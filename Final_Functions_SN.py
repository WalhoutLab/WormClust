#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Importing external modules
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
import


# In[3]:


# Converts a square value matrix to pair-value combinations 
#(example:converts an n*n gene correlation matrix to 
#3 column table with Gene A-Gene B- value respectively )
def StackedDF_SN(df,Index1,Index2,ValueName):
    x=df.unstack()
    x.index.rename([Index1,Index2], inplace=True)
    x = x.to_frame(ValueName).reset_index()
    x.dropna(inplace=True)
    x = x[(x.iloc[:,0])!=(x.iloc[:,1])]
    x['check_string'] = x.apply(lambda row: ''.join(sorted([row[Index1], row[Index2]])), axis=1)
    x.drop_duplicates('check_string',keep='first',inplace=True)
    x.drop(columns=['check_string'],inplace=True)
    return x


# In[4]:


def wb_to_gene_SN(matrix):
    mapper_df=pd.read_csv("/data/nandas/MEFIT/predicted/mapper_final.csv", header='infer',index_col=0)
    wb_to_gene = {};
    for wb in mapper_df.index:
        wb_to_gene[wb] = str(mapper_df.loc[wb]['GeneID'])
    matrix=matrix.rename(index=wb_to_gene,columns=wb_to_gene)
    return matrix

def PreRank_SN(gene, outdir,metabolic_corr_df):
    Combined=metabolic_corr_df[gene];
    Combined.dropna(inplace=True)
#     print(Mean)
    rnk=Combined.sort_values(ascending=False)
    #print("Rank: {}".format(rnk))    
    pre_res = gp.prerank(rnk=rnk, gene_sets=pathway_filename4, processes=4,min_size=2, outdir=outdir, format='png', weighted_score_type=1,verbose=True)
    return pre_res

def _is_regulated_pathway_SN(pre_res, gene):
#     print('Hello There: {}'.format(pre_res));
#     print('Shivani Here: {}'.format(pre_res.res2d))
    pathway_pre_res = pre_res.res2d.sort_values(by=['fdr'],ascending=True);
#     is_regulated_pathway = pathway_pre_res.es >= 0.70 and pathway_pre_res.fdr <= 0.05
    is_regulated_pathway=False;
    for index in pathway_pre_res.index:
        val = pathway_pre_res.loc[index];
#         print(val)
        if(val.fdr <= 0.05 and val.nes > 0 and val.es >0 and val.nes!=np.inf):
            is_regulated_pathway = True;
#             print("Changed IsRegulatedValue to True");
#             print("FDR: {}---NES: {}---: ES:{}".format(val.fdr, val.nes, val.es))
            break;
    print("OutsideLoop: {}".format(is_regulated_pathway));
    return is_regulated_pathway;

def PlotEnrichment_SN(pre_res,gene, outdir):
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
    
def PlotGSEA_SN(pre_res, gene, outdir):
    terms = pre_res.res2d.sort_values(by=['fdr'],ascending=True)
    term_index=terms.index
    if terms.nes[0]>0:
        fig=gseaplot(rank_metric=pre_res.ranking, term=term_index[0], **pre_res.results[term_index[0]])
#     ax2=gseaplot(rank_metric=pre_res.ranking, term=term_index[0], **pre_res.results[term_index[1]])
#     ax3=gseaplot(rank_metric=pre_res.ranking, term=term_index[0], **pre_res.results[term_index[2]])
#     ax4=gseaplot(rank_metric=pre_res.ranking, term=term_index[0], **pre_res.results[term_index[3]])
    
    else:
        fig=gseaplot(rank_metric=pre_res.ranking, term=term_index[1], **pre_res.results[term_index[1]])
    plt.savefig('{}/{}_gsea.png'.format(outdir,gene))
    
    
def ConvertPairsToMatrix_SN(bayesian_metabol_df):
    bayesian_metabol_df.set_axis(['Gene1','Gene2','weight'], axis=1,inplace=True)
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

def SeqToGene_SN(output_df):
    All_metabolic_genes_df=pd.read_excel("/data/nandas/Transcription/All_Worm_Metabolism_Genes.xlsx")
    All_metabolic_genes_df_1=All_metabolic_genes_df[['Gene Name','Sequence Name']]
    All_metabolic_genes_df_1.set_index('Sequence Name', inplace=True)
    temp_dict = {};
    for seq in All_metabolic_genes_df_1.index:
        temp_dict[seq] = str(All_metabolic_genes_df_1.loc[seq]['Gene Name']);
    output_df.rename(columns=temp_dict,index=temp_dict,inplace=True)
    return output_df
    


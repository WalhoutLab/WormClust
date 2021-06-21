#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


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

# Converts Wormbase IDs to gene IDs
def wb_to_gene_SN(matrix):
    mapper_df=pd.read_csv("/data/nandas/MEFIT/predicted/mapper_final.csv", header='infer',index_col=0)
    wb_to_gene = {};
    for wb in mapper_df.index:
        wb_to_gene[wb] = str(mapper_df.loc[wb]['GeneID'])
    matrix=matrix.rename(index=wb_to_gene,columns=wb_to_gene)
    return matrix

# Converts Gene IDs to Wormbase IDs
def gene_to_wb_SN(matrix):
    mapper_df=pd.read_csv("/data/nandas/MEFIT/predicted/mapper_final.csv", header='infer',index_col=1)
    gene_to_wb = {};
    for gene in mapper_df.index:
        gene_to_wb[gene] = str(mapper_df.loc[gene]['Wormbase ID'])
    matrix=matrix.rename(index=gene_to_wb,columns=gene_to_wb)
    return matrix

# Using the preranked file for gene set enrichment analysis (GSEA) [Individual genes]
def PreRank_SN(gene, outdir,metabolic_corr_df):
    Combined=metabolic_corr_df[gene];
    Combined.dropna(inplace=True)
#     print(Mean)
    rnk=Combined.sort_values(ascending=False)
    #print("Rank: {}".format(rnk))    
    pre_res = gp.prerank(rnk=rnk, gene_sets=pathway_filename4, processes=4,min_size=2, outdir=outdir, format='png', weighted_score_type=1,verbose=True)
    return pre_res

# Determine if the pathway is co-regulated based on pathway self-enrichment with FDR<=0.05
def _is_regulated_pathway_SN(pre_res, gene):
    pathway_pre_res = pre_res.res2d.sort_values(by=['fdr'],ascending=True);
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

# Plot the enrichment statistics
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

#Plot GSEA result for individual pathways
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
    
# Convert pairwise data to square matrix   
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

# Convert Sequence ID to Gene Symbol
def SeqToGene_SN(output_df):
    All_metabolic_genes_df=pd.read_excel("/data/nandas/Transcription/All_Worm_Metabolism_Genes.xlsx")
    All_metabolic_genes_df_1=All_metabolic_genes_df[['Gene Name','Sequence Name']]
    All_metabolic_genes_df_1.set_index('Sequence Name', inplace=True)
    temp_dict = {};
    for seq in All_metabolic_genes_df_1.index:
        temp_dict[seq] = str(All_metabolic_genes_df_1.loc[seq]['Gene Name']);
    output_df.rename(columns=temp_dict,index=temp_dict,inplace=True)
    return output_df

# Merge two dataframes together based on two columns
def MergeDFs_SN(df1,df2,par1,par2):
    result=pd.merge(df1, df2, on=[par1, par2])
    result = result[(result.iloc[:,0])!=(result.iloc[:,1])]
    result.reset_index(inplace=True)
#     result=result[~(result.Coflux.isna())]
    result['check_string'] = result.apply(lambda row: ''.join(sorted([row[par1], row[par2]])), axis=1)
    result.drop_duplicates('check_string',inplace=True)
    return result;


# In[3]:


def ExtractPathwayExpression(expressionfile,pathway,pathwayfile):
    genes = pd.DataFrame(pathwayfile.loc[pathway]).dropna();
    genes= genes[pathway]
    genes=list(genes)
    pathway_matrix=expressionfile.loc[genes]
    pathway_matrix=wb_to_gene_SN(pathway_matrix)
    Clustermap(pathway_matrix,pathway)
    return pathway_matrix,clusterpathway_matrix

def Clustermap(pathway_matrix,pathway):
    clusterpathway_matrix=pathway_matrix.replace(np.nan,0,inplace=True)
    sns.clustermap(pathway_matrix,figsize=(50, 25),method='average',cbar_kws={'label':'Mean_Correlation'},
               row_cluster=True,
                yticklabels=True,xticklabels=True,cmap="vlag")
    plt.savefig("ExpressionClusterMap_{}".format(pathway))
    return clusterpathway_matrix

# Extract gene pairs of specific pathways
def Filtering_Df_SN(par1,par2,dataFrame,PathwayFile,pathway):
    print("Inside Filtering_Df: {}".format(pathway));
    genes = pd.DataFrame(PathwayFile.loc[pathway]).dropna();
    genes= genes[pathway]
    genes=list(genes)
    print(genes)
    dataFrame2 = dataFrame[dataFrame[par1].isin(genes)];
    dataFrame2 = dataFrame2[dataFrame2[par2].isin(genes)];
    return dataFrame2

## Extract datasets of specific Pathways
def ExtractDatasets_SN(pathway):
    filename="PathwayGeneCorrs_{}.csv".format(pathway)
    print(filename)
    dataset=pd.read_csv(filename)
    dataset.set_index(['Gene1'],inplace=True)
    dataset=wb_to_gene(dataset)
    dataset.reset_index(inplace=True)
    dataset.set_index(['Gene2'],inplace=True)
    dataset=wb_to_gene(dataset)
    dataset.reset_index(inplace=True)
    dataset.set_index(['Gene1', 'Gene2'],inplace=True)
#     print(dataset)
    dataset.drop(columns=['Unnamed: 0'],inplace=True)
    #Replacing NA with zero
    dataset.replace(np.nan,0,inplace=True)
    sns.clustermap(dataset,figsize=(35, 20),method='average',cbar_kws={'label':'Coexpression'},
               row_cluster=True,
                yticklabels=True,xticklabels=True,cmap="vlag") 
    return dataset


# In[4]:


## Assigning coexpression to gene relations
def GeneRelation_to_Coexpression(RelationDf,CoexpressionDf):
    for i in RelationDf.index:
        print(i)
        val = RelationDf.iloc[i]
        #print(val.Gene1)
        #print(val.Gene2)
    #     if (i==100):
    #         break;


        combined_value_forward = CoexpressionDf[(CoexpressionDf.Gene1 == val.Gene1) & (CoexpressionDf.Gene2 == val.Gene2)]
        combined_value_backward = CoexpressionDf[(CoexpressionDf.Gene1 == val.Gene2) & (CoexpressionDf.Gene2 == val.Gene1)]
        if(combined_value_forward.empty == False):
            combined_value = combined_value_forward;
        elif (combined_value_backward.empty == False):
            combined_value = combined_value_backward;
        else:
            continue; 
        #print(combined_value)
        #print("Weight --- {}".format(float(combined_value.weight)));
        RelationDf.at[i,'Correlation'] = float(combined_value.weight);
    return RelationDf

def DropDuplicateColumns(x,Index1,Index2):
    x.dropna(inplace=True)
    x = x[(x.iloc[:,0])!=(x.iloc[:,1])]
    x['check_string'] = x.apply(lambda row: ''.join(sorted([row[Index1], row[Index2]])), axis=1)
    x.drop_duplicates('check_string',keep='first',inplace=True)
    x.drop(columns=['check_string'],inplace=True)
    return x

## Sort gene pairs in such a way that the first value in pair (index1) is smaller than second value in the pair(index2)
def SortGenePairs(Unsorteddf,index1,index2):
    Sorteddf=pd.DataFrame([]);
    for i in range(Unsorteddf.shape[0]):
        print (i)
        if (Unsorteddf.loc[i,index1])>(Unsorteddf.loc[i,index2]):
            Sorteddf.at[i,index1]=Unsorteddf.loc[i,index2];
            Sorteddf.at[i,index2]=Unsorteddf.loc[i,index1];
        else:
            Sorteddf.at[i,index1]=Unsorteddf.loc[i,index1];
            Sorteddf.at[i,index2]=Unsorteddf.loc[i,index2];            
    return Sorteddf 

def GeneRelationParser(GeneID,genes):
    fp = open('GeneRelations.csv', 'w');
    count = 0;
    for rxn in genes.index:
        count = count + 1;
        print(count);
        GENES = genes.loc[rxn][GeneID]; 
        #print(GENES)
        expression = ',';
        expression_name = 'OR';
        #print(geneExp)
        geneExp_list = GENES.split(expression)
        gene_number = len(geneExp_list);
        for ii in range(gene_number):
            for iii in range(ii+1, gene_number):
                line = '{},{},{},{}\n'.format(rxn, geneExp_list[ii].strip(), geneExp_list[iii].strip(), expression_name);
                #print(line)
                fp.write(line);

    fp.close();   


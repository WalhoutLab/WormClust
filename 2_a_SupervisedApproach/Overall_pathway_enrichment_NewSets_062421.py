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
import fnmatch


# In[2]:


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

# Converts Wormbase IDs to gene IDs
def wb_to_gene_SN(matrix):
    mapper_df=pd.read_csv("/data/nandas/mapper_062421.csv", header='infer',index_col=0)
    wb_to_gene = {};
    for wb in mapper_df.index:
        wb_to_gene[wb] = str(mapper_df.loc[wb]['GeneID'])
    matrix=matrix.rename(index=wb_to_gene,columns=wb_to_gene)
    return matrix


# In[3]:


get_ipython().system(u' pip install matplotlib.pyplot.gca')


# ## Setting base directory

# In[4]:


Base_dir='/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/'
output_dir='/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/OverallPathwayEnrichment062421/'
os.chdir(Base_dir)


# In[5]:


MetabolicPairs=pd.read_csv("/data/nandas/Combined_coexp/MetabolicCorrPairs_062321.dat",
                           sep='\t',header=None)


# In[6]:


# !mkdir /data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/OverallPathwayEnrichment062421/


# In[7]:


# metabolic_corr_df=ConvertPairsToMatrix_SN(MetabolicPairs)


# ## Reading required files: GeneSets(gmt), PathwayToGenes and Gene Correlations 

# In[8]:


pathway_filename1 = '/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/Genesets_NAME_090320_LEVEL_1.gmt';
pathway_filename2 = '/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/Genesets_NAME_090320_LEVEL_2.gmt';
pathway_filename3 = '/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/Genesets_NAME_090320_LEVEL_3.gmt';
pathway_filename4 = '/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/Genesets_NAME_090320_LEVEL_4.gmt';
# metabolic_corr_df=pd.read_csv("/data/nandas/Combined_coexp/Sleipnir/Final_data_080620/UMN/MetabolicCorrMatrix_083120.csv",index_col=0,header='infer')
Pathway_df1=pd.read_csv(pathway_filename1,index_col=0,sep='\t')
Pathway_df2=pd.read_csv(pathway_filename2,index_col=0,sep='\t')
Pathway_df3=pd.read_csv(pathway_filename3,index_col=0,sep='\t')
Pathway_df4=pd.read_csv(pathway_filename4,index_col=0,sep='\t')


# In[9]:


# metabolic_corr_df.to_csv("MetabolicCorrMatrix_021921.csv")


# In[10]:


metabolic_corr_df=pd.read_csv("/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/OrphanGenes/MetabolicCorrMatrix062321.csv",
                              index_col=0)


# In[11]:


## To ignore self-correlations, remove diagonals of correlation matrix
np.fill_diagonal(metabolic_corr_df.values,np.nan)


# In[12]:


metabolic_corr_df.min().min()


# In[13]:


metabolic_corr_df.shape


# In[14]:


metabolic_corr_df=wb_to_gene_SN(metabolic_corr_df)


# In[15]:


metabolic_corr_df.to_csv("MetabolicCorrMatrix_GeneID_062321.csv")


# In[16]:


metabolic_corr_df=(metabolic_corr_df+1)/2


# In[17]:


metabolic_corr_df.min().min()


# In[18]:


metabolic_corr_df['haly-1'].sort_values(ascending=False)[0:20]


# ### Setting default coregulated state of pathway

# In[19]:


Pathway_df1['IsRegulated'] = False

Pathway_df2['IsRegulated'] = False
Pathway_df3['IsRegulated'] = False
Pathway_df4['IsRegulated'] = False


# In[20]:


metabolic_corr_df=metabolic_corr_df[~metabolic_corr_df.index.duplicated(keep='first')]


# ## PreRank Gene set enrichment analyses for custom pathway annotations

# In[21]:


def wb_to_gene_SN(matrix):
    mapper_df=pd.read_csv("/data/nandas/mapper_062421.csv", header='infer',index_col=0)
    wb_to_gene = {};
    for wb in mapper_df.index:
        wb_to_gene[wb] = str(mapper_df.loc[wb]['GeneID'])
    matrix=matrix.rename(index=wb_to_gene,columns=wb_to_gene)
    return matrix

def gene_to_wb(matrix):
    mapper_df=pd.read_csv("/data/nandas/mapper_062421.csv", header='infer',index_col=1)
    gene_to_wb = {};
    for gene in mapper_df.index:
        gene_to_wb[gene] = str(mapper_df.loc[gene]['Wormbase ID'])
    matrix=matrix.rename(index=gene_to_wb,columns=gene_to_wb)
    return matrix

def SeqToGene(output_df):
    All_metabolic_genes_df=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_121120.csv")
    All_metabolic_genes_df_1=All_metabolic_genes_df[['Gene Name','Sequence Name']]
    All_metabolic_genes_df_1.set_index('Sequence Name', inplace=True)
    temp_dict = {};
    for seq in All_metabolic_genes_df_1.index:
        temp_dict[seq] = str(All_metabolic_genes_df_1.loc[seq]['Gene Name']);
    output_df.rename(columns=temp_dict,index=temp_dict,inplace=True)
    return output_df

def PreRank(genes, outdir,gene_sets):
#     print("Genes: {}".format(genes));
    print("Length of genes:{}".format(len(genes)))
    genes=pd.DataFrame(genes)
    genes.set_index([0],inplace=True)
    genes=gene_to_wb(genes)
    genes=wb_to_gene_SN(genes)
    genes=SeqToGene(genes)
    print(genes.index)
    intersection_list = list(set(metabolic_corr_df.index).intersection(set(genes.index)))
    missing_genes=list(set(genes).difference(set(intersection_list)))
#     print("IntersectionList: {}".format(intersection_list));
#     print("Length of intersection list:{}".format(len(intersection_list)))
    intersection_list
    print('Missing genes:{}\n{}'.format(len(missing_genes),missing_genes))
    if(len(missing_genes) == len(genes)):
        return;
    Combined=metabolic_corr_df[intersection_list];
    Mean=np.exp(np.log(Combined.prod(axis=1))/Combined.notna().sum(1))
#    print("Mean before scaling:{}".format(Mean))
    Mean=(Mean*2)-1
#    print("Mean after scaling to lie between -1 and +1:{}".format(Mean))
#     Mean.dropna(inplace=True)
    rnk=Mean.sort_values(ascending=False)
    plt.rcParams["font.family"] = "Arial"
#     print("Rank: {}".format(rnk))    
    pre_res = gp.prerank(rnk=rnk, gene_sets=gene_sets, processes=4,min_size=2, outdir=outdir, format='svg', weighted_score_type=1,verbose=True)
    return pre_res

def _is_regulated_pathway_(pre_res, pathway):
#     print('Hello There: {}'.format(pre_res));
#     print('Shivani Here: {}'.format(pre_res.res2d))
    if(pathway not in pre_res.res2d.index):
        return "NaN"
    else:
        pathway_pre_res = pre_res.res2d.loc[pathway];
#     is_regulated_pathway = pathway_pre_res.es >= 0.70 and pathway_pre_res.fdr <= 0.05
        is_regulated_pathway =  (pathway_pre_res.fdr <= 0.05) and (pathway_pre_res.nes>0) and (pathway_pre_res.nes!=np.inf) and (pathway_pre_res.es>0)
    return is_regulated_pathway;

def PlotEnrichment(pre_res,pathway, outdir):
    Sorted_values=pre_res.res2d.sort_values(ascending=False,by=['nes'])[0:20]
    fig = plt.figure(figsize=(8,15))
    df = pd.DataFrame({'Normalized Enrichment Score': Sorted_values.nes,
                   'p-value': Sorted_values.pval,'FDR':Sorted_values.fdr}, index=Sorted_values.index)
    ax = df.plot.barh(rot=0)
#     plt.gca.invert_yaxis()
#     plt.legend(loc='best', bbox_to_anchor=(1, 1))
    plt.savefig("{}/{}_plot.svg".format(outdir, pathway))
    plt.show()
    
def PlotGSEA(pre_res, pathway, outdir):
    terms = pre_res.res2d.sort_values(by=['fdr'],ascending=True).index
    fig=gseaplot(rank_metric=pre_res.ranking, term=pathway, **pre_res.results[pathway],ofname='{}/{}_gsea.svg'.format(outdir,pathway))
    plt.rcParams["font.family"] = "Arial"
    plt.savefig("{}/{}_gsea.svg".format(outdir,pathway),dpi=300)
    
    


# In[22]:


metabolic_corr_df=gene_to_wb(metabolic_corr_df)
metabolic_corr_df=wb_to_gene_SN(metabolic_corr_df)
metabolic_corr_df=SeqToGene(metabolic_corr_df)


# In[23]:


Pathway_df1_withoutIsRegulated = Pathway_df1.drop(['IsRegulated'], axis=1);
Pathway_df2_withoutIsRegulated = Pathway_df2.drop(['IsRegulated'], axis=1);
Pathway_df3_withoutIsRegulated = Pathway_df3.drop(['IsRegulated'], axis=1);
Pathway_df4_withoutIsRegulated = Pathway_df4.drop(['IsRegulated'], axis=1);


# In[37]:


# New_df1 = pd.DataFrame([])
# for pathway in Pathway_df1.index:
#     #pathway = 'GLY_CLEAVAGE_SYSTEM';
#     print(pathway)
# #     pathway = 'ALA_ASP_AND_GLU_METABOLISM';
#     genes = list(Pathway_df1_withoutIsRegulated.loc[pathway].dropna());
#     pre_res = PreRank(genes, pathway,gene_sets=pathway_filename1);
#     if(pre_res is None):
#         continue; 
#     Pathway_df1.at[pathway, 'IsRegulated'] = _is_regulated_pathway_(pre_res, pathway);
#     print("{} is regulated:{}".format(pathway,_is_regulated_pathway_(pre_res, pathway)))
#     PlotEnrichment(pre_res, pathway, outdir=pathway)
#     if(pathway in pre_res.res2d.index):
#         PlotGSEA(pre_res, pathway,pathway)
#         gsea_result_df=pre_res.res2d.loc[pathway];
#         New_df1=New_df1.append(gsea_result_df)
# # Pathway_df1.to_csv("Pathway_Regulation_status1.csv")
# New_df1.to_csv("Final_pathway1_gsea.csv")




# In[38]:


# New_df2 = pd.DataFrame([])
# for pathway in Pathway_df2.index:
#     #pathway = 'GLY_CLEAVAGE_SYSTEM';
#     print(pathway)
# #     pathway = 'ALA_ASP_AND_GLU_METABOLISM';
#     genes = list(Pathway_df2_withoutIsRegulated.loc[pathway].dropna());
#     pre_res = PreRank(genes, pathway,gene_sets=pathway_filename2);
#     if(pre_res is None):
#         continue; 
#     Pathway_df2.at[pathway, 'IsRegulated'] = _is_regulated_pathway_(pre_res, pathway);
#     print("{} is regulated:{}".format(pathway,_is_regulated_pathway_(pre_res, pathway)))
#     PlotEnrichment(pre_res, pathway, outdir=pathway)
#     if(pathway in pre_res.res2d.index):
#         PlotGSEA(pre_res, pathway,pathway)
#         gsea_result_df=pre_res.res2d.loc[pathway];
#         New_df2=New_df2.append(gsea_result_df)
# # Pathway_df2.to_csv("Pathway_Regulation_status2.csv")
# New_df2.to_csv("Final_pathway2_gsea.csv")


# In[39]:


# New_df3 = pd.DataFrame([])
# for pathway in Pathway_df3.index:
#     #pathway = 'GLY_CLEAVAGE_SYSTEM';
#     print(pathway)
# #     pathway = 'ALA_ASP_AND_GLU_METABOLISM';
#     genes = list(Pathway_df3_withoutIsRegulated.loc[pathway].dropna());
#     pre_res = PreRank(genes, pathway,gene_sets=pathway_filename3);
#     if(pre_res is None):
#         continue; 
#     Pathway_df3.at[pathway, 'IsRegulated'] = _is_regulated_pathway_(pre_res, pathway);
#     print("{} is regulated:{}".format(pathway,_is_regulated_pathway_(pre_res, pathway)))
#     PlotEnrichment(pre_res, pathway, outdir=pathway)
#     if(pathway in pre_res.res2d.index):
#         PlotGSEA(pre_res, pathway,pathway)
#         gsea_result_df=pre_res.res2d.loc[pathway];
#         New_df3=New_df3.append(gsea_result_df)
# # Pathway_df3.to_csv("Pathway_Regulation_status3.csv")
# New_df3.to_csv("Final_pathway3_gsea.csv")


# In[ ]:


New_df4 = pd.DataFrame([])
Final_gsea=pd.DataFrame([])
# count=1
for pathway in Pathway_df4.index:
#     pathway = 'ALA';
#     count=count+1
    print(pathway)
#     pathway = 'ALA_ASP_AND_GLU_METABOLISM';
    genes = list(Pathway_df4_withoutIsRegulated.loc[pathway].dropna());
    file_exist = False;
    for file in os.listdir(output_dir):
        if fnmatch.fnmatch(file, "Pathway_{}.csv".format(pathway)):
            print("File: {} found, skipping!!!".format(file))
            file_exist = True;
    if(not file_exist):
        pre_res = PreRank(genes=genes, outdir="{}/{}".format(output_dir,pathway),gene_sets=pathway_filename4);
        if(pre_res is None):
            continue; 
        Pathway_df4.at[pathway, 'IsRegulated'] = _is_regulated_pathway_(pre_res, pathway);
        print("{} is regulated:{}".format(pathway,_is_regulated_pathway_(pre_res, pathway)))
        PlotEnrichment(pre_res, pathway, outdir=pathway)
        if(pathway in pre_res.res2d.index):
            PlotGSEA(pre_res, pathway,pathway)
            gsea_result_df=pre_res.res2d.loc[pathway];
            gsea_result_df=pd.DataFrame(gsea_result_df)
            gsea_result_df.to_csv("{}/Pathway_self_enrichment_{}.csv".format(output_dir,pathway))
            gsea_result_df=gsea_result_df.transpose()
#             gsea_result_df['Pathway_Main']=pathway
            New_df4=New_df4.append(gsea_result_df)
#             New_df4.at[pathway,'Pathway_main']=pathway
#             New_df4.to_csv("{}/Pathway_{}.csv".format(output_dir,pathway))
#         if count>2:
#             break;
Pathway_df4.to_csv("Pathway_Regulation_status4.csv")
New_df4.to_csv("{}/Final_pathway4_gsea.csv".format(output_dir))


# In[123]:


metabolic_corr_df.shape


# In[ ]:


Pathway_df4=pd.read_csv("Pathway_Regulation_status4.csv",index_col=0)


# In[ ]:


New_df4=pd.read_csv("{}/Final_pathway4_gsea.csv".format(output_dir),index_col=0)


# In[ ]:


New_df4.index=New_df4.index.str.replace("_"," ")


# In[ ]:


New_df4.index=New_df4.index.str.title()


# In[ ]:


New_df4.to_csv("Final_pathway4_gsea.csv".format(output_dir))


# In[ ]:


New_df4.loc['His']


# In[80]:


Pathway_df4=Pathway_df4['IsRegulated']


# In[81]:


Pathway_df4=pd.DataFrame(Pathway_df4)


# In[82]:


Regulated=Pathway_df4[Pathway_df4['IsRegulated']==True]


# In[83]:


NotRegulated=Pathway_df4[Pathway_df4['IsRegulated']!=True]


# In[107]:


Enrichment=pd.read_csv("{}/Final_pathway4_gsea.csv".format(output_dir),index_col=0)


# In[84]:


Regulated.shape


# In[108]:


Enrichment_result=New_df4[['es','nes','fdr']]


# In[109]:


Enrichment_result=Enrichment_result[Enrichment_result.fdr<=0.05]


# In[110]:


Enrichment_result.sort_values(by=['fdr'],inplace=True,ascending=False)


# In[112]:


Enrichment_result.rename(columns={'es':'Enrichment Score','nes':'Normalized Enrichment Score','fdr':'FDR'},inplace=True)


# In[118]:


ax=plt.figure(figsize=(40,180))
plt.rcParams["font.family"] = "Arial"
# Enrichment_result=pd.DataFrame({'es':'Enrichment Score','nes':'Normalized Enrichment Score','fdr':'FDR'},index=Enrichment_result.index)
Enrichment_result.plot.barh()
plt.savefig("EnrichmentPlot.svg",dpi=300,format='svg')
plt.show()


# In[117]:


sns.barplot(Enrichment_result)


# In[32]:


def Pathway_Enrichment_Of_Regulated_Genes(Base_dir,Pathway_df4,level):
    path =Base_dir
    print(path)
    dfs = []
    for x in Pathway_df4.index:
        files='{}/Pathway_{}.csv'.format(path,x)
        df=pd.read_csv(files)
        df['Pathway_main']=x;
#         df['Class']=RegulatedMetabolic.loc[x]['Class']
        df.set_index(['Pathway_main'],inplace=True)
        df.to_csv(files)
        print(files)
        dfs.append(pd.read_csv(files))

    # # Concatenate all data into one DataFrame

    All_TFs = pd.concat(dfs, ignore_index=False)
    # # Filtering out negative or not applicable enrichment¶
    All_TFs=All_TFs[All_TFs.nes!=np.inf]
    All_TFs=All_TFs[All_TFs.nes>=0]
    All_TFs.to_csv("Combined_TF_Pathway_Enrichement_{}.csv".format(level))
    return dfs,All_TFs


# In[33]:


df4,All_TFs4=Pathway_Enrichment_Of_Regulated_Genes(Base_dir=".",Pathway_df4=Pathway_df4,level='Level_4')


# In[ ]:


New_df4=pd.read_csv("Combined_TF_Pathway_Enrichement_Level_4.csv",index_col=0)


# In[ ]:


# New_df4.reset_index(inplace=True)
Significant_Pathways4=New_df4[New_df4.fdr<=0.05]
Significant_Pathways4=Significant_Pathways4[Significant_Pathways4.nes>2]
Significant_Pathways4=Significant_Pathways4[['Pathway_main','Term','fdr']]


# In[ ]:


# Significant_Pathways4['fdr']=1-(Significant_Pathways4.fdr)


# In[ ]:


Significant_Pathways4.sort_values(by=['fdr'],ascending=True)


# In[ ]:


def ConvertPairsToMatrix(bayesian_metabol_df):
    a = np.unique(bayesian_metabol_df['Pathway_main'])
    b = np.unique(bayesian_metabol_df['Term'])
#     c = np.union1d(a,b);
    data = np.zeros((len(a), len(b)));
    output_df = pd.DataFrame(data, index=a, columns=b)
    for values in bayesian_metabol_df.values: 
        output_df.loc[values[0]][values[1]] = values[2];
#         output_df[values[1]][values[0]]=values[2];
#     np.fill_diagonal(output_df.values,1)
    return output_df


# In[ ]:


Significant_PathwaysMatrix4=ConvertPairsToMatrix(Significant_Pathways4)
Significant_PathwaysMatrix4=Significant_PathwaysMatrix4.transpose()


# In[ ]:


len(np.unique(Significant_Pathways4['Pathway_main']))


# In[ ]:


Significant_PathwaysMatrix4=Significant_PathwaysMatrix4.loc[Significant_PathwaysMatrix4.sum(axis=1)!=0]
#     Significant_PathwaysMatrix2=Significant_PathwaysMatrix2.loc[(Significant_PathwaysMatrix2.sum(axis=1)!=0).index]
# #     Significant_PathwaysMatrix2=Significant_PathwaysMatrix2.loc[Significant_PathwaysMatrix2.index.str.contains('nhr')==True]
# #     Significant_PathwaysMatrix2=Significant_PathwaysMatrix2[(Significant_PathwaysMatrix2.sum()!=0).index]
Significant_PathwaysMatrix4=Significant_PathwaysMatrix4.transpose()
# #     print(Significant_PathwaysMatrix2.shape)
#     Significant_PathwaysMatrix2=Significant_PathwaysMatrix2.loc[(Significant_PathwaysMatrix2.sum()!=0).index]
# #     Significant_PathwaysMatrix2.drop(index=Categories,inplace=True)
# #     print(Significant_PathwaysMatrix2.shape)
# # # Pathways1.set_index(['Gene'],inplace=True)
sns.clustermap(Significant_PathwaysMatrix4,figsize=(28, 28),method='average',cbar_kws={'label':'FDR'},col_cluster=True,
                  yticklabels=True,xticklabels=True)
plt.savefig("PathwayCluster_Leve4.png")


# In[ ]:


Significant_Pathways4


# In[ ]:


Significant_PathwaysMatrix4=0.05-Significant_PathwaysMatrix4


# In[ ]:


Significant_PathwaysMatrix4


# In[ ]:





# In[ ]:


status1=Pathway_df1['IsRegulated']
status1=pd.DataFrame(status1)
status2=Pathway_df1['IsRegulated']
status2=pd.DataFrame(status2)
status3=Pathway_df3['IsRegulated']
status3=pd.DataFrame(status3)
status4=Pathway_df4['IsRegulated']
status4=pd.DataFrame(status4)


# In[ ]:


status1.to_csv("Pathway_Regulation_status1.csv")
status2.to_csv("Pathway_Regulation_status2.csv")
status3.to_csv("Pathway_Regulation_status3.csv")
status4.to_csv("Pathway_Regulation_status4.csv")


# In[ ]:


status1=pd.read_csv("Pathway_Regulation_status1.csv",index_col=0)
status2=pd.read_csv("Pathway_Regulation_status2.csv",index_col=0)
status3=pd.read_csv("Pathway_Regulation_status3.csv",index_col=0)
status4=pd.read_csv("Pathway_Regulation_status4.csv",index_col=0)


# In[ ]:


status1.IsRegulated.groupby(status1.IsRegulated).count()
status2.IsRegulated.groupby(status2.IsRegulated).count()
status3.IsRegulated.groupby(status3.IsRegulated).count()
status4.IsRegulated.groupby(status4.IsRegulated).count()


# In[ ]:


#Status_1
my_labels='Non-Regulated','Regulated'
sums1 = status1.IsRegulated.groupby(status1.IsRegulated).count()
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

#Status_2
my_labels='Non-Regulated','Regulated'
sums2 = status2.IsRegulated.groupby(status2.IsRegulated).count()
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

#Status_3
my_labels='Non-Regulated','Regulated'
sums3 = status3.IsRegulated.groupby(status3.IsRegulated).count()
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct



# In[ ]:


#Status_4
my_labels='Non-Regulated','Regulated'
sums4 = status4.IsRegulated.groupby(status4.IsRegulated).count()
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct


# In[ ]:


sums4


# In[ ]:


# Level_1
fig, ax = plt.subplots(figsize=(12,10))
size = 0.3
plt.pie(sums1,labels=my_labels,autopct=make_autopct(sums1))
plt.savefig("Overall_result_piechart_level_1.png")
plt.show()


# In[ ]:


# Level_2
fig, ax = plt.subplots(figsize=(12,10))
size = 0.3
plt.pie(sums2,labels=my_labels,autopct=make_autopct(sums2))
plt.savefig("Overall_result_piechart_level_2.png")


# In[ ]:


# Level 3
fig, ax = plt.subplots(figsize=(12,10))
size = 0.3
plt.pie(sums3,labels=my_labels,autopct=make_autopct(sums3))
plt.savefig("Overall_result_piechart_level_3.png")


# In[ ]:


# Level 4
fig, ax = plt.subplots(figsize=(12,10))
size = 0.3
plt.pie(sums4,labels=my_labels,autopct=make_autopct(sums4))
plt.savefig("Overall_result_piechart_level_4.png")


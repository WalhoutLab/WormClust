
# coding: utf-8

# In[1]:


## Importing modules
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
from gseapy.parser import Biomart
import os
import numpy as np
import seaborn as sns
from gseapy.plot import gseaplot


# ## Setting base directory

# In[2]:


Base_dir='/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/OverallPathwayEnrichmentmean033122'
os.chdir(Base_dir)


# In[3]:


# ! mkdir /data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/OverallPathwayEnrichmentmean033122


# ## Reading required files

# In[4]:


pathway_filename = '/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/Genesets_NAME_090320_LEVEL_4.gmt';
metabolic_corr_df=pd.read_csv("/data/nandas/Combined_coexp/Sleipnir/Final_data_080620/UMN/MetabolicCorrMatrix_083120.csv",index_col=0,header='infer')
Pathway_df=pd.read_csv(pathway_filename,index_col=0,sep='\t')


# In[5]:


# ## PreRank Gene set enrichment analyses for custom pathway annotations

# In[60]:
def wb_to_gene(matrix):
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", 
                          header='infer',index_col=1)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    wb_to_gene = {};
    for wb in mapper_df.index:
        wb_to_gene[wb] = str(mapper_df.loc[wb]['GeneName']);
    matrix=matrix.rename(index=wb_to_gene,columns=wb_to_gene)
    return matrix

def gene_to_wb(matrix):
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv",
                          header='infer',index_col=2)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    gene_to_wb = {};
    for gene in mapper_df.index:
        gene_to_wb[gene] = str(mapper_df.loc[gene]['WormBaseID']);
    matrix=matrix.rename(index=gene_to_wb,columns=gene_to_wb)
    return matrix

def SeqToWB(output_df):
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv",
                          header='infer',index_col=3)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    Seq_to_Wb = {};
    mapper_df=mapper_df[mapper_df.index!=np.nan]
    for seq in mapper_df.index:
        Seq_to_Wb[seq] = str(mapper_df.loc[seq]['WormBaseID']);
    matrix=matrix.rename(index=Seq_to_Wb,columns=Seq_to_Wb)
    return matrix

def SeqToGene(matrix):
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", 
                          header='infer',index_col=3)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    Seq_to_Gene = {};
    mapper_df=mapper_df[mapper_df.index!=np.nan]
    for seq in mapper_df.index:
        Seq_to_Gene[seq] = str(mapper_df.loc[seq]['GeneName']);
    matrix=matrix.rename(index=Seq_to_Gene,columns=Seq_to_Gene)
    return matrix

def GeneToSeq(matrix):
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", 
                          header='infer',
                          index_col=2)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    Gene_to_Seq = {};
    mapper_df=mapper_df[mapper_df.index!=np.nan]
    for gene in mapper_df.index:
        Gene_to_Seq[gene] = str(mapper_df.loc[gene]['SequenceID']);
    matrix=matrix.rename(index=Gene_to_Seq,columns=Gene_to_Seq)
    return matrix

def PreRank(genes, outdir,gene_sets):
#     print("Genes: {}".format(genes));
    print("Length of genes:{}".format(len(genes)))
    genes=pd.DataFrame(genes)
    genes.set_index([0],inplace=True)
    genes=SeqToGene(genes)
    genes=list(genes.index)
    intersection_list = list(set(metabolic_corr_df.index).intersection(set(genes)))
#     print("intersection_list:{}".format(intersection_list))
    missing_genes=list(set(genes).difference(set(intersection_list)))
#     print("IntersectionList: {}".format(intersection_list));
#     print("Length of intersection list:{}".format(len(intersection_list)))
#     print('Missing genes:{}\n{}'.format(len(missing_genes),missing_genes))
    if(len(missing_genes) == len(genes)):
        return;
    Combined=metabolic_corr_df[intersection_list];
    Mean=Combined.mean(axis=1)
#    print("Mean before scaling:{}".format(Mean))
#     Mean=(Mean*2)-1
#     print("Mean after scaling to lie between -1 and +1:{}".format(Mean))
    Mean.dropna(inplace=True)
    rnk=Mean.sort_values(ascending=False)
    plt.rcParams["font.family"] = "Arial"
#     print("Rank: {}".format(rnk))    
    pre_res = gp.prerank(rnk=rnk, gene_sets=gene_sets, processes=4,min_size=2, outdir=outdir, format='svg', 
                         weighted_score_type=1,verbose=True)
    plt.close()
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
    Sorted_values=pre_res.res2d.sort_values(ascending=False,by=['nes'])[0:40]
    fig = plt.figure(figsize=(8,15))
    df = pd.DataFrame({'Enrichment Score': Sorted_values.es,
                   'p-value': Sorted_values.pval,'FDR':Sorted_values.fdr}, index=Sorted_values.index)
    ax = df.plot.barh(rot=0)
    plt.legend(loc='best', bbox_to_anchor=(1, 1))
    plt.rcParams["font.family"] = "Arial"
    plt.savefig("{}/{}_plot.svg".format(outdir, pathway))
    plt.show()
    plt.close()
    
def PlotGSEA(pre_res, pathway, outdir,term):
    terms = pre_res.res2d.sort_values(by=['es'],ascending=False).index
#     print(terms[17])
    print("term is: {}".format(term))
    fig=gseaplot(rank_metric=pre_res.ranking,term=term, **pre_res.results[term],ofname='{}/{}_gsea.png'.format(outdir,term))
    plt.show()
#     plt.close()
    


# In[6]:


metabolic_corr_df=SeqToGene(metabolic_corr_df)


# In[7]:


metabolic_corr_df=metabolic_corr_df[~metabolic_corr_df.index.duplicated(keep='first')]


# In[8]:


# metabolic_corr_df=SeqToGene(metabolic_corr_df)


# In[9]:


metabolic_corr_df=wb_to_gene(metabolic_corr_df)


# In[10]:


# metabolic_corr_df=(metabolic_corr_df+1)/2


# In[11]:


metabolic_corr_df.min().min()


# In[12]:


# ### Setting default coregulated state of pathway
Pathway_df['IsRegulated'] = False
Pathway_df


# In[13]:


Pathway_df[Pathway_df.index.str.startswith("P")]


# In[14]:


np.fill_diagonal(metabolic_corr_df.values,np.nan)


# In[15]:


metabolic_corr_df=(metabolic_corr_df*2)-1


# In[16]:


metabolic_corr_df.min().min()


# In[17]:


Pathway_df[0:40]


# In[18]:


# Pathway_df_withoutIsRegulated = Pathway_df.drop(['IsRegulated'], axis=1);
# New_df = pd.DataFrame([])
# # i=0
# # j=0
# for pathway in Pathway_df.index:
# #     for i in range(0,6):
# #         print(i)
# #         for j in range (0,11):
# #     print(j)
# #     pathway = 'PANTOTHENATE_AND_COA_BIOSYNTHESIS';
#     print(pathway)
# #     pathway = 'OTHER';
#     genes = list(Pathway_df_withoutIsRegulated.loc[pathway].dropna());
#     pre_res = PreRank(genes, pathway,gene_sets=pathway_filename);
#     if(pre_res is None):
#         continue; 
#     Pathway_df.at[pathway, 'IsRegulated'] = _is_regulated_pathway_(pre_res, pathway);
#     print("{} is regulated:{}".format(pathway,_is_regulated_pathway_(pre_res, pathway)))
#     PlotEnrichment(pre_res, pathway, outdir=pathway)
#     if(pathway in pre_res.res2d.index):
# #                 fig, axes = plt.subplots(nrows= 5, ncols=10)
#         PlotGSEA(pre_res, pathway,pathway,term=pathway)
#         plt.show()
#         gsea_result_df=pre_res.res2d.loc[pathway];
#         New_df=New_df.append(gsea_result_df)
#     #     break;
# # Pathway_df.to_csv("Pathway_Regulation_status_112821.csv")
# New_df.to_csv("Final_pathway_gsea_033122.csv")


# In[19]:


Pathway_df.index


# In[21]:


Pathway_df_withoutIsRegulated = Pathway_df.drop(['IsRegulated'], axis=1);
New_df = pd.DataFrame([])
for pathway in Pathway_df.index:
    pathway = 'HIS';
    print(pathway)
#     pathway = 'OTHER';
    genes = list(Pathway_df_withoutIsRegulated.loc[pathway].dropna());
    pre_res = PreRank(genes, pathway,gene_sets=pathway_filename);
    if(pre_res is None):
        continue; 
    Pathway_df.at[pathway, 'IsRegulated'] = _is_regulated_pathway_(pre_res, pathway);
    print("{} is regulated:{}".format(pathway,_is_regulated_pathway_(pre_res, pathway)))
    PlotEnrichment(pre_res, pathway, outdir=pathway)
    if(pathway in pre_res.res2d.index):
        PlotGSEA(pre_res=pre_res,outdir=".",pathway=pathway,term=pathway)
        plt.show()
        plt.close()
        gsea_result_df=pre_res.res2d.loc[pathway];
        gsea_result_df['Pathway']=pathway
        gsea_result_df.to_csv("Pathway_self_enrichment_{}.csv".format(pathway))
        print(gsea_result_df)
        New_df=New_df.append(gsea_result_df)
    break;
# Pathway_df.to_csv("Pathway_Regulation_status_112821.csv")
# New_df.to_csv("PropionateShunt_Enrichment.csv")


# In[27]:


New_df.loc['NICOTINATE_AND_NICOTINAMIDE_METABOLISM']


# In[ ]:


New_df[0:40]


# In[ ]:


# New_df=New_df[New_df.es>0]


# In[ ]:


# New_df=New_df[New_df.fdr<=0.05]


# In[ ]:


# New_df.to_csv("Ketone_body_metabolism.csv")


# In[ ]:


New_df.to_csv("GSEAResult_052722.csv")


# In[ ]:


# Pathway_df.loc['PROPIONATE_SHUNT']


# In[ ]:


# Pathway_df=pd.read_csv("Pathway_Regulation_status.csv")


# In[24]:


New_df=pd.read_csv("GSEAResult_052722.csv",index_col=0)


# In[ ]:


FDR=New_df[New_df.fdr>0.05]


# In[ ]:


FDR.sort_values(by=['fdr'])


# In[ ]:


FDR.set_index(['Unnamed: 0'],inplace=True)


# In[ ]:


FDR


# In[ ]:


# FDR.drop(index=['CHITIN_BREAKDOWN','UGT_ENZYME'],inplace=True)


# In[ ]:


FDR.to_csv("OverallPathwayEnrichment_032822.csv")


# In[ ]:


FDR.shape


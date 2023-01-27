#!/usr/bin/env python
# coding: utf-8

# ### Importing modules

# In[ ]:





# In[94]:


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
  
import numpy as np
import pandas as pd
import re
#Import libraries
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
get_ipython().magic(u'matplotlib inline')
from pylab import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
Base_dir='/data/nandas/Transcription/KimDevTime_071620'
os.chdir(Base_dir)
from scipy.signal import find_peaks
# from CatExp import *
import statsmodels as st
import statsmodels.api as sm
import seaborn as sns


# In[ ]:





# In[2]:


mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',
                      index_col=0)
mapper_df=mapper_df.loc[mapper_df.index.dropna()]


# In[3]:


mapper_df


# ## Functions

# In[4]:


def wb_to_gene(matrix):
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=1)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    wb_to_gene = {};
    for wb in mapper_df.index:
        wb_to_gene[wb] = str(mapper_df.loc[wb]['GeneName']);
    matrix=matrix.rename(index=wb_to_gene,columns=wb_to_gene)
    return matrix

def gene_to_wb(matrix):
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=2)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    gene_to_wb = {};
    for gene in mapper_df.index:
        gene_to_wb[gene] = str(mapper_df.loc[gene]['WormBaseID']);
    matrix=matrix.rename(index=gene_to_wb,columns=gene_to_wb)
    return matrix

def SeqToWB(matrix):
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=3)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    Seq_to_Wb = {};
    mapper_df=mapper_df[mapper_df.index!=np.nan]
    for seq in mapper_df.index:
        Seq_to_Wb[seq] = str(mapper_df.loc[seq]['WormBaseID']);
    matrix=matrix.rename(index=Seq_to_Wb,columns=Seq_to_Wb)
    return matrix

def SeqToGene(matrix):
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=3)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    Seq_to_Gene = {};
    mapper_df=mapper_df[mapper_df.index!=np.nan]
    for seq in mapper_df.index:
        Seq_to_Gene[seq] = str(mapper_df.loc[seq]['GeneName']);
    matrix=matrix.rename(index=Seq_to_Gene,columns=Seq_to_Gene)
    return matrix
    


# In[5]:


Allgenes_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv")


# In[6]:


Allgenes_df.columns


# ## Reading regulated metabolic and non-metabolic genes

# In[7]:


Kim_dir='/data/nandas/Transcription/KimDevTime_071620/'
Tissue_dir='/data/nandas/Combined_coexp/Part_1_TranscriptionallyRegulatedGenes/Tissue/'
AllDatasets_dir='/data/nandas/Transcription/AllDatasets_061922/'


# ## Combining VS, CV and fraction of CV from development, tissue and all datasets

# In[8]:


iCEL1314_Kim=pd.read_csv(Kim_dir+"iCEL1314_KimDevTime_Bin.csv",index_col=0)


# In[9]:


iCEL1314_Tissue=pd.read_csv(Tissue_dir+"CVTissue_iCEL1314.csv",index_col=0)


# In[10]:


iCEL1314_Compendium=pd.read_csv(AllDatasets_dir+"CVCompendium_iCEL1314genes.csv",index_col=0)


# In[11]:


MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",index_col=0)
ClassA=MetabolicClasses[MetabolicClasses.Class=='A']
# Combined_iCEL1314['Class']=MetabolicClasses['Class']


# In[12]:


ClassA.Class.replace("A","iCEL1314 gene",inplace=True)


# In[13]:


Combined_iCEL1314=pd.DataFrame([])


# In[14]:


Combined_iCEL1314['Class']=ClassA['Class']


# In[15]:


Combined_iCEL1314['Development Variation Score']=iCEL1314_Kim['Variation Score']
Combined_iCEL1314['Development Bin']=iCEL1314_Kim['Development Bin']


# In[16]:


iCEL1314_Tissue


# In[17]:


# for gene in iCEL1314_Tissue.index:
#     if (iCEL1314_Tissue.loc[gene]['Tissue bin'])=='Lowly expressed':
#         iCEL1314_Tissue.at[gene, 'Bin']='Lowly expressed'
        


# In[18]:


iCEL1314_Tissue


# In[19]:


Combined_iCEL1314['Tissue Bin']=iCEL1314_Tissue['Bin']


# In[20]:


Combined_iCEL1314['Tissue CV']=iCEL1314_Tissue['CoefVar']


# In[21]:


Combined_iCEL1314.rename(columns={'Variation Score':"Development variation score (VS)"},inplace=True)


# In[22]:


iCEL1314_Compendium.columns


# In[23]:


Combined_iCEL1314['Number of datasets with CV>=0.75']=iCEL1314_Compendium['Number of datasets with CV>=0.75']


# In[24]:


Combined_iCEL1314['Fraction of datasets with CV<0.3']=iCEL1314_Compendium['CVLowFraction']


# In[25]:


Combined_iCEL1314['Compendium bin']=iCEL1314_Compendium['Bin']


# In[26]:


GeneSets=pd.read_excel("/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/Gene2Pathways_090320.xlsx",
                       header='infer',index_col=0)


# In[27]:


Combined_iCEL1314['WormBase ID']=Combined_iCEL1314.index


# In[28]:


Combined_iCEL1314=wb_to_gene(Combined_iCEL1314)


# In[29]:



Combined_iCEL1314.reset_index(inplace=True)
Combined_iCEL1314.set_index(['WormBase ID'],inplace=True)


# In[30]:


Combined_iCEL1314['Level 1']=GeneSets['LEVEL 1']
Combined_iCEL1314['Level 2']=GeneSets['LEVEL 2']
Combined_iCEL1314['Level 3']=GeneSets['LEVEL 3']
Combined_iCEL1314['Level 4']=GeneSets['LEVEL 4']


# In[31]:


Combined_iCEL1314.to_csv("CombinedVS_CV_fraction_021022.csv")


# In[32]:


Combined_iCEL1314


# In[33]:


ReversibilityDict=pd.read_pickle("Drxn2rev.pkl")


# In[34]:


ReversibilityDict


# ## Plotting VS versus CV for all the iCEL1314 genes

# In[35]:


VS=pd.read_csv(Kim_dir+"Binning_021222.csv",index_col=0)


# In[36]:


VS


# In[37]:


metabolicCV=list(set(MetabolicClasses.index).intersection(set(VS.index)))


# In[38]:


MetabolicVS=VS.loc[metabolicCV]


# In[39]:


MetabolicVS=MetabolicVS[MetabolicVS['Variation Score']!='Lowly expressed']


# In[40]:


MetabolicVS.sort_values(ascending=False,by=['Variation Score'])


# In[41]:


CoefVar=pd.read_csv(Tissue_dir+"CVTissue_AllGenes.csv",index_col=0)


# In[42]:


MetabolicVS['CoefVar']=CoefVar['CoefVar']


# In[43]:


CoefVar


# In[44]:


MetabolicVS['GeneID']=CoefVar['gene_id']


# In[45]:


MetabolicVS=MetabolicVS[['GeneID','Variation Score','CoefVar']]


# In[46]:


MetabolicVS.sort_values(ascending=False,by=['Variation Score'],inplace=True)


# In[47]:


# MetabolicCoefVar=CoefVar.loc[metabolicCV]


# In[48]:


# MetabolicCoefVar=MetabolicCoefVar[['gene_id','CoefVar']]


# In[49]:


# MetabolicCoefVar=MetabolicCoefVar[MetabolicCoefVar.CoefVar!=nan]


# In[50]:


MetabolicVS.dropna(inplace=True)


# In[51]:


MetabolicVS=MetabolicVS[MetabolicVS['Variation Score']!='Lowly expressed']


# In[52]:


MetabolicVS=MetabolicVS[MetabolicVS['CoefVar']!='Lowly expressed']


# In[53]:


MetabolicVS


# In[54]:


MetabolicVS.corr(method='spearman')


# In[55]:


MetabolicVS['Class']=MetabolicClasses['Class']


# In[56]:


ClassA=MetabolicVS[MetabolicVS.Class=='A']


# In[57]:


ClassA['Median']=ClassA.median(axis=1)


# In[58]:


ClassA


# In[59]:


CoefVarClassA=ClassA.sort_values(by='Median',ascending=False)


# In[60]:


CoefVarClassA=CoefVarClassA[0:20]


# In[61]:


TopClassA=ClassA[0:15]


# In[62]:


TopClassA=TopClassA.append(CoefVarClassA)


# In[63]:


TopClassA=TopClassA[~TopClassA.index.duplicated(keep='first')]


# In[64]:


TopClassA.shape


# In[65]:


annotations=TopClassA.index


# In[66]:


annotations


# In[67]:


annotations=list(annotations)


# In[68]:


plt.figure(figsize=(24,20))
MetabolicVS.plot.scatter(x='Variation Score',y='CoefVar',c='#56B4E9')
plt.axvline(x=0.169,color='black',linestyle='--')
plt.axhline(0.75,color='black',linestyle='--')


for count,gene in enumerate(annotations):
    print(count,gene)
    X=MetabolicVS.loc[gene]['Variation Score']
    Y=MetabolicVS.loc[gene]['CoefVar']
    label=MetabolicVS.loc[gene]['GeneID']
    plt.text(X,Y,label)

plt.savefig("MetabolicVSCoefVarScatterPlot.svg",dpi=600)


# In[69]:


MetabolicVS
from scipy.stats import pearsonr
from scipy.stats import spearmanr
MetabolicVS.dropna(inplace=True)
spearmanr((MetabolicVS['Variation Score']),(MetabolicVS['CoefVar']))


# In[70]:


MetabolicVS


# ## Dividing metabolic genes into quadrants

# In[71]:


for gene in MetabolicVS.index:
    CV=MetabolicVS.loc[gene]['CoefVar']
    VS=MetabolicVS.loc[gene]['Variation Score']
    if VS>=0.169 and CV>=0.75:
        MetabolicVS.at[gene,'Quadrant']=4
    elif VS>=0.169 and CV<0.75:
        MetabolicVS.at[gene,'Quadrant']=2
    elif VS<0.169 and CV>=0.75:
        MetabolicVS.at[gene,'Quadrant']=3
    else:
        MetabolicVS.at[gene,'Quadrant']=1
        
        
        
        


# In[72]:


Quadrant4=MetabolicVS[MetabolicVS.Quadrant==4]
Quadrant3=MetabolicVS[MetabolicVS.Quadrant==3]
Quadrant2=MetabolicVS[MetabolicVS.Quadrant==2]
Quadrant1=MetabolicVS[MetabolicVS.Quadrant==1]


# In[73]:


Quadrant1.to_csv("Quadrant1.csv")
Quadrant2.to_csv("Quadrant2.csv")
Quadrant3.to_csv("Quadrant3.csv")
Quadrant4.to_csv("Quadrant4.csv")


# In[74]:


ClassA=MetabolicClasses[MetabolicClasses.Class=='A']


# In[75]:


classAQ1=list(set(ClassA.index).intersection(set(Quadrant1.index)))
ClassAQuadrant1=Quadrant1.loc[classAQ1]
ClassAQuadrant1.to_csv("ClassAQuadrant1.csv")
classAQ2=list(set(ClassA.index).intersection(set(Quadrant2.index)))
ClassAQuadrant2=Quadrant2.loc[classAQ2]
ClassAQuadrant2.to_csv("ClassAQuadrant2.csv")
classAQ3=list(set(ClassA.index).intersection(set(Quadrant3.index)))
ClassAQuadrant3=Quadrant3.loc[classAQ3]
ClassAQuadrant3.to_csv("ClassAQuadrant3.csv")


# In[76]:


classAQ4=list(set(ClassA.index).intersection(set(Quadrant4.index)))
ClassAQuadrant4=Quadrant4.loc[classAQ4]


# In[77]:


ClassAQuadrant4.to_csv("ClassAQuadrant4.csv")


# In[78]:


ClassAQuadrant3


# In[79]:


MetabolicVS.to_csv("QuadrantVSCV.csv")


# In[ ]:





# ## Regulated Metabolic genes

# In[80]:


Regulated_Kim=pd.read_csv(Kim_dir+"RegulatedMetabolic_Kim.csv",index_col=0)
Regulated_Tissue=pd.read_csv(Tissue_dir+"RegulatedMetabolic_Tissue.csv",index_col=0)
Regulated_AllDatasets=pd.read_csv(AllDatasets_dir+"Regulated_AllDatasets.csv",index_col=0)


# In[81]:


Regulated_Kim=gene_to_wb(Regulated_Kim)


# In[82]:


Regulated_AllDatasets['Class']=MetabolicClasses['Class']


# In[ ]:





# In[83]:


ClassARegulatedAllDatasets=Regulated_AllDatasets[Regulated_AllDatasets.Class=='A']


# In[84]:


ClassARegulatedAllDatasets.to_csv("ClassARegulatedAllDatasets.csv")


# In[85]:


KimTissue=list(set(Regulated_Kim.index).intersection(set(Regulated_Tissue.index)))


# In[86]:


TissueAllDataset=list(set(Regulated_AllDatasets.index).intersection(set(Regulated_Tissue.index)))


# In[87]:


KimAllDatasets=list(set(Regulated_AllDatasets.index).intersection(set(Regulated_Kim.index)))


# In[88]:


TissueAll=list(set(KimTissue).intersection(set(Regulated_AllDatasets.index)))


# In[89]:


# AllDatasetsGenes=list(set(RegulatedMetabolicAllDatasets.index).difference(set(Regulated_Kim_wb.index)))
# AllDatasetsGenes=list(set(AllDatasetsGenes).difference(set(Regulated_Tissue_wb.index)))


# In[90]:


# AllDatasetsGenes=pd.DataFrame(AllDatasetsGenes)
# AllDatasetsGenes.to_csv("Regulated_OnlyDatasetsExcludingKimTissue.csv")


# ## Finding enrichment of regulated metabolic genes only in all datasets

# In[91]:


def FindWormPathsEnrichment(filename,title):
    WormPathsResult=pd.read_csv(filename,sep='\t')
    WormPathsResult.set_index(['Category'],inplace=True)
# SignificantPathwayOnlyDatasets=FindWormPathsEnrichment(WormPathsResult=)
    pvalse=WormPathsResult['p_enrichment']
    WormPathsResult['FDR-corrected p-enrichment']=(st.stats.multitest.fdrcorrection(pvalse, 
                                                                                    alpha=0.05, 
                                                                                    method='indep', 
                                                                                    is_sorted=False))[1]
    PathwayEnrichmentRegulated=WormPathsResult
#     PathwayEnrichmentRegulated['FDR-corrected p-depletion']=(st.stats.multitest.fdrcorrection(pvalsd, alpha=0.05, method='indep', is_sorted=False))[1]
    SignificantPathwayEnrichment= PathwayEnrichmentRegulated[PathwayEnrichmentRegulated['FDR-corrected p-enrichment']<=0.05]
    SignificantPathwayEnrichment=SignificantPathwayEnrichment[['Enrichment score (n_Hits/n_Genes)',
                                                               'FDR-corrected p-enrichment']]
    SignificantPathwayEnrichment['-log10(BH FDR corrected p-value)']=-(np.log10(SignificantPathwayEnrichment['FDR-corrected p-enrichment']))
    g = SignificantPathwayEnrichment.reset_index()
    survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
    n = g['-log10(BH FDR corrected p-value)']

    norm = plt.Normalize(g['Enrichment score (n_Hits/n_Genes)'].min(), g['Enrichment score (n_Hits/n_Genes)'].max())
    sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
    sm.set_array([])

    ax = sns.barplot(x='-log10(BH FDR corrected p-value)', y='Category', 
                     hue='Enrichment score (n_Hits/n_Genes)', palette='Oranges', dodge=False,data=g)

    ax.set_ylabel('WormPaths Pathway/Category ')

    ax.get_legend().remove()
    from matplotlib import rcParams
    rcParams['font.family'] = 'Arial'
    plt.tight_layout()
    ax.figure.colorbar(sm)
    plt.savefig("WormPathsEnrichment_{}.svg".format(title),dpi=300)
    return SignificantPathwayEnrichment

def FindWormPathsDepletion(filename,title):
    WormPathsResult=pd.read_csv(filename,sep='\t')
    WormPathsResult.set_index(['Category'],inplace=True)
# SignificantPathwayOnlyDatasets=FindWormPathsEnrichment(WormPathsResult=)
    pvalse=WormPathsResult['p_depletion']
    WormPathsResult['FDR-corrected p-depletion']=(st.stats.multitest.fdrcorrection(pvalse, 
                                                                                    alpha=0.05, 
                                                                                    method='indep', 
                                                                                    is_sorted=False))[1]
    PathwayEnrichmentRegulated=WormPathsResult
#     PathwayEnrichmentRegulated['FDR-corrected p-depletion']=(st.stats.multitest.fdrcorrection(pvalsd, alpha=0.05, method='indep', is_sorted=False))[1]
    SignificantPathwayEnrichment= PathwayEnrichmentRegulated[PathwayEnrichmentRegulated['FDR-corrected p-depletion']<=0.05]
    SignificantPathwayEnrichment=SignificantPathwayEnrichment[['Enrichment score (n_Hits/n_Genes)',
                                                               'FDR-corrected p-depletion']]
    SignificantPathwayEnrichment['-log10(BH FDR corrected p-value)']=-(np.log10(SignificantPathwayEnrichment['FDR-corrected p-depletion']))
    g = SignificantPathwayEnrichment.reset_index()
    survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
    n = g['-log10(BH FDR corrected p-value)']

    norm = plt.Normalize(g['Enrichment score (n_Hits/n_Genes)'].min(), g['Enrichment score (n_Hits/n_Genes)'].max())
    sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
    sm.set_array([])

    ax = sns.barplot(x='-log10(BH FDR corrected p-value)', y='Category', 
                     hue='Enrichment score (n_Hits/n_Genes)', palette='Blues', dodge=False,data=g)

    ax.set_ylabel('WormPaths Pathway/Category ')

    ax.get_legend().remove()
    from matplotlib import rcParams
    rcParams['font.family'] = 'Arial'
    plt.tight_layout()
    ax.figure.colorbar(sm)
    plt.savefig("WormPathsDepletion_{}.svg".format(title),dpi=300)
    return SignificantPathwayEnrichment


# In[95]:


RegulatedDatasetsEnrichment=pd.read_csv("RegulatedOnlyCompendium.txt",sep='\t')
# SignificantPathwayOnlyDatasets=FindWormPathsEnrichment(WormPathsResult=)


# In[96]:


SignificantPathwayOnlyDatasets=FindWormPathsDepletion(filename='RegulatedOnlyCompendium.txt',
                                                      title='DepletedPathways_RegulatedOnlyCompendium')


# In[100]:


RegulatedDatasetsEnrichment.set_index(['Category'],inplace=True)
# SignificantPathwayOnlyDatasets=FindWormPathsEnrichment(WormPathsResult=)
pvalse=RegulatedDatasetsEnrichment['p_enrichment']
RegulatedDatasetsEnrichment['FDR-corrected p-enrichment']=(st.stats.multitest.fdrcorrection(pvalse, 
                                                                                    alpha=0.05, 
                                                                                    method='indep', 
                                                                                    is_sorted=False))[1]


# In[102]:


# RegulatedDatasetsEnrichment.set_index(['Category'],inplace=True)
# SignificantPathwayOnlyDatasets=FindWormPathsEnrichment(WormPathsResult=)
pvalse=RegulatedDatasetsEnrichment['p_depletion']
RegulatedDatasetsEnrichment['FDR-corrected p-depletion']=(st.stats.multitest.fdrcorrection(pvalse, 
                                                                                    alpha=0.05, 
                                                                                    method='indep', 
                                                                                    is_sorted=False))[1]


# In[103]:


RegulatedDatasetsEnrichment.sort_values(by='FDR-corrected p-enrichment')


# In[104]:


len(TissueAll)


# In[105]:


len(KimAllDatasets)


# In[106]:


##MetabolicDiet
MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",index_col=0)
metabolicdiet=list(set(MetabolicClasses.index).intersection(set(Regulated_AllDatasets.index)))
RegulatedMetabolicAllDatasets=Regulated_AllDatasets.loc[metabolicdiet]


# In[107]:


RegulatedMetabolicAllDatasets


# In[108]:


Regulated_Kim_wb=gene_to_wb(Regulated_Kim)
Regulated_Tissue_wb=gene_to_wb(Regulated_Tissue)
Regulated_AllDatasets=gene_to_wb(Regulated_AllDatasets)


# In[109]:


Regulated_Kim_wb.to_csv("RegulatedMetabolic_Kim_wb.csv")
Regulated_Tissue_wb.to_csv("Regulated_Tissue_wb.csv")
Regulated_AllDatasets.to_csv("Regulated_AllDatasets_wb.csv")


# In[110]:


Regulated_Kim_wb[~Regulated_Kim_wb.index.str.contains("WB")]


# In[111]:


# Regulated_AllDatasets.sort_values(by=['HighCVNumber'])


# In[112]:


Regulated_Tissue_wb


# In[113]:


## Calculating total regulated metabolic
RegulatedMetabolic=pd.concat([Regulated_Kim_wb,Regulated_Tissue_wb,Regulated_AllDatasets])


# In[114]:


RegulatedMetabolic


# In[115]:


RegulatedMetabolic=RegulatedMetabolic[~RegulatedMetabolic.index.duplicated(keep='first')]


# In[116]:


RegulatedMetabolic.drop(columns=['Variation Score','Gene Cluster','CLASS','HighCVNumber','Bin'],inplace=True)


# In[118]:


RegulatedMetabolic['Class']=MetabolicClasses['Class']


# In[119]:


RegulatedMetabolic.to_csv("RegulatedMetabolic.csv")


# In[120]:


ClassARegulatedMetabolic=RegulatedMetabolic[(RegulatedMetabolic.Class=='A')]


# In[121]:


ClassARegulatedMetabolic.to_csv("ClassARegulatedMetabolic.csv")


# In[122]:


ClassARegulatedMetabolic


# ## Finding common in development, tissue and all datasets

# In[123]:


KimTissue=list(set(Regulated_Kim_wb.index).intersection(set(Regulated_Tissue_wb.index)))


# In[124]:


KimAllDatasets=list(set(Regulated_Kim_wb.index).intersection(set(Regulated_AllDatasets.index)))


# In[125]:


TissueAllDatasets=list(set(Regulated_Tissue_wb.index).intersection(set(Regulated_AllDatasets.index)))


# In[126]:


len(TissueAllDatasets)


# In[127]:


len(list(set(TissueAllDatasets).intersection(set(KimAllDatasets))))


# In[128]:


# Regulated_Kim=gene_to_wb(Regulated_Kim)


# ## Low expressed genes

# In[129]:


LowExpKim=pd.read_csv(Kim_dir+"LowExp_Kim.csv",index_col=0)


# In[130]:


LowExpTissue=pd.read_csv(Tissue_dir+"LowExpGenes.csv",index_col=0)


# In[131]:


leg=list(set(LowExpKim.index).intersection(set(LowExpTissue.index)))
LowExpGenes=LowExpKim.loc[leg]


# In[132]:


LowExpGenes


# In[133]:


LowExpGenes.to_csv("TotalLowExpressedGenes_KimTissue.csv")


# In[134]:


LowExpGenes.to_csv("AllLowExpGenes_KimTissue.csv")


# In[135]:


MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",index_col=0)


# In[136]:


legmetabolic=list(set(LowExpGenes.index).intersection(set(MetabolicClasses.index)))
LowExpGenesMetabolic=LowExpGenes.loc[legmetabolic]


# In[137]:


LowExpGenesMetabolic['Class']=MetabolicClasses['Class']


# In[138]:


LowExpGenesMetabolic[LowExpGenesMetabolic.Class!='A']


# In[139]:


len(list(set(RegulatedMetabolic.index).difference(set(LowExpGenesMetabolic.index))))


# In[140]:


len(list(set(Regulated_AllDatasets.index).difference(set(LowExpGenesMetabolic.index))))


# ### Finding metabolic low expressed genes

# In[141]:


mleg=list(set(MetabolicClasses.index).intersection(set(LowExpGenes.index)))
MetabolicLowExpGenes=LowExpGenes.loc[mleg]


# In[142]:


MetabolicLowExpGenes=MetabolicLowExpGenes[~MetabolicLowExpGenes.index.duplicated(keep='first')]


# In[143]:


MetabolicLowExpGenes.to_csv("MetabolicLowExpGenes_KimTissue.csv")


# ### Finding non-metabolic low expressed genes

# In[144]:


nmleg=list(set(LowExpGenes.index).difference(set(MetabolicClasses.index)))


# In[145]:


NonMetabolicLowExpGenes=LowExpGenes.loc[nmleg]


# In[146]:


NonMetabolicLowExpGenes.to_csv("NonMetabolicLowExpGenes_KimTissue.csv")


# In[147]:


NonMetabolicLowExpGenes


# ## Assigning Class to Kim Regulated

# In[148]:


Regulated_Kim_wb['Class']=Regulated_Kim_wb['CLASS']
Regulated_Kim_wb.drop(columns=['CLASS'],inplace=True)


# ## Reading Bin1 or invariant genes

# #### Metabolic

# In[149]:


Bin1Kim=pd.read_csv(Kim_dir+"Bin1Metabolic_Kim.csv",index_col=0)


# In[150]:


Bin1Kim=gene_to_wb(Bin1Kim)
Bin1Kim=SeqToWB(Bin1Kim)


# In[151]:


Bin1Tissue=pd.read_csv(Tissue_dir+"Bin1_Tissue.csv",index_col=0)


# In[152]:


Bin1AllDatasets=pd.read_csv(AllDatasets_dir+"Bin1_AllDatasets.csv",index_col=0)


# In[153]:


Bin1Tissue.shape


# In[154]:


Bin1DevTissue_Metabolic=list(set(Bin1Kim.index).intersection(set(Bin1Tissue.index)))
Bin1DevTissue_Metabolic=Bin1Kim.loc[Bin1DevTissue_Metabolic]


# In[155]:


# Bin1DevTissue_Metabolic=list(set(Bin1Kim.index).intersection(set(Bin1Tissue.index)))


# In[ ]:


Bin1Kim


# In[ ]:


Bin1AllDatasets


# In[ ]:


Bin1DevTissue_AllDatasets_Metabolic=list(set(Bin1DevTissue_Metabolic.index).intersection(set(Bin1AllDatasets.index)))


# In[ ]:


Bin1DevTissue_Metabolic.index


# In[ ]:


Bin1DevTissue_AllDatasets_Metabolic=Bin1AllDatasets.loc[Bin1DevTissue_AllDatasets_Metabolic]


# In[ ]:


Bin1DevTissue_AllDatasets_Metabolic


# #### Non-metabolic

# In[ ]:


Bin1Kim_NonMetabolic=pd.read_csv(Kim_dir+"Bin1NonMetabolic_Kim.csv",index_col=1)
Bin1Kim_NonMetabolic=gene_to_wb(Bin1Kim_NonMetabolic)
Bin1Kim_NonMetabolic=SeqToWB(Bin1Kim_NonMetabolic)


# In[ ]:


Bin1Tissue_NonMetabolic=pd.read_csv(Tissue_dir+"Bin1NonMetabolic_Tissue.csv",index_col=0)


# In[ ]:


Bin1Tissue_NonMetabolic


# In[ ]:


Bin1AllDatasets_NonMetabolic=pd.read_csv(AllDatasets_dir+"Bin1_NonMetabolic.csv",index_col=0)


# In[ ]:


Bin1AllDatasets_NonMetabolic


# In[ ]:


Bin1Kim_NonMetabolic


# In[ ]:


Bin1NonMetabolic=list(set(Bin1Tissue_NonMetabolic.index).intersection(set(Bin1Kim_NonMetabolic.index)))
Bin1NonMetabolic=Bin1Tissue_NonMetabolic.loc[Bin1NonMetabolic]


# In[ ]:


bin1NonMetabolic=list(set(Bin1AllDatasets_NonMetabolic.index).intersection(set(Bin1NonMetabolic.index)))


# In[ ]:


bin1NonMetabolic


# ### No Bin1 metabolic genes in combined data

# In[ ]:


Regulated_Tissue=pd.DataFrame(Regulated_Tissue)


# In[ ]:


# Regulated_Kim=gene_to_wb(Regulated_Kim)
# Regulated_Temperature=gene_to_wb(Regulated_Temperature)
# Regulated_Tissue=gene_to_wb(Regulated_Tissue)


# In[ ]:


Regulated_Kim_wb.shape


# In[ ]:





# In[ ]:


Regulated_Tissue_wb.shape


# In[ ]:


KimTissue=list(set(Regulated_Kim_wb.index).intersection(set(Regulated_Tissue_wb.index)))


# In[ ]:


KimTissue


# ## Calculating Total regulated metabolic with time, tissue and diet

# In[ ]:


Regulated_Kim_wb=Regulated_Kim_wb['Class']


# In[ ]:


Regulated_Tissue_wb.shape


# In[ ]:


# TotalMetabolicKim['Class']=TotalMetabolicKim['CLASS']


# In[ ]:


# TotalMetabolicKim.drop(columns=['CLASS'],inplace=True)


# In[ ]:


Regulated_Kim_wb=pd.DataFrame(Regulated_Kim_wb)


# In[ ]:


Regulated_Kim_wb


# In[ ]:


RegulatedMetabolic=pd.concat([Regulated_Kim_wb,Regulated_Tissue_wb,Regulated_AllDatasets])


# In[ ]:


RegulatedMetabolic=RegulatedMetabolic[~RegulatedMetabolic.index.duplicated(keep='first')]


# In[ ]:


RegulatedMetabolic.to_csv("RegulatedMetabolic_KimTissue_AllDatasets.csv")


# In[ ]:


RegulatedMetabolic['Class']=MetabolicClasses['Class']


# In[ ]:


RegulatedMetabolic


# In[ ]:


RegulatedMetabolic[(RegulatedMetabolic.Class=='A')]


# In[ ]:


MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",index_col=0)


# In[ ]:


MetabolicClasses=MetabolicClasses[~MetabolicClasses.index.duplicated(keep='first')]


# In[ ]:


ClassA=MetabolicClasses[MetabolicClasses.Class=='A']


# In[ ]:


ClassB=MetabolicClasses[MetabolicClasses.Class=='B']


# In[ ]:


ClassD=MetabolicClasses[MetabolicClasses.Class=='D']


# In[ ]:


ClassA


# In[ ]:


ClassC=MetabolicClasses[MetabolicClasses.Class=='C']


# ## Calculating all C.elegans genes in time, tissue and diet

# In[ ]:


AllGenes_Regulated_Kim=pd.read_csv(Kim_dir+"RegulatedAllGenes.csv",index_col=0)
AllGenes_Regulated_Tissue=pd.read_csv(Tissue_dir+"RegulatedTissueAllGenes.csv",index_col=0)
NonMetabolic_Regulated_AllDatasets=pd.read_csv(AllDatasets_dir+"Regulated_NonMetabolic.csv",index_col=0)


# In[ ]:


AllGenes_Regulated_Tissue


# In[ ]:


AllGenes_Regulated_Kim_wb=gene_to_wb(AllGenes_Regulated_Kim)
AllGenes_Regulated_Tissue_wb=gene_to_wb(AllGenes_Regulated_Tissue)


# In[ ]:


AllGenes_Regulated_Kim_wb=pd.DataFrame(AllGenes_Regulated_Kim_wb.index)


# In[ ]:


AllGenes_Regulated_Kim_wb.set_index(['WormBaseID'],inplace=True)


# In[ ]:


AllGenes_Regulated_Kim_wb[~AllGenes_Regulated_Kim_wb.index.str.contains('WB')]


# In[ ]:


# AllGenes_Regulated_Kim_wb.drop(columns=['Variation Score','Gene Cluster','Bin'],inplace=True)


# In[ ]:


AllGenes_Regulated_Tissue_wb


# In[ ]:


# AllGenes_Regulated_Kim_wb.set_index(['GENEID'],inplace=True)


# In[ ]:


AllGenes_Regulated_Tissue_wb=pd.DataFrame(AllGenes_Regulated_Tissue_wb.index)
AllGenes_Regulated_Tissue_wb.set_index(['gene_id'],inplace=True)


# In[ ]:


AllGenes_Regulated_Tissue_wb


# In[ ]:


AllGenes_KimTissue=list(set(AllGenes_Regulated_Kim_wb.index).intersection(set(AllGenes_Regulated_Tissue_wb.index)))
# AllGenes_KimAll=list(set(AllGenes_Regulated_.index).intersection(set(AllGenes_Regulated_Kim_wb.index)))
# AllGenes_TissueDiet=list(set(AllGenes_Regulated_Tissue_wb.index).intersection(set(AllGenes_Regulated_Diet.index)))


# In[ ]:


AllGenes_TotalKim=pd.read_csv(Kim_dir+"TotalAllGenes_Kim.csv",index_col=1)
# AllGenes_TotalMetabolicKim=gene_to_wb(AllGenes_TotalMetabolicKim)
# AllGenes_TotalMetabolicKim=SeqToWB(AllGenes_TotalMetabolicKim)
AllGenes_TotalTissue=pd.read_csv(Tissue_dir+"TotalGenesTissue.csv",index_col=0)
# AllGenes_TotalMetabolicTissue=gene_to_wb(AllGenes_TotalMetabolicTissue)
# AllGenes_TotalMetabolicTissue=SeqToWB(AllGenes_TotalMetabolicTissue)
# AllGenes_TotalMetabolicDiet=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_121120.csv",index_col=0)
# AllGenes_TotalMetabolicDiet=gene_to_wb(TotalMetabolicDiet)
AllGenes_TotalMetabolicAllDatasets=pd.read_csv(AllDatasets_dir+"Regulated_NonMetabolic.csv",index_col=0)


# In[ ]:


AllGenes_TotalMetabolicAllDatasets


# In[ ]:


# AllGenes_TotalMetabolicKim.set_index(['WormBaseID'],inplace=True)


# In[ ]:


# AllGenes_TotalMetabolicTissue.set_index(['Gene'],inplace=True)


# In[ ]:


# AllGenes_TotalMetabolicTissue.drop(columns=[ 'Neurons', 'Gonad', 'Hypodermis', 'Pharynx',
#        'Body_wall_muscle', 'Glia', 'Intestine'],inplace=True)


# In[ ]:


# AllGenes_TotalMetabolicKim.set_index(['WormBaseID'],inplace=True)


# ## Calculating total genes from three datasets

# In[ ]:


AllGenes_total=pd.concat([AllGenes_TotalKim,AllGenes_TotalTissue,
                          NonMetabolic_Regulated_AllDatasets,AllGenes_TotalMetabolicAllDatasets])


# In[ ]:


AllGenes_total


# In[ ]:


AllGenes_total=AllGenes_total[~AllGenes_total.index.duplicated(keep='first')]


# In[ ]:


magtm=list(set(AllGenes_total.index).intersection(set(MetabolicClasses.index)))
AllGenes_total_metabolic=AllGenes_total.loc[magtm]


# In[ ]:


AllGenes_total_metabolic


# In[ ]:


notregulated=list(set(AllGenes_total_metabolic.index).difference(RegulatedMetabolic.index))
NotRegulated=AllGenes_total_metabolic.loc[notregulated]


# In[ ]:


NotRegulated['Class']=MetabolicClasses['Class']


# In[ ]:


NotRegulated.to_csv("NonRegulated_Metabolic_121822.csv")


# In[ ]:


nmagtm=list(set(AllGenes_total.index).difference(set(MetabolicClasses.index)))
NonMetabolicGenes_total=AllGenes_total.loc[nmagtm]
NonMetabolicGenes_total.to_csv("NonMetabolicGenes_total.csv")


# In[ ]:


NonMetabolicGenes_total


# In[ ]:


AllGenes_Regulated=pd.concat([AllGenes_Regulated_Kim_wb,AllGenes_Regulated_Tissue_wb,NonMetabolic_Regulated_AllDatasets])


# In[ ]:


AllGenes_Regulated


# In[ ]:


AllGenes_Regulated=AllGenes_Regulated[~AllGenes_Regulated.index.duplicated(keep='first')]


# In[ ]:


# AllGenes_Regulated[~AllGenes_Regulated.index.str.contains("WB")]


# In[ ]:


## Finding non-metabolic regulated genes
nmag=list(set(AllGenes_Regulated.index).difference(set(MetabolicClasses.index)))
NonMetabolicGenesRegulated=AllGenes_Regulated.loc[nmag]


# In[ ]:


NonMetabolicGenesRegulated


# In[ ]:


AllGenes_total=AllGenes_total[~AllGenes_total.index.duplicated(keep='first')]


# In[ ]:


agm=list(set(MetabolicClasses.index).intersection(set(AllGenes_Regulated.index)))


# In[ ]:


AllGenes_total


# In[ ]:





# ## Finding common in tissue, time and diet

# In[ ]:


Regulated_Kim_wb=Regulated_Kim_wb[~Regulated_Kim_wb.index.duplicated(keep='first')]
Regulated_Tissue_wb=Regulated_Tissue_wb[~Regulated_Tissue_wb.index.duplicated(keep='first')]


# In[ ]:


KimTissue=list(set(Regulated_Kim_wb.index).intersection(set(Regulated_Tissue_wb.index)))
len(KimTissue)


# In[ ]:


KimDatasets=list(set(Regulated_Kim_wb.index).intersection(set(RegulatedMetabolicAllDatasets.index)))


# In[ ]:


len(KimDatasets)


# In[ ]:


TissueDatasets=list(set(Regulated_Tissue_wb.index).intersection(Regulated_AllDatasets.index))


# In[ ]:


len(TissueDatasets)


# In[ ]:


diff=list(set(Regulated_AllDatasets.index).difference(set(TissueDatasets)))


# In[ ]:


diff=list(set(diff).difference(set(KimDatasets)))


# In[ ]:


KimTissueDatasets=list(set(KimTissue).intersection(set(KimDatasets)))


# In[ ]:


len(KimTissueDatasets)


# In[ ]:


fig=plt.figure(figsize=(4,3))
venn3(subsets = (578, 109, 350, 45,177,132,400),set_labels = ('Tissue', 'Compendium','Development'))
plt.rcParams["font.family"] = "Arial"
plt.savefig("KimTissueCompendiumVenn.svg",dpi=300)
plt.show()


# In[ ]:


Bin2KimNonMetabolic=pd.read_csv(Kim_dir+"Bin2NonMetabolic_Kim.csv",index_col=0)
# NonRegulatedKim=gene_to_wb(NonRegulatedKim)
Bin2TissueNonMetabolic=pd.read_csv(Tissue_dir+"Bin2NonMetabolic_Tissue.csv",index_col=0)
# NonRegulatedTissue=gene_to_wb(NonRegulatedTissue)
# NonRegulatedTemperature=pd.read_csv(Temperature_dir+"Zless3Temperature.csv",index_col=0)
# NonRegulatedTemperature=gene_to_wb(NonRegulatedTemperature)


# In[ ]:


# NonRegulatedKim.rename(index={'nol-2':'nsun-1'},inplace=True)


# In[ ]:


# NonRegulatedKim=gene_to_wb(NonRegulatedKim)


# In[ ]:


# NonRegulatedKim=SeqToWB(NonRegulatedKim)


# In[ ]:


# NonRegulatedKim[~NonRegulatedKim.index.str.contains('WB')]


# In[ ]:


# NonRegulatedTissue[~NonRegulatedTissue.index.str.contains('WB')]


# In[ ]:


# inter2=list(set(Bin1Kim.index).intersection(set(Bin1Temperature.index)))


# In[ ]:


# inter3=list(set(inter).intersection(set(inter2)))


# In[ ]:


# NonRegulated=pd.concat([NonRegulatedKim,NonRegulatedTissue,NonRegulatedTemperature])


# In[ ]:


# NonRegulated=NonRegulated[~NonRegulated.index.duplicated(keep='first')]


# In[ ]:


TotalMetabolicKim=pd.read_csv(Kim_dir+"TotalMetabolic_Kim.csv",index_col=0)


# In[ ]:


TotalMetabolicTissue=pd.read_csv(Tissue_dir+"TotalMetabolic_Tissue.csv",index_col=0)


# In[ ]:


# TotalMetabolicTissue.drop(columns=['Unnamed: 0'],inplace=True)


# In[ ]:


TotalMetabolicTissue


# In[ ]:


TotalMetabolicTissue['Class']=MetabolicClasses['Class']


# In[ ]:


TotalMetabolicTissue


# In[ ]:


TotalMetabolic=pd.concat([TotalMetabolicKim,TotalMetabolicTissue,TotalMetabolic])


# In[ ]:


TotalMetabolic=TotalMetabolic[~TotalMetabolic.index.duplicated(keep='first')]


# In[ ]:


# NonRegulated=NonRegulated[~NonRegulated.index.duplicated(keep='first')]


# In[ ]:


RegulatedMetabolic=RegulatedMetabolic[~RegulatedMetabolic.index.duplicated(keep='first')]


# In[ ]:


TotalMetabolic


# In[ ]:


RegulatedMetabolic


# In[ ]:


nrm=list(set(TotalMetabolic.index).difference(set(RegulatedMetabolic.index)))
NonRegulatedMetabolic=TotalMetabolic.loc[nrm]


# In[ ]:


NonRegulatedMetabolic.to_csv("NonRegulatedMetabolic_KimTissueDiet.csv")


# In[ ]:


RegulatedMetabolic.to_csv("RegulatedMetabolic_KimTissueDiet.csv")


# In[ ]:


RegulatedMetabolic


# In[ ]:


list1=set(TotalMetabolic.index)
list2=set(RegulatedMetabolic.index)


# In[ ]:


difference=list1.difference(list2)


# In[ ]:


len(difference)


# In[ ]:


Non_Regulated_diff=TotalMetabolic.loc[difference]
Non_Regulated_diff.to_csv("Non_Regulated_diff.csv")


# In[ ]:


Non_Regulated_diff


# In[ ]:


(~(RegulatedMetabolic.Class=='A')).sum()


# In[ ]:


(~(TotalMetabolic.Class=='A')).sum()


# In[ ]:


# NonRegulated=wb_to_gene(NonRegulated)


# In[ ]:


RegulatedMetabolic.to_csv("RegulatedMetabolic_gene.csv")


# In[ ]:


Regulated_Diet[~(Regulated_Diet.Class=='A')]


# In[ ]:


RegulatedMetabolic=gene_to_wb(RegulatedMetabolic)


# In[ ]:


RegulatedMetabolic.to_csv("RegulatedMetabolic_wb.csv")


# In[ ]:


RegulatedMetabolic.replace(np.nan,0,inplace=True)


# In[ ]:


RegulatedMetabolic


# In[ ]:


RegulatedClass=RegulatedMetabolic[RegulatedMetabolic.Class!=0]


# In[ ]:


Bin1DevTissue


# In[ ]:


RegulatedMetabolic_Final=RegulatedClass


# In[ ]:


RegulatedMetabolic_Final


# In[ ]:


MetabolicClasses=pd.read_csv("MetabolicClasses_SN_121020.csv",index_col=0)
MetabolicClasses=SeqToGene(MetabolicClasses)
MetabolicClasses=gene_to_wb(MetabolicClasses)


# In[ ]:


MetabolicClasses=MetabolicClasses[~MetabolicClasses.index.duplicated(keep='first')]


# In[ ]:


TotalMetabolic['Class']=MetabolicClasses['Class']


# In[ ]:


TotalMetabolic


# In[ ]:


## Finding transcriptionally regulated dev genes
classes=['A','B','C','D']
TotalClassA=(TotalMetabolic.Class=='A').sum()
RegulatedClassA=(RegulatedMetabolic_Final.Class=='A').sum()
TotalClassB=(TotalMetabolic.Class=='B').sum()
RegulatedClassB=(RegulatedMetabolic_Final.Class=='B').sum()
TotalClassC=(TotalMetabolic.Class=='C').sum()
RegulatedClassC=(RegulatedMetabolic_Final.Class=='C').sum()
TotalClassD=(TotalMetabolic.Class=='D').sum()
RegulatedClassD=(RegulatedMetabolic_Final.Class=='D').sum()
regulated_vals =[];
non_regulated_vals =[];
for val in classes:
    regulated_vals.append((RegulatedMetabolic_Final.Class==val).sum());
    non_regulated_vals=[(TotalClassA-RegulatedClassA),
                        (TotalClassB-RegulatedClassB),
                        (TotalClassC-RegulatedClassC),
                        (TotalClassD-RegulatedClassD)];
print(regulated_vals)
print(non_regulated_vals)


# In[ ]:


classes = list(classes)
print(classes)


# In[ ]:


RegulatedMetabolic_Final=RegulatedMetabolic_Final[~RegulatedMetabolic_Final.duplicated(keep='first')]
TotalMetabolic=TotalMetabolic[~TotalMetabolic.duplicated(keep='first')]


# In[ ]:


TotalMetabolic


# In[ ]:


fig, ax = plt.subplots(figsize=(12,10))

size = 0.3
vals = np.array([regulated_vals, non_regulated_vals])
outer_labels=["Regulated", "Invariant"]
inner_labels = classes
print(inner_labels)
print(vals)
explode = (0.1, 0)
cmap = plt.get_cmap("tab20c")
outer_colors = ['#54CB73','#F16718']
inner_colors = ['#F0EFC0','#92E7DA','#E7C1E5','cyan','#F0EFC0','#92E7DA','#E7C1E5','cyan']
inner_explode=(0.1,0.1,0.1,0.1,0,0,0,0)
ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w',linewidth=3),shadow=False,autopct='%1.1f%%',
       textprops={'size': 'larger','fontweight':'bold'},
      pctdistance=0.85,labeldistance=1.05,startangle=90,counterclock=False)

ax.pie(vals.flatten(),radius=1-size, colors= inner_colors, labeldistance=0.5,
       wedgeprops=dict(width=0.5, edgecolor='w',linewidth='3'),pctdistance=0.90,startangle=90,counterclock=False)
plt.legend(loc='upper right',bbox_to_anchor=(1.22, 0.85),labels=outer_labels+inner_labels)
ax.set_title("Pie chart showing percentage of transcriptionally regulated metabolic genes",pad=3,fontweight="bold",fontsize=12)

plt.savefig("DEVTissueTemperatureRegulationNoLabel.svg")
plt.tight_layout()
plt.show()


# In[ ]:


fig, ax = plt.subplots(figsize=(12,10))
size = 0.3
vals = np.array([regulated_vals, non_regulated_vals])
outer_labels=["Regulated", "Non Regulated"]
inner_labels = classes
print(inner_labels)
explode = (0.1, 0)
cmap = plt.get_cmap("tab20c")
outer_colors = ['#54CB73','#F16718']
inner_colors = ['#F0EFC0','#92E7DA','#E7C1E5','cyan','#F0EFC0','#92E7DA','#E7C1E5','cyan']
inner_explode=(0.1,0.1,0.1,0.1,0,0,0,0)
ax.pie(vals.sum(axis=1), radius=1, labels=outer_labels, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w',linewidth=3),autopct='%1.1f%%',shadow=False,
       textprops={'size': 'larger','fontweight':'bold'},
      pctdistance=0.85,labeldistance=1.05,startangle=90,counterclock=False)

ax.pie(vals.flatten(),radius=1-size, colors= inner_colors, labeldistance=0.5,autopct='%1.1f%%',
       wedgeprops=dict(width=0.5, edgecolor='w',linewidth='3'),pctdistance=0.90,startangle=90,counterclock=False)
plt.legend(loc='upper right',bbox_to_anchor=(1.02, 0.85),labels=outer_labels+inner_labels)
ax.set_title("Pie chart showing percentage of transcriptionally regulated metabolic genes",pad=3,fontweight="bold",fontsize=12)

plt.savefig("DEVTissueTemperatureRegulation.svg")
plt.tight_layout()
plt.show()


# In[ ]:


Bin1


# In[ ]:


## Calculating non-metabolic percentage regulation
AllgenesKim=pd.read_csv(Kim_dir+"HighModLowExpData.csv",index_col=0)


# In[ ]:


AllgenesTissue=pd.read_csv(Tissue_dir+"TotalGenesTissue.csv",index_col=0)


# In[ ]:


AllgenesDiet=pd.read_csv(Diet_dir+"RegulatedDietAllGenes.csv",index_col=0)


# In[ ]:


AllGenes_total=pd.DataFrame(AllGenes_total.index)


# In[ ]:


AllGenes_total.set_index([0],inplace=True)


# In[ ]:


AllGenes_Regulated


# In[ ]:


nmt=list(set(AllGenes_total.index).difference(MetabolicClasses.index))
NonMetabolicTotal=AllGenes_total.loc[nmt]


# In[ ]:


NonMetabolicTotal.to_csv("NonMetabolicTotal.csv")


# In[ ]:


nmr=list(set(AllGenes_Regulated.index).difference(set(MetabolicClasses.index)))
NonMetabolicRegulated=AllGenes_Regulated.loc[nmr]


# In[ ]:


NonMetabolicRegulated.to_csv("NonMetabolicRegulated.csv")


# ## Enrichment analysis

# ### Phenotype depletion analysis

# In[ ]:


# Phenotype Enrichment Analysis
PhenotypeDepleted=pd.read_csv("DepletedPhenotypes_Bin3_FC_=0.01.csv",index_col=0)


# In[ ]:


PhenotypeDepleted.sort_values(by=['-log10(BH FDR-corrected p-value)'],ascending=False,inplace=True)


# In[ ]:


import seaborn as sns
g = PhenotypeDepleted.reset_index()
survival_rates = g['Enrichment fold change'].mean()
n = g['-log10(BH FDR-corrected p-value)']

norm = plt.Normalize(g['Enrichment fold change'].min(), g['Enrichment fold change'].max())
sm = plt.cm.ScalarMappable(cmap="Blues", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR-corrected p-value)', y='WormBase phenotype name', 
                 hue='Enrichment fold change', palette='Blues', 
                 dodge=False,data=g)

ax.set_ylabel('WormBase phenotype term')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormBasePhenotypeDepletion.svg",dpi=300)


# In[ ]:


g


# ### Phenotye enrichment analysis

# In[ ]:


import seaborn as sns
PhenotypeEnriched=pd.read_csv("EnrichedPhenotypes_Bin3.xlsx - P_=0.01.csv",index_col=0)
g = PhenotypeEnriched.reset_index()
survival_rates = g['Enrichment fold change'].mean()
n = g['-log10(BH FDR corrected p-value)']

norm = plt.Normalize(g['Enrichment fold change'].min(), g['Enrichment fold change'].max())
sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR corrected p-value)', y='Phenotype Name', 
                 hue='Enrichment fold change', palette='Oranges', 
                 dodge=False,data=g)

ax.set_ylabel('Phenotype Name')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormBasePhenotypeEnrichment.svg",dpi=300)


# ### Pathway enrichment and depletion analysis with WormPaths

# In[ ]:


PathwayEnrichmentRegulated=pd.read_csv("PEA_StatsTable.txt",sep='\t')


# In[ ]:


PathwayEnrichmentRegulated.set_index(['Category'],inplace=True)


# In[ ]:


PathwayEnrichmentRegulated=PathwayEnrichmentRegulated[~PathwayEnrichmentRegulated.index.duplicated(keep='first')]


# In[ ]:


pvalse=PathwayEnrichmentRegulated.p_enrichment
PathwayEnrichmentRegulated['FDR-corrected p-enrichment']=(st.stats.multitest.fdrcorrection(pvalse, alpha=0.05, method='indep', 
                                                                                           is_sorted=False))[1]


# In[ ]:


pvalsd=PathwayEnrichmentRegulated.p_depletion
PathwayEnrichmentRegulated['FDR-corrected p-depletion']=(st.stats.multitest.fdrcorrection(pvalsd, alpha=0.05, method='indep', is_sorted=False))[1]


# In[ ]:


PathwayEnrichmentRegulated.sort_values(by=['FDR-corrected p-depletion'])[0:20]


# In[ ]:


PathwayEnrichmentRegulated.columns


# In[ ]:





# ### WormPaths significant pathway enrichment analysis (p<=0.01)

# In[ ]:


SignificantPathwayEnrichment= PathwayEnrichmentRegulated[PathwayEnrichmentRegulated['FDR-corrected p-enrichment']<=0.01]


# In[ ]:


SignificantPathwayEnrichment=SignificantPathwayEnrichment[['Enrichment score (n_Hits/n_Genes)',
                                                           'FDR-corrected p-enrichment']]


# In[ ]:


SignificantPathwayEnrichment['-log10(BH FDR corrected p-value)']=-(np.log10(SignificantPathwayEnrichment['FDR-corrected p-enrichment']))


# In[ ]:


SignificantPathwayEnrichment


# In[ ]:


g = SignificantPathwayEnrichment.reset_index()
survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
n = g['-log10(BH FDR corrected p-value)']

norm = plt.Normalize(g['Enrichment score (n_Hits/n_Genes)'].min(), g['Enrichment score (n_Hits/n_Genes)'].max())
sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR corrected p-value)', y='Category', 
                 hue='Enrichment score (n_Hits/n_Genes)', palette='Oranges', 
                 dodge=False,data=g)

ax.set_ylabel('WormPaths Pathway/Category ')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormPathsEnrichment.svg",dpi=300)


# ### WormPaths significant pathway depletion analysis

# In[ ]:


SignificantPathwayDepletion= PathwayEnrichmentRegulated[PathwayEnrichmentRegulated['FDR-corrected p-depletion']<=0.01]


# In[ ]:


SignificantPathwayDepletion=SignificantPathwayDepletion[['Enrichment score (n_Hits/n_Genes)',
                                                           'FDR-corrected p-depletion']]


# In[ ]:


SignificantPathwayDepletion['-log10(BH FDR corrected p-value)']=-(np.log10(SignificantPathwayDepletion['FDR-corrected p-depletion']))


# In[ ]:


SignificantPathwayDepletion.sort_values(by=['-log10(BH FDR corrected p-value)'],ascending=False,inplace=True)


# In[ ]:


SignificantPathwayDepletion


# In[ ]:


g = SignificantPathwayDepletion.reset_index()
survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
n = g['-log10(BH FDR corrected p-value)']

norm = plt.Normalize(g['Enrichment score (n_Hits/n_Genes)'].min(), g['Enrichment score (n_Hits/n_Genes)'].max())
sm = plt.cm.ScalarMappable(cmap="Blues", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR corrected p-value)', y='Category', 
                 hue='Enrichment score (n_Hits/n_Genes)', palette='Blues', 
                 dodge=False,data=g)

ax.set_ylabel('WormPaths Pathway/Category ')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormPathsDepletion.svg",dpi=300)


# ## Enrichment and depletion analysis of invariant or bin1 genes

# In[ ]:


Bin1Pathway=pd.read_csv("PEA_StatsTable_NonRegulated.txt",sep='\t')


# In[ ]:


Bin1Pathway.set_index(['Category'],inplace=True)


# In[ ]:


Bin1Pathway


# In[ ]:


Bin1Pathway=Bin1Pathway[~Bin1Pathway.index.duplicated(keep='first')]


# In[ ]:


pvalse=Bin1Pathway['p_enrichment']
Bin1Pathway['FDR-corrected p-enrichment']=(st.stats.multitest.fdrcorrection(pvalse, alpha=0.05, method='indep', is_sorted=False))[1]


# In[ ]:





# In[ ]:


pvalsd=Bin1Pathway['p_depletion']
Bin1Pathway['FDR-corrected p-depletion']=(st.stats.multitest.fdrcorrection(pvalsd, alpha=0.05, 
                                                                              method='indep', is_sorted=False))[1]


# In[ ]:


Bin1Pathway


# In[ ]:


SignificantPathwayDepletion_NonRegulated= Bin1Pathway[Bin1Pathway['FDR-corrected p-depletion']<=0.01]


# In[ ]:


SignificantPathwayEnrichment_NonRegulated= Bin1Pathway[Bin1Pathway['FDR-corrected p-enrichment']<=0.01]


# In[ ]:


SignificantPathwayDepletion_NonRegulated


# In[ ]:


SignificantPathwayDepletion_NonRegulated['-log10(BH FDR corrected p-value)']=-(np.log10(SignificantPathwayDepletion_NonRegulated['FDR-corrected p-depletion']))


# In[ ]:


SignificantPathwayEnrichment_NonRegulated['-log10(BH FDR corrected p-value)']=-(np.log10(SignificantPathwayEnrichment_NonRegulated['FDR-corrected p-enrichment']))


# In[ ]:


g = SignificantPathwayEnrichment_NonRegulated.reset_index()
survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
n = g['-log10(BH FDR corrected p-value)']

norm = plt.Normalize(g['Enrichment score (n_Hits/n_Genes)'].min(), g['Enrichment score (n_Hits/n_Genes)'].max())
sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR corrected p-value)', y='Category', 
                 hue='Enrichment score (n_Hits/n_Genes)', palette='Oranges', 
                 dodge=False,data=g)

ax.set_ylabel('WormPaths Pathway/Category ')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormPathsEnrichment_NonRegulated.svg",dpi=300)


# In[ ]:


g = SignificantPathwayDepletion_NonRegulated.reset_index()
g.sort_values(by=['-log10(BH FDR corrected p-value)'],ascending=False,inplace=True)
survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
n = g['-log10(BH FDR corrected p-value)']

norm = plt.Normalize(g['Enrichment score (n_Hits/n_Genes)'].min(), g['Enrichment score (n_Hits/n_Genes)'].max())
sm = plt.cm.ScalarMappable(cmap="Blues", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR corrected p-value)', y='Category', 
                 hue='Enrichment score (n_Hits/n_Genes)', palette='Blues', 
                 dodge=False,data=g)

ax.set_ylabel('WormPaths Pathway/Category ')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormPathsDepletion_NonRegulated.svg",dpi=300)


# In[ ]:


Bin1Phenotype=pd.read_csv("Bin1PhenotypeEnrichment_102721.xlsx - p_=0.01.csv",index_col=0)


# In[ ]:


Bin1Phenotype


# In[ ]:


Bin1Phenotype=Bin1Phenotype[Bin1Phenotype['BH FDR corrected P-value']<=0.01]


# In[ ]:


Bin1Phenotype['-log10(BH FDR corrected p-value)']=-(np.log10(Bin1Phenotype['BH FDR corrected P-value']))


# In[ ]:


g = Bin1Phenotype.reset_index()
g.sort_values(by=['-log10(BH FDR corrected p-value)'],ascending=False,inplace=True)
survival_rates = g['Enrichment fold change (Observed/ Expected)'].mean()
n = g['-log10(BH FDR corrected p-value)']

norm = plt.Normalize(g['Enrichment fold change (Observed/ Expected)'].min(), 
                     g['Enrichment fold change (Observed/ Expected)'].max())
sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR corrected p-value)', y='Phenotype Name', 
                 hue='Enrichment fold change (Observed/ Expected)', palette='Oranges', 
                 dodge=False,data=g)

ax.set_ylabel('WormBase Phenotype term ')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormBase_Phenotype_Enrichment_Bin1.svg",dpi=300)


# ## Non metabolic genes

# In[ ]:


NMRegulated=pd.read_csv("NonMetabolicEnrichment_FDR_0.01.csv")


# In[ ]:


NMRegulated.columns


# In[ ]:


g = NMRegulated.reset_index()
g.sort_values(by=['-log10(BH-FDR corrected p-value)'],ascending=False,inplace=True)
survival_rates = g['Enrichment fold change'].mean()
n = g['-log10(BH-FDR corrected p-value)']

norm = plt.Normalize(g['Enrichment fold change'].min(), 
                     g['Enrichment fold change'].max())
sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH-FDR corrected p-value)', y='Phenotype Name', 
                 hue='Enrichment fold change', palette='Oranges', 
                 dodge=False,data=g)

ax.set_ylabel('WormBase Phenotype term ')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormBase_Phenotype_Enrichment_NonMetabolicRegulated.svg",dpi=300)


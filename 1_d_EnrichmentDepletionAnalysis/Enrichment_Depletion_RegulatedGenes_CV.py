
# coding: utf-8

# In[19]:


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
get_ipython().run_line_magic('matplotlib', 'inline')
from pylab import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
Base_dir='/data/nandas/Transcription/CombinedTimeSpaceConditions/'
os.chdir(Base_dir)
from scipy.signal import find_peaks
# from CatExp import *
import statsmodels as st
# import statsmodels.api as sm


# In[86]:


Tissue_dir='/data/nandas/Combined_coexp/Part_1_TranscriptionallyRegulatedGenes/Tissue/'
Kim_dir='/data/nandas/Transcription/KimDevTime_071620/'


# In[3]:


mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',
                      index_col=0)
mapper_df=mapper_df.loc[mapper_df.index.dropna()]


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
    


# ## Reading development and tissue binning files

# In[89]:


VS=pd.read_csv(Kim_dir+"Binning_021222.csv",index_col=0)


# In[90]:


CoefVar=pd.read_csv(Tissue_dir+"CVTissue_AllGenes.csv",index_col=0)


# In[91]:


intersect=list(set(VS.index).intersection(set(CoefVar.index)))


# In[92]:


CommonVS=VS.loc[intersect]


# In[93]:


CommonCoefVar=CoefVar.loc[intersect]


# In[94]:


VS


# In[95]:


CommonVS=CommonVS[['Variation Score','Bin']]
CommonCoefVar=CommonCoefVar[['CoefVar','Bin']]


# In[96]:


CommonVS['DevelopmentBin']=CommonVS['Bin']


# In[97]:


CommonCoefVar['TissueBin']=CommonCoefVar['Bin']


# In[98]:


CommonVSCoefVar=CommonCoefVar[['CoefVar','TissueBin']]


# In[99]:


CommonVSCoefVar[['Variation Score','DevelopmentBin']]=CommonVS[['Variation Score','DevelopmentBin']]


# In[100]:


CommonVSCoefVar


# In[101]:


CommonVS=pd.DataFrame(CommonVS)
CommonCoefVar=pd.DataFrame(CommonCoefVar)


# In[102]:


# CommonVSCoefVar['DevelopmentBin']=str(CommonVSCoefVar['DevelopmentBin'])


# In[103]:


CommonVSCoefVar


# In[104]:


LowExpKim=pd.read_csv(Kim_dir+"LowExp_Kim.csv",index_col=0)


# In[105]:


LowExpKim


# In[106]:


list(set(LowExpKim.index).intersection(set(CommonVSCoefVar.index)))


# In[107]:


(CommonVSCoefVar['DevelopmentBin']).replace(3,'Highly variant',inplace=True)
(CommonVSCoefVar['DevelopmentBin']).replace(2,'Moderately variant',inplace=True)
(CommonVSCoefVar['DevelopmentBin']).replace(1,'Invariant',inplace=True)


# In[108]:


np.unique(CommonVSCoefVar.TissueBin)


# In[109]:


CommonVSCoefVar[CommonVSCoefVar.TissueBin=='Lowly expressed']


# In[110]:


for genes in CommonVSCoefVar.index:
#     print(genes)
    if genes in LowExpKim.index:
        print(genes)
        CommonVSCoefVar.at[genes,'DevelopmentBin']='Lowly expressed'
        print(CommonVSCoefVar.loc[genes,'DevelopmentBin'])
#     elif (CommonVSCoefVar.loc[genes,'DevelopmentBin'])==1:
#         CommonVSCoefVar.at[genes,'DevelopmentBin']='Invariant'
#     elif (CommonVSCoefVar.loc[genes,'DevelopmentBin'])==3:
#         print(type((CommonVSCoefVar.loc[genes,'DevelopmentBin'])))
#         CommonVSCoefVar.at[genes,'DevelopmentBin']='Highly variant'
#     elif (CommonVSCoefVar.loc[genes,'DevelopmentBin'])==2:
#         CommonVSCoefVar.at[genes,'DevelopmentBin']='Moderately variant'


# In[111]:


np.unique(CommonVSCoefVar.DevelopmentBin)


# ## Removing genes that are lowly expressed in both development and tissues

# In[112]:


count=0
for genes in CommonVSCoefVar.index:
    if (CommonVSCoefVar.loc[genes,'DevelopmentBin'])=='Lowly expressed' :
        if  (CommonVSCoefVar.loc[genes,'TissueBin'])=='Lowly expressed':
            print(genes)
            count=count+1
            print(count)
            CommonVSCoefVar.drop(index=genes,inplace=True)


# In[113]:


CommonVSCoefVar


# ## Quadrant enrichment of all genes

# In[114]:


for gene in CommonVSCoefVar.index:
    CV=CommonVSCoefVar.loc[gene]['TissueBin']
    VS=CommonVSCoefVar.loc[gene]['DevelopmentBin']
    if VS=='Highly variant' and CV=='Highly variant':
        CommonVSCoefVar.at[gene,'Quadrant']=4
    elif VS=='Highly variant' and CV!='Highly variant':
        CommonVSCoefVar.at[gene,'Quadrant']=2
    elif VS!='Highly variant' and CV=='Highly variant':
        CommonVSCoefVar.at[gene,'Quadrant']=3
    else:
        CommonVSCoefVar.at[gene,'Quadrant']=1


# In[115]:


MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",index_col=0)


# In[116]:


CommonVSCoefVar


# In[117]:


CommonVSCoefVar['Class']=MetabolicClasses['Class']


# In[118]:


CommonVSCoefVar.Class.sort_values(ascending=False)


# In[119]:


CommonVSCoefVar


# In[120]:


AllGenesQuadrant1=CommonVSCoefVar[CommonVSCoefVar.Quadrant==1]
AllGenesQuadrant2=CommonVSCoefVar[CommonVSCoefVar.Quadrant==2]
AllGenesQuadrant3=CommonVSCoefVar[CommonVSCoefVar.Quadrant==3]
AllGenesQuadrant4=CommonVSCoefVar[CommonVSCoefVar.Quadrant==4]


# In[121]:


totalmetabolic=list(set(CommonVSCoefVar.index).intersection(set(MetabolicClasses.index)))


# In[122]:


MetabolicCommonVSCoefVar=CommonVSCoefVar.loc[totalmetabolic]


# In[123]:


MetabolicCommonVSCoefVar['WormBase ID']=MetabolicCommonVSCoefVar.index


# In[124]:


MetabolicCommonVSCoefVar=wb_to_gene(MetabolicCommonVSCoefVar)


# In[125]:


MetabolicCommonVSCoefVar.to_csv("QuadrantCategoriesMetabolic.csv")


# In[126]:


metabolicq1=list(set(MetabolicClasses.index).intersection(set(AllGenesQuadrant1.index)))
MetabolicGenesQuadrant1=AllGenesQuadrant1.loc[metabolicq1]


# In[127]:


metabolicq2=list(set(MetabolicClasses.index).intersection(set(AllGenesQuadrant2.index)))
MetabolicGenesQuadrant2=AllGenesQuadrant2.loc[metabolicq2]
metabolicq3=list(set(MetabolicClasses.index).intersection(set(AllGenesQuadrant3.index)))
MetabolicGenesQuadrant3=AllGenesQuadrant3.loc[metabolicq3]
metabolicq4=list(set(MetabolicClasses.index).intersection(set(AllGenesQuadrant4.index)))
MetabolicGenesQuadrant4=AllGenesQuadrant4.loc[metabolicq4]


# In[128]:


len(AllGenesQuadrant2.index)


# In[129]:


PercentageQ1=(len(metabolicq1)*100)/len(AllGenesQuadrant1.index)
PercentageQ2=(len(metabolicq2)*100)/len(AllGenesQuadrant2.index)
PercentageQ3=(len(metabolicq3)*100)/len(AllGenesQuadrant3.index)
PercentageQ4=(len(metabolicq4)*100)/len(AllGenesQuadrant4.index)


# In[130]:


PercentageQ4


# In[131]:


def hypergeom(x,M,n,N):
    #M is the population size (previously N)
    #n is the number of successes in the population (previously K)
    #N is the sample size (previously n)
    #X is still the number of drawn “successes”.
    from scipy.stats import hypergeom
    pval = hypergeom.sf(x-1, M, n, N)
    return pval


# In[132]:


len(totalmetabolic)/len(CommonVSCoefVar.index)


# In[133]:


len(totalmetabolic)


# In[134]:


pvalq1=hypergeom(x=len(metabolicq1),N=len(AllGenesQuadrant1.index),n=len(totalmetabolic),
                 M=len((CommonVSCoefVar.index)))


# In[135]:


pvalq2=hypergeom(x=len(metabolicq2),N=len(AllGenesQuadrant2.index),n=len(totalmetabolic),
                 M=len((CommonVSCoefVar.index)))


# In[136]:


pvalq3=hypergeom(x=len(metabolicq3),N=len(AllGenesQuadrant3.index),n=len(totalmetabolic),
                 M=len((CommonVSCoefVar.index)))


# In[137]:


pvalq4=hypergeom(x=len(metabolicq4),N=len(AllGenesQuadrant4.index),n=len(totalmetabolic),
                 M=len((CommonVSCoefVar.index)))


# ## Extracting metabolic VS and CV from development and tissue dataset

# In[138]:


MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",index_col=0)


# In[139]:


CommonVS


# In[140]:


MetabolicVS=list(set(CommonVSCoefVar.index).intersection(set(MetabolicClasses.index)))
MetabolicVS=MetabolicClasses.loc[MetabolicVS]


# In[141]:


MetabolicVS ['Variation Score']=CommonVS ['Variation Score']
MetabolicVS ['DevelopmentBin']= CommonVS ['DevelopmentBin']


# In[142]:


MetabolicCoefVar=list(set(CommonCoefVar.index).intersection(set(MetabolicClasses.index)))
MetabolicCoefVar=MetabolicClasses.loc[MetabolicCoefVar]


# In[143]:


MetabolicCoefVar ['CoefVar']=CommonCoefVar.CoefVar
MetabolicCoefVar['TissueBin']=CommonCoefVar['TissueBin']


# In[144]:


MetabolicVS


# In[145]:


for gene in MetabolicCoefVar.index:
    CV=MetabolicCoefVar.loc[gene]['TissueBin']
    VS=MetabolicVS.loc[gene]['DevelopmentBin']
    if VS=='Highly variant' and CV=='Highly variant':
        MetabolicCoefVar.at[gene,'Quadrant']=4
    elif VS=='Highly variant' and CV!='Highly variant':
        MetabolicCoefVar.at[gene,'Quadrant']=2
    elif VS!='Highly variant' and CV=='Highly variant':
        MetabolicCoefVar.at[gene,'Quadrant']=3
    else:
        MetabolicCoefVar.at[gene,'Quadrant']=1


# In[146]:


# MetabolicCoefVar ['TissueBin']=CoefVar['Tissue bin']


# In[147]:


MetabolicCoefVar


# In[148]:


Lowlyexpressed=MetabolicCoefVar[MetabolicCoefVar.TissueBin=='Lowly expressed']


# In[149]:


# MetabolicCoefVar.TissueBin)


# In[150]:


# MetabolicCoefVar['TissueBin']=CoefVar['Bin']


# In[151]:


MetabolicCoefVar


# In[152]:


for index in MetabolicCoefVar.index:
    if index in Lowlyexpressed.index:
        MetabolicCoefVar.at[index,'TissueBin']="Lowlyexpressed"


# In[153]:


MetabolicCoefVar['VS']=MetabolicVS['Variation Score']


# In[168]:


MetabolicCoefVar['DevelopmentBin']=MetabolicVS['DevelopmentBin']


# In[169]:


MetabolicCoefVar


# ## Dividing metabolic genes into quadrants

# In[156]:


## Common lowly expressed


# In[157]:


LowExpKim


# In[158]:


MetabolicLowExpKim=list(set(MetabolicClasses.index).intersection(set(LowExpKim.index)))


# In[159]:


MetabolicLowExpKim


# In[161]:


MetabolicCoefVar['DevelopmentBin']


# In[71]:


MetabolicLowExpKim


# In[166]:


MetabolicCoefVar['DevelopmentBin']=str(MetabolicCoefVar['DevelopmentBin'])


# In[167]:


MetabolicCoefVar


# In[72]:


# (MetabolicCoefVar.loc[genes,'DevelopmentBin'])==2


# In[170]:


# MetabolicCoefVar['DevelopmentBin']=str(MetabolicCoefVar['DevelopmentBin'])
for genes in MetabolicCoefVar.index:
    print(genes)
    if genes in MetabolicLowExpKim:
        print(genes)
        MetabolicCoefVar.at[genes,'Development Bin']='Lowlyexpressed'
        print(MetabolicCoefVar.loc[genes,'Development Bin'])
    elif (MetabolicCoefVar.loc[genes,'DevelopmentBin'])==1:
        MetabolicCoefVar.at[genes,'Development Bin']='Invariant'
    elif (MetabolicCoefVar.loc[genes,'DevelopmentBin'])==3:
        MetabolicCoefVar.at[genes,'Development Bin']='Highly variant'
    elif (MetabolicCoefVar.loc[genes,'DevelopmentBin'])==2:
        MetabolicCoefVar.at[genes,'Development Bin']='Moderately variant'
        


# In[171]:


MetabolicCoefVar


# In[172]:


MetabolicCoefVar.drop(columns=['DevelopmentBin'],inplace=True)


# In[173]:


MetabolicCoefVar.columns


# ## Removing genes that are lowly expressed in both development and tissues

# In[174]:


count=0
for genes in MetabolicCoefVar.index:
    if (MetabolicCoefVar.loc[genes,'Development Bin'])=='Lowlyexpressed' :
        if  (MetabolicCoefVar.loc[genes,'TissueBin'])=='Lowlyexpressed':
            print(genes)
            count=count+1
            print(count)
            MetabolicCoefVar.drop(index=genes,inplace=True)
        


# In[175]:


MetabolicCoefVar


# In[176]:


from scipy.stats import pearsonr
from scipy.stats import spearmanr
FilteredMetabolicVS=MetabolicCoefVar.dropna(inplace=False)
spearmanr((FilteredMetabolicVS['VS']),(FilteredMetabolicVS['CoefVar']))


# In[177]:


pearsonr((FilteredMetabolicVS['VS']),(FilteredMetabolicVS['CoefVar']))


# ## Quadrant distinction

# In[178]:


for gene in MetabolicCoefVar.index:
    CV=MetabolicCoefVar.loc[gene]['TissueBin']
    VS=MetabolicCoefVar.loc[gene]['Development Bin']
    if VS=='Highly variant' and CV=='Highly variant':
        MetabolicCoefVar.at[gene,'Quadrant']=4
    elif VS=='Highly variant' and CV!='Highly variant':
        MetabolicCoefVar.at[gene,'Quadrant']=2
    elif VS!='Highly variant' and CV=='Highly variant':
        MetabolicCoefVar.at[gene,'Quadrant']=3
    else:
        MetabolicCoefVar.at[gene,'Quadrant']=1


# In[107]:


# type(MetabolicCoefVar.loc['WBGene00009981']['Development Bin'])


# In[108]:


FilteredMetabolicVS=(MetabolicCoefVar[['Development Bin','TissueBin']]).dropna()


# In[109]:


FilteredMetabolicVS


# In[ ]:


# FilteredMetabolicVS=MetabolicCoefVar[MetabolicCoefVar['Development Bin']!=np.nan]
# FilteredMetabolicVS=FilteredMetabolicVS[FilteredMetabolicVS['Development Bin']!=np.NaN]
# FilteredMetabolicVS=FilteredMetabolicVS[FilteredMetabolicVS['Development Bin']!='nan']
# FilteredMetabolicVS=FilteredMetabolicVS[FilteredMetabolicVS['TissueBin']!=np.NaN]
# FilteredMetabolicVS=FilteredMetabolicVS[FilteredMetabolicVS['TissueBin']!=np.nan]
# FilteredMetabolicVS=FilteredMetabolicVS[FilteredMetabolicVS['TissueBin']!='nan']


# In[110]:


for gene in FilteredMetabolicVS.index:
    CV=FilteredMetabolicVS.loc[gene]['TissueBin']
    VS=FilteredMetabolicVS.loc[gene]['Development Bin']
    if VS=='Highly variant' and CV=='Highly variant':
        FilteredMetabolicVS.at[gene,'Quadrant']=4
    elif VS=='Highly variant' and CV!='Highly variant':
        FilteredMetabolicVS.at[gene,'Quadrant']=2
    elif VS!='Highly variant' and CV=='Highly variant':
        FilteredMetabolicVS.at[gene,'Quadrant']=3
    else:
        FilteredMetabolicVS.at[gene,'Quadrant']=1


# In[111]:


FilteredMetabolicVS


# In[112]:


Quadrant4=FilteredMetabolicVS[FilteredMetabolicVS.Quadrant==4]
Quadrant3=FilteredMetabolicVS[FilteredMetabolicVS.Quadrant==3]
Quadrant2=FilteredMetabolicVS[FilteredMetabolicVS.Quadrant==2]
Quadrant1=FilteredMetabolicVS[FilteredMetabolicVS.Quadrant==1]


# In[113]:


Quadrant1.to_csv("Quadrant1.csv")
Quadrant2.to_csv("Quadrant2.csv")
Quadrant3.to_csv("Quadrant3.csv")
Quadrant4.to_csv("Quadrant4.csv")


# ## Dividing only iCEL1314 genes into quadrants

# In[114]:


ClassA=MetabolicClasses[MetabolicClasses.Class=='A']
classAQ1=list(set(ClassA.index).intersection(set(Quadrant1.index)))
ClassAQuadrant1=Quadrant1.loc[classAQ1]
ClassAQuadrant1.to_csv("ClassAQuadrant1.csv")
classAQ2=list(set(ClassA.index).intersection(set(Quadrant2.index)))
ClassAQuadrant2=Quadrant2.loc[classAQ2]
ClassAQuadrant2.to_csv("ClassAQuadrant2.csv")
classAQ3=list(set(ClassA.index).intersection(set(Quadrant3.index)))
ClassAQuadrant3=Quadrant3.loc[classAQ3]
ClassAQuadrant3.to_csv("ClassAQuadrant3.csv")
classAQ4=list(set(ClassA.index).intersection(set(Quadrant4.index)))
ClassAQuadrant4=Quadrant4.loc[classAQ4]
ClassAQuadrant4.to_csv("ClassAQuadrant4.csv")


# In[115]:


MetabolicCoefVar['Quadrant']=FilteredMetabolicVS.Quadrant
MetabolicCoefVar.to_csv("QuadrantVSCV.csv")


# In[116]:




RegulatedGenes=pd.read_csv("Differentiallyenriched_Regulated_p0.01.csv",index_col=0)


# In[117]:


RegulatedGenes.columns


# ## Total Class A Regulated genes patwhay enrichment

# In[118]:


EnrichedWormPathsClassARegulated=pd.read_csv("TotalClassAMetabolic_PEAStats.txt",sep='\t',index_col=1)


# In[119]:


EnrichedWormPathsClassARegulated=EnrichedWormPathsClassARegulated[~EnrichedWormPathsClassARegulated.index.duplicated(keep='first')]


# In[122]:


import statsmodels as st
# import statsmodels.api as sm
pvalse=EnrichedWormPathsClassARegulated.p_enrichment
EnrichedWormPathsClassARegulated['FDR-corrected p-enrichment']=(st.stats.multitest.fdrcorrection(pvalse, alpha=0.05, method='indep', 
                                                                                           is_sorted=False))[1]


# In[123]:


SignificantPathwayEnrichment= EnrichedWormPathsClassARegulated[EnrichedWormPathsClassARegulated['p_enrichment']<=0.3]


# In[124]:


SignificantPathwayEnrichment.index


# In[125]:


SignificantPathwayEnrichment=SignificantPathwayEnrichment.loc[['LIPIDS','STEROID METABOLISM','BIOSYNTHESIS OF BILE ACID-LIKE MOLECULES',
       'PEROXISOMAL FATTY ACID DEGRADATION','FATTY ACID DEGRADATION',
       'SPHINGOLIPID METABOLISM','GLYCEROPHOSPHOLIPID METABOLISM','ETHER LIPID METABOLISM','FATTY ACID BIOSYNTHESIS OTHER',
                                                               'FATTY ACID BIOSYNTHESIS','MITOCHONDRIAL FATTY ACID DEGRADATION']]


# In[5]:


def SignificantPathwayEnrichmentF(ClassAQuadrant1Pathway):
    PathwayEnrichmentRegulated1=ClassAQuadrant1Pathway[~ClassAQuadrant1Pathway.index.duplicated(keep='first')]
    pvalse=PathwayEnrichmentRegulated1.p_enrichment
    PathwayEnrichmentRegulated1['FDR-corrected p-enrichment']=(st.stats.multitest.fdrcorrection(pvalse, alpha=0.05, method='indep', 
                                                                                           is_sorted=False))[1]
    SignificantPathwayEnrichment1=PathwayEnrichmentRegulated1[['Enrichment score (n_Hits/n_Genes)',
                                                           'p_enrichment','FDR-corrected p-enrichment']]
    SignificantPathwayEnrichment1['-log10(p-enrichment)']=-(np.log10(SignificantPathwayEnrichment1['p_enrichment']))
    
    return SignificantPathwayEnrichment1

def PlotSignificantPathwayEnrichment(SignificantPathwayEnrichment1,title,cutoff):
    import seaborn as sns
    SignificantPathwayEnrichment1=SignificantPathwayEnrichment1[(SignificantPathwayEnrichment1['p_enrichment']<=cutoff)]
    g = SignificantPathwayEnrichment1.reset_index()
    survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
    n = g['-log10(p-enrichment)']

    norm = plt.Normalize(0,1)
    sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
    sm.set_array([])

    ax = sns.barplot(x='-log10(p-enrichment)', y='Category', 
                 hue='Enrichment score (n_Hits/n_Genes)', palette='Oranges', 
                 dodge=False,data=g)

    ax.set_ylabel('WormPaths Pathway/Category ')
    plt.xlim(0,16)

    ax.get_legend().remove()
    from matplotlib import rcParams
    rcParams['font.family'] = 'Arial'
    plt.tight_layout()
    ax.figure.colorbar(sm)
    plt.savefig("WormPathsEnrichment_{}.svg".format(title),dpi=300)


# In[6]:


SignificantPathwayEnrichment


# In[128]:


SignificantPathwayEnrichment=SignificantPathwayEnrichmentF(SignificantPathwayEnrichment)


# In[129]:


PlotSignificantPathwayEnrichment(SignificantPathwayEnrichment,title='ClassAOverallRegulated',cutoff=0.05)


# ## Pathways enriched in 25% of the regulated

# In[160]:


EnrichedPathwaysModeratelyVariant=pd.read_csv("ModeratelyVariant_PathwayEnrich.csv",index_col=0)


# In[134]:


Enri


# In[156]:


EnrichedPathwaysModeratelyVariant=EnrichedPathwaysModeratelyVariant[EnrichedPathwaysModeratelyVariant.index=='LEVEL 4']


# In[145]:


EnrichedPathwaysModeratelyVariant.sort_values(ascending=False,inplace=True,by=['-log10(p_enrichment)'])


# In[161]:


EnrichedPathwaysModeratelyVariant=EnrichedPathwaysModeratelyVariant[EnrichedPathwaysModeratelyVariant.p_enrichment<=0.05]


# In[155]:


EnrichedPathwaysModeratelyVariant


# In[162]:


fig=plt.figure(figsize=(10,10))
import seaborn as sns
g = EnrichedPathwaysModeratelyVariant.reset_index()
survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
n = g['-log10(p_enrichment)']
# cMAP=sns.diverging_palette(220, 20, as_cmap=True)

norm = plt.Normalize(g['Enrichment score (n_Hits/n_Genes)'].min(), g['Enrichment score (n_Hits/n_Genes)'].max())
sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(p_enrichment)', y='Category', 
                 hue='Enrichment score (n_Hits/n_Genes)', palette="Oranges", 
                 dodge=False,data=g)

ax.set_ylabel('Category')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("ModeratelyVariantPathwayEnrichment.svg",dpi=300)


# ## Phenotypes

# In[ ]:


EnrichedPhenotypesModeratelyvariant=pd.read_csv("ModeratelyVariantGenesPhenotypeEnrichment.csv",index_col=0)


# In[ ]:


EnrichedPhenotypesModeratelyvariant.columns


# In[ ]:


EnrichedPhenotypesModeratelyvariant.sort_values(ascending=False,inplace=True,by=['-log10(BH-FDR Corrected FDR)'])


# In[ ]:


fig=plt.figure(figsize=(10,10))
import seaborn as sns
g = EnrichedPhenotypesModeratelyvariant.reset_index()
survival_rates = g['Enrichment fold change'].mean()
n = g['-log10(BH-FDR Corrected FDR)']
# cMAP=sns.diverging_palette(220, 20, as_cmap=True)

norm = plt.Normalize(g['Enrichment fold change'].min(), g['Enrichment fold change'].max())
sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH-FDR Corrected FDR)', y='Phenotype Name', 
                 hue='Enrichment fold change', palette="Oranges", 
                 dodge=False,data=g)

ax.set_ylabel('Phenotype Name')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("ModeratelyVariantPhenotypeEnrichment.svg",dpi=300)


# In[ ]:


DepletedPhenotypes=RegulatedGenes[RegulatedGenes['Enrichment fold change']<1]


# In[ ]:


EnrichedPhenotypes=RegulatedGenes[RegulatedGenes['Enrichment fold change']>1]


# In[ ]:


EnrichedPhenotypes.columns


# In[ ]:


EnrichedPhenotypes.sort_values(ascending=False,inplace=True,by=['-log10(BH FDR-corrected p-value)'])


# In[ ]:


DepletedPhenotypes.sort_values(by=['-log10(BH FDR-corrected p-value)'],ascending=False)


# In[ ]:


EnrichedDepletedPhenotypes=pd.concat([EnrichedPhenotypes,DepletedPhenotypes])


# In[ ]:


EnrichedDepletedPhenotypes.sort_values(ascending=False,inplace=True,by=['-log10(BH FDR-corrected p-value)'])


# In[ ]:


EnrichedDepletedPhenotypes


# In[ ]:


fig=plt.figure(figsize=(10,10))
import seaborn as sns
g = EnrichedPhenotypes.reset_index()
survival_rates = g['Enrichment fold change'].mean()
n = g['-log10(BH FDR-corrected p-value)']
# cMAP=sns.diverging_palette(220, 20, as_cmap=True)

norm = plt.Normalize(g['Enrichment fold change'].min(), g['Enrichment fold change'].max())
sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR-corrected p-value)', y='Phenotype Name', 
                 hue='Enrichment fold change', palette="Oranges", 
                 dodge=False,data=g)

ax.set_ylabel('Phenotype Name')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormBasePhenotypeEnrichment.svg",dpi=300)


# In[ ]:


DepletedPhenotypes.sort_values(ascending=False,inplace=True,by=['-log10(BH FDR-corrected p-value)'])


# In[ ]:


DepletedPhenotypes['Enrichment fold change']=1-DepletedPhenotypes['Enrichment fold change']


# In[ ]:


fig=plt.figure(figsize=(10,10))
import seaborn as sns
g = DepletedPhenotypes.reset_index()
survival_rates = g['Enrichment fold change'].mean()
n = g['-log10(BH FDR-corrected p-value)']
# cMAP=sns.diverging_palette(220, 20, as_cmap=True)

norm = plt.Normalize(g['Enrichment fold change'].min(), g['Enrichment fold change'].max())
sm = plt.cm.ScalarMappable(cmap="Blues", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR-corrected p-value)', y='Phenotype Name', 
                 hue='Enrichment fold change', palette="Blues", 
                 dodge=False,data=g)

ax.set_ylabel('Phenotype Name')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormBasePhenotypeDepleted.svg",dpi=300)


# In[ ]:


DepletedPhenotypes['-log10(BH FDR-corrected p-value)']=-1*(DepletedPhenotypes['-log10(BH FDR-corrected p-value)'])


# In[ ]:


fig=plt.figure(figsize=(10,10))
import seaborn as sns
g = EnrichedDepletedPhenotypes.reset_index()
survival_rates = g['Enrichment fold change'].mean()
n = g['-log10(BH FDR-corrected p-value)']
# cMAP=sns.diverging_palette(220, 20, as_cmap=True)

norm = plt.Normalize(g['Enrichment fold change'].min(), g['Enrichment fold change'].max())
sm = plt.cm.ScalarMappable(cmap="vlag", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(BH FDR-corrected p-value)', y='Phenotype Name', 
                 hue='Enrichment fold change', palette="vlag", 
                 dodge=False,data=g)

ax.set_ylabel('Phenotype Name')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormBasePhenotypeEnrichmentDepletion.svg",dpi=300)


# ## Quadrant-wise enrichment analysis

# ### WormPaths

# In[7]:


get_ipython().system('ls')


# In[13]:


ClassAQuadrant1Pathway=pd.read_excel("ClassAQuadrant1PEA.xlsx",index_col=1)


# In[14]:


ClassAQuadrant1Pathway


# In[15]:


ClassAQuadrant2Pathway=pd.read_excel("ClassAQuadrant2PEA.xlsx",index_col=1)
ClassAQuadrant3Pathway=pd.read_excel("ClassAQuadrant3PEA.xlsx",index_col=1)
ClassAQuadrant4Pathway=pd.read_excel("ClassAQuadrant4PEA.xlsx",index_col=1)


# In[16]:


ClassAQuadrant1Pathway


# In[25]:


def SignificantPathwayEnrichment(ClassAQuadrant1Pathway):
    PathwayEnrichmentRegulated1=ClassAQuadrant1Pathway[~ClassAQuadrant1Pathway.index.duplicated(keep='first')]
    pvalse=PathwayEnrichmentRegulated1.p_enrichment
    PathwayEnrichmentRegulated1['FDR-corrected p-enrichment']=(st.stats.multitest.fdrcorrection(pvalse, alpha=0.05, method='indep', 
                                                                                           is_sorted=False))[1]
    SignificantPathwayEnrichment1=PathwayEnrichmentRegulated1[['Enrichment score (n_Hits/n_Genes)',
                                                           'p_enrichment','FDR-corrected p-enrichment']]
    SignificantPathwayEnrichment1['-log10(p-enrichment)']=-(np.log10(SignificantPathwayEnrichment1['p_enrichment']))
    
    return SignificantPathwayEnrichment1

def PlotSignificantPathwayEnrichment(SignificantPathwayEnrichment1,title):
    import seaborn as sns
    SignificantPathwayEnrichment1=SignificantPathwayEnrichment1[(SignificantPathwayEnrichment1['p_enrichment']<=0.05)]
    g = SignificantPathwayEnrichment1.reset_index()
    survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
    n = g['-log10(p-enrichment)']

    norm = plt.Normalize(0,1)
    sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
    sm.set_array([])

    ax = sns.barplot(x='-log10(p-enrichment)', y='Category', 
                 hue='Enrichment score (n_Hits/n_Genes)', palette='Oranges', 
                 dodge=False,data=g)

    ax.set_ylabel('WormPaths Pathway/Category ')
#     plt.xlim(0,16)

    ax.get_legend().remove()
    from matplotlib import rcParams
    rcParams['font.family'] = 'Arial'
    plt.tight_layout()
    ax.figure.colorbar(sm)
    plt.savefig("WormPathsEnrichment_{}.svg".format(title),dpi=300)

                                                               


# In[26]:


SignificantPathwayEnrichment1=SignificantPathwayEnrichment(ClassAQuadrant1Pathway)
PlotSignificantPathwayEnrichment(SignificantPathwayEnrichment1,title='ClassAQuadrant1_0.01')


# In[27]:


SignificantPathwayEnrichment2=SignificantPathwayEnrichment(ClassAQuadrant2Pathway)
PlotSignificantPathwayEnrichment(SignificantPathwayEnrichment2,title='ClassAQuadrant2_0.01')


# In[ ]:


SignificantPathwayEnrichment2


# In[28]:


PlotSignificantPathwayEnrichment(SignificantPathwayEnrichment1,title='ClassAQuadrant1_0.01')


# In[29]:


SignificantPathwayEnrichment3=SignificantPathwayEnrichment(ClassAQuadrant3Pathway)
PlotSignificantPathwayEnrichment(SignificantPathwayEnrichment3,title='ClassAQuadrant3_0.03')


# In[30]:


SignificantPathwayEnrichment4=SignificantPathwayEnrichment(ClassAQuadrant4Pathway)
PlotSignificantPathwayEnrichment(SignificantPathwayEnrichment4,title='ClassAQuadrant4_0.03')


# In[ ]:


SignificantPathwayEnrichment4[0:50]


# In[ ]:


SignificantPathwayEnrichment3[]


# ## Regulated genes

# In[ ]:


MetabolicRegulatedGenes=pd.read_csv("RegulatedMetabolic_KimTissue_AllDatasets.csv",index_col=0)


# In[ ]:


MetabolicRegulatedGenes


# ## Checking which regulated metabolic genes have stress-response variants

# In[ ]:


StressResponse=pd.read_csv("StressRespone_genes_direct_and_inferred_for_WBPhenotype_0000067.txt",index_col=0,
                           sep='\t',header=None)


# In[ ]:


MetabolicStressResponse=list(set(StressResponse.index).intersection(set(MetabolicRegulatedGenes.index)))


# In[ ]:


MetabolicStressResponse=MetabolicRegulatedGenes.loc[MetabolicStressResponse]


# In[ ]:


MetabolicStressResponse


# In[ ]:


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


# In[ ]:


MetabolicStressResponse=wb_to_gene(MetabolicStressResponse)


# In[ ]:


MetabolicStressResponse=gene_to_wb(MetabolicStressResponse)


# ## Reading CV of all genes across all datasets

# In[ ]:


CV_alldatasets=pd.read_csv("/data/nandas/Transcription/AllDatasets_012422/CV_allgenes.csv",index_col=0)


# In[ ]:


CV_alldatasets


# In[ ]:


MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",index_col=0)
metabolicCV=list(set(MetabolicClasses.index).intersection(set(CV_alldatasets.index)))


# In[ ]:


metabolicCV


# In[ ]:


MetabolicCV=CV_alldatasets.loc[metabolicCV]


# In[ ]:


# MetabolicCV=MetabolicCV['Max_CV']


# In[ ]:


metabolicCVall=list(set(MetabolicRegulatedGenes.index).intersection(set(CV_alldatasets.index)))
MetabolicRegulated_CV_alldatasets=CV_alldatasets.loc[metabolicCVall]


# In[ ]:


MetabolicRegulated_CV_alldatasets.replace(np.NaN,-6,inplace=True)


# In[ ]:


MetabolicRegulated_CV_alldatasets


# In[ ]:


data=MetabolicRegulated_CV_alldatasets.transpose()
data=wb_to_gene(data)


# In[ ]:


data=data.loc['Max_CV']


# In[ ]:


data=pd.DataFrame(data)


# In[ ]:


data.sort_values(ascending=False,by=['Max_CV'],inplace=True)


# In[ ]:


data[data.Max_CV>=0.75]


# In[ ]:


data


# ## Adding variation score and Tissue CV to list of regulated metabolic genes
# 

# In[ ]:


VS=pd.read_csv(Kim_dir+"BinningVS.csv",index_col=0)


# In[ ]:


CoefVar=pd.read_csv(Tissue_dir+"CoefVarHighModTissue.csv",
                    index_col=0)


# In[ ]:


MetabolicCV['Variation Score']=VS['Variation Score']


# In[ ]:


MetabolicCV['CV across tissues']=CoefVar['CoefVar']


# In[ ]:


MetabolicCV=MetabolicCV[['Max_CV','Variation Score','CV across tissues']]


# In[ ]:


df=MetabolicCV
normalized_df=(df-df.mean())/df.std()
normalized_df.replace(np.NaN,-6,inplace=True)
normalized_df.rename(columns={'Max_CV':'Max CV across datasets'},inplace=True)


# In[ ]:


normalized_df


# In[ ]:


# normalized_df=normalized_df.transpose()
mrg=list(set(MetabolicRegulatedGenes.index).intersection(set(normalized_df.index)))
MetabolicRegulatedGenes_CV=normalized_df.loc[mrg]


# In[ ]:


sns.clustermap(MetabolicRegulatedGenes_CV,mask=MetabolicRegulatedGenes_CV==-6,
               vmin=-2,vmax=2,cmap='vlag',figsize=(5,30))
plt.tight_layout()


# In[ ]:


# Invariantgenes=normalized_df.loc[InvariantCV.index]


# In[ ]:


Ivg=list(set())


# In[ ]:


normalized_df=normalized_df.transpose()


# In[ ]:


normalized_df.min().min()


# In[ ]:


MetabolicRegulatedGenes


# In[ ]:


sns.clustermap(normalized_df,mask=normalized_df==-6,
               vmin=-2,vmax=2,cmap='vlag',row_cluster=False,figsize=(30,5))
plt.tight_layout()


# In[ ]:


VS


# In[ ]:


CoefVar


# In[ ]:


CoefVar.columns


# In[ ]:


data=gene_to_wb(data)


# In[ ]:


data['Variation Score']=VS['Variation Score']


# In[ ]:


data['CV across tissues']=CoefVar['CoefVar']


# In[ ]:


data.rename(columns={'Max_CV':'Max CV across datasets'},inplace=True)


# In[ ]:


data.replace(np.NaN,-6,inplace=True)


# In[ ]:


data=data.transpose()


# In[ ]:


data=wb_to_gene(data)


# In[ ]:


sns.clustermap(data,mask=data==-6,
               vmin=0,vmax=2,cmap='YlOrRd',row_cluster=False,figsize=(50,10))
plt.tight_layout()


# ## Plotting invariant genes

# In[ ]:


Invariantgenes=pd.read_csv("/data/nandas/Transcription/AllDatasets_012422/Bin1_AllDatasets.csv",index_col=0)


# In[ ]:


invariantcv=list(set(Invariantgenes.index).intersection(set(CV_alldatasets.index)))


# In[ ]:


InvariantCV=CV_alldatasets.loc[invariantcv]


# In[ ]:


InvariantCV.replace(np.NaN,-6,inplace=True)


# In[ ]:


data2=InvariantCV.transpose()
data2=wb_to_gene(data2)


# In[ ]:


sns.clustermap(data2,mask=data2==-6,
               vmin=0,vmax=0.75,cmap='YlOrRd')
plt.tight_layout()


# In[ ]:


data2.index


# ## Lipid enrichment
# 

# In[ ]:


PathwayEnrichmentRegulated=pd.read_csv("Lipid_metabolism_PEA.xlsx - lipid.csv",index_col=0)


# In[ ]:


PathwayEnrichmentRegulated


# In[ ]:


import statsmodels as st
import statsmodels.api as sm
PathwayEnrichmentRegulated=pd.read_csv("Lipid_metabolism_PEA.xlsx - lipid.csv",index_col=0)
PathwayEnrichmentRegulated.set_index(['Category'],inplace=True)
PathwayEnrichmentRegulated=PathwayEnrichmentRegulated[~PathwayEnrichmentRegulated.index.duplicated(keep='first')]
pvalse=PathwayEnrichmentRegulated.p_enrichment
PathwayEnrichmentRegulated['FDR-corrected p-enrichment']=(st.stats.multitest.fdrcorrection(pvalse, alpha=0.05, method='indep', 
                                                                                           is_sorted=False))[1]
SignificantPathwayEnrichment=PathwayEnrichmentRegulated[['Enrichment score (n_Hits/n_Genes)',
                                                           'p_enrichment']]
SignificantPathwayEnrichment['-log10(p-enrichment)']=-(np.log10(SignificantPathwayEnrichment['p_enrichment']))


# In[ ]:


SignificantPathwayEnrichment=SignificantPathwayEnrichment[SignificantPathwayEnrichment.p_enrichment<=0.25]


# In[ ]:


import seaborn as sns
g = SignificantPathwayEnrichment.reset_index()
survival_rates = g['Enrichment score (n_Hits/n_Genes)'].mean()
n = g['-log10(p-enrichment)']

norm = plt.Normalize(g['Enrichment score (n_Hits/n_Genes)'].min(), g['Enrichment score (n_Hits/n_Genes)'].max())
sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
sm.set_array([])

ax = sns.barplot(x='-log10(p-enrichment)', y='Category', 
                 hue='Enrichment score (n_Hits/n_Genes)', palette='Oranges', 
                 dodge=False,data=g)

ax.set_ylabel('WormPaths Pathway/Category ')

ax.get_legend().remove()
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
plt.tight_layout()
ax.figure.colorbar(sm)
plt.savefig("WormPathsEnrichment_lipids_alldatasets.svg",dpi=300)


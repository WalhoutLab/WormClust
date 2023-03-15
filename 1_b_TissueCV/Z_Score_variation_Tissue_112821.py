#!/usr/bin/env python
# coding: utf-8

# ## Importing modules

# In[226]:


import pandas as pd
# import gseapy as gp
import matplotlib.pyplot as plt
# from gseapy.parser import Biomart
import os
import numpy as np
import seaborn as sns
# from gseapy.plot import gseaplot
from scipy import stats as st
import fnmatch
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# ## Setting base directory

# In[227]:


Base_dir='/data/nandas/Combined_coexp/Part_1_TranscriptionallyRegulatedGenes/Tissue/'
os.chdir(Base_dir)


# In[228]:


mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=0)


# In[229]:


mapper_df


# ## Functions

# In[230]:


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

def CalculateZScore(df):
    from scipy import stats
    Zscore=stats.zscore(df.values,axis=1,nan_policy='omit')
    ZScore_df=pd.DataFrame(Zscore,index=df.index,columns=df.columns)
    return ZScore_df

## Calculate coefficient of variation
def CalculateCoefVar(ZScore_df):
    Std=np.nanstd(ZScore_df.values,axis=1)
    Mean=ZScore_df.mean(axis=1,skipna=True)
    CoefVar=pd.DataFrame([])
#     from scipy.stats import variation 
#     Variation=variation(ZScore_df.values, axis = 1)
    CoefVar['CoefVar']=(Std/Mean)
    CoefVar.index=ZScore_df.index
    return CoefVar

def Log2Transform(df,title):
    #df=df.drop(columns=['GENEID','GENENAME'])
    df=df.stack()
    df=df.loc[~(df==0)]
    df=pd.DataFrame(df)
    df=df.applymap(np.log2)
    ## Plot the values of all cells in the dataframe to understand the distribution
    hist=df.hist(grid=False,bins=100,color='skyblue')
#     plt.axvline(1.3,color='r')
#     plt.axvline(df.mode()[1],color='g')
    plt.xlabel("log2(TPM)")
    plt.ylabel("Frequency")
    plt.title("{}".format(title))
    plt.savefig("{}.svg".format(title))
    return df,hist

def PlotCoefVar(CoefVar):
    fig, ax = plt.subplots(figsize=(7,5))
    AllGenes=CoefVar.CoefVar.hist(ax=ax,color='green',bins=100,label='All genes',alpha=0.5)
    MetabolicGenes=CoefVar.loc[intersect2].hist(ax=ax,color='midnightblue',bins=100,label='Metabolic genes',alpha=0.5)  
    ax.grid(False)
    plt.title("Distribution of Coefficient of variation")
    ax.set_xlabel('Coefficient of Variation')
    ax.set_ylabel('Number of genes')
#     for xc in xcoords:
#         ax.axvline(xc,color='blue',linestyle='--')
    plt.legend(loc='best')
    plt.savefig("CoefVar.svg", dpi=300)


# ## Reading correlation files

# In[231]:


Tissue_data=pd.read_csv("/data/nandas/Combined_coexp/Pathway_enrichment/Gmean/Tissue7/normalized_7tissue_metabolic.dat",
                        sep='\t')


# In[232]:


Tissue_exp=pd.read_csv("Tissue_WB_7_exprsssion.csv",index_col=0)


# In[233]:


Tissue_exp[~Tissue_exp.index.str.contains('WB')]


# In[234]:


# Tissue_exp=gene_to_wb(Tissue_exp)  


# In[235]:


Tissue_cat=pd.read_csv("Tcat_tissues_allgenes_08122020.tsv",sep='\t',header='infer',index_col=0)


# In[236]:


# Tissue_exp=SeqToWB(Tissue_exp)


# In[237]:


# Tissue_cat=wb_to_gene(Tissue_cat)
# Tissue_exp=wb_to_gene(Tissue_exp)
# Tissue_data=wb_to_gene(Tissue_data)


# ## Filtering for only live, protein coding genes

# In[238]:


MasterList=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv",index_col=0)
MasterList.set_index(['WormBaseID'],inplace=True)


# In[239]:


MasterList=MasterList[MasterList.Type=='protein_coding_gene']
MasterList=MasterList[MasterList.Status=='Live']


# In[240]:


MasterList[~MasterList.index.str.contains('WB')]


# In[241]:


masterlist=list(set(MasterList.index).intersection(set(Tissue_exp.index)))
x=Tissue_exp.loc[masterlist]


# In[242]:


diff=list(set(Tissue_exp.index).difference(set(MasterList.index)))


# In[243]:


len(diff)


# In[244]:



x.to_csv("TotalGenesTissue.csv")


# In[245]:


x


# In[246]:


## Filtering for live protein-coding genes
Tissue_cat=Tissue_cat.loc[x.index]


# ##### Tissue_exp=Tissue_exp.loc[x.index]

# In[247]:


Tissue_cat


# In[248]:


# Tissue_cat=SeqToGene(Tissue_cat)
# Tissue_exp=SeqToGene(Tissue_exp)
# Tissue_data=SeqToGene(Tissue_data)


# In[249]:


# Tissue_cat=gene_to_wb(Tissue_cat)


# In[250]:


HighExp_genes=Tissue_cat[(Tissue_cat.values=='High').sum(axis=1)>=1]


# In[251]:


ModerateExp_genes=Tissue_cat[(Tissue_cat.values=='Moderate').sum(axis=1)>=1]


# In[252]:


Tissue_exp


# In[253]:


set(HighExp_genes.index.str.strip())


# ## Finding genes that show either high or moderate expression in at least one tissue

# In[254]:


HighMod= list(set(HighExp_genes.index.str.strip()).union(set(ModerateExp_genes.index.str.strip())))


# In[255]:


# HighMod=Tissue_exp.loc[HighMod]
HighMod=[x for x in HighMod if str(x) != 'nan']


# In[256]:


HighMod=Tissue_exp.loc[HighMod]


# In[257]:


HighMod=HighMod[~HighMod.index.duplicated(keep='first')]


# In[258]:


HighExp_genes


# In[259]:


leg=list(set(x.index).difference(set(HighMod.index)))
LowExpGenes=Tissue_exp.loc[leg]


# In[260]:


LowExpGenes.to_csv("LowExpGenes.csv")


# In[261]:


LowExpGenes


# In[262]:


# HighMod=SeqToGene(HighMod)


# In[263]:


HighMod=HighMod[~HighMod.index.duplicated(keep='first')]


# In[264]:


HighMod.to_csv("Tissue_HighMod.csv")


# In[265]:


HighMod.index=HighMod.index.str.strip()


# In[266]:


MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",index_col=0)


# In[267]:


MetabolicClasses=MetabolicClasses[~MetabolicClasses.index.duplicated(keep='first')]


# In[ ]:





# In[268]:


# MetabolicClasses=gene_to_wb(MetabolicClasses)


# In[269]:


# MetabolicClasses.to_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv")


# In[270]:


# MetabolicClasses['WormBaseID']='WORMBASE'


# In[271]:


# MetabolicClasses=gene_to_wb(MetabolicClasses)


# In[272]:


# MetabolicClasses.reset_index(inplace=True)
# MetabolicClasses.set_index(['Sequence Name'],inplace=True)


# In[273]:


# MetabolicClasses=SeqToWB(MetabolicClasses)


# In[274]:


# MetabolicClasses.reset_index(inplace=True)


# In[275]:


# MetabolicClasses


# In[276]:


# # MetabolicClasses.reset_index(inplace=True)
# # df.reset_index(inplace=True)
# for i in range(MetabolicClasses.shape[0]):
#     if "WB" in str(MetabolicClasses.iloc[i]['Gene Name']):
#          if "\n" not in str(MetabolicClasses.iloc[i]['Gene Name']):
#             print ("Checking in index: {}:{}".format(i,MetabolicClasses.iloc[i]['Gene Name']))
#             MetabolicClasses.at[i,'WormBaseID']=MetabolicClasses.iloc[i]['Gene Name']
#     if "WB" in str(MetabolicClasses.iloc[i]['Sequence Name']):
#         if "\n" not in str(MetabolicClasses.iloc[i]['Sequence Name']):
#             print ("Checking in Sequence name: {}:{}".format(i,MetabolicClasses.iloc[i]['Sequence Name']))
#             MetabolicClasses.at[i,'WormBaseID']=MetabolicClasses.iloc[i]['Sequence Name']


# In[277]:


# MetabolicClasses.set_index(['WormBaseID'],inplace=True)


# In[278]:


# MetabolicClasses.drop(columns=['Sequence Name','Gene Name'],inplace=True)


# In[279]:


MetabolicClasses.loc[['WBGene00000481']]


# In[280]:


MetabolicClasses=MetabolicClasses[~MetabolicClasses.index.duplicated(keep='first')]


# In[281]:


Tissue_exp.index=Tissue_exp.index.str.strip()


# In[282]:


x[x.index.str.contains('WB')]


# In[283]:


MetabolicClasses.index=MetabolicClasses.index.str.strip()


# In[284]:


MetabolicClasses[~MetabolicClasses.index.str.contains('WB')]


# In[285]:


diff=list(set(MetabolicClasses.index).difference(set(Tissue_exp.index)))


# In[286]:


# HighMod=wb_to_gene(HighMod)
# HighMod=SeqToGene(HighMod)


# In[287]:


len(diff)


# ## Extracting metabolic genes

# In[288]:


intersect=list(set(x.index).intersection(set(MetabolicClasses.index)))


# In[289]:


len(intersect)


# In[290]:


MetabolicTissueHighModLow=Tissue_exp.loc[intersect]


# In[291]:


MetabolicTissueHighModLow=MetabolicTissueHighModLow[~MetabolicTissueHighModLow.index.duplicated(keep='first')]


# In[292]:


# MetabolicTissueHighModLow.to_csv("/data/nandas/Combined_coexp/Part_1_TranscriptionallyRegulatedGenes/Tissue/MetabolicHighModLow.csv")


# In[293]:


# MetabolicTissueHighModLow=pd.read_csv("/data/nandas/Combined_coexp/Part_1_TranscriptionallyRegulatedGenes/Tissue/MetabolicHighModLow.csv",
#                                       index_col=0)


# In[294]:


CoefVar_HighModLow=CalculateCoefVar(MetabolicTissueHighModLow)


# In[295]:


CoefVar_HighModLow.to_csv("CoefVarHighModLow.csv")


# In[296]:


MetabolicTissueHighModLow.to_csv("/data/nandas/Resolve_OR_genes/CaoTissue/MetabolicTissueExp.csv")


# In[297]:


# MetabolicTissueHighModLow.to_csv("/data/nandas/Combined_coexp/Part_1_TranscriptionallyRegulatedGenes/Tissue/MetabolicHighModLow.csv")


# In[298]:


intersect2=list(set(MetabolicClasses.index).intersection(set(HighMod.index)))


# In[299]:


len(intersect2)


# In[300]:


MetabolicHighMod=HighMod.loc[intersect2]


# In[301]:


MetabolicHighMod.to_csv("MetabolicHighMod.csv")


# In[302]:


LogTransformed,hist=Log2Transform(Tissue_exp,
                                  title="Distribution_of_log2(TPM)_of_all_the_genes")


# In[303]:


LogTransformed_HighMod,hist_HighMod=Log2Transform(HighMod,"Distribution of log2(TPM) of genes that are\n atleast moderately or highly expressed in one tissue")


# In[304]:


HighMod


# ## Calculate CV

# In[305]:


CoefVar_HighMod=CalculateCoefVar(HighMod)


# In[306]:


CoefVar_HighMod['Class']=MetabolicClasses['Class']


# In[307]:


CoefVar_HighMod


# In[308]:


LowExpGenesBin=LowExpGenes[['Neurons']]


# In[309]:


LowExpGenesBin['Tissue bin']='Lowly expressed'


# In[310]:


LowExpGenesBin.drop(columns=['Neurons'],inplace=True)


# In[311]:


LowExpGenesBin['Class']=MetabolicClasses['Class']


# In[312]:


LowExpGenesBin


# In[313]:


CoefVar_HighMod_Bin=CoefVar_HighMod.append(LowExpGenesBin)


# In[314]:


# np.unique(CoefVar_HighMod_Bin.Class)


# In[315]:


GeneSets=pd.read_excel("/data/nandas/Combined_coexp/Pathway_enrichment/NewSets_090420/Gene2Pathways_090320.xlsx",
                       header='infer',index_col=0)


# In[316]:


CoefVar_HighMod_Bin['LEVEL 1']=GeneSets['LEVEL 1']
CoefVar_HighMod_Bin['LEVEL 2']=GeneSets['LEVEL 2']
CoefVar_HighMod_Bin['LEVEL 3']=GeneSets['LEVEL 3']
CoefVar_HighMod_Bin['LEVEL 4']=GeneSets['LEVEL 4']


# In[317]:


count=0
CoefVar_HighMod_Bin.Bin=np.NaN
for r in CoefVar_HighMod_Bin.index:
    print(count)
    count=count+1
    x = CoefVar_HighMod_Bin.loc[r]['CoefVar']
    if x < 0.3:
        CoefVar_HighMod_Bin.loc[r,'Bin'] = 'Invariant'
    elif x >= 0.75:
        CoefVar_HighMod_Bin.loc[r,'Bin'] = 'Highly variant'
    elif (x>=0.3) and (x<0.75):
        CoefVar_HighMod_Bin.loc[r,'Bin'] = 'Moderately variant'


# In[318]:


CoefVar_HighMod_Bin.Class.replace('A','iCEL1314 gene',inplace=True)
CoefVar_HighMod_Bin.Class.replace('B','Other metabolic gene',inplace=True)
CoefVar_HighMod_Bin.Class.replace('C','Other metabolic gene',inplace=True)
CoefVar_HighMod_Bin.Class.replace('D','Other metabolic gene',inplace=True)


# In[319]:


r='WBGene00000480'
CoefVar_HighMod_Bin.loc[r]['Tissue bin']


# In[320]:


for r in CoefVar_HighMod_Bin.index:
    print(r)
    x = CoefVar_HighMod_Bin.loc[r]['Tissue bin']
    print(x)
    if x=='Lowly expressed':
        print("Low expressed: {}".format(x))
        CoefVar_HighMod_Bin.at[r,'Bin']='Lowly expressed'


# In[321]:


for index in CoefVar_HighMod_Bin.index:
    print(index)
#     print((HighModCV.loc[index]['Class'])
    Class=(CoefVar_HighMod_Bin.loc[index]['Class']) 
    A=(Class=='iCEL1314 gene')
    print ("Class A: {}".format(A))
    B=(Class=='Other metabolic gene')
#     C=(Class=='Lowly expressed')
    print ("Class B: {}".format(B))
    
    if not(A or B ):
        print("non-metabolic")
        CoefVar_HighMod_Bin.at[index,['Class']]='Non-metabolic gene'
    else:
        continue


# In[322]:


np.unique(CoefVar_HighMod_Bin.Class)


# In[323]:


CoefVar_HighMod_Bin['WormBaseID']=CoefVar_HighMod_Bin.index


# In[324]:


CoefVar_HighMod_Bin=wb_to_gene(CoefVar_HighMod_Bin)


# In[325]:


CoefVar_HighMod_Bin.drop(columns=['Tissue bin'],inplace=True)


# In[326]:


CoefVar_HighMod_Bin.reset_index(inplace=True)
CoefVar_HighMod_Bin.set_index(['WormBaseID'],inplace=True)


# In[327]:


np.unique(CoefVar_HighMod_Bin.Class)


# In[328]:


np.unique(CoefVar_HighMod_Bin.Bin)


# In[329]:


CoefVar_HighMod_Bin.to_csv("CVTissue_AllGenes.csv")


# In[330]:


CoefVar_HighMod_Bin


# In[331]:


iCEL1314genes=CoefVar_HighMod_Bin[CoefVar_HighMod_Bin.Class=='iCEL1314 gene']


# In[332]:


iCEL1314genes.to_csv("CVTissue_iCEL1314.csv")


# In[333]:


iCEL1314genes


# In[334]:


CoefVar_HighMod.to_csv("CoefVarHighModTissue.csv")


# In[335]:


CoefVar_HighMod['Class']=MetabolicClasses['Class']
# CoefVar_HighMod_Bin.index=CoefVar_HighMod.index


# In[336]:


PlotCoefVar(CoefVar_HighMod)


# In[337]:


CoefVar_HighMod


# In[338]:


CoefVar_HighMod.sort_values(by=['CoefVar'])


# In[339]:


metaboliccoefvar=list(set(MetabolicClasses.index).intersection(set(CoefVar_HighMod.index)))
MetabolicCoefVarHighMod=CoefVar_HighMod.loc[metaboliccoefvar]


# In[340]:


MetabolicClasses


# In[341]:


MetabolicCoefVarHighMod.sort_values(by=['CoefVar'],inplace=True)


# In[342]:


MetabolicCoefVarHighMod


# In[343]:


xcoords=[0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2]
fig, ax = plt.subplots(figsize=(7,5))
AllGenes=CoefVar_HighMod.CoefVar.hist(ax=ax,color='green',bins=100,label='All genes',alpha=0.5)
MetabolicGenes=CoefVar_HighMod.loc[intersect2].hist(ax=ax,color='midnightblue',bins=100,label='Metabolic genes',alpha=0.3)  
ax.grid(False)
plt.title("Distribution of Coefficient of variation")
ax.set_xticks(np.arange(0, 2.6, step=0.3))
ax.set_xlabel('Coefficient of Variation')
ax.set_ylabel('Number of genes')
for xc in xcoords:
    ax.axvline(xc,color='blue',linestyle='--')
plt.legend(loc='best')
plt.savefig("CoefVar_differentCV.svg", dpi=300)


# In[344]:


CoefVar_HighMod


# In[345]:


# CoefVar_HighMod=CoefVar_HighMod[CoefVar_HighMod]


# ## Defining genes as regulated using CV=0.75

# In[346]:


CoefVar_1=CoefVar_HighMod[CoefVar_HighMod.CoefVar>=0.75]


# In[347]:


CoefVar_1


# In[348]:


CoefVar_0_5=CoefVar_HighMod[CoefVar_HighMod.CoefVar<=0.3]
CoefVar_0_75=CoefVar_HighMod[CoefVar_HighMod.CoefVar>0.3]
CoefVar_0_75=CoefVar_0_75[CoefVar_0_75.CoefVar<0.75]


# In[448]:


intersect2=list(set(MetabolicClasses.index).intersection(set(CoefVar_0_75.index)))
metabolicCoefVar_0_75=CoefVar_0_75.loc[intersect2]


# In[451]:


intersect3=list(set(MetabolicClasses.index).intersection(set(CoefVar_0_5.index)))
metabolicCoefVar_0_5=CoefVar_0_5.loc[intersect3]


# In[459]:


len(metabolicCoefVar_0_5[metabolicCoefVar_0_5.Class=='A'])


# In[349]:


# Bin3=CoefVar_HighMod[CoefVar_HighMod.CoefVar>=1]


# In[350]:


CoefVar_test=CoefVar_HighMod[CoefVar_HighMod.CoefVar>=0.3]


# In[351]:


intersect=list(set(MetabolicClasses.index).intersection(set(CoefVar_test.index)))
CoefVar_test.loc[intersect]


# In[352]:


LowExpressed=(Tissue_exp.shape[0])-((CoefVar_1.shape[0])+(CoefVar_0_75.shape[0])+(CoefVar_0_5.shape[0]))


# In[353]:


CoefVar_1


# In[354]:


CoefVar_1.to_csv("RegulatedTissueAllGenes.csv")


# In[355]:


nmbin3=list(set(CoefVar_1.index).difference(set(MetabolicClasses.index)))


# In[356]:


NonMetabolicRegulated=CoefVar_1.loc[nmbin3]
NonMetabolicRegulated.to_csv("NonMetabolicRegulated_Tissue.csv")


# In[357]:


nmbin2=list(set(CoefVar_0_75.index).difference(set(MetabolicClasses.index)))
Bin2NonMetabolic=CoefVar_0_75.loc[nmbin2]
Bin2NonMetabolic.to_csv("Bin2NonMetabolic_Tissue.csv")


# In[358]:


nmbin1=list(set(CoefVar_0_5.index).difference(set(MetabolicClasses.index)))
Bin1NonMetabolic=CoefVar_0_5.loc[nmbin1]
Bin1NonMetabolic.to_csv("Bin1NonMetabolic_Tissue.csv")


# In[359]:


fig, ax = plt.subplots(figsize=(12,10))

size = 0.3
vals = [(CoefVar_1.shape[0]),(CoefVar_0_75.shape[0]),(CoefVar_0_5.shape[0]),LowExpressed]
outer_labels=["Bin 3", "Bin 2", "Bin 1","Low Expressed"]
print(vals)

explode = (0.1, 0)
cmap = plt.get_cmap("tab20c")
outer_colors = ['#54CB73','#F16718','#ff0000','#b3b3cc']
ax.pie(vals, radius=1, labels=outer_labels, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w',linewidth=3),autopct='%1.1f%%',shadow=False,
       textprops={'size': 'larger','fontweight':'bold'},
      pctdistance=0.85,labeldistance=1.05,startangle=90,counterclock=False)
plt.title("Percentage of transcriptionally regulated genes across tissues",fontweight="bold",fontsize=12)


# In[360]:


Tissue_exp.to_csv("TissueHighModLow.csv")


# In[361]:


# fig, ax = plt.subplots(figsize=(12,10))

# size = 0.3
# vals = [CoefVar_1.shape[0],(x.shape[0]-CoefVar_1.shape[0])]
# outer_labels=["Regulated", "Non Regulated"]
# print(vals)

# explode = (0.1, 0)
# cmap = plt.get_cmap("tab20c")
# outer_colors = ['#54CB73','#F16718']
# ax.pie(vals, radius=1, labels=outer_labels, colors=outer_colors,
#        wedgeprops=dict(width=size, edgecolor='w',linewidth=3),autopct='%1.1f%%',shadow=False,
#        textprops={'size': 'larger','fontweight':'bold'},
#       pctdistance=0.85,labeldistance=1.05,startangle=90,counterclock=False)
# plt.title("Percentage of transcriptionally regulated genes across tissues",fontweight="bold",fontsize=12)


# ## Finding metabolic genes as regulated

# In[362]:


Metabolic_HighMod_Class=MetabolicClasses.loc[intersect2]


# In[363]:


Metabolic_HighMod_Class=Metabolic_HighMod_Class.Class


# In[364]:


Metabolic_HighMod_Class


# In[365]:


Metabolic_HighMod_Class=pd.DataFrame(Metabolic_HighMod_Class)


# In[366]:


MetabolicCoefVar=CoefVar_HighMod.loc[intersect2]


# In[367]:


MetabolicCoefVar=MetabolicCoefVar[~MetabolicCoefVar.index.duplicated(keep='first')]


# In[368]:


MetabolicCoefVar.shape


# In[369]:


# Metabolic_HighMod_Class=SeqToGene(Metabolic_HighMod_Class)


# In[370]:


# MetabolicCoefVar=SeqToGene(MetabolicCoefVar)


# In[371]:


x=list(set(MetabolicCoefVar.index).intersection(Metabolic_HighMod_Class.index))


# In[372]:


y=MetabolicCoefVar.loc[x]


# In[373]:


y=y[~y.index.duplicated(keep='first')]


# In[374]:




Metabolic_HighMod_Class['CV']=y['CoefVar']


# In[375]:


# Metabolic_HighMod_Class=wb_to_gene(Metabolic_HighMod_Class)


# In[376]:


Metabolic_HighMod_Class[Metabolic_HighMod_Class.CV>=0.75].sort_values(by=['CV'])


# In[ ]:





# In[377]:


# CoefVar_1=wb_to_gene(CoefVar_1)


# ## Calculating transcriptionally regulated for metabolic genes only

# In[378]:


intersect3=list(set(CoefVar_1.index).intersection(set(intersect2)))
RegulatedMetabolic=Metabolic_HighMod_Class.loc[intersect3]


# In[379]:


len(intersect3)


# In[380]:


RegulatedMetabolic['Class'].to_csv("RegulatedMetabolic_Tissue.csv")


# In[445]:


RegulatedMetabolic[~(RegulatedMetabolic.Class=='A')]


# In[382]:


MetabolicClasses=MetabolicClasses[~MetabolicClasses.index.duplicated(keep='first')]


# In[383]:


MetabolicTissueHighModLow['Class']=MetabolicClasses['Class']


# In[384]:


MetabolicTissueHighModLow[~(MetabolicTissueHighModLow.Class=='A')]


# In[385]:


MetabolicTissueHighModLow['Class'].to_csv("TotalMetabolic_Tissue.csv")


# In[386]:


MetabolicTissueHighModLow


# In[387]:


RegulatedMetabolic_number=RegulatedMetabolic.shape[0]/(MetabolicTissueHighModLow.shape[0])


# In[388]:


RegulatedMetabolic_number


# In[389]:


Bin1=Metabolic_HighMod_Class[Metabolic_HighMod_Class.CV<0.3]


# In[390]:


Bin1


# In[391]:


Bin2=Metabolic_HighMod_Class[Metabolic_HighMod_Class.CV>=0.3]
Bin2=Bin2[Bin2.CV<0.75]


# In[392]:


Bin1.to_csv("Bin1_Tissue.csv")


# In[393]:


Bin2.to_csv("Bin2_Tissue.csv")


# In[ ]:





# In[394]:


NonRegulated_Metabolic=Metabolic_HighMod_Class.drop(RegulatedMetabolic.index)


# In[395]:


NonRegulated_Metabolic['Class'].to_csv("NonRegulated_Tissue.csv")


# In[396]:


(RegulatedMetabolic.Class=='D').sum()


# In[397]:


(MetabolicTissueHighModLow.Class=='D').sum()


# In[398]:


MetabolicHighMod


# In[399]:


metaboliclowexpressed=list(set(MetabolicTissueHighModLow.index).difference(set(MetabolicHighMod.index)))


# In[400]:


metaboliclowexpressed=pd.DataFrame(metaboliclowexpressed)


# In[462]:


metaboliclowexpressed[metaboliclowexpressed.Class!='A']


# In[401]:


metaboliclowexpressed.set_index([0],inplace=True)


# In[402]:


metaboliclowexpressed.to_csv("LowExpressedTissue.csv")


# In[403]:


## Extracting CV values of propionate shunt genes
prpn=['WBGene00016943','WBGene00001155','WBGene00017301','WBGene00012608','WBGene00000114']
prpn_Class_CV=Metabolic_HighMod_Class.loc[prpn]


# In[404]:


prpn_Class_CV


# In[405]:


prpn_Class_CV=wb_to_gene(prpn_Class_CV)


# In[406]:


mle=list(set(MetabolicClasses.index).intersection(set(metaboliclowexpressed.index)))


# In[407]:


mle=MetabolicClasses.loc[mle]


# In[408]:


metaboliclowexpressed['Class']=mle['Class']


# In[409]:


metaboliclowexpressed


# In[410]:


# (RegulatedMetabolic.Class=='D').sum()
# classes = set(RegulatedMetabolic.Class.values)
classes=['A','B','C','D']
regulated_vals =[];
low_expressed_vals=[]
non_regulated_vals =[];
for val in classes:
    regulated_vals.append((RegulatedMetabolic.Class==val).sum());
    non_regulated_vals.append((NonRegulated_Metabolic.Class==val).sum());
    low_expressed_vals.append((metaboliclowexpressed.Class==val).sum())
print(regulated_vals)
print(non_regulated_vals)
print(low_expressed_vals)


# In[411]:


classes = list(classes)
print(classes)


# In[412]:


MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",index_col=0)


# In[413]:


MetabolicClasses=SeqToGene(MetabolicClasses)
MetabolicClasses=gene_to_wb(MetabolicClasses)


# In[414]:


CoefVar_HighMod.sort_values(by=['CoefVar'],inplace=True)


# In[415]:


CoefVar_HighMod[CoefVar_HighMod.CoefVar>=0.15]


# In[416]:


CoefVar_HighMod.loc['WBGene00010809']


# In[417]:


Tissue_exp.loc['WBGene00004453']


# In[418]:


MetabolicClassesClassA=list(set(MetabolicClasses[MetabolicClasses.Class=='A'].index).intersection(set(MetabolicCoefVarHighMod.index)))


# In[419]:


x=(MetabolicCoefVarHighMod.loc[MetabolicClassesClassA]).sort_values(by=['CoefVar'])


# In[420]:


x[x.CoefVar>=1.05]


# In[421]:


MetabolicCoefVarHighMod[MetabolicCoefVarHighMod.CoefVar>=0.3]


# In[422]:


MetabolicCoefVarHighMod.sort_values(by=['CoefVar'],inplace=True)


# In[423]:


fig, ax = plt.subplots(figsize=(12,10))

size = 0.3
vals = [(RegulatedMetabolic.shape[0]),(Bin2.shape[0]),(Bin1.shape[0]),metaboliclowexpressed.shape[0]]
outer_labels=["Bin 3", "Bin 2", "Bin 1","Low Expressed"]
print(vals)
explode = (0.1, 0)
cmap = plt.get_cmap("tab20c")
outer_colors = ['#54CB73','#F16718','#ff0000','#b3b3cc']
ax.pie(vals, radius=1, labels=outer_labels, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w',linewidth=3),autopct='%1.1f%%',shadow=False,
       textprops={'size': 'larger','fontweight':'bold'},
      pctdistance=0.85,labeldistance=1.05,startangle=90,counterclock=False)
# plt.title("Percentage of transcriptionally regulated genes across tissues",fontweight="bold",fontsize=12)
plt.savefig("PieChartTissueRegulation.svg")


# In[424]:


fig, ax = plt.subplots(figsize=(12,10))

size = 0.3
vals = [(RegulatedMetabolic.shape[0]),(Bin2.shape[0]),(Bin1.shape[0]),metaboliclowexpressed.shape[0]]
# outer_labels=["Bin 3", "Bin 2", "Bin 1","Low Expressed"]
print(vals)
explode = (0.1, 0)
cmap = plt.get_cmap("tab20c")
outer_colors = ['#54CB73','#F16718','#ff0000','#b3b3cc']
ax.pie(vals, radius=1,  colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w',linewidth=3),shadow=False,
       textprops={'size': 'larger','fontweight':'bold'},
      pctdistance=0.85,labeldistance=1.05,startangle=90,counterclock=False)
# plt.title("Percentage of transcriptionally regulated genes across tissues",fontweight="bold",fontsize=12)
plt.savefig("PieChartTissueRegulation.svg")


# In[425]:


fig, ax = plt.subplots(figsize=(12,10))
size = 0.3
vals = np.array([regulated_vals])
outer_labels=["Regulated"]
inner_labels = classes
print(inner_labels)
print(vals)
explode = (0.1, 0)
cmap = plt.get_cmap("tab20c")
outer_colors = ['#54CB73','#F16718']
inner_colors = ['#F0EFC0','#92E7DA','#E7C1E5','cyan','#F0EFC0','#92E7DA','#E7C1E5','cyan']
inner_explode=(0.1,0.1,0.1,0.1,0,0,0,0)
ax.pie(vals.sum(axis=1), radius=1, labels=outer_labels, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w',linewidth=3),shadow=False,
       textprops={'size': 'larger','fontweight':'bold'},
      pctdistance=0.85,labeldistance=1.05,startangle=90,counterclock=False)

ax.pie(vals.flatten(),radius=1-size, colors= inner_colors, labeldistance=0.5,autopct='%1.1f%%',
       wedgeprops=dict(width=0.5, edgecolor='w',linewidth='3'),pctdistance=0.90,startangle=90,counterclock=False)
plt.legend(loc='upper right',bbox_to_anchor=(1.02, 0.85),labels=outer_labels+inner_labels)
ax.set_title("Pie chart showing percentage of transcriptionally regulated metabolic genes",pad=3,fontweight="bold",fontsize=12)
plt.tight_layout()
plt.savefig("TranscriptionalRegulationTissueNoLabels.svg")
plt.show()


# In[426]:


MetabolicTissueHighModLow=SeqToGene(MetabolicTissueHighModLow)
MetabolicTissueHighModLow=gene_to_wb(MetabolicTissueHighModLow)


# In[427]:


MetabolicTissueHighModLow=MetabolicTissueHighModLow[~MetabolicTissueHighModLow.index.duplicated(keep='first')]
MetabolicClasses=MetabolicClasses[~MetabolicClasses.index.duplicated(keep='first')]
MetabolicTissueHighModLow['Class']=MetabolicClasses['Class']


# In[428]:


MetabolicTissueHighModLow[MetabolicTissueHighModLow.Class=='D']


# In[429]:


fig, ax = plt.subplots(figsize=(12,10))
size = 0.3
vals = np.array([regulated_vals])
outer_labels=["Regulated"]
inner_labels = classes
print(inner_labels)
explode = (0.1, 0)
cmap = plt.get_cmap("tab20c")
outer_colors = ['#54CB73','#F16718']
inner_colors = ['#F0EFC0','#92E7DA','#E7C1E5','cyan','#F0EFC0','#92E7DA','#E7C1E5','cyan']
inner_explode=(0.1,0.1,0.1,0.1,0,0,0,0)
ax.pie(vals.sum(axis=1), radius=1, labels=outer_labels, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w',linewidth=3),shadow=False,
       textprops={'size': 'larger','fontweight':'bold'},
      pctdistance=0.85,labeldistance=1.05,startangle=90,counterclock=False)

ax.pie(vals.flatten(),radius=1-size, colors= inner_colors, labeldistance=0.5,
       wedgeprops=dict(width=0.5, edgecolor='w',linewidth='3'),pctdistance=0.90,startangle=30,counterclock=False)
plt.legend(loc='upper right',bbox_to_anchor=(1.02, 0.85),labels=outer_labels+inner_labels)
# ax.set_title("Pie chart showing percentage of transcriptionally regulated metabolic genes",pad=3,fontweight="bold",fontsize=12)
plt.tight_layout()
plt.savefig("TranscriptionalRegulationTissueNoLabels.svg")
plt.show()


# In[430]:


fig, ax = plt.subplots(figsize=(12,10))

size = 0.3
vals = np.array([regulated_vals, non_regulated_vals])
outer_labels=["Regulated", "Non Regulated"]
inner_labels = classes
print(vals)
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
plt.tight_layout()
plt.savefig("TranscriptionalRegulationTissue.svg")
plt.show()


# In[431]:


RegulatedMetabolic=gene_to_wb(RegulatedMetabolic)
NonRegulated_Metabolic=gene_to_wb(NonRegulated_Metabolic)
MetabolicClasses


# In[432]:


NonRegulated_Metabolic


# In[433]:


RegulatedMetabolic


# In[434]:


RegulatedMetabolic.to_csv("RegulatedMetabolicTissue.csv")
NonRegulated_Metabolic.to_csv("NonRegulatedTissue.csv")


# In[435]:


ClassARegulatedMetabolic=RegulatedMetabolic[RegulatedMetabolic.Class=='A']


# In[436]:


OtherClassesRegulated=RegulatedMetabolic[~(RegulatedMetabolic.Class=='A')]


# In[437]:


TotalMetabolic=list(set(MetabolicClasses.index).intersection(Tissue_exp.index))


# In[438]:


TotalMetabolic=MetabolicClasses.loc[TotalMetabolic]


# In[439]:


ClassAMetabolic=TotalMetabolic[TotalMetabolic.Class=='A']


# In[440]:


OtherClassesMetabolic=TotalMetabolic[~(TotalMetabolic.Class=='A')]


# In[441]:


OtherClassesMetabolic


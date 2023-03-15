#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import re
from pylab import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
import fnmatch
# from CatExp import *
from sklearn import preprocessing


# In[2]:


# ! mkdir /data/nandas/Transcription/AllDatasets_061922/


# In[3]:


Base_dir='/data/nandas/Transcription/AllDatasets_061922/'
os.chdir(Base_dir)


# In[4]:


#Calculate maximum CV from all the datasets
def MaxCV(x):
    x['Max_CV']=x.max(axis=1)
    return x

# Calculate maximum Z score across all the datasets
def Max_Z(x):
    x['Max_Z']=x.max(axis=1)
    return x


# In[5]:


# !mkdir /data/nandas/Transcription/AllDatasets_012422/


# In[6]:


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

def CalculateCoefVar(ZScore_df,count):
    Std=ZScore_df.std(axis=1,skipna=True)
    Mean=ZScore_df.mean(axis=1,skipna=True)
    CoefVar=pd.DataFrame([])
#     from scipy.stats import variation 
#     Variation=variation(ZScore_df.values, axis = 1)
    CoefVar['CoefVar_{}'.format(count)]=(Std/Mean)
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
    plt.savefig("{}.png".format(title))
    return df,hist

def StackedDF(df):
    x=df.unstack()
    x.index.rename(['Condition', 'Gene'], inplace=True)
    x = x.to_frame('Log2RPKM').reset_index()
    x.dropna(inplace=True)
    x = x[x.Condition!=x.Gene]
    x['check_string'] = x.apply(lambda row: ''.join(sorted([row['Condition'], row['Gene']])), axis=1)
    x.drop_duplicates('check_string',keep='first',inplace=True)
    x.drop(columns=['check_string'],inplace=True)
    return x

def normal(x,mu,sigma,N):
    '''Returns a normal distribution scaled with <N> as the population size.'''
    deltax=x[2]-x[1];
    coef=deltax*N*(1./(sqrt(2.*pi*sigma**2)));
    return coef*exp(-(x-mu)**2/2/sigma**2)

def superimpose_bimodal(x,mu1,sigma1,mu2,sigma2,Np):
    '''Superimposition function for bi-modal curve fit. Returns a bimodal distribution with two subpopulations,
    one defined by <mu1> and <sigma1> and scaled with global Ndetected-<Np> serving 
    as the subpopulation size, and the other defined by 
    <mu2> and <sigma2> and scaled with <Np>. <Np> is the subpopulation size for curve 2.'''    
    return normal(x,mu1,sigma1,Ndetected-Np)+normal(x,mu2,sigma2,Np);
    
def bimodal(datao,expectedStats=[0,1.5,5.,2.5],Nbin=100,figsize=(12,6),xlims=None,showPlot=True):
    ''' 
    Fits the superimposition of two Gaussian curves to a histogram of data in <datao> with the number of bins indicated by <Nbin>.
    
    <expectedStats> is an estimate of the mean (mu) and standard deviation (sigma) of the two curves to be fitted (as [mu1,sigma1,mu2,sigma2]). 
    If <showPlot> is True, the histogram and fit are plotted. The cfunction should be run iteratively to change <expectedStats> in case of misfit.
    <figsize> is the width and height of figure to be plotted in the respective order.
    <xlims> is the range of the values (bins) used in the histogram. If not provided (default), this range is automatically calculated.
    
    Returns a dictionary of fitted parameters and statistics:
        -mu and sigma are as defined above
        -N indicates subpopulation size for the corresponding curves
        -n is number of bins
        -RSS is Residual Sum of Squares
        -Rsquared is R2 of the fit
        
    USAGE
    
        Stats=bimodal(data,expectedStats=[0,1.5,5.,2.5],Nbin=100,figsize=(12,6),xlims=None,showPlot=True);
    '''
    global Ndetected
    #Static Input
    sumcolor='black';
    bgcolor='black';
    poscolor='black';
    bgline=':';
    posline='--';

    #Data processing and initial plotting
    V=np.array(datao); 
    V=V[np.nonzero(V)];
    data=np.log2(V);
    if xlims is None:
        txlim=(-(int(abs(min(data))/5.)+1)*5.,(int(max(data)/5.)+1)*5.);
    else:
        txlim=(xlims[0],xlims[1]);
    figure(figsize=figsize);
    y,x,_=hist(data,Nbin,range=txlim,alpha=.3,label='data');
    x=(x[1:]+x[:-1])/2 # for len(x)==len(y)
    #Fitting
    Ndetected=len(data);
    expected=tuple(expectedStats+[Ndetected/2.]);
    params,_=curve_fit(superimpose_bimodal,x,y,expected);
    yhat=superimpose_bimodal(x,*params);
    resid=y-yhat;
    ss_res=np.sum(resid**2);
    rss=ss_res/len(y);
    ss_tot=np.sum((y-np.mean(y))**2);
    Rsquared=1.-(ss_res/ss_tot);
    Dparams={'mu1':params[0],'sigma1':params[1],'mu2':params[2],'sigma2':params[3],             'N1':Ndetected-params[4],'N2':params[4],'Rsquared':Rsquared,'RSS':rss,'n':len(y)};            
    #Plot the fits
    plot(x,superimpose_bimodal(x,*params),color=sumcolor,lw=3,label='superimposed');
    plot(x,normal(x,params[0],params[1],Ndetected-params[-1]),color=bgcolor,lw=3,linestyle=bgline,label='curve 1');
    plot(x,normal(x,*params[2:]),color=poscolor,lw=3,linestyle=posline,label='curve 2');
    xlim(txlim[0],txlim[1]);
    legend();
    if showPlot:
        show();
    else:
        close();
    return Dparams

def categorize_absCutoff(T,cutoffs,excludeCols=[]):
    '''
    Categorizes genes (rows) of expression table <T> according to <cutoffs>, list-like object with [rare cutoff,low cutoff,high cutoff].
    
    Genes with expression levels less than rare cutoff are categorized as Rare.
    Other genes with expression levels less than low cutoff are categorized as Low.
    Genes with expression levels greater than high cutoff are categorized as High.
    All other genes are categorized as Moderate.
    
    <excludeCols> list indicates which columns of <T> should not be categorized.
    
    USAGE
    
        Tcat=categorize_absCutoff(T,cutoffs,excludeCols=[])
    '''
    Lcols=[];
    for col in T.columns:
        if not col in excludeCols:
            Lcols.append(col);
    Tcat=pd.DataFrame(index=T.index,columns=Lcols);
    for i in Tcat.index:
        for col in Tcat.columns:
            val=T[col][i];
            if val<=0:
                Tcat[col][i]='Rare';
            elif val<cutoffs[0]:
                Tcat[col][i]='Rare';
            elif val<cutoffs[1]:
                Tcat[col][i]='Low';
            elif val>cutoffs[2]:
                 Tcat[col][i]='High';
            else:
                Tcat[col][i]='Moderate';       
    return Tcat
    
def stackedCat(Tcat):
    '''Plots a stacked bar graph of categories (Rare, Low, Moderate, High) for each column of <Tcat>.'''
    h,m,l,r=[],[],[],[];
    N=len(Tcat.columns);
    ind = np.arange(N); 
    width = 0.8;
    
    for col in Tcat.columns:
        s=pd.value_counts(Tcat[col]);
        h.append(s['High']);
        m.append(s['Moderate']);
        l.append(s['Low']);
        r.append(s['Rare']);

    pr = plt.bar(ind, r, color='red',width=width);
    pl = plt.bar(ind, l, color='orange',bottom=r,width=width);
    pm = plt.bar(ind, m, color='grey',bottom=[r[i]+l[i] for i in range(N)],width=width);
    ph = plt.bar(ind, h, color='green',bottom=[r[i]+l[i]+m[i] for i in range(N)],width=width);

    plt.ylabel('Number of genes');
    plt.xticks([i+0.5 for i in ind], Tcat.columns);
    plt.xticks(rotation=90, fontsize=8);
    plt.legend((pr[0],pl[0],pm[0],ph[0]), ('Rare', 'Low','Moderate','High'));

    plt.show()   

def relativeExp(Texpo,Tcato,tao_rel,fc_mid=1.5,fc_end=4.,fcHigh=None,fcLow=None):
    '''
    Recategorizes some Moderate genes in <Tcato>, a table of categorized genes (rows) in tissues or conditions (columns) as Low or High based on relative expression levels.
    
    Expression levels are provided as <Texpo>, a table that matches <Tcato> in rows and columns, but has expression levels instead of categories.
    
    A heuristic method is used such that the expression profile of each gene (row) is first obtained by sorting the row in <Texpo> from low to high expression.
    Then fold changes from one tissue or condition to the next is monitored and significant jumps are tracked.
    For an increase to be considered significant, the fold change (FC) should be greater than a threshold and the larger number should be greater than <tao_rel>.
    If the increase is  greater than an FC threshold (FC for low), then the lower value and all below are labeled Low.
    If the increase is  greater than an FC threshold (FC for high), the higher value and all above are labeled High.
    Variable (High and Low) labeling of the same value in different steps result in no categorization for that value.
    FC for low and FC for high values for every increment (for N columns, there are N-1 increments) can be provided as lists (<fcLow> and <fcHigh>).
    If low vs high thresholds are to be the same and the same threshold is going to be used for all middle increments and another threshold for terminal increments, then <fc_mid> and <fc_end> thresholds can be provided instead, as single thresholds for the respective increment sets.  
    If the same threshold is to be used throughout, then just enter this value for both <fc_mid> and <fc_low>.
    If <fcHigh> or <fcLow> list is provided, <fc_mid> and <fc_end> will be automatically null for the thresholds determined by the list.
    
    In the end, a Moderate gene in a tissue can be labeled as Low only if its expression level is less than <tao_rel> and labeled as High only if its expression level is higher than <tao_rel>. 
    
    USAGE
    
        Tcat_final=relativeExp(Texp,Tcat,tao_rel,fc_mid=1.5,fc_end=4.,fcHigh=None,fcLow=None)
 
    '''
    Tcatf=Tcato.copy();
    Texp=Texpo.drop('ave',axis=1);

    n=len(Texp.columns);
    if not fcHigh:
        fcHigh=[fc_end]+(n-2)*[fc_mid];
    elif not len(fcHigh)==n-1:
        raise Exception('fcHigh vector must of size n-1, where n is the number of conditions');
    if not fcLow:
        fcLow=(n-2)*[fc_mid]+[fc_end];
    elif not len(fcLow)==n-1:
        raise Exception('fcLow vector must of size n-1, where n is the number of conditions');
    for gene in Texp.index:
        row=Texp.loc[gene,Texp.columns].copy(); ###
        row=row.sort_values(); ###
        S=[set() for i in range(n)];
        fc=[0];
        for k in range(1,n):
            if row[k-1]:
                fc.append(row[k]/row[k-1]);
            else:
                fc.append(inf);       
        for k in range(1,n):
            if row[k]>tao_rel:
                #check relatively high
                if fc[k]>fcHigh[k-1]:
                    for i in range(0,k):
                        S[i].add('not high');
                    for i in range(k,n):
                        S[i].add('high');   
                #check relatively low
                if fc[k]>fcLow[k-1]:
                    for i in range(0,k):
                        S[i].add('low');
                    for i in range(k,n):
                        S[i].add('not low');                   
        for i in range(n):
            cond=row.index[i];
            if Tcatf[cond][gene]=='Moderate':
                if row[i]>tao_rel:
                    if S[i].issubset(['high','not low']) and 'high' in S[i]:
                        Tcatf [cond][gene]='High';
                elif S[i].issubset(['low','not high']) and 'low' in S[i]:
                    Tcatf [cond][gene]='Low';    
    return Tcatf

def plotCatExp(gene,Texp,Tcato,Tcatf):
    '''
    Plots a bar chart that shows ascending tissue (or condition) profile of <gene> based on expression level table <Texp> with genes as rows and tissues (coditions) as columns.
    Categorization tables from categorize_absCutoff() function (<Tcato>) and relativeExp() function (<Tcatf>) are to be provided.
    
    Bars are colored according to category ({'High':'green','Moderate':'gray','Low':'orange','Rare':'red'}).
    Categories that differ between <Tcato> and <Tcatf> are indicated by hatches.
    
    USAGE
        
        plotCatExp(gene,Texp,Tcato,Tcatf)
    '''
    width=0.5;
    dcolor={'High':'green','Moderate':'gray','Low':'orange','Rare':'red'};

    x=Texp.loc[gene,Texp.columns];
    x=x.sort_values();  
    ind=np.arange(len(x));      
    rects = plt.bar(ind, x.values, width=width);
    plt.xticks([i+width/2 for i in ind],x.index);
    plt.xticks(rotation=90, fontsize=8);
    for i in ind:
        cond=x.index[i];
        rects[i].set_facecolor(dcolor[Tcatf[cond][gene]]);
        if not Tcatf[cond][gene]==Tcato[cond][gene]:
            rects[i].set_hatch('/');
    plt.show(); 
    
def PlotDistributionofLogData(df,title):
    #df=df.drop(columns=['GENEID','GENENAME'])
    df=df.stack()
#     df=df.loc[~(df==0)]
    df=pd.DataFrame(df)
    ## Plot the values of all cells in the dataframe to understand the distribution
    hist=df.hist(grid=False,bins=500,color='skyblue')
#     plt.axvline(1.3,color='r')
#     plt.axvline(df.mode()[1],color='g')
    plt.xlabel("log2 transformed expression")
    plt.ylabel("Frequency")
    plt.title("{}".format(title))
    plt.savefig("{}.png".format(title))
    df2=2**df
    hist2=df2.hist(grid=False,color='skyblue')
    plt.savefig("{}_antilog.png".format(title))
    return df,hist

# def PlotDistributionofUnLogData(df,title):
#     #df=df.drop(columns=['GENEID','GENENAME'])
#     df=(2**df)
#     print(df)
#     ## Plot the values of all cells in the dataframe to understand the distribution
#     hist=df.hist(grid=False,bins=500,color='skyblue')
# #     plt.axvline(1.3,color='r')
# #     plt.axvline(df.mode()[1],color='g')
#     plt.xlabel("expression")
#     plt.ylabel("Frequency")
#     plt.title("{}".format(title))
#     plt.savefig("{}_antilog.png".format(title))
#     return df,hist



# In[7]:


# Reading metabolic classes
MetabolicClasses=pd.read_csv("/data/nandas/MetabolicClasses_August_SN_090221.csv",
                             index_col=0)


# In[8]:


MetabolicClasses=SeqToGene(MetabolicClasses)
MetabolicClasses=gene_to_wb(MetabolicClasses)
MetabolicClasses=MetabolicClasses[~MetabolicClasses.index.duplicated(keep='first')]


# In[9]:


# `!ls


# In[10]:


## Calculating Combined CV
count=0
x=pd.DataFrame([])
Dataset_dir='/data/nandas/Combined_coexp/Compendium/BatchCorrectedFiles061822'
for files in os.listdir(Dataset_dir):
    if fnmatch.fnmatch(files, "WB*_a.ce*csv"):
        count=count+1
        print(files)
        df=pd.read_csv("{}/{}".format(Dataset_dir,files),index_col=0,header='infer',low_memory=False)
#         df.drop(columns=['IDENTIFIER', 'GWEIGHT'],index=['EWEIGHT'],inplace=True)
        
        df.replace("\\N",np.nan,inplace=True)
        df.replace("Inf",np.nan,inplace=True)
        df.replace(np.inf,np.nan,inplace=True)
        
#         print(df.sum(axis=1).sort_values(ascending=False))
        df=df.astype(float)
        l = df.values #returns a numpy array
        min_max_scaler = preprocessing.MinMaxScaler()
        x_scaled = min_max_scaler.fit_transform(l)
        df = pd.DataFrame(x_scaled,index=df.index,columns=df.columns)
        df=df[((df<0).sum(axis=1))!=(df.shape[1])]
#         df[df<0]=np.nan
#         print(df.isna().sum(axis=1).sort_values())
        #print(df)
        print(np.unique(df.values))
        
        cv=CalculateCoefVar(df,count=files.split(".csv")[0])
        print(cv)
        cv.to_csv("CV_{}.csv".format(files.split(".csv")[0]))
        x=pd.concat([x,cv],axis=1)
        print(cv.shape)
    
    elif fnmatch.fnmatch(files, "WB*_b.ce*csv"):
        count=count+1
        print(files)
        df=pd.read_csv("{}/{}".format(Dataset_dir,files),index_col=0,header='infer',low_memory=False)
#         df.drop(columns=['IDENTIFIER', 'GWEIGHT'],index=['EWEIGHT'],inplace=True)
        
        df.replace("\\N",np.nan,inplace=True)
        df.replace("Inf",np.nan,inplace=True)
        df.replace(np.inf,np.nan,inplace=True)
        
#         print(df.sum(axis=1).sort_values(ascending=False))
        df=df.astype(float)
        l = df.values #returns a numpy array
        min_max_scaler = preprocessing.MinMaxScaler()
        x_scaled = min_max_scaler.fit_transform(l)
        df = pd.DataFrame(x_scaled,index=df.index,columns=df.columns)
        df=df[((df<0).sum(axis=1))!=(df.shape[1])]
#         df[df<0]=np.nan
#         print(df.isna().sum(axis=1).sort_values())
        #print(df)
        print(np.unique(df.values))
        
        cv=CalculateCoefVar(df,count=files.split(".csv")[0])
        print(cv)
        cv.to_csv("CV_{}.csv".format(files.split(".csv")[0]))
        x=pd.concat([x,cv],axis=1)
        print(cv.shape)
        
    elif fnmatch.fnmatch(files, "WB*ce*csv"):
        count=count+1
        print(files)
        df=pd.read_csv("{}/{}".format(Dataset_dir,files),sep='\t',index_col=0,header='infer',low_memory=False)
#         df.drop(columns=['IDENTIFIER', 'GWEIGHT'],index=['EWEIGHT'],inplace=True)
        
        df.replace("\\N",np.nan,inplace=True)
        df.replace("Inf",np.nan,inplace=True)
        df.replace(np.inf,np.nan,inplace=True)
        
#         print(df.sum(axis=1).sort_values(ascending=False))
        df=df.astype(float)
        l = df.values #returns a numpy array
        min_max_scaler = preprocessing.MinMaxScaler()
        x_scaled = min_max_scaler.fit_transform(l)
        df = pd.DataFrame(x_scaled,index=df.index,columns=df.columns)
        df=df[((df<0).sum(axis=1))!=(df.shape[1])]
#         df[df<0]=np.nan
#         print(df.isna().sum(axis=1).sort_values())
        #print(df)
        print(np.unique(df.values))
        
        cv=CalculateCoefVar(df,count=files.split(".csv")[0])
        print(cv)
        cv.to_csv("CV_{}.csv".format(files.split(".csv")[0]))
        x=pd.concat([x,cv],axis=1)
        print(cv.shape)
        


# In[11]:


# !ls -lhrt


# In[12]:


# files.split(".csv")[0]


# ## Selecting only live and protein-coding genes

# In[13]:


mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=1)


# In[14]:


mapper_df=mapper_df[mapper_df.Status=='Live']


# In[15]:


mapper_df=mapper_df[mapper_df.Type=='protein_coding_gene']


# In[16]:


MetabolicClasses


# In[17]:


# mapper_df.loc['WBGene00000989']


# In[18]:


xmapper=list(set(x.index).intersection(set(mapper_df.index)))
x=x.loc[xmapper]


# In[19]:


x


# In[20]:


x.to_csv("TotalAllGenes_AllDatasets.csv")


# In[21]:


missing=list(set(mapper_df.index).difference(set(x.index)))


# In[22]:


missing=pd.DataFrame(missing)


# In[23]:


missing.set_index([0],inplace=True)


# In[24]:


missing.to_csv("LowExpGenes.csv")


# In[25]:


# mapper_df.loc['WBGene00086558']


# In[26]:


# MetabolicClasses=SeqToGene(MetabolicClasses)
# MetabolicClasses=gene_to_wb(MetabolicClasses)


# In[27]:


MetabolicClasses=wb_to_gene(MetabolicClasses)


# In[28]:


# mapper_df=wb_to_gene(mapper_df)


# In[29]:


common=list(set(MetabolicClasses.index).intersection(set(mapper_df.index)))


# In[30]:


common


# In[31]:


mle=(set(common).difference(set(MetabolicX.index)))


# In[ ]:


#Dataset1
# Temp_df1=pd.read_csv("WBPaper00041697.ce.rs.csv",sep='\t',index_col=0,header='infer',low_memory=False)
# Temp_df1.replace("\\N",np.nan,inplace=True)
# Temp_df1=Temp_df1.astype(float)
# l = Temp_df1.values #returns a numpy array
# PlotDistributionofLogData(df=Temp_df1,title='Distribution of expression values in temperature dataset 1')


# In[ ]:


# #Dataset2
# Temp_df2=pd.read_csv("WBPaper00039792.ce.mr.csv",sep='\t',index_col=0,header='infer',low_memory=False)
# Temp_df2.replace("\\N",np.nan,inplace=True)
# Temp_df2=Temp_df2.astype(float)
# l = Temp_df2.values #returns a numpy array
# x,hist=PlotDistributionofLogData(df=Temp_df2,title='Distribution of expression values in temperature dataset 2')
# # df=(2**Temp_df2)
# # df=df.stack()
# # df=pd.DataFrame(df)
# # df=df[df!=0]
# # df.hist(bins=300, color='skyblue')


# In[ ]:


# x.values


# In[ ]:


# y=(2 ** x)


# In[ ]:


# plt.hist(y.values)


# In[ ]:


# #Dataset3
# Temp_df3=pd.read_csv("WBPaper00037147.ce.mr.csv",sep='\t',index_col=0,header='infer',low_memory=False)
# Temp_df3.replace("\\N",np.nan,inplace=True)
# Temp_df3=Temp_df3.astype(float)
# l = Temp_df3.values #returns a numpy array
# PlotDistributionofLogData(df=Temp_df3,title='Distribution of expression values in temperature dataset 3')


# In[ ]:


# #Dataset4
# Temp_df4=pd.read_csv("WBPaper00035227.ce.mr.csv",sep='\t',index_col=0,header='infer',low_memory=False)
# Temp_df4.replace("\\N",np.nan,inplace=True)
# Temp_df4=Temp_df4.astype(float)
# l = Temp_df4.values #returns a numpy array
# PlotDistributionofLogData(df=Temp_df4,title='Distribution of expression values in temperature dataset 4')


# In[ ]:


# #Dataset5
# Temp_df5=pd.read_csv("WBPaper00046212.ce.mr.csv",sep='\t',index_col=0,header='infer',low_memory=False)
# Temp_df5.replace("\\N",np.nan,inplace=True)
# Temp_df5=Temp_df5.astype(float)
# l = Temp_df5.values #returns a numpy array
# PlotDistributionofLogData(df=Temp_df5,title='Distribution of expression values in temperature dataset 5')


# In[ ]:


# #Dataset6
# Temp_df6=pd.read_csv("WBPaper00027104.ce.mr.csv",sep='\t',index_col=0,header='infer',low_memory=False)
# Temp_df6.replace("\\N",np.nan,inplace=True)
# Temp_df6=Temp_df6.astype(float)
# l = Temp_df6.values #returns a numpy array
# PlotDistributionofLogData(df=Temp_df6,title='Distribution of expression values in temperature dataset 6')


# In[ ]:


# #Dataset7
# Temp_df7=pd.read_csv("WBPaper00049736.ce.mr.csv",sep='\t',index_col=0,header='infer',low_memory=False)
# Temp_df7.replace("\\N",np.nan,inplace=True)
# Temp_df7=Temp_df7.astype(float)
# l = Temp_df7.values #returns a numpy array
# PlotDistributionofLogData(df=Temp_df7,title='Distribution of expression values in temperature dataset 7')


# In[ ]:


# #Dataset8
# Temp_df8=pd.read_csv("WBPaper00037695.ce.mr.csv",sep='\t',index_col=0,header='infer',low_memory=False)
# Temp_df8.replace("\\N",np.nan,inplace=True)
# Temp_df8=Temp_df8.astype(float)
# l = Temp_df8.values #returns a numpy array
# PlotDistributionofLogData(df=Temp_df8,title='Distribution of expression values in temperature dataset 8')


# In[ ]:


# count=0
# x=pd.DataFrame([])
# for files in os.listdir('.'):
#     if fnmatch.fnmatch(files, "WB*csv"):
#         count=count+1
#         print(files)
#         df=pd.read_csv(files,sep='\t',index_col=0,header='infer',low_memory=False)
# #         df.drop(columns=['IDENTIFIER', 'GWEIGHT'],index=['EWEIGHT'],inplace=True)
        
#         df.replace("\\N",np.nan,inplace=True)
        
# #         print(df.sum(axis=1).sort_values(ascending=False))
#         df=df.astype(float)
#         if count==6:
#             break;


# In[ ]:


# Stacked_df=StackedDF(df)


# In[ ]:


# Stacked_df.Log2RPKM.hist(bins=100)


# In[ ]:


# LowExpGenes=pd.read_csv("/data/nandas/Transcription/KimDevTime_071620/AllLowExpGenes_KimTissue.csv",index_col=0)


# ## Filtering out Low exp genes

# In[ ]:


# LowExpGenes


# In[ ]:


# Calculating max CV across the dataframe
# Max=MaxCV(x)


# In[ ]:


# legCV=list(set(LowExpGenes.index).intersection(set(x.index)))


# In[ ]:


# LowExpGenes=x.loc[legCV]


# In[ ]:


x


# ## excluding low expressed genes

# In[47]:


# HMlegCV=list(set(x.index).difference(set(LowExpGenes.index)))


# In[32]:


HighModCV=x


# In[34]:


HighModCV.to_csv("CV_allgenes_alldatasets.csv")


# In[25]:


HighModCVHigh=pd.DataFrame((HighModCV>=0.75).sum(axis=1))


# In[26]:


LowVariantCV=pd.DataFrame((HighModCV<0.3).sum(axis=1))


# In[27]:


LowVariantCV.set_axis(['Number of datasets with CV<0.3'],axis=1,inplace=True)


# In[28]:


LowVariantCV['NumberOfValues']=HighModCV.notna().sum(axis=1).sort_values()


# In[29]:


LowVariantCV['CVLowFraction']=(LowVariantCV['Number of datasets with CV<0.3']/LowVariantCV.NumberOfValues)


# In[30]:


HighModCV


# In[31]:


HighModCVHigh.set_axis(["Number of datasets with CV>=0.75"],axis=1,inplace=True)


# In[32]:


HighModCV["Number of datasets with CV>=0.75"]=HighModCVHigh["Number of datasets with CV>=0.75"]


# In[33]:


HighModCV['Number of datasets with CV<0.3']=LowVariantCV['Number of datasets with CV<0.3']


# In[34]:


HighModCV.sort_values(ascending=False,by=["Number of datasets with CV>=0.75"],inplace=True)


# In[35]:


HighModCV['CVLowFraction']=LowVariantCV['CVLowFraction']


# In[36]:


HighModCV


# In[37]:


for gene in HighModCV.index:
    if (HighModCV.loc[gene]["Number of datasets with CV>=0.75"])>2:
        HighModCV.at[gene,'Bin']='Highly variant'
    elif (HighModCV.loc[gene]['CVLowFraction'])>=0.95:
        HighModCV.at[gene,'Bin']='Invariant'
    else:
        HighModCV.at[gene,'Bin']='Moderately variant'
        


# In[38]:


HighModCV[(HighModCV.Bin=='Highly variant')]


# In[39]:


HighModCV['Class']=MetabolicClasses['Class']


# In[40]:


HighModCV.Class.replace('A','iCEL1314 gene',inplace=True)


# In[41]:


HighModCV.Class.replace('B','Other metabolic gene',inplace=True)
HighModCV.Class.replace('C','Other metabolic gene',inplace=True)
HighModCV.Class.replace('D','Other metabolic gene',inplace=True)


# In[42]:


HighModCV.Class.replace('np.nan','Non-metabolic gene',inplace=True)


# In[43]:


HighModCV


# In[44]:


for index in HighModCV.index:
    print(index)
#     print((HighModCV.loc[index]['Class'])
    Class=(HighModCV.loc[index]['Class']) 
    A=(Class=='iCEL1314gene')
    B=(Class=='Other metabolic gene')
    
    if not(A or B):
        print("non-metabolic")
        HighModCV.at[index,['Class']]='Non-metabolic gene'


# In[45]:


HighModCV['WormBase ID']=HighModCV.index


# In[46]:


HighModCV=wb_to_gene(HighModCV)


# In[47]:


HighModCV.reset_index(inplace=True)


# In[48]:


HighModCV.rename(columns={'index':'Gene ID'},inplace=True)


# In[49]:


HighModCV.set_index(['WormBase ID'],inplace=True)


# In[50]:


HighModCV.to_csv("CVAllDatasets_AllGenes_062722.csv")


# In[51]:


MetabolicClasses=SeqToWB(MetabolicClasses)
MetabolicClasses=gene_to_wb(MetabolicClasses)


# In[52]:


metabolic=list(set(MetabolicClasses.index).intersection(set(HighModCV.index)))


# In[53]:


HighModCV


# In[56]:


MetabolicHighModCV=HighModCV.loc[metabolic]


# In[57]:


MetabolicHighModCV['Class']=MetabolicClasses['Class']


# In[58]:


MetabolicHighModCV.columns


# In[79]:


# mle=list(set(LowExpGenes.index).intersection(set(MetabolicClasses.index)))
# MetabolicLowExpGenes=LowExpGenes.loc[mle]


# In[80]:


# MetabolicLowExpGenes['Class']=MetabolicClasses['Class']


# In[81]:


# MetabolicLowExpGenes['Number of datasets with CV>=0.75']='Lowly expressed'
# MetabolicLowExpGenes['Number of datasets with CV<0.3']='Lowly expressed'
# MetabolicLowExpGenes['CVLowFraction']='Lowly expressed'
# MetabolicLowExpGenes['Bin']='Lowly expressed'


# In[82]:


# MetabolicLowExpGenes['WormBase ID']=MetabolicLowExpGenes.index
# MetabolicLowExpGenes=wb_to_gene(MetabolicLowExpGenes)


# In[83]:


# MetabolicLowExpGenes.reset_index(inplace=True)
# MetabolicLowExpGenes.set_index(['WormBase ID'],inplace=True)


# In[59]:


MetabolicHighModCV.sort_values(by=['Number of datasets with CV>=0.75'],ascending=False,inplace=True)


# In[60]:


MetabolicHighModCVBin=MetabolicHighModCV


# In[58]:


# MetabolicLowExpGenes


# In[87]:


# MetabolicLowExpGenes['Gene ID']=MetabolicLowExpGenes['index']


# In[61]:


MetabolicHighModCVBin.Class.replace("A","iCEL1314 gene",inplace=True)
MetabolicHighModCVBin.Class.replace("B","Other metabolic gene",inplace=True)
MetabolicHighModCVBin.Class.replace("C","Other metabolic gene",inplace=True)
MetabolicHighModCVBin.Class.replace("D","Other metabolic gene",inplace=True)


# In[62]:


MetabolicHighModCVBin.to_csv("CVCompendium_iCEL1314genes.csv")


# In[63]:


MetabolicHighModCVBin


# In[64]:


MetabolicHighModCV.to_csv("MetabolicCVAllDatasets_AllGenes_061922.csv")


# In[65]:


MetabolicHighModCV[MetabolicHighModCV['Number of datasets with CV>=0.75']>=3]


# In[66]:


MetabolicHighModCV[MetabolicHighModCV['CVLowFraction']>=.95]


# ## Metabolic

# In[81]:


# MetabolicX


# In[72]:


x


# In[73]:


metabolicx=list(set(MetabolicClasses.index).intersection(set(x.index)))


# In[74]:


MetabolicX=x.loc[metabolicx]


# In[75]:


MetabolicX.to_csv("TotalMetabolic_AllDatasets.csv")


# In[76]:


MetabolicX


# In[81]:


MetabolicX['Class']=MetabolicClasses['Class']
# ClassA=MetabolicX[MetabolicX.Class=='A']


# In[82]:


np.unique(MetabolicX.Class)


# In[98]:


# mlegCV=list(set(LowExpGenes.index).intersection(set(MetabolicX.index)))
# MetabolicLowExpGenes=MetabolicX.loc[mlegCV]


# In[105]:


# MetabolicLowExpGenes['Class']=MetabolicClasses['Class']


# In[107]:


# MetabolicLowExpGenes[MetabolicLowExpGenes.Class!='A']


# In[98]:


MetabolicX


# In[83]:


Regulated=MetabolicX[MetabolicX['Number of datasets with CV>=0.75']>2]


# In[100]:


# CVHigh


# In[109]:


# nonlowexpressed=list(set(CVHigh.index).difference(set(MetabolicLowExpGenes.index)))
# CVHigh=CVHigh.loc[nonlowexpressed]


# In[110]:


# CVHigh.set_axis(['HighCVNumber'],axis=1,inplace=True)


# In[111]:


# Regulated=CVHigh[CVHigh.HighCVNumber>2]


# In[84]:


Regulated['Class']=MetabolicClasses['Class']


# In[86]:


Regulated[Regulated.Class=='A']


# In[47]:


# Regulated=Regulated[~Regulated.index.str.contains("WB")]


# In[48]:


# Regulated=wb_to_gene(Regulated)


# In[49]:


# Regulated=gene_to_wb(Regulated)


# In[87]:


Regulated.to_csv("Regulated_AllDatasets.csv")


# In[88]:


x.to_csv("CV_allgenes.csv")


# In[52]:


# CVHi


# ## Non-metabolic

# In[89]:


nonmetabolicx=list(set(x.index).difference(set(MetabolicClasses.index)))


# In[90]:


Totalnm=list(set(mapper_df.index).difference(set(MetabolicClasses.index)))


# In[91]:


len(Totalnm)


# In[94]:


NonMetabolicX=x.loc[nonmetabolicx]


# In[95]:


NonMetabolicX


# In[339]:




# legnonmetabolic=list(set(NonMetabolicX.index).intersection(set(LowExpGenes.index)))


# In[340]:


# LowExpressedNonMetabolic=NonMetabolicX.loc[legnonmetabolic]


# In[388]:


# LowExpressedNonMetabolic


# In[342]:


# LowExpressedNonMetabolic.to_csv("LowExpressedNonMetabolic.csv")


# In[389]:


# legexcludingnonmetabolic=list(set(NonMetabolicX.index).difference(set(LowExpGenes.index)))


# In[96]:


NonMetabolicX=NonMetabolicX


# In[97]:


NonMetabolicX=NonMetabolicX[~NonMetabolicX.index.duplicated(keep='first')]


# In[98]:


NonMetabolicX=gene_to_wb(NonMetabolicX)


# In[99]:


NonMetabolicX


# In[100]:


NonMetabolicX[NonMetabolicX.index.str.contains('WB')]


# In[117]:


# NonMetabolicCVHigh=pd.DataFrame((NonMetabolicX>=0.75).sum(axis=1))


# In[101]:


RegulatedNonMetabolic=NonMetabolicX[NonMetabolicX['Number of datasets with CV>=0.75']>2]


# In[102]:


RegulatedNonMetabolic.to_csv("Regulated_NonMetabolic.csv")


# In[103]:


RegulatedNonMetabolic


# In[396]:


# NonMetabolicCVLow=pd.DataFrame((NonMetabolicX<0.3).sum(axis=1))


# In[397]:


# NonMetabolicCVLow['NumberOfValues']=NonMetabolicX.notna().sum(axis=1).sort_values()


# In[104]:


Invariant=NonMetabolicX[NonMetabolicX['CVLowFraction']>=0.95]


# In[398]:


# NonMetabolicCVLow['CVLowFraction']=(NonMetabolicCVLow[0]/NonMetabolicCVLow.NumberOfValues)


# In[399]:


# NonMetabolicCVLow=NonMetabolicCVLow[NonMetabolicCVLow.CVLowFraction>=0.95]


# In[105]:


Invariant.to_csv("Bin1_NonMetabolic.csv")


# In[106]:


Invariant


# In[401]:


# MetabolicX=MetabolicX.loc[nonlowexpressed]


# In[403]:


# NonMetabolicX


# In[116]:


# CVLow=pd.DataFrame((MetabolicX<0.3).sum(axis=1))


# In[117]:


# CVLow=CVLow.loc[nonlowexpressed]


# In[118]:


# CVLow['NumberOfValues']=MetabolicX.notna().sum(axis=1).sort_values()


# In[119]:


# CVLow['CVLowFraction']=(CVLow[0]/CVLow.NumberOfValues)


# In[120]:


# CVLow.sort_values(ascending=False,by=['NumberOfValues'])


# In[121]:


# CVLow=CVLow[CVLow.CVLowFraction>=0.95]


# In[107]:


CVLow=Invariant


# In[108]:


CVLow['Class']=MetabolicClasses.Class


# In[109]:


CVLow['Class'].to_csv("Bin1_AllDatasets.csv")


# In[112]:


HighModCVHigh['Class']=MetabolicClasses['Class']


# In[113]:


HighModCVHigh


# In[114]:


MetabolicX=wb_to_gene(MetabolicX)


# In[115]:


MetabolicX=MetabolicX[~MetabolicX.index.str.contains("WB")]


# In[116]:


MetabolicX=gene_to_wb(MetabolicX)


# In[117]:


MetabolicX.to_csv("TotalMetabolic_AllDatasets.csv")


# In[118]:


ClassACVHigh=HighModCVHigh[HighModCVHigh.Class=='A']


# In[119]:


ClassACVHigh


# In[120]:


ClassACVHigh['Class']=MetabolicClasses['Class']


# In[121]:


ClassACVHigh=ClassACVHigh[ClassACVHigh.Class=='A']


# In[122]:


ClassACVHigh=ClassACVHigh[ClassACVHigh['Number of datasets with CV>=0.75']>2]


# In[123]:


ClassACVHigh


# In[149]:


MetabolicX


# In[150]:


MetabolicCVLow=MetabolicX[MetabolicX.CVLowFraction>=0.95]


# In[152]:


ClassAMetabolicCVLow=MetabolicCVLow[MetabolicCVLow.Class=='A']
OtherMetabolicCVLow=MetabolicCVLow[MetabolicCVLow.Class!='A']


# In[154]:


OtherMetabolicCVLow


# In[155]:


ModerateCV=MetabolicX[(MetabolicX['CVLowFraction'])<0.95]


# In[156]:


ModerateCV=ModerateCV[(ModerateCV['Number of datasets with CV>=0.75'])<=2]


# In[161]:


# MetabolicX[(MetabolicX['CVLowFraction'])>=0.95]


# In[158]:


ClassAModerateCV=ModerateCV[ModerateCV.Class=='A']
OtherMetabolicModerateCV=ModerateCV[ModerateCV.Class!='A']


# In[163]:


OtherMetabolicModerateCV


# In[164]:


MetabolicX['Class']=MetabolicClasses['Class']


# In[165]:


# ClassACVLow=CVLow[CVLow.Classes=='A']


# In[168]:


HighModCVHigh


# In[170]:


le=list(set(mapper_df.index).difference(set(HighModCV.index)))


# In[172]:


MetabolicLowExp=list(set(MetabolicClasses.index).intersection(set(le)))


# In[175]:


ClassOfMetabolicLowExp=MetabolicClasses.loc[MetabolicLowExp]


# In[176]:


ClassOfMetabolicLowExp


# In[127]:


# Plotting coefficient of variation
def PlotCoefVar(CoefVar,MetabolicClasses):
    MetabolicClasses=gene_to_wb(MetabolicClasses)
    fig, ax = plt.subplots(figsize=(6,5))
    AllGenes=CoefVar.Max_CV.hist(ax=ax,color='green',label='All genes',alpha=0.5)
    metabolic=list(set(MetabolicClasses.index).intersection(set(CoefVar.index)))
    MetabolicGenes=CoefVar.loc[metabolic].hist(ax=ax,color='midnightblue',bins=100,label='Metabolic genes',alpha=0.5)  
    ax.grid(False)
    plt.title("Distribution of Coefficient of variation")
    ax.set_xlabel('Coefficient of Variation')
    ax.set_ylabel('Number of genes')
#     for xc in xcoords:
#         ax.axvline(xc,color='blue',linestyle='--')
    plt.legend(loc='best')
    plt.savefig("CoefVar.png", dpi=300)
    plt.show()


# In[143]:


# Observing maximum CV of each gene
# Max.Max_CV.sort_values()[0:60]


# In[128]:


# Converting gene IDs to WormBase IDs for metabolic classes
MetabolicClasses=gene_to_wb(MetabolicClasses)


# In[129]:


# Finding intersection of metabolioc genes and genes in temperature datasets
metabolic=list(set(MetabolicClasses.index).intersection(set(Max.index)))
Max_Metabolic_CV=Max.loc[metabolic]


# In[62]:



# Max_metabolic=Max.loc[metabolic]


# In[384]:


Max_Metabolic_CV


# In[195]:


# Extracting only the Maximum CV of each gene
Max_CoefVar=Max['Max_CV']


# In[196]:


Max_Metabolic=Max_CoefVar.loc[metabolic]


# In[197]:


Max_Metabolic=pd.DataFrame(Max_Metabolic)


# In[198]:


# Extracting classes of genes only available in temperature dataset
Classes=MetabolicClasses.loc[metabolic]
Classes=Classes[~Classes.index.duplicated(keep='first')]


# In[199]:



Classes


# In[69]:


# Max_Metabolic.at['Class']=Classes['Class']


# In[200]:


Max_Metabolic=pd.DataFrame(Max_Metabolic)


# In[201]:


Max_Metabolic


# In[202]:


# Finding bin 1 genes that show least fluctuation in response to temperature stimulus (Max_CV<0.25)
Bin1=Max_Metabolic[Max_Metabolic.Max_CV<0.3]
Bin1.dropna(inplace=True)


# In[203]:


Bin1.to_csv("Bin1_AllDatasets.csv")


# In[204]:


Bin1


# In[205]:


Regulated_All=Max_Metabolic[Max_Metabolic['Max_CV']>=0.75]
Regulated_All.dropna(inplace=True)


# In[206]:


Regulated_All['Classes']=Classes['Class']


# In[385]:


Regulated_All


# In[78]:





# In[79]:


Bin2=Max_Metabolic[Max_Metabolic.Max_CV>=0.3]
Bin2=Bin2[Bin2<0.75]
Bin2.dropna(inplace=True)
Bin2.to_csv("Bin2_AllDatasets.csv")


# In[80]:


Bin1_gene=wb_to_gene(Bin1)


# In[81]:


Bin1


# In[82]:


Non_Regulated_Temperature=Max_Metabolic[Max_Metabolic['Max_CV']<1]


# In[83]:


Non_Regulated_Temperature['Classes']=Classes['Class']


# In[84]:


# Regulated_classes=Regulated_Temperature['Classes']
# Regulated_classes.to_csv("Regulated_Temperature.csv")
Non_Regulated_Temperatureclasses=Non_Regulated_Temperature['Classes']
Non_Regulated_Temperatureclasses.to_csv("Non_Regulated_Temperature.csv")


# In[85]:


Max_CoefVar=pd.DataFrame(Max_CoefVar)


# In[86]:


Max_CoefVar


# In[87]:


CV_nonzero=Max_CoefVar[Max_CoefVar!=0]


# In[88]:


Max_CoefVar.hist()


# In[89]:


PlotCoefVar(CoefVar=Max_CoefVar,MetabolicClasses=MetabolicClasses)


# In[90]:


Regulated_Temperature.shape


# In[ ]:


Bin2.dropna()


# In[ ]:


fig, ax = plt.subplots(figsize=(12,10))

size = 0.3
vals = [(Regulated_Temperature.shape[0]),(Bin2.shape[0]),(Bin1.shape[0])]
outer_labels=["Bin 3", "Bin 2", "Bin 1"]
print(vals)

explode = (0.1, 0)
cmap = plt.get_cmap("tab20c")
outer_colors = ['#54CB73','#F16718','#ff0000']
ax.pie(vals, radius=1, labels=outer_labels, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w',linewidth=3),autopct='%1.1f%%',shadow=False,
       textprops={'size': 'larger','fontweight':'bold'},
      pctdistance=0.85,labeldistance=1.05,startangle=90,counterclock=False)
plt.title("Percentage of transcriptionally regulated genes across temperatures",fontweight="bold",fontsize=12)


# In[ ]:


# ## Calculating highest Z
# count=0
# x=pd.DataFrame([])
# Maxz=pd.DataFrame([])
# for files in os.listdir('.'):
#     if fnmatch.fnmatch(files, "z_*pcl"):
#         count=count+1
#         print(files)
#         df=pd.read_csv(files,sep='\t',index_col=0,header='infer',low_memory=False)
#         Maxz=pd.DataFrame([])
#         df.drop(columns=['IDENTIFIER', 'GWEIGHT'],index=['EWEIGHT'],inplace=True)

#         df.replace("\\N",np.nan,inplace=True)
#         print(df.shape)
#         #         print(df.sum(axis=1).sort_values(ascending=False))
#         df=df.astype(float)
#         df=df.abs()
# #         print(df)
#         #         df[df<0]=np.nan
#         #         print(df.isna().sum(axis=1).sort_values())
#         #print(df)
#         MaxZ=df.max(skipna=True,axis=1)
#         Maxz['MaxZ_{}'.format(count)]=MaxZ
#         Maxz.index=df.index
# #         print(Maxz)
#         Maxz.to_csv("MaxZ_{}.csv".format(count))
#         x=pd.concat([x,Maxz],axis=1)
# #         print(Maxz.shape)


# In[ ]:





# In[ ]:


# Max=Max_Z(x)


# In[ ]:


# Max.Max_Z.hist(bins=20,grid=False)
# plt.xlabel("Max Z Score across temperature")
# plt.ylabel("Number of genes")


# In[ ]:


MetabolicClasses=pd.read_csv("/data/nandas/Combined_coexp/Part_1_TranscriptionallyRegulatedGenes/Tissue/MetabolicClasses_August_SN_082320.csv",index_col=0)


# In[ ]:


MetabolicClasses.shape


# In[ ]:


MetabolicClasses=gene_to_wb(MetabolicClasses)


# In[ ]:


# intersect=list(set(MetabolicClasses.index).intersection(set(Max.index)))


# In[ ]:


# Max_Metabolic=Max.loc[intersect]


# In[ ]:


# Max_Metabolic.Max_Z.hist(bins=100)


# In[ ]:


# Max_Metabolic_Zless2=Max_Metabolic[Max_Metabolic.Max_Z<=3]
# Max_Metabolic_Zmore2=Max_Metabolic[Max_Metabolic.Max_Z>3]


# In[ ]:


# Max_Metabolic_Zmore2.to_csv("Regulated_Temperature.csv")


# In[ ]:


# Max_Metabolic_Zless2.to_csv("Zless3Temperature.csv")


# In[ ]:


# Max_Metabolic['Class'].to_csv("TotalMetabolicTemperature.csv")


# In[ ]:


# Max_Metabolic


# In[ ]:


# df=pd.read_csv('z_normalized_WBPaper00037147.ce.mr.pcl',sep='\t',index_col=0,header='infer',low_memory=False)
# Maxz=pd.DataFrame([])
# df.drop(columns=['IDENTIFIER', 'GWEIGHT'],index=['EWEIGHT'],inplace=True)

# df.replace("\\N",np.nan,inplace=True)
# print(df)
# #         print(df.sum(axis=1).sort_values(ascending=False))
# df=df.astype(float)
# df=df.abs()
# print(df)
# #         df[df<0]=np.nan
# #         print(df.isna().sum(axis=1).sort_values())
# #print(df)
# MaxZ=df.max(skipna=True,axis=1)
# Maxz['MaxZ']=MaxZ
# Maxz.index=df.index
# print(Maxz)
# Maxz.to_csv("MaxZ_{}.csv".format(count))
# x=pd.concat([x,MaxZ],axis=1)
# print(Maxz.shape)


# In[ ]:



df=StackedDF(df)


# In[ ]:


# df[df<0]=np.nan


# In[ ]:


# df


# In[ ]:


df.hist(bins=100)


# In[ ]:


# y=df['Log2RPKM'][df.Log2RPKM>4]


# In[ ]:


# y.std()


# In[ ]:


x = df.values #returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
df = pd.DataFrame(x_scaled,index=df.index,columns=df.columns)


# In[ ]:


# df


# In[ ]:


# df.loc['WBGene00012942'].sort_values(ascending=False)


# In[ ]:


Std=df.std(axis=1,skipna=True)
Mean=df.mean(axis=1,skipna=True)


# In[ ]:


CoefVar=pd.DataFrame([])
#     from scipy.stats import variation 
#     Variation=variation(ZScore_df.values, axis = 1)
CoefVar['CoefVar_{}'.format(count)]=(Std/Mean)
CoefVar.index=df.index
CoefVar['Std']=Std
CoefVar['Mean']=Mean


# In[ ]:


CoefVar.CoefVar_7.sort_values()


#!/usr/bin/env python
# coding: utf-8

# ## Importing packages

# In[1]:


import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import to_hex
from matplotlib import gridspec
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from multiprocessing import Pool
# import xlrd
from dynamicTreeCut import cutreeHybrid
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
import scipy.cluster.hierarchy as spc
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import fnmatch
import graphviz
#import pygraphviz

# from plotnine import *


# ## Setting base directory

# In[2]:


Base_dir='/data/nandas/Combined_coexp_TFplusMetabolic/Pathway_centric/'
os.chdir(Base_dir)


# ## Reading files and matrices

# ### Reading combined coexpression files

# In[3]:


# output_df=pd.read_csv("../pearson_imputed_combined_total_z_normalised.dat",header=None,sep='\t')
# output_df_z=pd.read_csv("/data/nandas/Resolve_OR_genes/zpearson_imputed_combined_total_z_normalised.dat",
#                         header=None,sep='\t')


# In[4]:


output_matrix=pd.read_csv('/data/nandas/Combined_coexp_TFplusMetabolic/TF_metabol_matrix_genesymbol.csv',index_col=0,header='infer')


# ### Reading Gene to pathway file

# In[5]:


genes_df=pd.read_excel("PATHWAYS AND CATEGORIES APRIL 17 2020.xlsx",sheet_name='Gene2Pathway')


# ### Reading Pathway to Gene file

# In[6]:


pathway_df=pd.read_excel("PATHWAYS AND CATEGORIES APRIL 17 2020.xlsx",sheet_name='Pathway2Gene')


# ### Reading the list of Transciption factors

# In[7]:


TF=pd.read_csv("/data/nandas/Resolve_OR_genes/TF.csv",header=None,index_col=0)
TF.drop(index=['WBGene00021924','WBGene00001155'],inplace=True)


# ### Reading only metabolic genes

# In[8]:


metabolic_genes=pd.read_csv("/data/nandas/MEFIT/Combined/z_normalised/combined_imputed_total_corr_matrix.csv",header='infer',index_col=0)


# In[9]:


metabolic_genes_list=metabolic_genes.index


# ### Convert empty spaces if any to NaNs

# In[10]:


output_matrix=output_matrix.reindex()


# ### Check if there are any NaNs
# 

# In[11]:


output_matrix.columns[output_matrix.isnull().any()]


# In[12]:


output_matrix.columns[output_matrix.isnull().any()].tolist()


# In[13]:


# output_matrix = output_matrix[output_matrix.index.duplicated(keep='first')]


# In[14]:


missing=output_matrix[output_matrix.isnull()==True]


# ### Check any missing_zero values

# In[15]:


def missing_zero_values_table(df):
        zero_val = (df == 0.00).astype(int).sum(axis=0)
        mis_val = df.isnull().sum()
        mis_val_percent = 100 * df.isnull().sum() / len(df)
        mz_table = pd.concat([zero_val, mis_val, mis_val_percent], axis=1)
        mz_table = mz_table.rename(
        columns = {0 : 'Zero Values', 1 : 'Missing Values', 2 : '% of Total Values'})
        mz_table['Total Zero Missing Values'] = mz_table['Zero Values'] + mz_table['Missing Values']
        mz_table['% Total Zero Missing Values'] = 100 * mz_table['Total Zero Missing Values'] / len(df)
        mz_table['Data Type'] = df.dtypes
        mz_table = mz_table[
            mz_table.iloc[:,1] != 0].sort_values(
        '% of Total Values', ascending=False).round(1)
        print ("Your selected dataframe has " + str(df.shape[1]) + " columns and " + str(df.shape[0]) + " Rows.\n"      
            "There are " + str(mz_table.shape[0]) +
              " columns that have missing values.")
#         mz_table.to_excel('D:/sampledata/missing_and_zero_values.xlsx', freeze_panes=(1,0), index = False)
        return mz_table


# In[16]:


missing_zero_values_table(output_matrix)


# ### Check if there are duplicate indices

# In[17]:


output_matrix[output_matrix.index.duplicated()]


# ## Functions

# In[18]:


# Convert WBIDs to Gene symbol
def wb_to_gene(matrix):
    mapper_df=pd.read_csv("/data/nandas/MEFIT/predicted/mapper_final.csv", header='infer',index_col=0)
    wb_to_gene = {};
    for wb in mapper_df.index:
        wb_to_gene[wb] = mapper_df.loc[wb]['GeneID'];
    matrix=matrix.rename(index=wb_to_gene,columns=wb_to_gene)
    return matrix
# Return flat leaves for each node in linkage matrix
def get_node_leafs(Z):
    """
    For each node in a linkage matrix Z return the flat leafs
    :param Z: linkage matrix returned by scipy.cluster.hierarchy.linkage
    :return:
    """
    n = len(Z)
    leaf_nodes = [[x, ] for x in range((n + 1))]
    for i, node in enumerate(Z):
        id1, id2 = list(map(int, node[:2]))
        if np.all(node[:2] <= n):  # If all indices are leaves
            leaf_nodes.append([id1, id2])
        else:
            leaf_nodes.append([])
            if id1 <= n:
                leaf_nodes[-1].append(id1)
            else:
                leaf_nodes[-1].extend(leaf_nodes[id1])
            if id2 <= n:
                leaf_nodes[-1].append(id2)
            else:
                leaf_nodes[-1].extend(leaf_nodes[id2])
    return leaf_nodes

def find_merge(l1, Z):
    for j, node in enumerate(Z):
        if l1 in list(map(int, node[:2])):
            l2 = int([l for l in node[:2] if l !=l1][0])
            l3 = j+len(Z)+1
        d = node[2]
    return l1, l2, l3, d

def calculate_dist_matrix(matrix):
    print("---Calculating Dist Matrix----");
    dist=np.ones(matrix.shape)-matrix
    link = spc.linkage(squareform(dist), method='average')
    clusters = cutreeHybrid(link, squareform(dist), deepSplit=4,minClusterSize=5)
    return clusters, link

# Using networkx display the structure of every cluster
def display_the_gene_in_respective_cluster_or_subtree(matrix, gene_list, folder_name):
    clusters, link = calculate_dist_matrix(matrix)
    node_leaves=get_node_leafs(link)
    print("Here")
    for gene_name in gene_list:
        print(gene_name)
        if (gene_name in matrix.index):
            index = list(matrix.index).index(gene_name);
            cluster_label = clusters['labels'][index]

            print("outside : {} ".format(matrix.shape))
            ##Saving cluster gene list in a separate file for each gene
            get_cluster_gene_list(clusters, cluster_label, matrix, gene_name, folder_name);

            cluster1 = [x for x in range(matrix.shape[0]) if (clusters['labels']==cluster_label)[x]]
            G = nx.DiGraph()

            # Appending the graph
            for leaf in cluster1:
                l1, l2, l3,d = find_merge(leaf, link)
                if l3 not in G.nodes:
                    G.add_node(l3)
                G.add_node(l1)
                G.add_edge(l3, l1)
                if not len(set(cluster1).intersection(node_leaves[l2])):
                    continue
                else:
                    cluster1.append(l3)

            # Plotting the graph
            fig = plt.figure(figsize=(12,12))
            nx.nx_agraph.write_dot(G,'test.dot')

            # same layout using matplotlib with no labels
            plt.title('{}'.format(folder_name))
            pos=graphviz_layout(G, prog='dot')
            labels = dict([(i, gn) for i, gn in enumerate(matrix.index) if i in cluster1])
            text=nx.draw(G, pos, with_labels=False, arrows=False,node_size=500,node_color='#62CFB7')
            text=nx.draw_networkx_labels(G,pos,labels,font_size=14,font_weight='bold')
#             for _,t in text.items():
#                 t.set_rotation('vertical')
            plt.savefig("{}{}.png".format(folder_name, gene_name))
            plt.close()
            
#Get the list of genes in respective clusters    
def get_cluster_gene_list(clusters,cluster_label, matrix, gene_name, folder_name):
#     print("inside : {} ".format(matrix.shape));

    indices = [i for i, x in enumerate(clusters['labels']) if x == cluster_label]
    cluster_gene_list = matrix.index[indices];
    fp = open('{}{}_cluster_{}_gene_list.txt'.format(folder_name, gene_name, cluster_label), 'w')
#     print(cluster_gene_list)
    for gene in cluster_gene_list:
        fp.write(gene + "\n");
        if gene in TF.index:
            print ("file_name:{}_TF_{}".format(fp,gene))
    fp.close();


# In[19]:


# Converting Combined Coexpression matrix with WBIDs to Gene Symbols
output_matrix = wb_to_gene(output_matrix);
TF=wb_to_gene(TF)


# In[20]:


count = 0
#for tf in TF_final.index:
#     count = count + 1;
#     if(count == 5):
#        break;
for index in pathway_df.index:
    gene_list = pathway_df.loc[index].Name
    pathway = pathway_df.loc[index]['Categories/Pathways']
    pathway = pathway.replace('/', '_')
    pathway = pathway.replace(' - ', '_')
    pathway = pathway.replace(', ', '_')
    pathway = pathway.replace(' ', '_')
    print("Pathway is {}".format(pathway))

# for i in range(0,862):
#     tf = TF_final.index[i];
    gene_list = gene_list.split(';')
    print("Size of gene list is {}".format(len(gene_list)));
    print("\n------Gene list----\n")
    print(gene_list)
    drop_list = metabolic_genes.drop(index=gene_list,errors='ignore')
    drop_list = drop_list.index
    print("Size of drop list is {}".format(len(drop_list)))
    matrix = output_matrix.drop(index=drop_list, columns=drop_list, errors='ignore');
    file_exist = False;
    folder_name = '{}/'.format(pathway)
    print(folder_name)
    print(matrix.shape)

    if(os.path.exists(folder_name)==False):
        os.mkdir(folder_name);
#    for file in os.listdir(folder_name):
#         if fnmatch.fnmatch(file, "{}_cluster_*".format(gene)):
#             print("File: {} found, skipping!!!".format(file))
#             file_exist = True;
#             break;
#         if(not file_exist):
    display_the_gene_in_respective_cluster_or_subtree(matrix, gene_list, folder_name)

    
    


#!/usr/bin/env python
# coding: utf-8

# In[14]:


import os
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.colors import to_hex
from matplotlib import gridspec
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from multiprocessing import Pool
import mygene
# import xlrd
from goatools import go_search
from dynamicTreeCut import cutreeHybrid
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
from biokit.viz import corrplot
from biokit import corrplot as cp
import scipy.cluster.hierarchy as spc
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import fnmatch
# from plotnine import *


# In[15]:


# !pip install plotnine

# In[16]:


Base_dir='/data/nandas/Combined_coexp_TFplusMetabolic/Pathway_centric/'
os.chdir(Base_dir)

# In[17]:


output_df=pd.read_csv("../pearson_imputed_combined_total_z_normalised.dat",header=None,sep='\t')
output_df_z=pd.read_csv("/data/nandas/Resolve_OR_genes/zpearson_imputed_combined_total_z_normalised.dat",
                        header=None,sep='\t')
output_matrix=pd.read_csv('../pearson_matrix.csv',index_col=0,header='infer')

# In[18]:


# output_matrix = output_matrix[output_matrix.index.duplicated(keep='first')]

# ## Functions

# In[19]:


output_matrix

# In[40]:


def wb_to_gene(matrix):
    mapper_df=pd.read_csv("/data/nandas/MEFIT/predicted/mapper_final.csv", header='infer',index_col=0)
    wb_to_gene = {};
    for wb in mapper_df.index:
        wb_to_gene[wb] = mapper_df.loc[wb]['GeneID'];
    matrix=matrix.rename(index=wb_to_gene,columns=wb_to_gene)
    return matrix
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
    dist=np.ones(matrix.shape)-matrix
    link = spc.linkage(squareform(dist), method='average')
    print("--------Link----------")
    print(link)
    print("----------squareForm(dist)---------")
    print(dist)
    clusters = cutreeHybrid(link, squareform(dist), deepSplit=4,minClusterSize=6)
    return clusters, link
    
def display_the_gene_in_respective_cluster_or_subtree(matrix, gene_list, folder_name):
    clusters, link = calculate_dist_matrix(matrix)
    node_leaves=get_node_leafs(link)
    print("Here")
    for gene_name in gene_list:
        print(gene_name)
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
        plt.title('draw_networkx')
        pos=graphviz_layout(G, prog='dot')
        labels = dict([(i, gn) for i, gn in enumerate(matrix.index) if i in cluster1])
        text=nx.draw(G, pos, with_labels=False, arrows=False,node_size=500,node_color='r')
        text=nx.draw_networkx_labels(G,pos,labels,font_size=10,font_weight='bold',rotation='vertical')
        for _,t in text.items():
            t.set_rotation('vertical')
        plt.savefig("{}{}.png".format(folder_name, gene_name))
        plt.show()
    
def get_cluster_gene_list(clusters,cluster_label, matrix, gene_name, folder_name):
#     print("inside : {} ".format(matrix.shape));

    indices = [i for i, x in enumerate(clusters['labels']) if x == cluster_label]
    cluster_gene_list = matrix.index[indices];
    fp = open('{}{}_cluster_{}_gene_list.txt'.format(foler_name, gene_name, cluster_label), 'w')
#     print(cluster_gene_list)
    for gene in cluster_gene_list:
        fp.write(gene + "\n");
    fp.close();

# In[21]:


pathway_df=pd.read_excel("PATHWAYS AND CATEGORIES APRIL 17 2020.xlsx",sheet_name='Pathway2Gene')

# In[22]:


genes_df=pd.read_excel("PATHWAYS AND CATEGORIES APRIL 17 2020.xlsx",sheet_name='Gene2Pathway')

# In[23]:


genes_df


# In[24]:


pathway_df

# In[25]:


output_matrix = wb_to_gene(output_matrix);


# In[35]:


genes_df = genes_df.set_index('Gene')

# In[37]:


count = 0
#for tf in TF_final.index:
#     count = count + 1;
#     if(count == 5):
#        break;
for index in pathway_df.index:
    gene_list = pathway_df.loc[index].Name
    pathway = pathway_df.loc[index]['Categories/Pathways']
    print(pathway)
# for i in range(0,862):
#     tf = TF_final.index[i];
    gene_list = gene_list.split(';')
    print(gene_list)
    
    drop_list = genes_df.drop(index=gene_list)
    

    drop_list = drop_list.index
    matrix = output_matrix.drop(index=drop_list, columns=drop_list, errors='ignore');
    matrix.reindex;
    file_exist = False;
    folder_name = './{}/'.format(pathway)
    print(folder_name)

    if(os.path.exists(folder_name)==False):
        os.mkdir(folder_name);
    for file in os.listdir(folder_name):
        if fnmatch.fnmatch(file, "{}_cluster_*".format(gene)):
            print("File: {} found, skipping!!!".format(file))
            file_exist = True;
            break;
    if(not file_exist):
        display_the_gene_in_respective_cluster_or_subtree(matrix, gene_list, folder_name)
    break;
    

# In[ ]:



# In[ ]:


pathway_df=pathway_df[['Categories/Pathways','WBID']]

# In[ ]:


pathway_df

# In[ ]:


TF=pd.read_csv("/data/nandas/Resolve_OR_genes/TF.csv",header=None,index_col=0)
TF.drop(index=['WBGene00021924'],inplace=True)


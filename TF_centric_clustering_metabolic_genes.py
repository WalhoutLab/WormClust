#!/usr/bin/env python
# coding: utf-8

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
import mygene
# import xlrd

from dynamicTreeCut import cutreeHybrid
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
#from biokit.viz import corrplot
#from biokit import corrplot as cp
import scipy.cluster.hierarchy as spc
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import fnmatch


# In[2]:


Base_dir='/data/nandas/Combined_coexp_TFplusMetabolic/TF_centric_042920/'
os.chdir(Base_dir)


# In[3]:


TF_metabol_pairs=pd.read_csv("/data/nandas/Combined_coexp_TFplusMetabolic/pearson_imputed_combined_total_z_normalised.dat",header=None,sep='\t')


# In[18]:


def convertGenePairToGeneMatrix(output_df):
    output_df.set_axis(['Gene1','Gene2','weight'], axis=1,inplace=True)
    a = np.unique(output_df['Gene1'])
    b = np.unique(output_df['Gene2'])
    c = np.union1d(a,b);
    data = np.zeros((len(c), len(c)));
    output_matrix = pd.DataFrame(data, index=c, columns=c)
    for values in output_df.values: 
        output_matrix[values[0]][values[1]] = values[2];
        output_matrix[values[1]][values[0]]=values[2];
    np.fill_diagonal(output_matrix.values,1)
    return output_matrix

# Convert WBIDs to Gene symbol
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
    clusters = cutreeHybrid(link, squareform(dist), deepSplit=4,minClusterSize=5)
    return clusters, link
    
def display_the_gene_in_respective_cluster_or_subtree(matrix, gene_name):
    clusters, link = calculate_dist_matrix(matrix)
    node_leaves=get_node_leafs(link)
    
    index = list(matrix.index).index(gene_name);
    cluster_label = clusters['labels'][index]
    
    print("outside : {} ".format(matrix.shape))
    ##Saving cluster gene list in a separate file for each gene
    get_cluster_gene_list(clusters, cluster_label, matrix, gene_name);

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
    plt.title('{}'.format(gene_name))
    pos=graphviz_layout(G, prog='dot')
    labels = dict([(i, gn) for i, gn in enumerate(matrix.index) if i in cluster1])
    text=nx.draw(G, pos, with_labels=False, arrows=False,node_size=500,node_color='#62CFB7')
    text=nx.draw_networkx_labels(G,pos,labels,font_size=14,font_weight='bold')
#     for _,t in text.items():
#         t.set_rotation('vertical')
    plt.savefig("{}.png".format(gene_name))
    plt.show()
    
def get_cluster_gene_list(clusters,cluster_label, matrix, gene_name):
    print("inside : {} ".format(matrix.shape));

    indices = [i for i, x in enumerate(clusters['labels']) if x == cluster_label]
    cluster_gene_list = matrix.index[indices];
    fp = open('TF_clusters/{}_cluster_{}_gene_list.txt'.format(gene_name, cluster_label), 'w')
    print(cluster_gene_list)
    for gene in cluster_gene_list:
        fp.write(gene + "\n");
    fp.close();


# In[5]:


TF_metabol_matrix=convertGenePairToGeneMatrix(TF_metabol_pairs)


# In[ ]:


np.fill_diagonal(TF_metabol_matrix.values,1)


# In[ ]:


TF_metabol_matrix.to_csv("TF_metabol_matrix.csv")


# In[ ]:


TF_metabol_matrix=wb_to_gene(TF_metabol_matrix)


# In[28]:


TF_metabol_matrix.to_csv("/data/nandas/Combined_coexp_TFplusMetabolic/TF_metabol_matrix_genesymbol.csv")


# ## Reading the combined coexpression matrix and list of TFs

# In[29]:


TF_metabol_matrix=pd.read_csv("/data/nandas/Combined_coexp_TFplusMetabolic/TF_metabol_matrix_genesymbol.csv",index_col=0)


# In[30]:


TF_metabol_matrix.shape


# In[32]:


TF=pd.read_csv("/data/nandas/Resolve_OR_genes/TF.csv",header=None,index_col=0)
TF.drop(index=['WBGene00021924'],inplace=True)


# In[33]:


TF=wb_to_gene(TF)


# In[34]:


# dropped ech-6 from the list of TFs
TF.drop(index=['ech-6'],inplace=True)


# In[35]:


# Get the list of TFs whose coexpression values are present in the combined matrix
intersected_list = list(set(TF_metabol_matrix.index).intersection(set(TF.index)))
print(len(intersected_list))

TF_final = TF.loc[intersected_list]
# TF_final = TF_final.drop_duplicates()
TF_final.reset_index(inplace=True)
TF_final.drop_duplicates(inplace=True)
TF_final.set_index(0, inplace=True)
TF_final.shape


# In[43]:


set(TF_metabol_matrix.columns).symmetric_difference(TF_metabol_matrix.index)


# In[48]:


TF_metabol_matrix.rename(columns={'gei-3.1':'gei-3'},inplace=True)


# In[56]:


TF_metabol_matrix.drop_duplicates(keep='first',inplace=True)
TF_metabol_matrix.transpose().drop_duplicates(keep='first',inplace=True)


# In[60]:


count = 0
#for tf in TF_final.index:
#     count = count + 1;
#     if(count == 5):
#        break;
for tf in TF_final.index:
# for i in range(0,862):
#     tf = TF_final.index[i];
    print("tfvalue: {}".format(tf))
    drop_list = TF_final.drop(index=tf)
    drop_list = drop_list.index
    print(len(drop_list))
    print(TF_metabol_matrix.shape)
    matrix = TF_metabol_matrix.drop(index=drop_list, columns=drop_list);
    print(matrix.shape)
    file_exist = False;
    for file in os.listdir('./TF_clusters/'):
        if fnmatch.fnmatch(file, "{}_cluster_*".format(tf)):
            print("File: {} found, skipping!!!".format(file))
            file_exist = True;
    if(not file_exist):
        display_the_gene_in_respective_cluster_or_subtree(matrix, tf)


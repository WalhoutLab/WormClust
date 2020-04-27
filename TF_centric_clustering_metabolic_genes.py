#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

# In[2]:


Base_dir='/data/nandas/Combined_coexp_TFplusMetabolic/'
os.chdir(Base_dir)

# ## Reading combined coexpression of metabolic genes & TFs

# In[3]:


output_df=pd.read_csv("pearson_imputed_combined_total_z_normalised.dat",header=None,sep='\t')
output_df_z=pd.read_csv("/data/nandas/Resolve_OR_genes/zpearson_imputed_combined_total_z_normalised.dat",header=None,sep='\t')

# In[4]:


output_df

# ## Creating mapper dictionary WBID to gene_symbol

# In[5]:


mapper_df=pd.read_csv("/data/nandas/MEFIT/predicted/mapper_final.csv", header='infer',index_col=0)
mapper_df

wb_to_gene = {};
for wb in mapper_df.index:
    wb_to_gene[wb] = mapper_df.loc[wb]['GeneID'];
wb_to_gene

# ## Reading the list of Transcription Factors

# In[6]:


TF=pd.read_csv("/data/nandas/Resolve_OR_genes/TF.csv",header=None,index_col=0)

# In[7]:


TF.drop(index=['WBGene00021924'],inplace=True)

# In[8]:


TF.rename(index=wb_to_gene,inplace=True)

# In[9]:


output_df.set_axis(['Gene1','Gene2','weight'], axis=1,inplace=True)

# In[10]:


output_matrix=pd.read_csv("pearson_matrix.csv",header='infer',index_col=0)

# In[11]:


output_matrix.rename(columns=wb_to_gene,index=wb_to_gene,inplace=True)

# In[12]:


output_matrix['hphd-1'].sort_values(ascending=False)

# In[13]:


intersected_list = list(set(output_matrix.index).intersection(set(TF.index)))
print(len(intersected_list))

TF_final = TF.loc[intersected_list]
# TF_final = TF_final.drop_duplicates()
TF_final.reset_index(inplace=True)

# In[14]:


TF_final.drop_duplicates(inplace=True)

# In[15]:


TF_final

# In[16]:


TF_final.set_index(0, inplace=True)

# In[17]:


# def clustering(matrix):
#     print("Starting clustering ")
#     edist=np.ones(matrix.shape)-matrix
#     dist=edist
#     link = spc.linkage(squareform(dist), method='average')
#     clusters = cutreeHybrid(link, squareform(dist), deepSplit=3,minClusterSize=7)
#     #print(link)
#     #print(clusters)
#     print("Ending")
#     return clusters;
    
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
    clusters = cutreeHybrid(link, squareform(dist), deepSplit=4,minClusterSize=7)
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
    plt.title('draw_networkx')
    pos=graphviz_layout(G, prog='dot')
    labels = dict([(i, gn) for i, gn in enumerate(matrix.index) if i in cluster1])
    text=nx.draw(G, pos, with_labels=False, arrows=False,node_size=500,node_color='r')
    text=nx.draw_networkx_labels(G,pos,labels,font_size=10,font_weight='bold',rotation='vertical')
    for _,t in text.items():
        t.set_rotation('vertical')
    plt.savefig("{}.png".format(gene_name))
    plt.show()
    
def get_cluster_gene_list(clusters,cluster_label, matrix, gene_name):
    print("inside : {} ".format(matrix.shape));

    indices = [i for i, x in enumerate(clusters['labels']) if x == cluster_label]
    cluster_gene_list = matrix.index[indices];
    fp = open('TF_clusters/genes_list/{}_cluster_{}_gene_list.txt'.format(gene_name, cluster_label), 'w')
    print(cluster_gene_list)
    for gene in cluster_gene_list:
        fp.write(gene + "\n");
    fp.close();

# In[25]:




# In[24]:


# output_matrix

clusters, link = calculate_dist_matrix(output_matrix);

# In[19]:


# gene_name ='nhr-90'
# index = list(output_matrix.index).index(gene_name);
# cluster_label = clusters['labels'][index]
# get_cluster_gene_list(clusters, cluster_label, output_matrix, gene_name)
# for i in range(0,200):
#     tf = TF_final.index[i];
#     print("{}---{}".format(i,tf));


# In[20]:


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
    matrix = output_matrix.drop(index=drop_list, columns=drop_list);
    file_exist = False;
    for file in os.listdir('./TF_clusters/genes_list/'):
        if fnmatch.fnmatch(file, "{}_cluster_*".format(tf)):
            print("File: {} found, skipping!!!".format(file))
            file_exist = True;
    if(not file_exist):
        display_the_gene_in_respective_cluster_or_subtree(matrix, tf)
    

# In[21]:


# default = "#808080"   # Unclustered gray
# n_clusters = len(np.unique(clusters['labels']))

# colors = plt.cm.jet(np.linspace(0, 1, n_clusters)) # Create colors from colormap

# leaf_color_dict = {}

# for i, l in enumerate(clusters['labels']):
#     if l:
#         leaf_color_dict[i] = to_hex(colors[l-1])
#     else:
#         leaf_color_dict[i] = default
# # notes:
# # * rows in Z correspond to "inverted U" links that connect clusters
# # * rows are ordered by increasing distance
# # * if the colors of the connected clusters match, use that color for link
# link_cols = {}
# for i, i12 in enumerate(link[:,:2].astype(int)):
#     c1, c2 = (link_cols[x] if x > len(link) else leaf_color_dict[x] for x in i12)
#     link_cols[i+1+len(link)] = c1 if np.all(c1 == c2) else default

# In[22]:


# ## Plotting dendogram
# fig = plt.figure(figsize=(15, 8))
# gs  = gridspec.GridSpec(1, 2, width_ratios=[.5, 1])

# ax1 = plt.subplot(gs[0])
# ax2 = plt.subplot(gs[1])
# ax1.axis('off')
# dend=dendrogram(link,ax=ax1, orientation='left', no_labels=True, link_color_func=lambda x: link_cols[x])

# heatmap = ax2.pcolor(dist.values[dend['leaves']].T[dend['leaves']].T)
# ax1.yaxis.tick_right()

# #for s, si in zip(onecc, onecc_ndx): 
#   #  s=s.strip()
#   #  for i, l in enumerate(dend['leaves']):
        
#      #   if l == si:
#      #       if s == 'dao-3':
#          #       ax1.scatter(0.4,i*10,s=80, marker='>', color='blue', edgecolor='blue', linewidth=3, zorder=10)
#         #    elif s == 'metr-1':
#                # ax1.scatter(0.4,i*10,s=80, marker='>', color='midnightblue', edgecolor='midnightblue', linewidth=3, zorder=10)
#                 #ax.text(-1,i*10,'${}$'.format(s), color='Black', fontsize=28, zorder=10)
#           #  elif s == 'sams-5' or s == 'dhfr-1':
#             #    ax1.scatter(0.4,i*10,s=80, marker='>', color='red', edgecolor='red', linewidth=3, zorder=10)
#     #        elif s == 'sams-1'or s=='ahcy-1':
#      #           ax1.scatter(0.4,i*10,s=80, marker='>', color='purple', edgecolor='purple', linewidth=3, zorder=10)
#          #   else:
#        #         print(i,s)
#            #     ax1.scatter(0.4,i*10,s=100, marker='>', color='cyan', edgecolor='cyan', linewidth=3, zorder=10)
#                 #ax.text(-1,i*10,'${}$'.format(s), color='Black', fontsize=28, zorder=10)
                
                
# #Place textbox with cluster labels
# #for i, color in enumerate(colors):
#     #txt = ax1.text(0.0, 1-(0.1*i), 'cluster {}'.format(i+1), transform=ax.transAxes, 
#      #            **{'color': color, 'fontsize': 26})
#     #txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
# ax1.set_yticklabels([''])
# ax2.set_xticklabels([''])
# ax2.set_yticklabels([''])
# cbar = plt.colorbar(heatmap)
# cbar.set_label('Distance', fontsize=24)
# cbar.set_ticks(np.arange(0,1.3,0.1))
# cbar.set_ticklabels(['{:.1f}'.format(x) for x in np.arange(0,1.3,0.1)])
# cbar.ax.tick_params(labelsize=20)
# #plt.tight_layout()

# ax1.set_ylabel("Genes",fontsize=28)
# ax1.set_xlabel("Distance",fontsize=28)
# #xticklabels=['']*len(dend['leaves'])
# #ax.set_xticklabels(xticklabels,fontsize=28)
# gs.update(wspace=0.001, hspace=0.)
# #plt.tight_layout()
# plt.savefig('Combinedcoexp_TF_metabolic_032920.png')

# In[23]:



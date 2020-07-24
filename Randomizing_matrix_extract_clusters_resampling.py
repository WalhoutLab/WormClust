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
#from goatools import go_search
from dynamicTreeCut import cutreeHybrid
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
from biokit.viz import corrplot
from biokit import corrplot as cp
import scipy.cluster.hierarchy as spc
import networkx as nx
import seaborn as sns
import random
import fnmatch
import sys


# In[7]:


Base_dir='/data/nandas/Combined_coexp/RandomizeClusterExtraction/Resample/'
os.chdir(Base_dir)


# In[6]:


# !mkdir Resample


# In[8]:


output_df=pd.read_csv("/data/nandas/Combined_coexp/Pathway_enrichment/normalized_metaboliccorr_df.csv",header='infer',index_col=0)


# In[48]:


np.fill_diagonal(output_df.values,1)


# In[38]:


x=convertGenePairToGeneMatrix(x)


# In[44]:


x['acdh-1']['alh-8']


# In[50]:


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

def Resample(x):
    x=output_df.unstack()
    x.index.rename(['GeneA', 'GeneB'], inplace=True)
    x = x.to_frame('Correlation').reset_index()
    x.dropna(inplace=True)
    x = x[x.GeneA!=x.GeneB]
    y=x.Correlation.sample(n=x.shape[0],replace=True)
    y=pd.DataFrame(y);
    y.reset_index(inplace=True)
    x['Correlation']=y.Correlation
    x=convertGenePairToGeneMatrix(x)
    return x

def shuffle(metabolic_corr_df):
    a = metabolic_corr_df.index
    # print(a)
    # print(len(a))
    indices = random.sample(range(len(a)), len(a));
    metabolic_corr_df['NewCol']=indices
    metabolic_corr_df['NewIndices']=metabolic_corr_df.iloc[indices].index
    shuffled_df=metabolic_corr_df.set_index('NewIndices')
    shuffled_df.drop(columns=['NewCol'],inplace=True)
    # shuffled_df = metabolic_corr_df.iloc[indices,indices]
    # return shuffled_df
    shuffled_df=shuffled_df.transpose()
    shuffled_df['NewCol']=indices
    shuffled_df['NewIndices']=metabolic_corr_df.iloc[indices].index
    shuffled_df.set_index('NewIndices',inplace=True)
    shuffled_df.drop(columns=['NewCol'],inplace=True)
    return shuffled_df

def Clustering(output_df,method, deepSplit,MinClustSize):
    dist=np.ones(output_df.shape)-output_df
    link = spc.linkage(squareform(dist), method=method)
    clusters = cutreeHybrid(link, squareform(dist), deepSplit=deepSplit,minClusterSize = MinClustSize)
    return dist,link,clusters

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

def display_the_gene_in_respective_cluster_or_subtree(matrix, gene_name,link,clusters):
#     clusters, link = calculate_dist_matrix(matrix)
    node_leaves=get_node_leafs(link)
    
    index = list(matrix.index).index(gene_name);
    cluster_label = clusters['labels'][index]
    
#     print("outside : {} ".format(matrix.shape))
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
    plt.title('Cluster_{}'.format(cluster1))
    pos=graphviz_layout(G, prog='dot')
    labels = dict([(i, gn) for i, gn in enumerate(matrix.index) if i in cluster1])
    text=nx.draw(G, pos, with_labels=False, arrows=False,node_size=500,node_color='r')
    text=nx.draw_networkx_labels(G,pos,labels,font_size=10,font_weight='bold',rotation='vertical')
    for _,t in text.items():
        t.set_rotation('vertical')
    plt.savefig("{}.png".format(gene_name))
    plt.show()
    
def get_cluster_gene_list(clusters,cluster_label, matrix, gene_name):
#     print("inside : {} ".format(matrix.shape));

    indices = [i for i, x in enumerate(clusters['labels']) if x == cluster_label]
    cluster_gene_list = matrix.index[indices];
    
    file_exist = False;
    for file in os.listdir('./'):
        if fnmatch.fnmatch(file, "genes_list_{}/{}_cluster_{}_gene_list.txt".format(n,gene_name, cluster_label)):
            print("File: {} found, skipping!!!".format(file))
            file_exist = True;
    if(not file_exist):
        fp = open('genes_list_{}/{}_cluster_{}_gene_list.txt'.format(n,gene_name, cluster_label), 'w')
    #     print(cluster_gene_list)
        for gene in cluster_gene_list:
            fp.write(gene + "\n");
        fp.close();

def GetClusters(clusters):
    #np.unique(clusters['labels']).to_csv("Cluster_labels_{}".format(n))
    for cluster_id in np.unique(clusters['labels']): 
    #print("Running for Cluster: {}".format(cluster_id));
        indices = [i for i, x in enumerate(clusters['labels']) if x == cluster_id]
        print("Number of genes in cluster {} : {}".format(cluster_id, len(indices)))
        genes_list = output_df.index[indices];
        gene_filename = "ClusterDistribution_gene_{}/cluster_gene_list_{}.csv".format(n,str(cluster_id));
        fw_gene_list = open(gene_filename, 'w');
        for gene in genes_list:
            line = gene + "\n";
            fw_gene_list.write(line);
        fw_gene_list.close();
        #print(genes_list);
        file_name = "ClusterDistribution_gene_{}/cluster_{}.csv".format(n,str(cluster_id));
        fw = open(file_name, 'w');
        for i in range(len(genes_list)-1):
            for j in range(i+1, len(genes_list)):
                #print("{}--{}".format(i,j));
                if(genes_list[i] in output_df.index and genes_list[j] in output_df.index):
                    line = "{}\t{}\t{}\n".format(genes_list[i], genes_list[j], output_df[genes_list[i]][genes_list[j]])
                    fw.write(line);
        fw.close();
        
def PlotClusteredHeatMap(clusters,link, dist):
    default = "#808080"   # Unclustered gray
    n_clusters = len(np.unique(clusters['labels']))

    colors = plt.cm.jet(np.linspace(0, 1, n_clusters)) # Create colors from colormap

    leaf_color_dict = {}

    for i, l in enumerate(clusters['labels']):
        if l:
            leaf_color_dict[i] = to_hex(colors[l-1])
        else:
            leaf_color_dict[i] = default
    # notes:
    # * rows in Z correspond to "inverted U" links that connect clusters
    # * rows are ordered by increasing distance
    # * if the colors of the connected clusters match, use that color for link
    link_cols = {}
    for i, i12 in enumerate(link[:,:2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(link) else leaf_color_dict[x] for x in i12)
        link_cols[i+1+len(link)] = c1 if np.all(c1 == c2) else default
    fig = plt.figure(figsize=(15, 8))
    gs  = gridspec.GridSpec(1, 2, width_ratios=[.5, 1])

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax1.axis('off')
    dend=dendrogram(link,ax=ax1, orientation='left', no_labels=True, link_color_func=lambda x: link_cols[x])

    heatmap = ax2.pcolor(dist.values[dend['leaves']].T[dend['leaves']].T)
    ax1.yaxis.tick_right()

    #for s, si in zip(onecc, onecc_ndx): 
      #  s=s.strip()
      #  for i, l in enumerate(dend['leaves']):

         #   if l == si:
         #       if s == 'dao-3':
             #       ax1.scatter(0.4,i*10,s=80, marker='>', color='blue', edgecolor='blue', linewidth=3, zorder=10)
            #    elif s == 'metr-1':
                   # ax1.scatter(0.4,i*10,s=80, marker='>', color='midnightblue', edgecolor='midnightblue', linewidth=3, zorder=10)
                    #ax.text(-1,i*10,'${}$'.format(s), color='Black', fontsize=28, zorder=10)
              #  elif s == 'sams-5' or s == 'dhfr-1':
                #    ax1.scatter(0.4,i*10,s=80, marker='>', color='red', edgecolor='red', linewidth=3, zorder=10)
        #        elif s == 'sams-1'or s=='ahcy-1':
         #           ax1.scatter(0.4,i*10,s=80, marker='>', color='purple', edgecolor='purple', linewidth=3, zorder=10)
             #   else:
           #         print(i,s)
               #     ax1.scatter(0.4,i*10,s=100, marker='>', color='cyan', edgecolor='cyan', linewidth=3, zorder=10)
                    #ax.text(-1,i*10,'${}$'.format(s), color='Black', fontsize=28, zorder=10)


    #Place textbox with cluster labels
    #for i, color in enumerate(colors):
        #txt = ax1.text(0.0, 1-(0.1*i), 'cluster {}'.format(i+1), transform=ax.transAxes, 
         #            **{'color': color, 'fontsize': 26})
        #txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
    ax1.set_yticklabels([''])
    ax2.set_xticklabels([''])
    ax2.set_yticklabels([''])
    cbar = plt.colorbar(heatmap)
    cbar.set_label('Distance', fontsize=24)
    cbar.set_ticks(np.arange(0,1.3,0.1))
    cbar.set_ticklabels(['{:.1f}'.format(x) for x in np.arange(0,1.3,0.1)])
    cbar.ax.tick_params(labelsize=20)
    #plt.tight_layout()

    ax1.set_ylabel("Genes",fontsize=28)
    ax1.set_xlabel("Distance",fontsize=28)
    #xticklabels=['']*len(dend['leaves'])
    #ax.set_xticklabels(xticklabels,fontsize=28)
    gs.update(wspace=0.001, hspace=0.)
    #plt.tight_layout()
    plt.savefig('Combined_coexp_cluster_sizemin6_random{}.png'.format(n))


# In[52]:


#Calculating for original matrix
dist,link,clusters = Clustering(output_df,method='average',deepSplit=4,MinClustSize=2);
n=0;
os.mkdir("genes_list_{}".format(n));
os.mkdir("ClusterDistribution_gene_{}".format(n))
shunt = ['acdh-1','hphd-1','alh-8','hach-1','ech-6']
GetClusters(clusters);
for gene_name in shunt:
    print(gene_name)
    index = list(output_df.index).index(gene_name);
    cluster_label = clusters['labels'][index]
    get_cluster_gene_list(clusters, cluster_label, output_df, gene_name)
    break;


# In[10]:


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
# plt.savefig('Combined_coexp_cluster_sizemin6.png')


# In[53]:


NoOfClusters_original=len(np.unique(clusters['labels']))


# In[54]:


NoOfClusters_original


# In[ ]:


# PlotClusteredHeatMap(clusters=clusters,link=link,dist=dist)


# In[ ]:


n = sys.argv[0];
Genes_In_Same_Cluster = pd.DataFrame([])
Genes_In_Same_Cluster2=pd.DataFrame([])
Number_of_clusters=pd.DataFrame([])
New_df = pd.DataFrame([])
New_df2 = pd.DataFrame([])
New_df3=pd.DataFrame([])
for n in range(1,1001):
    output_df=Resample(output_df)
    print(output_df);
#     os.mkdir("genes_list_{}".format(n));
#     os.mkdir("ClusterDistribution_gene_{}".format(n))
    dist,link,clusters = Clustering(output_df,method='average',deepSplit=3,MinClustSize=2);
    NumberOfClusters=len(np.unique(clusters['labels']))
    print("No. of clusters:{}".format(NumberOfClusters));
    shunt = ['acdh-1','hphd-1','alh-8','hach-1','ech-6']
    onecc=['sams-1','sams-3','sams-4','sams-5','metr-1','mtrr-1','ahcy-1']
    GetClusters(clusters)
    #Check for Shunt
    common_cluster_id=-1; 
    count = 1;
    print(n)
    in_same_cluster = True;
    for gene_name in shunt:
        print(gene_name)
        index = list(output_df.index).index(gene_name);
        cluster_label = clusters['labels'][index]
        if(count == 1):
            common_cluster_id = cluster_label;
        else:
            if(common_cluster_id != cluster_label):
                in_same_cluster = False;
                break;
        count = count + 1;
        print(cluster_label)
        #get_cluster_gene_list(clusters, cluster_label, output_df, gene_name)
    print("InSameCluster: {}".format(in_same_cluster));
    print("Count: {}".format(count));
#     print(Genes_In_Same_Cluster)
    Number_of_clusters.at[n,'Number_Of_Clusters']=NumberOfClusters;
    print("Dataframe of no. of clusters: {}".format(Number_of_clusters));
    Genes_In_Same_Cluster.at[n,'In_Same_Cluster']=in_same_cluster;
    print(Genes_In_Same_Cluster)
    New_df = New_df.append(Genes_In_Same_Cluster)
    New_df3=New_df3.append(Number_of_clusters);
    #Check for MET/SAM cycle
    common_cluster_id=-1; 
    count2 = 1;
    in_same_cluster = True;
    for gene_name in onecc:
        print(gene_name)
        index = list(output_df.index).index(gene_name);
        cluster_label = clusters['labels'][index]
        if(count2 == 1):
            common_cluster_id = cluster_label;
        else:
            if(common_cluster_id != cluster_label):
                in_same_cluster = False;
                break;
        count2 = count2 + 1;
        print(cluster_label)
        #get_cluster_gene_list(clusters, cluster_label, output_df, gene_name)
    print("InSameCluster: {}".format(in_same_cluster));
    print("Count: {}".format(count));
#     Genes_In_Same_Cluster2 = pd.DataFrame(index=['In_Same_Cluster'])
    Genes_In_Same_Cluster2.at[n,'In_Same_Cluster']=in_same_cluster;
    New_df2=New_df2.append(Genes_In_Same_Cluster2)
    

Number_of_clusters.to_csv("Number_of_clusters.csv")   
New_df.to_csv('Gene_shunt_cluster_check.csv')
New_df2.to_csv('Gene_MetSAM_cluster_check.csv')



# In[ ]:


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


# In[ ]:


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
# plt.savefig('Combined_coexp_cluster_sizemin6.png')


# In[ ]:


dist.max().max()


# In[ ]:


# output_df['CLUSTERS']=clusters['labels']
# shunt = ['acdh-1','hphd-1','alh-8','hach-1','ech-6']
# shunt_ndx = []
# try:
#     for s in shunt:
#         shunt_ndx.append(output_df.index.get_loc(s))
# except:
#     print('could not find gene {}'.format(s))
# shunt_cluster = np.unique(clusters['labels'][shunt_ndx])
# shunt_cluster
# shnt = np.zeros((len(shunt),len(shunt)))
# for i, gn1 in enumerate(shunt):
#     print('gene: {}\tCluster: {}'.format(gn1, output_df.loc[gn1,'CLUSTERS']))
#     for j, gn2 in enumerate(shunt):
#         shnt[(i,j)] = output_df.loc[gn1,gn2]

# shunt_corr=pd.DataFrame(shnt,index=shunt,columns=shunt)
# #shunt_corr=shunt_corr.rename(index={'F13D12.4':'alh-8','C55B7.4':'acdh-1','F09F7.4':'hach-1',
#  #                                   'T05G5.6':'ech-6','Y38F1A.6':'hphd-1','F27D9.5':'pcca-1','F52E4.1':'pccb-1',
#   #                                  'D2030.5':'mce-1','ZK1058.1':'mmcm-1','LLC1.3':'dld-1'},
#    #                          columns={'F13D12.4':'alh-8','C55B7.4':'acdh-1','F09F7.4':'hach-1',
#     #                                'T05G5.6':'ech-6','Y38F1A.6':'hphd-1','F27D9.5':'pcca-1','F52E4.1':'pccb-1',
#      #                               'D2030.5':'mce-1','ZK1058.1':'mmcm-1','LLC1.3':'dld-1'})
# d=cp.Corrplot(shunt_corr)
# d.plot(lower='text',upper='square')
# plt.show()
# # plt.savefig("Shunt_canonical_corr.png")


# In[ ]:


# onecc= ['acdh-1','pyk-1','pyk-2','pck-1','pck-2','pck-3','gpi-1','R05F9.6','Y43F4B.5','hxk-1','hxk-2','hxk-3','C50D2.7','C01B4.6'
#         ,'Y19D10A.16','aldo-2','aldo-1','fbp-1','pfk-1.2','pfk-1.1','tpi-1','gpd-4','gpd-3','gpd-2','gpd-1','pgk-1',
#        'ipgm-1','enol-1','pdha-1','pdhb-1','dlat-2','dlat-1', 'dld-1','ldh-1','alh-5','alh-2','alh-4',
#         'alh-11','alh-12','alh-1','alh-2','alh-9','acs-19','alh-5','H24K24.3','D2063.1','sodh-1','sodh-2'
#         ,'sams-1','sams-3','sams-4','sams-5','sars-1']
# onecc_ndx = []
# try:
#     for s in onecc:
#         onecc_ndx.append(output_df.index.get_loc(s))
# except:
#     print('could not find gene {}'.format(s))
# #shunt_cluster = np.unique(clusters['labels'][shunt_ndx])
# #shunt_cluster
# oncc = np.zeros((len(onecc),len(onecc)))
# for i, gn1 in enumerate(onecc):
#     print('gene: {}\tCluster: {}'.format(gn1, output_df.loc[gn1,'CLUSTERS']))
#     for j, gn2 in enumerate(onecc):
#         oncc[(i,j)] = output_df.loc[gn1,gn2]

# onecc_corr=pd.DataFrame(oncc,index=onecc,columns=onecc)
# #shunt_corr=shunt_corr.rename(index={'F13D12.4':'alh-8','C55B7.4':'acdh-1','F09F7.4':'hach-1',
#  #                                   'T05G5.6':'ech-6','Y38F1A.6':'hphd-1','F27D9.5':'pcca-1','F52E4.1':'pccb-1',
#   #                                  'D2030.5':'mce-1','ZK1058.1':'mmcm-1','LLC1.3':'dld-1'},
#    #                          columns={'F13D12.4':'alh-8','C55B7.4':'acdh-1','F09F7.4':'hach-1',
#     #                                'T05G5.6':'ech-6','Y38F1A.6':'hphd-1','F27D9.5':'pcca-1','F52E4.1':'pccb-1',
#      #                               'D2030.5':'mce-1','ZK1058.1':'mmcm-1','LLC1.3':'dld-1'})
# d=cp.Corrplot(onecc_corr)
# d.plot(lower='text',upper='square')
# plt.savefig("Onecc_corr.png")


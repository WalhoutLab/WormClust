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
import seaborn as sns

# In[2]:


Base_dir='/data/nandas/MEFIT/Combined/z_normalised/'
os.chdir(Base_dir)

# In[3]:


output_df=pd.read_csv("combined_imputed_total_corr_matrix.csv",header='infer',index_col=0)

# In[4]:


output_df

# In[5]:


dist=np.ones(output_df.shape)-output_df

# In[6]:


output_df

# In[7]:


link = spc.linkage(squareform(dist), method='average')
clusters = cutreeHybrid(link, squareform(dist), deepSplit=4,minClusterSize = 6)
# clusters=spc.fcluster(link,t=60,criterion='maxclust')

# In[8]:


output_df.shape

# In[9]:


sns.set(color_codes=True)

sns.clustermap(output_df, method="average",  row_linkage=link,figsize=(200, 150))


# In[10]:


np.unique(clusters['labels'])

# In[11]:


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

# In[12]:


dist.max().max()

# In[13]:


## Plotting dendogram
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
plt.savefig('Combined_coexp_cluster_sizemin5.png')

# In[14]:


len(np.unique(clusters['labels']))

# In[15]:


link

# In[16]:


output_df.shape

# In[17]:


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

# In[18]:


node_leaves=get_node_leafs(link)

# In[19]:


node_leaves[3800]

# In[20]:


np.unique((clusters['labels']))

# In[21]:


link.shape

# In[22]:


import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

# In[23]:


# cluster1 = [x for x in range(output_df.shape[0]) if (clusters['labels']==37)[x]]
# G = nx.DiGraph()

# In[24]:


cluster1 = [x for x in range(output_df.shape[0]) if (clusters['labels']==62)[x]]
G = nx.DiGraph()

# In[25]:


len(np.unique(clusters['labels']))

# In[26]:


def find_merge(l1, Z):
    for j, node in enumerate(Z):
        if l1 in list(map(int, node[:2])):
            l2 = int([l for l in node[:2] if l !=l1][0])
            l3 = j+len(Z)+1
        d = node[2]
    return l1, l2, l3, d
            
            

# In[27]:



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

# In[28]:


fig = plt.figure(figsize=(15,12))
nx.nx_agraph.write_dot(G,'test.dot')

# same layout using matplotlib with no labels
plt.title('Cluster_1233')
pos=graphviz_layout(G, prog='dot')
labels = dict([(i, gn) for i, gn in enumerate(output_df.index) if i in cluster1])
text=nx.draw(G, pos, with_labels=False, arrows=False,node_size=500,node_color='#62CFB7')
text=nx.draw_networkx_labels(G,pos,labels,font_size=14,font_weight='bold')
# for _,t in text.items():
#     t.set_rotation('vertical')
plt.show()

# In[29]:


G.edges

# In[30]:


l = [1,2,3]
for i in l:
    print(i)
    if i == 2:
        l.append(4)

# In[30]:


output_df['CLUSTERS']=clusters['labels']
shunt = ['acdh-1','hphd-1','alh-8','hach-1','ech-6']
shunt_ndx = []
try:
    for s in shunt:
        shunt_ndx.append(output_df.index.get_loc(s))
except:
    print('could not find gene {}'.format(s))
shunt_cluster = np.unique(clusters['labels'][shunt_ndx])
shunt_cluster
shnt = np.zeros((len(shunt),len(shunt)))
for i, gn1 in enumerate(shunt):
    print('gene: {}\tCluster: {}'.format(gn1, output_df.loc[gn1,'CLUSTERS']))
    for j, gn2 in enumerate(shunt):
        shnt[(i,j)] = output_df.loc[gn1,gn2]

shunt_corr=pd.DataFrame(shnt,index=shunt,columns=shunt)
#shunt_corr=shunt_corr.rename(index={'F13D12.4':'alh-8','C55B7.4':'acdh-1','F09F7.4':'hach-1',
 #                                   'T05G5.6':'ech-6','Y38F1A.6':'hphd-1','F27D9.5':'pcca-1','F52E4.1':'pccb-1',
  #                                  'D2030.5':'mce-1','ZK1058.1':'mmcm-1','LLC1.3':'dld-1'},
   #                          columns={'F13D12.4':'alh-8','C55B7.4':'acdh-1','F09F7.4':'hach-1',
    #                                'T05G5.6':'ech-6','Y38F1A.6':'hphd-1','F27D9.5':'pcca-1','F52E4.1':'pccb-1',
     #                               'D2030.5':'mce-1','ZK1058.1':'mmcm-1','LLC1.3':'dld-1'})
d=cp.Corrplot(shunt_corr)
d.plot(lower='text',upper='square')
plt.show()
# plt.savefig("Shunt_canonical_corr.png")

# In[31]:


for cluster_id in np.unique(clusters['labels']): 
    #print("Running for Cluster: {}".format(cluster_id));
    indices = [i for i, x in enumerate(clusters['labels']) if x == cluster_id]
    print("Number of genes in cluster {} : {}".format(cluster_id, len(indices)))
    genes_list = output_df.index[indices];
    gene_filename = "ClusterDistribution_gene/cluster_gene_list_{}.csv".format(str(cluster_id));
    fw_gene_list = open(gene_filename, 'w');
    for gene in genes_list:
        line = gene + "\n";
        fw_gene_list.write(line);
    fw_gene_list.close();
    #print(genes_list);
    file_name = "ClusterDistribution_gene/cluster_{}.csv".format(str(cluster_id));
    fw = open(file_name, 'w');
    for i in range(len(genes_list)-1):
        for j in range(i+1, len(genes_list)):
            #print("{}--{}".format(i,j));
            if(genes_list[i] in output_df.index and genes_list[j] in output_df.index):
                line = "{}\t{}\t{}\n".format(genes_list[i], genes_list[j], output_df[genes_list[i]][genes_list[j]])
                fw.write(line);
    fw.close();

# In[32]:


clusters['labels']==0

# In[33]:


onecc= ['acdh-1','pyk-1','pyk-2','pck-1','pck-2','pck-3','gpi-1','R05F9.6','Y43F4B.5','hxk-1','hxk-2','hxk-3','C50D2.7','C01B4.6'
        ,'Y19D10A.16','aldo-2','aldo-1','fbp-1','pfk-1.2','pfk-1.1','tpi-1','gpd-4','gpd-3','gpd-2','gpd-1','pgk-1',
       'ipgm-1','enol-1','pdha-1','pdhb-1','dlat-2','dlat-1', 'dld-1','ldh-1','alh-5','alh-2','alh-4',
        'alh-11','alh-12','alh-1','alh-2','alh-9','acs-19','alh-5','H24K24.3','D2063.1','sodh-1','sodh-2'
        ,'sams-1','sams-3','sams-4','sams-5','sars-1']
onecc_ndx = []
try:
    for s in onecc:
        onecc_ndx.append(output_df.index.get_loc(s))
except:
    print('could not find gene {}'.format(s))
#shunt_cluster = np.unique(clusters['labels'][shunt_ndx])
#shunt_cluster
oncc = np.zeros((len(onecc),len(onecc)))
for i, gn1 in enumerate(onecc):
    print('gene: {}\tCluster: {}'.format(gn1, output_df.loc[gn1,'CLUSTERS']))
    for j, gn2 in enumerate(onecc):
        oncc[(i,j)] = output_df.loc[gn1,gn2]

onecc_corr=pd.DataFrame(oncc,index=onecc,columns=onecc)
#shunt_corr=shunt_corr.rename(index={'F13D12.4':'alh-8','C55B7.4':'acdh-1','F09F7.4':'hach-1',
 #                                   'T05G5.6':'ech-6','Y38F1A.6':'hphd-1','F27D9.5':'pcca-1','F52E4.1':'pccb-1',
  #                                  'D2030.5':'mce-1','ZK1058.1':'mmcm-1','LLC1.3':'dld-1'},
   #                          columns={'F13D12.4':'alh-8','C55B7.4':'acdh-1','F09F7.4':'hach-1',
    #                                'T05G5.6':'ech-6','Y38F1A.6':'hphd-1','F27D9.5':'pcca-1','F52E4.1':'pccb-1',
     #                               'D2030.5':'mce-1','ZK1058.1':'mmcm-1','LLC1.3':'dld-1'})
d=cp.Corrplot(onecc_corr)
d.plot(lower='text',upper='square')
plt.savefig("Onecc_corr.png")

# In[ ]:





# coding: utf-8

# ## Import modules

# In[1]:


import os
import numpy as np
import pandas as pd
import psutil
from multiprocessing import Process, Queue
from scipy.stats import pearsonr
from goatools import go_search
from dynamicTreeCut import cutreeHybrid
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
from biokit.viz import corrplot
from biokit import corrplot as cp
import scipy.cluster.hierarchy as spc
import matplotlib
import matplotlib.pyplot as plt
from dynamicTreeCut import cutreeHybrid
from dynamicTreeCut import *
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
from biokit.viz import corrplot
from biokit import corrplot as cp
import scipy.cluster.hierarchy as spc
from matplotlib.colors import to_hex
from matplotlib import gridspec
import seaborn as sns
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.cm as cm
from scipy.stats import norm
import statsmodels as st
from matplotlib.patches import Rectangle
# import statsmodels.api as sm


# ## Setting base directory

# In[2]:


Base_dir='/data/nandas/FinalFileShivani/Unsupervised/'
os.chdir(Base_dir)


# ## Reading files

# ### Reading rxn X rxn coflux matrix

# In[3]:


Rxn_Flux_Matrix=pd.read_csv("RxnFluxMatrix.csv",index_col=0)


# ### Sanity check of matrix

# In[5]:


Rxn_Flux_Matrix.max().max()


# In[6]:


Rxn_Flux_Matrix.loc['RM04432'].sort_values(ascending=False)[0:10]


# ## Important functions

# In[7]:


def SymmetricRxnFluxMatrix(Rxn_flux_matrix):
    """Description: This function takes a reaction flux matrix as input 
    and ensures that each pair of elements (i, j) and (j, i) in the matrix are equal, 
    by setting both elements to the maximum of the two. 
    Arguments: Rxn_flux_matrix: A pandas DataFrame representing the reaction flux matrix. 
    Returns: The modified reaction flux matrix, where each element (i, j) is equal to (j, i). 
    """
    for i in Rxn_flux_matrix.index:
        for j in Rxn_flux_matrix.columns:
            if Rxn_flux_matrix.loc[i][j]< Rxn_flux_matrix.loc[j][i]:
                Rxn_flux_matrix.loc[i][j]=Rxn_flux_matrix.loc[j][i]
            else:
                Rxn_flux_matrix.loc[j][i]=Rxn_flux_matrix.loc[i][j]
    return Rxn_flux_matrix

def wb_to_gene(matrix):
    """Description: Converts WormBase IDs to gene names in a given matrix. 
    Arguments: matrix: A pandas DataFrame with indices and columns representing WormBase IDs. 
    Returns: The DataFrame with indices and columns renamed to gene names. gene_to_wb(matrix)"""
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=1)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    wb_to_gene = {};
    for wb in mapper_df.index:
        wb_to_gene[wb] = str(mapper_df.loc[wb]['GeneName']);
    matrix=matrix.rename(index=wb_to_gene,columns=wb_to_gene)
    return matrix

def gene_to_wb(matrix):
    """Description: Converts gene names to WormBase IDs in a given matrix. 
    Arguments: matrix: A pandas DataFrame with indices and columns representing gene names. 
    Returns: The DataFrame with indices and columns renamed to WormBase IDs."""
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=2)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    gene_to_wb = {};
    for gene in mapper_df.index:
        gene_to_wb[gene] = str(mapper_df.loc[gene]['WormBaseID']);
    matrix=matrix.rename(index=gene_to_wb,columns=gene_to_wb)
    return matrix

def SeqToWB(output_df):
    """Description: Converts sequence identifiers to gene names in a matrix. 
    Arguments: matrix: A pandas DataFrame with indices and columns representing sequence identifiers. 
    Returns: The DataFrame with indices and columns renamed to gene names. GeneToSeq"""
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=3)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    Seq_to_Wb = {};
    mapper_df=mapper_df[mapper_df.index!=np.nan]
    for seq in mapper_df.index:
        Seq_to_Wb[seq] = str(mapper_df.loc[seq]['WormBaseID']);
    matrix=matrix.rename(index=Seq_to_Wb,columns=Seq_to_Wb)
    return matrix

def SeqToGene(matrix):
    """Description: Converts sequence identifiers to gene names in a matrix. 
    Arguments: matrix: A pandas DataFrame with indices and columns representing sequence identifiers. 
    Returns: The DataFrame with indices and columns renamed to gene names."""
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',index_col=3)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    Seq_to_Gene = {};
    mapper_df=mapper_df[mapper_df.index!=np.nan]
    for seq in mapper_df.index:
        Seq_to_Gene[seq] = str(mapper_df.loc[seq]['GeneName']);
    matrix=matrix.rename(index=Seq_to_Gene,columns=Seq_to_Gene)
    return matrix

def GeneToSeq(matrix):
    """Description: Converts sequence identifiers to gene names in a matrix. 
    Arguments: matrix: A pandas DataFrame with indices and columns representing sequence identifiers. 
    Returns: The DataFrame with indices and columns renamed to gene names."""
    mapper_df=pd.read_csv("/data/nandas/WormBase_282/MasterProteinCodingGenesAnnotation_WS282.csv", header='infer',
                          index_col=2)
    mapper_df=mapper_df.loc[mapper_df.index.dropna()]
    Gene_to_Seq = {};
    mapper_df=mapper_df[mapper_df.index!=np.nan]
    for gene in mapper_df.index:
        Gene_to_Seq[gene] = str(mapper_df.loc[gene]['SequenceID']);
    matrix=matrix.rename(index=Gene_to_Seq,columns=Gene_to_Seq)
    return matrix

def Reaction2Genedf(modeltable):
    """Description: Processes a model table to extract reaction to gene mappings, filtering out certain values. 
    Arguments: modeltable: A pandas DataFrame representing the model table. 
    Returns: Three DataFrames: one with all reactions and genes, one with and-logic genes, and one with final genes."""
    Rxn2Gene_df = modeltable[['ID', 'GENES']]
    Rxn2Gene_df = Rxn2Gene_df.dropna()
    Rxn2Gene_df = Rxn2Gene_df.set_index('ID')
    Rxn2Gene_df = Rxn2Gene_df[Rxn2Gene_df.GENES != 'ND']
    Rxn2Gene_df = Rxn2Gene_df[Rxn2Gene_df.GENES != 'TBD']
    Rxn2Gene_df = Rxn2Gene_df[Rxn2Gene_df.GENES != 'Unknown']
    Rxn2Gene_final=Rxn2Gene_df[~(Rxn2Gene_df.GENES.str.contains("\|",regex=True))]
    Rxn2Gene_and=Rxn2Gene_final[Rxn2Gene_final.GENES.str.contains("\&",regex=True)]
    return Rxn2Gene_df,Rxn2Gene_and,Rxn2Gene_final

def GenerateG2RMap(Rxn2Gene_df):
    """Description: Generates a gene to reaction map from a reaction to gene DataFrame. 
    Arguments: Rxn2Gene_df: A pandas DataFrame mapping reactions to genes. 
    Returns: A dictionary mapping genes to a list of associated reactions."""
    G2RMap = {};

    for reaction in Rxn2Gene_df.index:
        reaction = reaction.strip();
        #print(reaction)
        genes = Rxn2Gene_df.loc[reaction]['GENES'];
        genes = genes.replace('(', '');
        genes = genes.replace(')', '');
        genes = genes.replace('[', '');
        genes = genes.replace(']', '');
        genes = genes.replace('|', ',');
        genes = genes.replace('&', ',');
        genes = genes.split(',');
        for gene in genes: 
            gene = gene.strip();
    #         print("---{}---".format(reaction));
            if(gene in G2RMap):
                reactionIndexList = G2RMap[gene];
                reactionIndexList.append(reaction);
                G2RMap[gene] = reactionIndexList;
            else:
                G2RMap[gene] = [reaction];
    return G2RMap
            
#     print(G2RMap)
                                                                  
def CalculateGeneCoflux(G2RMap,Rxn_flux_matrix):
    """Description: Calculates gene co-flux for each pair of genes based on their associated reaction fluxes. 
    Arguments: G2RMap: A dictionary mapping genes to reactions. 
    Rxn_flux_matrix: A matrix of reaction fluxes. 
    Returns: A DataFrame representing the gene co-flux matrix. """
    gene_dist_df = pd.DataFrame(np.zeros((len(G2RMap), len(G2RMap))), index=list(G2RMap.keys()), columns=list(G2RMap.keys()))
    np.fill_diagonal(gene_dist_df.values, 1)
    genes_list = list(G2RMap.keys());
    for i in range(len(genes_list)-1):
        for j in range(i+1, len(genes_list)): 
            gene_i = genes_list[i];
            gene_j = genes_list[j];
            reactions_list_i = G2RMap[gene_i];
            reactions_list_j = G2RMap[gene_j];
            dist_i_j = set();
            for k in reactions_list_i: 
                for l in reactions_list_j:
                    #print('---{}---{}---'.format(k,l))
                    dist_i_j.add(Rxn_flux_matrix[k][l]);
            #print(dist_i_j);
            dist_i_j = list(dist_i_j);
            dist_i_j.sort();
            max_dist=dist_i_j[-1]
    #         if(dist_i_j[0] == 0 and len(dist_i_j) > 1):
    #             min_dist = dist_i_j[1]
    #         else:
    #             min_dist = dist_i_j[0];
            #print(min_dist)
            gene_dist_df[gene_i][gene_j] = max_dist;
            gene_dist_df[gene_j][gene_i] = max_dist;
    np.fill_diagonal(gene_dist_df.values, 1)
    return gene_dist_df

def calculate_dist_matrix(matrix,deepSplit,minClusterSize):
    """Description: Calculates a distance matrix for clustering, along with linkage and cluster data. 
    Arguments: matrix: The matrix to calculate distances from. deepSplit: Parameter for hierarchical clustering. minClusterSize: Minimum cluster size for clustering. 
    Returns: A tuple containing the distance matrix, linkage data, and cluster information. """
    dist=np.ones(matrix.shape)-matrix
    link = spc.linkage(squareform(dist), method='average')
    clusters = cutreeHybrid(link, squareform(dist), deepSplit=deepSplit,minClusterSize=minClusterSize)
    return dist,link,clusters

def get_node_leafs(Z):
    """Description: Retrieves the leaf nodes for each node in a hierarchical clustering linkage matrix. 
    Arguments: Z: Linkage matrix from hierarchical clustering. 
    Returns: A list of leaf nodes for each node in the linkage matrix. 
    """
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
    """Description: Identifies the merge operation in a linkage matrix that includes a specified leaf.
    Arguments: l1: A leaf node. Z: Linkage matrix. 
    Returns: A tuple representing the merge operation, including the merged nodes and distance. """
    for j, node in enumerate(Z):
        if l1 in list(map(int, node[:2])):
            l2 = int([l for l in node[:2] if l !=l1][0])
            l3 = j+len(Z)+1
        d = node[2]
    return l1, l2, l3, d
    
def display_the_gene_in_respective_cluster_or_subtree(matrix, gene_name,deepSplit,minClusterSize):
    
    """Description: Displays a dendrogram of the cluster or subtree where the specified gene is located. 
    Arguments: matrix: The gene matrix. gene_name: Name of the gene to locate. deepSplit, minClusterSize: Parameters for hierarchical clustering. 
    Returns: A dictionary mapping node indices to gene names."""
    dist, link, clusters = calculate_dist_matrix(matrix,deepSplit,minClusterSize)
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
    text=nx.draw_networkx_labels(G,pos,labels,font_size=11,font_weight='bold')
#     for _,t in text.items():
#         t.set_rotation('vertical')
    plt.savefig("{}_{}_{}_Dendrogram.png".format(gene_name,deepSplit,minClusterSize))
    plt.show()
    return labels
    
    
    
def get_cluster_gene_list(clusters,cluster_label, matrix, gene_name):
    """### Function: `clustersheatmap`

    **Description:** 
    Displays a dendrogram of the cluster or subtree where the specified gene is located. 

    **Arguments:** 
    - `matrix`: The gene matrix.
    - `gene_name`: Name of the gene to locate.
    - `deepSplit`, `minClusterSize`: Parameters for hierarchical clustering.

    **Returns:** 
    A dictionary mapping node indices to gene names.
    """
    print("inside : {} ".format(matrix.shape));

    indices = [i for i, x in enumerate(clusters['labels']) if x == cluster_label]
    cluster_gene_list = matrix.index[indices];
    fp = open('{}_cluster_{}_gene_list.txt'.format(gene_name, cluster_label), 'w')
    print(cluster_gene_list)
    for gene in cluster_gene_list:
        fp.write(gene + "\n");
    fp.close();
    
def ClusterList(clusters,output_df):
    """Description:
    Generates and saves files listing genes and their pairwise relationships for each cluster.

    Arguments:

    clusters: Cluster object containing cluster labels.
    output_df: Dataframe with gene data.
    Returns:
    None. Creates files for each cluster containing gene lists and pairwise relationships.
    """
    for cluster_id in np.unique(clusters['labels']): 
        #print("Running for Cluster: {}".format(cluster_id));
        indices = [i for i, x in enumerate(clusters['labels']) if x == cluster_id]
        print("Number of genes in cluster {} : {}".format(cluster_id, len(indices)))
        genes_list = output_df.index[indices];
        gene_filename = "cluster_gene_list_{}.csv".format(str(cluster_id));
        fw_gene_list = open(gene_filename, 'w');
        for gene in genes_list:
            line = gene + "\n";
            fw_gene_list.write(line);
        fw_gene_list.close();
        
        file_name = "cluster_{}.csv".format(str(cluster_id));
        fw = open(file_name, 'w');
        for i in range(len(genes_list)-1):
            for j in range(i+1, len(genes_list)):
                #print("{}--{}".format(i,j));
                if(genes_list[i] in output_df.index and genes_list[j] in output_df.index):
                    line = "{}\t{}\t{}\n".format(genes_list[i], genes_list[j], output_df[genes_list[i]][genes_list[j]])
                    fw.write(line);
        fw.close();

def SilhoutteSample(dist,labels,metric):
    """Description:
    Calculates silhouette scores for each sample in a dataset.

    Arguments:

        dist: Dataset for silhouette score computation.
        labels: Cluster labels for the dataset.
        metric: Metric for silhouette score calculation.
    Returns:
        Silhouette score values for each sample.
    """
    silhouette_score_values=silhouette_samples(dist, labels, metric)
    return silhouette_score_values

def Silhoutte_Values(output_df,silhouette_score_values,labels):
    """Description:
    Adds silhouette scores and cluster labels to a dataframe.

    Arguments:

    output_df: Dataframe to be augmented.
    silhouette_score_values: Array of silhouette scores.
    labels: Cluster labels for each sample.
    Returns:
    Modified dataframe with added silhouette scores and cluster labels."""
    x=output_df
    x['Silhoutte_Score']=silhouette_score_values;
    x['Cluster']=labels
    x=x[['Silhoutte_Score','Cluster']]
    return x

def clustersheatmap(genename,deepSplit,minClusterSize,matrix):
    """Description:
    Generates and saves a heatmap visualization of gene clusters, highlighting a specified gene.

    Arguments:

    genename: Gene to be highlighted.
    deepSplit, minClusterSize, matrix: Parameters for clustering and heatmap generation.
    Returns:
    None. Saves heatmap visualization in SVG and PNG formats.
    """
    genedict=display_the_gene_in_respective_cluster_or_subtree(deepSplit=deepSplit,gene_name=genename,
                                                               matrix=matrix,
                                                               minClusterSize=minClusterSize)
    genelist=list()
    for key in genedict.keys():
#         print("key is {}".format(key))
        gene=genedict[key]
        print("gene is {}".format(gene))
        genelist.append(gene)
    print("gene list is {}".format(genelist))
    fig=plt.figure(figsize=(2,2))
    aagr_cluster=genelist
    aagr=matrix.loc[aagr_cluster][aagr_cluster]
    xticks=aagr_cluster
    yticks=aagr_cluster
#     yticks = ['' for y in yticks]
#     xticks = ['' for x in xticks]
#     annot_kws={'fontsize':10, 
#            'fontstyle':'italic'}
#     g=sns.clustermap(aagr,method='average',
#                      cbar_kws={'label':'Coflux * Coexpression'},
# #                      vmin=0,vmax=0.8,cmap='YlOrRd',annot=True,fmt='.2g')
#     g=sns.clustermap(aagr,method='average',
#                      cbar_kws={'label':'Coflux * Coexpression'},
#                      vmin=0,vmax=0.8,cmap='YlOrRd',annot=True,fmt='.2g')
    g=sns.clustermap(aagr,method='average',
                     cbar_kws={'label':'Coflux * Coexpression'},
                     vmin=0,vmax=0.8,cmap='YlOrRd')
#     Reordered_matrix=aagr.iloc[clustergrid.dendrogram_row.reordered_ind]
#     Reordered_matrix.set_index(['index'],inplace=True)
#     Reordered_matrix=Reordered_matrix.loc[:,Reordered_matrix.index]
    aagr.to_csv("{}_{}_{}_cluster.csv".format(genename,minClusterSize,deepSplit))
    rotation = 90 
    for i, ax in enumerate(g.fig.axes):   ## getting all axes of the fig object
         ax.set_xticklabels(ax.get_xticklabels(), rotation = rotation, fontsize=12)
    rotation = 0 
    for i, ax in enumerate(g.fig.axes):   ## getting all axes of the fig object
         ax.set_yticklabels(ax.get_yticklabels(), rotation = rotation, fontsize=12)
#     plt.yticks(rotation=45) 
#     plt.xticks(rotation=0) 
#     g.set_yticklabels(g.get_yticklabels(), rotation = 90, fontsize = 8)
    plt.rcParams['svg.fonttype']='none';
    plt.rcParams["font.family"] = "Arial"
    N=len(aagr_cluster)
    wanted_label = genename
    columns=aagr_cluster
# wanted_row = np.where(np.array(columns)[g.dendrogram_row.reordered_ind] == wanted_label)[0]
# wanted_col = np.where(np.array(columns)[g.dendrogram_col.reordered_ind] == wanted_label)[0]
    wanted_row = g.dendrogram_row.reordered_ind.index(columns.index(wanted_label))
    wanted_col = g.dendrogram_col.reordered_ind.index(columns.index(wanted_label))

    xywh_row = (0, wanted_row, N, 1)
    xywh_col = (wanted_col, 0, 1, N)
    for x, y, w, h in (xywh_row, xywh_col):
        g.ax_heatmap.add_patch(Rectangle((x, y), w, h, fill=False, edgecolor='green', lw=4, clip_on=False))
    g.ax_heatmap.set_xticklabels([r'$\it{' + ticklabel.get_text().replace('_', '\\ ') + '}$'
                        for ticklabel in g.ax_heatmap.get_xticklabels()])
    g.ax_heatmap.set_yticklabels([r'$\it{' + ticklabel.get_text().replace('_', '\\ ') + '}$'
                        for ticklabel in g.ax_heatmap.get_yticklabels()])
    plt.savefig("{}_{}_{}_cluster.svg".format(genename,minClusterSize,deepSplit),dpi=600)
    plt.savefig("{}_{}_{}_cluster.png".format(genename,minClusterSize,deepSplit),dpi=600)
    plt.show()
    
def ConvertPairsToMatrix_SN(bayesian_metabol_df):
    """Description:
    Converts a dataframe of gene pairs with weights into a symmetric matrix.

    Arguments:

    bayesian_metabol_df: Dataframe with gene pairs and weights.
    Returns:
    Symmetric matrix representing gene relationships.
    """
    bayesian_metabol_df.set_axis(['Gene1','Gene2','weight'], axis=1,inplace=True)
    a = np.unique(bayesian_metabol_df['Gene1'])
    b = np.unique(bayesian_metabol_df['Gene2'])
    c = np.union1d(a,b);
    data = np.zeros((len(c), len(c)));
    output_df = pd.DataFrame(data, index=c, columns=c)
    for values in bayesian_metabol_df.values: 
        output_df[values[0]][values[1]] = values[2];
        output_df[values[1]][values[0]]=values[2];
    np.fill_diagonal(output_df.values,1)
    return output_df


# ### Ensuring Symmetry of Reaction Flux Matrix
# 
# #### Description:
# - This code segment ensures the symmetry of a reaction flux matrix.
# 
# #### Code Explanation:
# 
# 1. **SymmetricRxnFluxMatrix Function:**
#    - The `SymmetricRxnFluxMatrix` function takes a reaction flux matrix as input.
#    - It ensures that each pair of elements (i, j) and (j, i) in the matrix are equal by setting both elements to the maximum of the two.
#    - The function iterates through the matrix elements and checks if `Rxn_flux_matrix.loc[i][j]` is less than `Rxn_flux_matrix.loc[j][i]`. If true, it sets both elements to `Rxn_flux_matrix.loc[j][i]`. Otherwise, it sets both elements to `Rxn_flux_matrix.loc[i][j]`.
#    - The modified reaction flux matrix is returned, where each element (i, j) is equal to (j, i).
# 
# #### Objective of This Code Segment:
# - This code segment ensures that the reaction flux matrix is symmetric, which is important for certain calculations and analyses in systems biology.
# 
# 
# 

# In[ ]:


Rxn_Flux_Matrix=SymmetricRxnFluxMatrix(Rxn_Flux_Matrix)


# In[ ]:


# Rxn_Flux_Matrix.to_csv("/data/nandas/Coflux_matrix/FluxRed_051121/RxnFlux_063021.csv")


# In[8]:


#!cp /data/nandas/Coflux_matrix/FluxRed_061722/modelTable_wbid.xlsx .


# ### Reading iCEL1314 model table

# In[9]:


modeltable=pd.read_excel("modelTable_wbid.xlsx")


# In[ ]:


# modeltable.loc['EX00081']


# ### Processing and Intersecting DataFrames
# 
# This section of the code involves processing a model table for gene and reaction data and identifying the intersection between two datasets. The code and its functionality are explained below:
# 
# ```python
# Rxn2Gene_df, Rxn2Gene_and, Rxn2Gene_final = Reaction2Genedf(modeltable)
# intersect = list(set(Rxn2Gene_df.index).intersection(set(Rxn_Flux_Matrix.index)))
# 
# len(intersect)
# Rxn2Gene_df = Rxn2Gene_df.loc[intersect]
# Rxn2Gene_df.shape
# ```
#  The main objective of this snippet is to align and filter the reaction and gene data based on common elements in two different datasets.
#  The final shape of Rxn2Gene_df indicates how many reactions (rows) and attributes (columns) are present in the dataset after filtering, setting the stage for further detailed analysis.
# 

# In[10]:


Rxn2Gene_df,Rxn2Gene_and,Rxn2Gene_final=Reaction2Genedf(modeltable)
intersect = list(set(Rxn2Gene_df.index).intersection(set(Rxn_Flux_Matrix.index)))

len(intersect)
Rxn2Gene_df = Rxn2Gene_df.loc[intersect];
Rxn2Gene_df.shape


# In[ ]:


G2RMap=GenerateG2RMap(Rxn2Gene_df)


# In[ ]:


Gene_Coflux=CalculateGeneCoflux(G2RMap,Rxn_flux_matrix=Rxn_Flux_Matrix)


# In[ ]:


# Gene_Coflux.to_csv("Gene_Coflux_062622.csv")


# In[13]:


## Read if already calculated
Gene_Coflux=pd.read_csv("Gene_Coflux_062622.csv",index_col=0)


# In[14]:


Gene_Coflux=wb_to_gene(Gene_Coflux)


# In[15]:


Gene_Coflux=SeqToGene(Gene_Coflux)


# In[16]:


Gene_Coflux.loc['haly-1']['cpin-1']


# In[17]:


Gene_Coflux


# ### Reading coexpression matrix

# In[18]:


MetabolicCoexp=pd.read_csv("MetabolicCorrMatrix_083120.csv",
                           header='infer',index_col=0)


# In[19]:


Gene_Coexp=MetabolicCoexp


# In[ ]:


# Gene_Coexp=(Gene_Coexp*2)-1


# In[20]:


Gene_Coexp=SeqToGene(Gene_Coexp)
Gene_Coexp=wb_to_gene(Gene_Coexp)


# ### Gene Data Intersection
# 
# The following code snippet focuses on finding the common elements (genes) between two datasets and updating the matrices based on this shared subset. Here is a step-by-step explanation of each part of the code:
# 
# ```python
# intersect2 = list(set(Gene_Coexp.index).intersection(set(Gene_Coflux.index)))
# FinalFluxMatrix = Gene_Coflux.loc[intersect2][intersect2]
# FinalCoexp = Gene_Coexp.loc[intersect2][intersect2]```
# 

# In[21]:


intersect2=list(set(Gene_Coexp.index).intersection(set(Gene_Coflux.index)))
FinalFluxMatrix=Gene_Coflux.loc[intersect2][intersect2]
FinalCoexp=Gene_Coexp.loc[intersect2][intersect2]


# ### Optional Step: Transforming `FinalCoexp` Data
# 
# #### Functionality:
# - This line of code transforms the data in `FinalCoexp`.
# - It scales the original data, which should range between 0 and 1, to a new range of -1 to +1.
# - The formula `(value * 2) - 1` is a standard method for scaling data from a [0, 1] range to a [-1, 1] range.
# 
# #### When to Use This Step:
# 
# ##### Optional Nature:
# - It's important to note that this step should only be applied if your original data is within a 0 to 1 range.
# - If your data already spans from -1 to +1, applying this transformation would be inappropriate as it could distort the data.
# 

# In[ ]:


# FinalFluxMatrix=wb_to_gene(FinalFluxMatrix)
# FinalCoexp=wb_to_gene(FinalCoexp)


# In[22]:


FinalCoexp=(FinalCoexp*2)-1


# ### Code Description
# 
# **Purpose:** 
# Modifies a given dataframe `FinalCoexp` by setting negative values to zero and then saves the updated dataframe to a CSV file.
# 
# **Process:**
# 1. Initializes a counter `count` to zero.
# 2. Iterates through each index (row) of the dataframe `FinalCoexp`.
# 3. Increments `count` by 1 for each row and prints the current count.
# 4. Within each row, iterates through each column.
# 5. Checks if the value at the current row and column is less than zero.
# 6. If the value is negative, it is set to zero.
# 7. After completing the iterations, saves the modified dataframe to a CSV file named "ClassAExpressionForProductMatrix_062722_2_5.csv".
# 
# 

# In[23]:


count=0
for i in FinalCoexp.index:
    count=count+1
    print(count)
    for j in FinalCoexp.columns:
#         print("i={},j={}".format(i,j))
        if (FinalCoexp[i][j])<0:
            FinalCoexp[i][j]=0
FinalCoexp.to_csv("ClassAExpressionForProductMatrix_062722_2_5.csv")


# In[ ]:


for i in FinalFluxMatrix.index:
    for j in FinalFluxMatrix.columns:
        print("i={},j={}".format(i,j))
        if (FinalFluxMatrix[i][j])<0:
            FinalFluxMatrix[i][j]=0
FinalFluxMatrix.to_csv("ClassACofluxForProductMatrix_100921.csv")


# In[ ]:


ProductMatrix=FinalCoexp*FinalFluxMatrix


# In[ ]:


ProductMatrix.to_csv("ProductMatrix_062822.csv")


# In[ ]:


ProductMatrix=pd.read_csv("ProductMatrix_062822.csv",index_col=0)


# In[ ]:


ProductMatrix=SymmetricRxnFluxMatrix(ProductMatrix)


# In[ ]:


ProductMatrix.replace(np.nan,0,inplace=True)
ProductMatrix.replace(np.inf,0,inplace=True)


# In[ ]:


np.fill_diagonal(ProductMatrix.values,1)
np.fill_diagonal(FinalCoexp.values,1)
np.fill_diagonal(FinalFluxMatrix.values,1)


# ## ### Function Call: `clustersheatmap`
# 
# **Description:** 
# This  describes a function call to `clustersheatmap`, which is used to generate and display a heatmap of gene clusters. The heatmap specifically focuses on the cluster or subtree where the specified gene  is located.
# 
# **Function Call Details:**
# ```python
# clustersheatmap(genename='ard-1', deepSplit=3, minClusterSize=6, matrix=ProductMatrix)
# 

# In[ ]:


clustersheatmap(genename='ard-1',deepSplit=3,minClusterSize=6,matrix=ProductMatrix)


# In[ ]:


clustersheatmap(genename='sams-3',deepSplit=3,minClusterSize=6,matrix=ProductMatrix)


# In[ ]:


clustersheatmap(genename='sams-3',deepSplit=2,minClusterSize=3,matrix=ProductMatrix)


# In[ ]:


clustersheatmap(genename='pmt-2',deepSplit=2,minClusterSize=3,matrix=ProductMatrix)


# In[ ]:


clustersheatmap(genename='pmt-1',deepSplit=2,minClusterSize=3,matrix=ProductMatrix)


# In[ ]:


clustersheatmap(genename='metr-1',deepSplit=2,minClusterSize=3,matrix=ProductMatrix)


# In[ ]:


clustersheatmap(genename='cbs-1',deepSplit=2,minClusterSize=3,matrix=ProductMatrix)


# In[ ]:


clustersheatmap(genename='cbs-2',deepSplit=2,minClusterSize=3,matrix=ProductMatrix)


# In[ ]:


clustersheatmap(genename='T13G4.4',deepSplit=2,minClusterSize=3,matrix=ProductMatrix)


# In[ ]:


clustersheatmap(genename='cbs-1',deepSplit=3,minClusterSize=6,matrix=ProductMatrix)


# In[ ]:


# def CombinedClusterList(dee)


# In[ ]:


ProductMatrix=gene_to_wb(ProductMatrix)
FinalCoexp=gene_to_wb(FinalCoexp)
FinalFluxMatrix=gene_to_wb(FinalFluxMatrix)


# In[ ]:


ProductMatrix.to_csv("ProductMatrix_WB_062722.csv")


# In[ ]:


ProductMatrix=pd.read_csv("ProductMatrix_WB_062722.csv",index_col=0)


# In[ ]:


# ProductMatrix=wb_to_gene(ProductMatrix)


# ## Clustering of Product Matrix

# ## Clustering Analysis Explanation
# 
# The clustering analysis in the Jupyter Notebook is carried out using a function `calculate_dist_matrix`, which is a key part of dynamic tree clustering. This function is used to generate a hierarchical clustering of data based on certain parameters. Below is an explanation of the code snippet and its significant arguments:
# 
# ```python
# dist, link, clusters = calculate_dist_matrix(deepSplit=3, matrix=ProductMatrix, minClusterSize=6)
# labels = clusters['labels']
# ClusterList(clusters=clusters, output_df=ProductMatrix)
# ```
# ### Arguments Explanation
# 
# #### `deepSplit`
# - **Type:** `int`
# - **Description:** This argument controls the depth of the tree in hierarchical clustering. A higher value leads to a deeper tree, potentially resulting in more distinct clusters.
# - **Usage in Code:** ### Understanding the `deepSplit` Parameter in Hierarchical Clustering
# 
# The `deepSplit` parameter is an integral part of hierarchical clustering algorithms, particularly when using dynamic tree cutting. It controls how deeply the tree is split while forming clusters. Here's an overview of the different `deepSplit` values and their implications:
# 
# ##### deepSplit = 0
# - **Description:** Minimal or no additional splitting, resulting in larger clusters.
# - **Use Case:** Best when over-segmentation is a concern and larger, more generalized clusters are preferred.
# 
# ##### deepSplit = 1
# - **Description:** Slightly more aggressive splitting compared to 0, revealing more structure without much granularity.
# - **Use Case:** Suitable for an initial understanding of the cluster groupings in the data.
# 
# ##### deepSplit = 2
# - **Description:** A balanced approach that allows for moderate splitting, unveiling distinct groups within the data.
# - **Use Case:** Optimal for datasets with some known or expected sub-groupings.
# 
# ##### deepSplit = 3
# - **Description:** More aggressive splitting, leading to smaller and more specific clusters.
# - **Use Case:** Ideal for datasets known to contain nuanced or subtle groupings.
# 
# ##### deepSplit = 4
# - **Description:** The most aggressive splitting within the typical allowed range, resulting in the finest clustering.
# - **Use Case:** Best for complex datasets with many potential subgroups, where detailed exploration is necessary.
# 
# ##### Key Notes:
# - The choice of `deepSplit` value depends on the dataset's nature and the specific objectives of the clustering analysis.
# - Higher `deepSplit` values increase granularity but may lead to over-segmentation.
# - Always consult the documentation of the specific hierarchical clustering tool or library for precise guidelines, as implementations of `deepSplit` can vary.
# 
# Selecting the appropriate `deepSplit` value is crucial for achieving meaningful clustering results, with the range typically being between 0 and 4.
# 
# 
# #### `ProductMatrix`
# - **Type:** `DataFrame`
# - **Description:** Represents the data to be clustered, containing relationships or distances between different data points.
# - **Usage in Code:** `ProductMatrix` is used as the input matrix for the clustering function, providing the necessary data structure for cluster analysis.
# 
# #### `minClusterSize`
# - **Type:** `int`
# - **Description:** Determines the minimum number of elements required to form a valid cluster.
# - **Usage in Code:** `minClusterSize=6` ensures that any formed cluster will contain at least six elements, avoiding overly small or insignificant clusters.
# 
# ### Function Output
# 
# - **`dist`**: A matrix representing the pairwise distances between elements in the dataset.
# - **`link`**: The linkage matrix showing how different clusters are linked in the hierarchical tree.
# - **`clusters`**: An object containing detailed information about the clusters formed, including their labels and sizes.
# - **`labels`**: An array extracted from `clusters['labels']`, assigning a unique cluster label to each data point in `ProductMatrix`.
# 
# Following the computation of clusters, the function `ClusterList` is invoked. It processes the clustering results and lists out the elements in each cluster, utilizing both the clustering information and the original data in `ProductMatrix`.
# 
# This method of clustering is particularly valuable for unsupervised learning scenarios where the underlying structure of the data is unknown and needs to be discovered through analysis.
# 

# In[ ]:


dist,link,clusters=calculate_dist_matrix(deepSplit=3,matrix=ProductMatrix,minClusterSize=6)
labels=clusters['labels']
ClusterList(clusters=clusters,output_df=ProductMatrix)


# ## Finding list of genes in each cluster

# In[ ]:


Clusters=pd.DataFrame([])
matrix=ProductMatrix
Clusters['Genes']=matrix.index
Clusters['Cluster']=clusters['labels']
Cluster_list={}
labels=np.unique(clusters['labels'])
for label in labels:
    print(label)
    Cluster_list[label]=list(Clusters[Clusters.Cluster==label]['Genes']);
Cluster_list=pd.DataFrame.from_dict(Cluster_list,orient='index')
Cluster_list.reset_index(inplace=True)
Cluster_list.drop(columns=['index'],inplace=True)


# ### Clustering Analysis Steps
# 
# #### Initialization of Parameters:
# - `deepSplit` and `minClusterSize` are set to 3 and 6, respectively. These values are crucial for defining the clustering behavior.
# - `deepSplit` controls the depth of the clustering tree, impacting how granular the clusters are.
# - `minClusterSize` specifies the minimum number of elements required for a cluster, ensuring that clusters are of a significant size.
# 
# #### Iterating Over Clusters:
# - The loop iterates over each index in `Cluster_list`, which contains information about clustered genes.
# - For each cluster, a new DataFrame, `FinalCluster`, is initialized to store combined data relevant to that cluster.
# 
# #### Processing Cluster Genes:
# - Genes within the current cluster are extracted and stored in the `genes` list.
# - Three sub-matrices – `ProductCluster`, `CoexpCluster`, and `CofluxCluster` – are created for each gene cluster. These matrices contain data from `ProductMatrix`, `FinalCoexp`, and `FinalFluxMatrix`, respectively.
# 
# #### Combining Data:
# - For each pair of genes in the cluster, the corresponding values from `CofluxCluster`, `CoexpCluster`, and `ProductCluster` are concatenated into a single string. This composite string is then stored in the `FinalCluster` DataFrame.
# - This process effectively merges information from different matrices, providing a consolidated view of data for each cluster.
# 
# #### Saving the Final Cluster Data:
# - Each `FinalCluster` DataFrame is saved as a CSV file. The naming convention for these files incorporates the cluster index, `minClusterSize`, and `deepSplit` values.
# - The files are saved in the specified directory, in this case, `/data/nandas/WormClust/product/`, facilitating easy access and organization of the cluster-specific data.
# 
# ### Objective of This Code Segment:
# - The primary aim is to synthesize and organize data from multiple sources into distinct cluster-focused datasets.
# - This approach allows for in-depth analysis of each gene cluster, utilizing a comprehensive dataset that encapsulates various aspects of gene behavior and interactions.
# 

# In[ ]:


deepSplit=3
minClusterSize=6
for index in Cluster_list.index:
    FinalCluster=pd.DataFrame([])
#     index=49
    print(index)
    genes=list((Cluster_list.loc[index]).dropna())
    print(genes)
    ProductCluster=ProductMatrix.loc[genes][genes]
    CoexpCluster=FinalCoexp.loc[genes][genes]
    CofluxCluster=FinalFluxMatrix.loc[genes][genes]
    for i in CoexpCluster.index:
        for j in CoexpCluster.columns:
#             print("{}{}".format(i,j))
#             print(j)
            FinalCluster.at[i,j]=str(CofluxCluster.loc[i][j])+"_"+str(CoexpCluster.loc[i][j])+"_"+str(ProductCluster.loc[i][j])
#     FinalCluster=gene_to_wb(FinalCluster)
    FinalCluster.to_csv("/data/nandas/WormClust/product/product{}_{}_{}.csv".format(index,minClusterSize,deepSplit))
#     break;
#     print(ProductCluster)
#     break;


# In[ ]:


ProductMatrix.replace(np.nan,0,inplace=True)
ProductMatrix.replace(np.inf,0,inplace=True)


# In[ ]:


np.fill_diagonal(ProductMatrix.values,1)


# In[ ]:


ProductMatrix=wb_to_gene(ProductMatrix)


# ### Clustering of product matrix

# In[ ]:


dist,link,clusters=calculate_dist_matrix(deepSplit=3,matrix=ProductMatrix,minClusterSize=6)


# In[ ]:


FinalCoexp.index=FinalCoexp.index.str.strip()


# In[ ]:


FinalCoexp=FinalCoexp.transpose()


# In[ ]:


FinalCoexp


# In[ ]:


labels=clusters['labels']


# In[ ]:


len(labels)


# In[ ]:


Number_of_clusters=np.unique(clusters['labels'])


# In[ ]:


ClusterList(clusters=clusters,output_df=ProductMatrix)


# In[ ]:


Number_of_clusters


# In[ ]:


clusters


# In[ ]:


(clusters['labels']==0)


# ### Cluster Analysis and Visualization
# 
# #### Creating Clusters DataFrame:
# - A new DataFrame `Clusters` is created to store cluster-related information.
# - The 'Genes' column in `Clusters` is populated with the gene names extracted from `ProductMatrix`'s index.
# - The 'Cluster' column in `Clusters` is assigned values from the 'labels' obtained from the clustering analysis.
# 
# #### Cluster List Initialization:
# - `Cluster_list` is initialized as an empty dictionary to store cluster information.
# - Unique cluster labels are extracted from the 'labels' array obtained from the clustering result.
# - For each unique label:
#   - A new entry in `Cluster_list` is created with the label as the key.
#   - The corresponding genes that belong to that cluster are extracted and stored as a list.
# 
# #### Converting Cluster List to DataFrame:
# - `Cluster_list` is converted to a DataFrame using `pd.DataFrame.from_dict()`, where each row represents a cluster.
# - The index is reset to create a more structured DataFrame.
# - A new 'labels' column is created to store cluster labels for reference.
# 
# #### Cluster Size Calculation:
# - For each cluster in `Cluster_list`, the size of the cluster (number of genes) is calculated.
# - The cluster sizes are stored in a new DataFrame `ClusterSize` with a column named 'Cluster_size'.
# - `ClusterSize` is sorted in ascending order based on cluster size.
# 
# #### Cluster Size Visualization:
# - A bar plot is created to visualize the distribution of cluster sizes.
# - The x-axis represents the cluster size, and the y-axis represents the count of clusters with that size.
# - The resulting plot provides insights into the distribution of genes across clusters.
# 
# #### Saving Cluster Size Data:
# - The cluster size information is saved as a CSV file named "ClusterSize_6_3.csv" for further analysis and reference.
# 
# ### Objective of This Code Segment:
# - This code segment performs cluster analysis on gene data and provides insights into the distribution and sizes of gene clusters.
# - The resulting cluster size visualization helps in understanding the composition of gene clusters in the dataset.
# 

# In[ ]:


Clusters=pd.DataFrame([])
Clusters['Genes']=ProductMatrix.index
Clusters['Cluster']=clusters['labels']
Cluster_list={}
labels=np.unique(clusters['labels'])
for label in labels:
    print(label)
    Cluster_list[label]=list(Clusters[Clusters.Cluster==label]['Genes']);
Cluster_list=pd.DataFrame.from_dict(Cluster_list,orient='index')
Cluster_list.reset_index(inplace=True)
for i in Cluster_list.index:
    print(i)
    Cluster_list.at[i,'labels']="{}".format(i)
Cluster_list.drop(columns=['index'],inplace=True)
ClusterSize=pd.DataFrame(Cluster_list.notna().sum(axis=1))
ClusterSize.set_axis(['Cluster_size'],inplace=True,axis=1)
ClusterSize.sort_values(by=['Cluster_size'],inplace=True)
ClusterSize['Cluster_size'].value_counts().plot(kind='bar')
# plt.xticks(np.arange(0,65,5))
ClusterSize.to_csv("ClusterSize_6_3.csv")


# ### Visualization of Cluster Sizes
# 
# #### Cluster Size Histogram:
# - A histogram of cluster sizes is created using `ClusterSize.hist()`.
# - The `grid` parameter is set to `False` to remove gridlines.
# - The `bins` parameter is set to `60` to control the number of bins in the histogram.
# 
# #### Customizing Plot Limits and Labels:
# - `plt.ylim(0, 28)` sets the y-axis limits to focus on the relevant frequency range.
# - `plt.xlim(0, 65)` sets the x-axis limits to show cluster sizes up to 65.
# - `plt.xlabel("Size of cluster")` sets the x-axis label to describe the data being plotted.
# - `plt.ylabel("Frequency")` sets the y-axis label to represent the count of clusters.
# - `plt.xticks(np.arange(0, 65, 5))` customizes the x-axis tick marks at intervals of 5.
# 
# #### Saving the Plot:
# - The histogram plot is saved in two formats:
#   - "ClusterSizePlot.svg" as an SVG file with a DPI of 300.
#   - "ClusterSizePlot.png" as a PNG file with a DPI of 300.
# 
# #### Saving Cluster List Data:
# - `Cluster_list` is saved as a GMT (Gene Matrix Transposed) file named "ProductMatrixClusterListFluxRed_6_3_062722.gmt".
# - `Cluster_list` is also saved as a CSV file named "ProductMatrixClusterListFluxRed_6_3_062722.csv" using a tab separator.
# 
# ### Objective of This Code Segment:
# - This code segment visualizes the distribution of cluster sizes using a histogram.
# - Customizations are applied to the plot, and the resulting histogram is saved in both SVG and PNG formats.
# - Additionally, the gene cluster list is saved in GMT and CSV formats for reference and further analysis.
# 

# In[ ]:


ClusterSize.hist(grid=False,bins=60)
plt.ylim(0,28)
plt.xlim(0,65)
plt.xlabel("Size of cluster")
plt.ylabel("Frequency")
plt.xticks(np.arange(0,65,5))
plt.savefig("ClusterSizePlot.svg",dpi=300)
plt.savefig("ClusterSizePlot.png",dpi=300)
Cluster_list.to_csv("ProductMatrixClusterListFluxRed_6_3_062722.gmt",sep='\t')
Cluster_list.to_csv("ProductMatrixClusterListFluxRed_6_3_062722.csv",sep='\t')


# In[ ]:


# ProductMatrix.drop(columns=['Silhoutte_Score', 'Cluster'],inplace=True)


# In[ ]:


# ProductMatrix.drop(columns=['Silhoutte_Score','Cluster'],inplace=True)


# ### Function: `clustersheatmap`
# 
# **Description:** 
# This function generates a heatmap for a specified gene and its associated cluster. It first identifies the cluster where the gene belongs and then creates a heatmap of the genes within that cluster.
# 
# **Parameters:**
# - `genename`: The name of the gene of interest.
# - `deepSplit`: Parameter for hierarchical clustering.
# - `minClusterSize`: Minimum size of clusters to be considered.
# - `matrix`: The gene matrix from which the clusters are derived.
# 
# **Functionality:**
# 1. Obtains a dictionary `genedict` mapping the genes to their respective clusters or subtrees using `display_the_gene_in_respective_cluster_or_subtree` function.
# 2. Extracts a list `genelist` of genes from the dictionary.
# 3. Prints the list of genes in the cluster associated with the specified gene.
# 4. Initializes a figure with specified dimensions.
# 5. Extracts a submatrix `aagr` from the main matrix, containing only the genes in the identified cluster.
# 6. Generates a clustermap (heatmap with hierarchical clustering) using seaborn's `clustermap` function. The heatmap displays the relationships between genes in the cluster with specific configurations (e.g., color map 'YlOrRd', value range from 0 to 0.8).
# 7. Saves the generated heatmap in both SVG and PNG formats with filenames based on the specified `genename`.
# 
# **Example Usage:**
# ```python
# clustersheatmap(genename='example_gene', deepSplit=3, minClusterSize=6, matrix=your_matrix)
# 

# In[ ]:



def clustersheatmap(genename,deepSplit,minClusterSize,matrix):
    genedict=display_the_gene_in_respective_cluster_or_subtree(deepSplit=deepSplit,gene_name=genename,
                                                               matrix=matrix,
                                                               minClusterSize=minClusterSize)
    genelist=list()
    for key in genedict.keys():
#         print("key is {}".format(key))
        gene=genedict[key]
        print("gene is {}".format(gene))
        genelist.append(gene)
    print("gene list is {}".format(genelist))
    fig=plt.figure(figsize=(2,2))
    aagr_cluster=genelist
    aagr=matrix.loc[aagr_cluster][aagr_cluster]
    sns.clustermap(aagr,method='average',cbar_kws={'label':'Coflux * Coexpression'},vmin=0,vmax=0.8,cmap='YlOrRd')
    plt.savefig("{}_cluster.svg".format(genename),dpi=600)
    plt.savefig("{}_cluster.png".format(genename),dpi=600)


# ## Silhouette scores of clusters

# ### Code Snippet Description
# 
# **Purpose:** 
# This snippet demonstrates the calculation and integration of silhouette scores into a given dataframe. Silhouette scores are used to assess the separation distance between the resulting clusters.
# 
# **Process:**
# 1. Calculates silhouette scores for each sample in a dataset using the `SilhoutteSample` function.
# 2. Incorporates these scores into the `ProductMatrix` dataframe using the `Silhoutte_Values` function.
# 
# 
# 

# In[ ]:


silhoutte_score_values=SilhoutteSample(dist=dist,labels=clusters['labels'],metric='precomputed')
Silhoutte_Val=Silhoutte_Values(output_df=ProductMatrix,labels=clusters['labels'],
                               silhouette_score_values=silhoutte_score_values)


# In[ ]:


# Silhoutte_Val=wb_to_gene(Silhoutte_Val)


# ### Find which cluster does propionate shunt fall in to test for positive control

# In[ ]:


Silhoutte_Val.loc['acdh-1']['Cluster']


# ### Code Snippet Description
# 
# **Purpose:** 
# This snippet calculates the average silhouette score for a given set of clusters, providing a measure of how distinct the clusters are within the dataset.
# 
# **Process:**
# 1. Determines the number of unique clusters in the dataset.
# 2. Computes the average silhouette score across all clusters.
# 

# In[ ]:


n_clusters=len(np.unique(clusters['labels']))
silhouette_avg = silhouette_score(dist, labels=clusters['labels'],metric='precomputed')


# 
# **Purpose:** 
# This snippet filters the dataframe `Silhoutte_Val` to include only the rows corresponding to the same cluster as the sample labeled 'acdh-1'. This is useful for analyzing data from a dynamically determined cluster based on a specific sample.
# 
# 

# In[ ]:


y=Silhoutte_Val[Silhoutte_Val.Cluster==(Silhoutte_Val.loc['acdh-1']['Cluster'])]


# In[ ]:


y.index


# In[ ]:


y=y.loc[['acdh-1','ech-6','hach-1', 'hphd-1','bckd-1A','acdh-2' , 'ard-1', 'B0250.5', 'acdh-9', 
        'bckd-1B', 'dbt-1', 'Y43F4A.4' , 'acdh-3']]


# 
# **Purpose:** 
# This snippet calculates the mean of the silhouette scores from a subset of the `Silhoutte_Val` dataframe. This subset (`y`) is defined by a specific cluster label, and the mean is computed for the silhouette scores within this cluster.
# 
# 

# In[ ]:


Mean_threshold=Silhoutte_Val.loc[y.index].Silhoutte_Score.mean()


# In[ ]:


Mean_threshold


# 
# **Purpose:** 
# This snippet creates a dataframe `Clusters` that contains the unique cluster labels from `clusters['labels']` and calculates the mean silhouette values for each cluster. It then computes the median of these mean silhouette values.
# 
# **Process:**
# 1. Initializes an empty dataframe `Clusters`.
# 2. Assigns unique cluster labels to `Clusters['labels']`.
# 3. Calculates the silhouette scores for each sample using `silhouette_samples`.
# 4. Iterates through each unique cluster label to calculate the mean silhouette value for that cluster and adds it to `Clusters`.
# 5. Calculates the median of the mean silhouette values across all clusters.
# 
# 

# In[ ]:


Clusters=pd.DataFrame([]);
Clusters['labels']=np.unique(clusters['labels']);
sample_silhouette_values = silhouette_samples(dist, labels=clusters['labels'], metric='precomputed')
for i in np.unique(clusters['labels']):
    x=np.mean(sample_silhouette_values[clusters['labels'] == i])
    Clusters.at[i,'Mean Silhoutte Values']=x
median=Clusters['Mean Silhoutte Values'].median()   


# In[ ]:


median


# 
# ### Silhouette Score Visualization
# 
# #### Scatter Plot of Silhouette Scores:
# - A scatter plot is generated using `plt.scatter()` to visualize silhouette scores (`y.Silhouette_Score`) for genes in the shunt cluster.
# - The `vmin` and `vmax` parameters are used to set the color scale for the scatter plot.
# - `vmin` is set to `0`, and `vmax` is set to `0.8`.
# - Silhouette scores are plotted on the y-axis, and gene indices (e.g., `y.index`) are plotted on the x-axis.
# 
# #### Horizontal Threshold Lines:
# - Two horizontal threshold lines are added using `plt.hlines()`:
#   - A green line represents the mean silhouette score of the shunt cluster (`Mean_threshold`) with a label.
#   - A red line represents the threshold of the mean silhouette score of selected clusters (`median`) with a label.
# 
# #### Legend and Labels:
# - A legend is displayed in the upper-left corner (`plt.legend(loc='upper left')`) to explain the threshold lines.
# - The x-axis is labeled as "Genes in shunt cluster" with a font size of 15.
# - The y-axis is labeled as "Silhouette Score" with a font size of 15.
# - The y-axis limits are set to range from 0 to 0.6 using `plt.ylim(0, 0.6)`.
# 
# #### Font Customization:
# - The font family is set to Arial (`plt.rcParams["font.family"] = "Arial"`).
# 
# #### Saving the Plot:
# - The scatter plot is saved as an SVG file named "Shunt_Genes_SS_{minClusterSize}_{deepSplit}.svg".
# - The DPI (dots per inch) is set to 600 for high resolution.
# - The plot is saved with a transparent background.
# 
# ### Objective of This Code Segment:
# - This code segment creates a scatter plot to visualize silhouette scores for genes in the shunt cluster.
# - Threshold lines are added to highlight specific silhouette score values.
# - The resulting plot is customized with labels, legends, and font settings.
# - It is saved as an SVG file for further analysis and documentation.
# 

# In[ ]:


plt.figure(figsize=(12,4))
plt.scatter(y.index,y.Silhoutte_Score,vmin=0,vmax=0.8)
plt.hlines(Mean_threshold,xmin=-0.5,xmax=12,color='green',label='Mean of Shunt cluster SS=0.17')
plt.hlines(median,xmin=-0.5,xmax=12,color='red',label='Threshold of mean SS of selected clusters=0.069')
plt.legend(loc='upper left')
plt.xlabel("Genes in shunt cluster",fontsize=15)
plt.ylabel("Silhouette Score",fontsize=15)
plt.ylim(0,0.6)
plt.rcParams["font.family"] = "Arial"
# plt.title("Plot showing Silhouette scores(SS) of cluster \ncontaining shunt genes",fontsize=18)
plt.savefig("Shunt_Genes_SS_{}_{}.svg".format(minClusterSize,deepSplit),dpi=600,transparent=True)


# ### Distribution of Mean Silhouette Scores
# 
# #### Histogram Plot:
# - A histogram plot is generated using `Silhoutte_Val.Silhoutte_Score.hist()` to visualize the distribution of mean silhouette scores.
# - The `bins` parameter is set to `100` for better granularity.
# - Grid lines are turned off using `grid=False`.
# 
# #### Labels and Vertical Threshold Line:
# - The x-axis is labeled as "Mean Silhouette Scores".
# - The y-axis is labeled as "No. of genes".
# - A vertical threshold line is added using `plt.axvline()` with `x=median` to represent the threshold of mean silhouette scores of clusters.
# - The line is colored red and labeled as 'Threshold of MSS of clusters=0.069'.
# 
# #### Font Customization:
# - The font family is set to Arial (`plt.rcParams["font.family"] = "Arial"`).
# 
# #### Title:
# - The plot is given a title: "Distribution of Silhouette scores of genes as per their placement in clusters".
# 
# #### Saving the Plot (Optional):
# - You can choose to save the plot using `plt.savefig()` if desired.
# 
# ### Objective of This Code Segment:
# - This code segment creates a histogram plot to visualize the distribution of mean silhouette scores for genes based on their placement in clusters.
# - It helps in understanding the spread of silhouette scores and identifies the threshold value for cluster quality.
# - The resulting plot is customized with labels, legends, and font settings.
# 

# In[ ]:


plt.figure(figsize=(5,4))
Silhoutte_Val.Silhoutte_Score.hist(bins=100,grid=False)
plt.xlabel("Mean Silhoutte Scores")
plt.ylabel("No. of genes")
plt.axvline(x=median,color='red',label='Threshold of MSS of clusters=0.069')
plt.rcParams["font.family"] = "Arial"
plt.title("Distribution of Silhoutte scores of genes \nas per their placement in clusters")


# ### Distribution of Mean Silhouette Scores for Clusters
# 
# #### Histogram Plot:
# - A histogram plot is generated using `plt.hist()` to visualize the distribution of mean silhouette scores for clusters.
# - The `density=True` parameter is used to normalize the histogram.
# - `bins=40` is set to define the number of bins for the histogram.
# 
# #### Vertical Threshold Line:
# - A vertical threshold line is added using `plt.axvline()` with `x=median` to represent the threshold of mean silhouette scores of selected clusters.
# - The line is colored red and labeled as 'Selected Clusters threshold (Median)=0.069'.
# 
# #### Labels and Font Customization:
# - The x-axis is labeled as "Mean Silhouette Score".
# - The y-axis is labeled as "Number of clusters".
# - The font family is set to Arial (`plt.rcParams["font.family"] = "Arial"`).
# 
# #### Saving the Plot (Optional):
# - You can choose to save the plot using `plt.savefig()` if desired.
# 
# ### Objective of This Code Segment:
# - This code segment creates a histogram plot to visualize the distribution of mean silhouette scores for clusters.
# - It helps in understanding the distribution of cluster quality based on their mean silhouette scores.
# - The resulting plot is customized with labels, legends, and font settings.
# 

# In[ ]:


fig=plt.figure(figsize=(5,4))
# x = np.linspace(-0.4,0.4,50)

# fitted_data = norm.pdf(x, mean, var)
# plt.axvline(mean,label='Mean =0.09986',color='darkgreen')
plt.axvline(median,label='Selected Clusters \nthreshold (Median)=0.069',color='red')
plt.hist(Clusters['Mean Silhoutte Values'], density=True,bins=40)
# plt.legend(loc='upper left')
plt.xlabel("Mean Silhoutte Score",fontsize=15)
plt.ylabel("Number of clusters",fontsize=15)
plt.rcParams["font.family"] = "Arial"
# plt.plot(x,fitted_data,'r-')
plt.savefig("MeanthresholdSS_{}_{}.svg".format(minClusterSize,deepSplit),dpi=300,transparent=True)


# ### Silhouette Plot for Cluster Evaluation
# 
# #### Description:
# - The code segment defines a function `PlotSilhouttePlot` to create a silhouette plot for evaluating clustering results.
# 
# #### Function Parameters:
# - `output_df`: A DataFrame containing the data to be clustered.
# - `labels`: Cluster labels assigned to data points.
# 
# #### Silhouette Plot Details:
# - The function creates a vertical figure with a specified size (`figsize=[10, 100]`).
# - It sets the x-axis limits to [-0.2, 0.4] and y-axis limits to [0, number of data points + (number of clusters + 1) * 10].
# - Silhouette scores for each sample are computed using `silhouette_samples()` and stored in `sample_silhouette_values`.
# - The silhouette plot is generated for each cluster, with clusters represented by filled areas in different colors.
# - Cluster numbers are labeled in the middle of the silhouette plots.
# - The vertical dashed line represents the average silhouette score (`ax.axvline(x=median, color="red", linestyle="--")`).
# - Y-axis labels and ticks are removed, and specific x-axis ticks are set.
# - The resulting silhouette plot is saved as "Silhoutte_plot.svg" and displayed.
# 
# #### Objective of This Code Segment:
# - This code segment creates a silhouette plot to visually assess the quality of clusters.
# - Silhouette plots help determine how well-separated clusters are and identify poorly clustered data points.
# - The average silhouette score threshold is indicated by the red dashed line.
# - It provides insights into the clustering effectiveness.
# 

# In[ ]:


def PlotSilhouttePlot(output_df,labels):    
    fig, ax = plt.subplots(figsize=[10, 100])
    # fig = plt.figure(figsize=(3, 100))
    ax.set_xlim([-0.2, 0.4])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax.set_ylim([0, output_df.shape[0] + (n_clusters + 1) * 10])
    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(dist, labels=labels, metric='precomputed')

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values =             sample_silhouette_values[labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title("The silhouette plot for the various clusters.")
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=median, color="red", linestyle="--")

    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks([-0.2,-0.1, 0, 0.2, 0.4,0.6])

    # 2nd Plot showing the actual clusters formed
    # colors = cm.nipy_spectral(labels.astype(float) / n_clusters)
    # ax2.scatter(dist[:, 0], dist[:, 1], marker='.', s=30, lw=0, alpha=0.7,
    #             c=colors, edgecolor='k')

    # # Labeling the clusters
    # centers = clusterer.cluster_centers_
    # # Draw white circles at cluster centers
    # ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
    #             c="white", alpha=1, s=200, edgecolor='k')

    # for i, c in enumerate(centers):
    #     ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
    #                 s=50, edgecolor='k')

    # ax2.set_title("The visualization of the clustered data.")
    # ax2.set_xlabel("Feature space for the 1st feature")
    # ax2.set_ylabel("Feature space for the 2nd feature")

    # plt.suptitle(("Silhouette analysis for Hierarchical clustering on sample data "
    #               "with n_clusters = %d" % n_clusters),
    #              fontsize=14, fontweight='bold')
    plt.savefig("Silhoutte_plot.svg")
    plt.show()
    
PlotSilhouttePlot(output_df=ProductMatrix,labels=clusters['labels'])


# ### Cluster Labeling and Selection
# 
# #### Description:
# - The code segment assigns formatted label names to clusters and selects clusters based on specific criteria.
# 
# #### Steps in the Code:
# 
# 1. **Assign Label Names to Clusters:**
#    - Each cluster index in the `Clusters` DataFrame is assigned a formatted label name of the form "Label_{}". The labels are stored in a new column called 'LabelName'.
#    - The first cluster (index 0) is labeled as "Unclustered_Cluster" for clarity.
# 
# 2. **Assign Label Names to Cluster List:**
#    - Similarly, each cluster index in the `Cluster_list` DataFrame is assigned a formatted label name.
#    - The first cluster (index 0) is labeled as "Unclustered_Cluster" for consistency.
# 
# 3. **Set Indexes:**
#    - The 'LabelName' column is set as the index for both the `Cluster_list` and `Clusters` DataFrames.
# 
# 4. **Select Clusters:**
#    - The code identifies the intersection of indexes between `Clusters` and `Cluster_list` and stores them as `selected_indices`.
# 
# 5. **Create Selected Clusters DataFrame:**
#    - A new DataFrame called `Selected_Clusters` is created by combining relevant columns from both `Cluster_list` and `Clusters`. This includes the 'Mean Silhoutte Values' and up to 38 genes associated with each cluster.
#    - Columns with integer names are renamed to 'Gene_X', where X represents the gene number.
#    - The index of `Selected_Clusters` is set to 'ClusterName', with labels of the form "Cluster_{}". The original 'LabelName' column is dropped.
# 
# 6. **Save Selected Clusters:**
#    - The `Selected_Clusters` DataFrame is saved as "ClusterList_6_3_021523.csv".
# 
# #### Objective of This Code Segment:
# - This code segment prepares and selects clusters for further analysis based on their silhouette scores and associated genes.
# - It assigns meaningful label names to clusters for easier reference.
# - The selected clusters are saved as a CSV file for future reference and analysis.
# 

# In[ ]:


for index in Clusters.index:
    # Assign a formatted label name to each cluster index
    Clusters.at[index, 'LabelName'] = "Label_{}".format(index)
for index in Clusters.index:
    if index==0:
        Clusters.at[index,'LabelName']="Unclustered_Cluster"
for index in Cluster_list.index:
#     print(index)
    Cluster_list.at[index,'LabelName']="Label_{}".format((index))
for index in Cluster_list.index:
    if index==0:
        Cluster_list.at[index,'LabelName']="Unclustered_Cluster"
Cluster_list.set_index(['LabelName'],inplace=True)
Clusters.set_index(['LabelName'],inplace=True)
selected_indices=list(set(Clusters.index).intersection(set(Cluster_list.index)))
Selected_Clusters=Cluster_list
Selected_Clusters['Mean Silhoutte Values']=Clusters['Mean Silhoutte Values']
Selected_Clusters.sort_values(by=['Mean Silhoutte Values'],inplace=True,ascending=False)
for column in Selected_Clusters.columns:
    print(column)
    if type(column)==int:
        print("int")
        Selected_Clusters.rename(columns={column:'Gene_{}'.format(column+1)},inplace=True)
Selected_Clusters=Selected_Clusters[['Mean Silhoutte Values','Gene_1', 'Gene_2', 'Gene_3', 'Gene_4', 'Gene_5', 'Gene_6', 'Gene_7',
       'Gene_8', 'Gene_9', 'Gene_10', 'Gene_11', 'Gene_12', 'Gene_13',
       'Gene_14', 'Gene_15', 'Gene_16', 'Gene_17', 'Gene_18', 'Gene_19', 'Gene_20', 'Gene_21', 'Gene_22', 'Gene_23',
       'Gene_24', 'Gene_25', 'Gene_26', 'Gene_27', 'Gene_28','Gene_29', 'Gene_30', 'Gene_31', 'Gene_32', 'Gene_33',
       'Gene_34', 'Gene_35', 'Gene_36', 'Gene_37', 'Gene_38']]   
Selected_Clusters.reset_index(inplace=True)
for index in Selected_Clusters.index:
#     print(index)
    Selected_Clusters.at[index,'ClusterName']="Cluster_{}".format((index+1))
Selected_Clusters.set_index(['ClusterName'],inplace=True)
Selected_Clusters.drop(columns=['LabelName'],inplace=True)
Selected_Clusters.to_csv("ClusterList_6_3_021523.csv")


# ### Histogram of Mean Silhouette Scores
# 
# #### Description:
# - This code segment generates a histogram to visualize the distribution of mean silhouette scores for clusters that are above a specified threshold.
# 
# #### Steps in the Code:
# 
# 1. **Create Histogram:**
#    - The code segment generates a histogram of the 'Mean Silhoutte Values' from the `Selected_Clusters` DataFrame.
#    - The `grid` parameter is set to `False` to remove grid lines from the plot.
# 
# 2. **Add Title and Labels:**
#    - The title of the histogram is set to 'Mean Silhoutte Scores of clusters that are above threshold'.
#    - A vertical line is added to indicate the threshold value, labeled as 'Threshold', and colored green.
#    - The x-axis label is set as "Mean Silhoutte Score", and the y-axis label as "No. of clusters".
# 
# 3. **Add Legend:**
#    - A legend is added to the plot, positioned at the best location, to label the threshold line.
# 
# #### Objective of This Code Segment:
# - This code segment provides a visualization of the distribution of mean silhouette scores for clusters that meet a specified threshold.
# - It helps users understand the distribution of cluster quality based on silhouette scores.
# 
# 

# In[ ]:


Selected_Clusters['Mean Silhoutte Values'].hist(grid=False)
plt.title('Mean Silhoutte Scores of clusters that are above threshold')
plt.axvline(x=median,label='Threshold',color='green')
plt.xlabel("Mean Silhoutte Score")
plt.legend(loc='best')
plt.ylabel("No. of clusters")


# ### Exporting Cluster and Selected Cluster Data
# 
# #### Description:
# - This code segment exports cluster-related data to CSV files.
# 
# #### Code Explanation:
# 
# 1. **Export Cluster List:**
#    - The code exports the `Cluster_list` DataFrame to a GMT format file named "ClusterSetsFluxRed_062722_6_3.gmt" using tab (`\t`) as the separator.
# 
# 2. **Export Selected Clusters:**
#    - The code exports the `Selected_Clusters` DataFrame to two different files:
#      - One in GMT format named "ProductMatrix_AllClusterSets_063022_6_3.gmt" using tab (`\t`) as the separator.
#      - Another in CSV format named "ProductMatrix_AllClusterSets_063022_6_3.csv".
# 
# #### Objective of This Code Segment:
# - This code segment allows you to export cluster-related data to files in the specified formats for further analysis or sharing with others.
# 
# You can copy and paste this markdown into your project's `README.md` file to document the export of cluster data.
# 

# In[ ]:


Cluster_list.to_csv("ClusterSetsFluxRed_062722_6_3.gmt",sep='\t')
Selected_Clusters.to_csv("ProductMatrix_AllClusterSets_063022_6_3.gmt",sep='\t')
Selected_Clusters.to_csv("ProductMatrix_AllClusterSets_063022_6_3.csv")


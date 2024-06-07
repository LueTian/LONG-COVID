# packages
import pandas as pd
import numpy as np
import os
import random
import scipy.stats as stats
import time
import itertools
import multiprocessing
from joblib import Parallel,delayed
import networkx as nx
random.seed(0)

#### Function 
# Compute Distance of the permutation sample
class PPI_Class:
    '''
    该类用于计算正负样本的 PPI distance 
    包含下列步骤：
    1.构建PPI LCC 网络
    2.构建 node-node-distance dict
    3.计算 group-group-distance 
    '''
    def __init__(self,net_df):
        self.net_df = net_df
        self.PPI_dict = {}

    
    def Net_Construct(self):
        G = nx.Graph()
        G.add_edges_from(np.array(self.net_df[['source','target']]))
        c = max(nx.connected_components(G), key=len)
        G = G.subgraph(c)
        self.net = G 
        # following codes for sampling
        node_degree_df = pd.DataFrame([list(i) for i in G.degree()],columns = ['node','degree'])
        node_degree_df = node_degree_df.sort_values('degree').reset_index(drop = True)     
        self.node_degree_df = node_degree_df

    def node_node_distance(self,node1,node2):
        '''
        :param G:
        :param node1: 
        :param node2:  
        :return: the length of the shortest path of node1 and node2
        '''
        G = self.net
        outcome = nx.shortest_path_length(G,source = node1,target=node2)
        return outcome


    def group_group_distance(self,input): 
        '''
        G_S: source condition proteins 
        G_T: target condition proteins

        outcome:
        update self.PPI_dict 
        Group Group Distance [three types] Dict Frame
        '''
        reference_dict = self.PPI_dict
        G_S = input[0]
        G_T = input[1]
        GG_dis = []
        # method one  the mean of min dis(G_S,each of G_T)
        for node1 in set(G_T):
            NN_dis = [] # the distance {each of G_S, A of G_T}
            for node2 in set(G_S):
                node_pair = sorted([node1,node2])
                node_pair_str = str(node_pair)
                if node_pair_str in reference_dict.keys():
                    distance = reference_dict[node_pair_str]
                else:
                    distance = self.node_node_distance(node1,node2)
                    reference_dict[node_pair_str] = distance
                NN_dis.append(distance)
            NG_value = np.min(NN_dis) # the min distance between the target node and the source disease
            GG_dis.append(NG_value)   # the distance between the target disease and the source disease 
        outcome_1 = np.mean(GG_dis)   # the mean of sum of min distance(each of G_T,G_S)

        GG_dis_2 = []
        for node1 in set(G_S):
            NN_dis = []
            for node2 in set(G_T):
                node_pair = sorted([node1,node2])
                node_pair_str = str(node_pair)
                if node_pair_str in reference_dict.keys():
                    distance = reference_dict[node_pair_str]
                else:
                    distance = self.node_node_distance(node1,node2)
                    reference_dict[node_pair_str] = distance
                NN_dis.append(distance)            
            NG_value = np.min(NN_dis)
            GG_dis.append(NG_value)
            GG_dis_2.append(NG_value)
        # method two the mean of sum of min dis(G_T,each of G_S)
        outcome_2 = np.mean(GG_dis_2)
        # method three the mean of sum of (min dis(G_T,each of G_S) + min dis(G_S,each of G_T))
        outcome_3 = np.mean(GG_dis)
        self.PPI_dict = reference_dict
        return [outcome_1,outcome_2,outcome_3]
    



#############################################################################
#############################################################################
#############################################################################
#### Data Load 
COVID_LCC = pd.read_csv(r"Data\Code_02_Outcome\COVID_LCC.csv")
ICD_LCC   = pd.read_csv(r"Data\Code_02_Outcome\ICD_LCC.csv")

ICD_LCC.shape

PPI_LCC   = pd.read_csv(r"Data\Code_02_Outcome\PPI_LCC.csv")
ICD_LCC_dict = {i:eval(j) for i,j in zip(ICD_LCC['ICDID'],ICD_LCC['geneId_LCC'])}

PPI_LCC   = PPI_LCC[['proteinA_entrezid','proteinB_entrezid']]
PPI_LCC.columns = ['source','target']

ppi_net = PPI_Class(PPI_LCC)
ppi_net.Net_Construct()
ppi_net.PPI_dict = ppi_net.PPI_dict = {i[0]:int(i[1]) for i in np.load('Data\Code_03_Outcome\PPI_dict.npy')}

file_path = "Data\\Code_04_Outcome"
file_path_save = "Data\\Code_05_Outcome"
file_ls = os.listdir(file_path)
for i in file_ls:
    df = pd.read_csv(os.path.join(file_path,i))
    df['source_protein'] = [eval(i) for i in df['source_protein']]
    df['target_protein'] = [eval(i) for i in df['target_protein']]
    df['source_protein_add'] = [sorted(set(i+COVID_LCC['EntrezID'].tolist())) for i in df['source_protein']]
    df[['outcome1','outcome2','outcome3']] = df[['source_protein','target_protein']].apply(ppi_net.group_group_distance,axis = 1,result_type='expand')
    df[['outcome1_add','outcome2_add','outcome3_add']] = df[['source_protein_add','target_protein']].apply(ppi_net.group_group_distance,axis = 1,result_type='expand')
    df.to_csv(os.path.join(file_path_save,i),index=False)

array = np.array(list(ppi_net.PPI_dict.items()))
np.save("Data\Code_03_Outcome\PPI_dict_update.npy",array)


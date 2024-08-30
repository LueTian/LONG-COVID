'''
compute ppi distance for positive sample and negative sample
'''
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
# compute ppi distance
class PPI_Class:
    '''
    the class used for compute positive and negative sample PPI distance
    1. construct PPI LCC network 
    2. construct node-node-distance dict 
    3. compute group-group distance
    '''
    
    def __init__(self,net_df):
        '''
        initial dataset
        '''
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

#### Data Load 
COVID_LCC = pd.read_csv(r"Data\Code_02_Outcome\COVID_LCC.csv")
ICD_LCC   = pd.read_csv(r"Data\Code_02_Outcome\ICD_LCC.csv")
PPI_LCC   = pd.read_csv(r"Data\Code_02_Outcome\PPI_LCC.csv")
PPI_LCC   = PPI_LCC[['proteinA_entrezid','proteinB_entrezid']]
PPI_LCC.columns = ['source','target']
PPI_LCC.shape 
len(set(PPI_LCC['source'].tolist()+PPI_LCC['target'].tolist()))

ICD_LCC_dict = {i:eval(j) for i,j in zip(ICD_LCC['ICDID'],ICD_LCC['geneId_LCC'])}

ICD_pair_T = pd.read_csv(r"Data\Code_02_Outcome\ICD_pair_T(2).csv")
ICD_pair_F = pd.read_csv(r"Data\Code_02_Outcome\ICD_pair_F(2).csv")

ICD_pair_F.shape

ICD_pair_T['source_portein'] = [ICD_LCC_dict[i] for i in ICD_pair_T['source']]
ICD_pair_T['target_portein'] = [ICD_LCC_dict[i] for i in ICD_pair_T['target']]
ICD_pair_F['source_portein'] = [ICD_LCC_dict[i] for i in ICD_pair_F['source']]
ICD_pair_F['target_portein'] = [ICD_LCC_dict[i] for i in ICD_pair_F['target']]

ppi_net = PPI_Class(PPI_LCC)
ppi_net.Net_Construct()
ppi_net.PPI_dict = np.load('Data\Code_04_Outcome\PPI_dict.npy', allow_pickle = True).item()

###【consider refer to Distance Computation Result】
ICD_pair_T_ref = pd.read_csv(r"Data\Code_03_Outcome\ICD_pair_T_add.csv")
ICD_pair_F_ref = pd.read_csv(r"Data\Code_03_Outcome\ICD_pair_F_add.csv")
def Ref_Dict(DF):
    dict = {}
    for idx in DF.index:
        edge = DF.loc[idx, 'edge']
        outcome1 = DF.loc[idx, 'outcome1']
        outcome2 = DF.loc[idx, 'outcome2']
        outcome3 = DF.loc[idx, 'outcome3']
        outcome1_add = DF.loc[idx, 'outcome1_add']
        outcome2_add = DF.loc[idx, 'outcome2_add']
        outcome3_add = DF.loc[idx, 'outcome3_add']
        dict[edge] = {'outcome1':outcome1, 'outcome2':outcome2, 'outcome3':outcome3, 'outcome1_add':outcome1_add, 'outcome2_add':outcome2_add, 'outcome3_add':outcome3_add}
    return dict 

Ref_dict_T = Ref_Dict(ICD_pair_T_ref)
Ref_dict_F = Ref_Dict(ICD_pair_F_ref)


### PPI distance before COVID infection
ICD_pair_T_Y = ICD_pair_T[ICD_pair_T['edge'].isin(Ref_dict_T.keys())]
ICD_pair_T_N = ICD_pair_T.drop(index = ICD_pair_T_Y.index)

ICD_pair_F_Y = ICD_pair_F[ICD_pair_F['edge'].isin(Ref_dict_F.keys())]
ICD_pair_F_N = ICD_pair_F.drop(index = ICD_pair_F_Y.index)


ICD_pair_T_N[['outcome1','outcome2','outcome3']] = ICD_pair_T_N[['source_portein','target_portein']].apply(ppi_net.group_group_distance,axis = 1,result_type='expand')
ICD_pair_F_N[['outcome1','outcome2','outcome3']] = ICD_pair_F_N[['source_portein','target_portein']].apply(ppi_net.group_group_distance,axis = 1,result_type='expand')

ICD_pair_T_Y['outcome1'] = ICD_pair_T_Y['edge'].apply(lambda x: Ref_dict_T[x]['outcome1'])
ICD_pair_F_Y['outcome1'] = ICD_pair_F_Y['edge'].apply(lambda x: Ref_dict_F[x]['outcome1'])
ICD_pair_T_Y['outcome2'] = ICD_pair_T_Y['edge'].apply(lambda x: Ref_dict_T[x]['outcome2'])
ICD_pair_F_Y['outcome2'] = ICD_pair_F_Y['edge'].apply(lambda x: Ref_dict_F[x]['outcome2'])
ICD_pair_T_Y['outcome3'] = ICD_pair_T_Y['edge'].apply(lambda x: Ref_dict_T[x]['outcome3'])
ICD_pair_F_Y['outcome3'] = ICD_pair_F_Y['edge'].apply(lambda x: Ref_dict_F[x]['outcome3'])

ICD_pair_T_Y['outcome1_add'] = ICD_pair_T_Y['edge'].apply(lambda x: Ref_dict_T[x]['outcome1_add'])
ICD_pair_F_Y['outcome1_add'] = ICD_pair_F_Y['edge'].apply(lambda x: Ref_dict_F[x]['outcome1_add'])
ICD_pair_T_Y['outcome2_add'] = ICD_pair_T_Y['edge'].apply(lambda x: Ref_dict_T[x]['outcome2_add'])
ICD_pair_F_Y['outcome2_add'] = ICD_pair_F_Y['edge'].apply(lambda x: Ref_dict_F[x]['outcome2_add'])
ICD_pair_T_Y['outcome3_add'] = ICD_pair_T_Y['edge'].apply(lambda x: Ref_dict_T[x]['outcome3_add'])
ICD_pair_F_Y['outcome3_add'] = ICD_pair_F_Y['edge'].apply(lambda x: Ref_dict_F[x]['outcome3_add'])

# ICD_pair_T.to_csv(r"Data\Code_03_Outcome\ICD_pair_T.csv",index = False)
# ICD_pair_F.to_csv(r"Data\Code_03_Outcome\ICD_pair_F.csv",index = False)

### PPI distance after COVID infection
ICD_pair_T_N['source_portein_add'] = [sorted(set(i+COVID_LCC['EntrezID'].tolist())) for i in ICD_pair_T_N['source_portein']]
ICD_pair_F_N['source_portein_add'] = [sorted(set(i+COVID_LCC['EntrezID'].tolist())) for i in ICD_pair_F_N['source_portein']]
ICD_pair_T_Y['source_portein_add'] = [sorted(set(i+COVID_LCC['EntrezID'].tolist())) for i in ICD_pair_T_Y['source_portein']]
ICD_pair_F_Y['source_portein_add'] = [sorted(set(i+COVID_LCC['EntrezID'].tolist())) for i in ICD_pair_F_Y['source_portein']]

ICD_pair_T_N[['outcome1_add','outcome2_add','outcome3_add']] = ICD_pair_T_N[['source_portein_add','target_portein']].apply(ppi_net.group_group_distance,axis = 1,result_type='expand')
ICD_pair_F_N[['outcome1_add','outcome2_add','outcome3_add']] = ICD_pair_F_N[['source_portein_add','target_portein']].apply(ppi_net.group_group_distance,axis = 1,result_type='expand')

ICD_pair_T = pd.concat([ICD_pair_T_Y,ICD_pair_T_N],axis = 0)
ICD_pair_F = pd.concat([ICD_pair_F_Y,ICD_pair_F_N],axis = 0)

ICD_pair_T.to_csv(r"Data\Code_03_Outcome\ICD_pair_T_add.csv",index = False)
ICD_pair_F.to_csv(r"Data\Code_03_Outcome\ICD_pair_F_add.csv",index = False)


### PPI distance before COVID infection (simulation)
permutation_sample = np.load("Data\\Code_04_Outcome\\Disease_Sample.npy",allow_pickle=True).item()
repeat_num = 1000
for idx in ICD_pair_T.index:
    permutation_df = []
    edge = ICD_pair_T.loc[idx, 'edge']
    source = ICD_pair_T.loc[idx, 'source']
    target = ICD_pair_T.loc[idx, 'target']
    permutation_df = [[edge,source,target,permutation_sample[source][i],permutation_sample[target][i]] for i in range(repeat_num)]
    permutation_df = pd.DataFrame(permutation_df, columns = ['edge', 'source', 'target', 'source_portein', 'target_portein'])
    permutation_df[['outcome1','outcome2','outcome3']] = permutation_df[['source_portein','target_portein']].apply(ppi_net.group_group_distance,axis = 1,result_type='expand')

    permutation_df['source_portein_add'] = [sorted(set(i+COVID_LCC['EntrezID'].tolist())) for i in permutation_df['source_portein']]
    permutation_df[['outcome1_add','outcome2_add','outcome3_add']] = permutation_df[['source_portein_add','target_portein']].apply(ppi_net.group_group_distance,axis = 1,result_type='expand')
    permutation_df.to_csv(f'Data/Code_04_Outcome/{source}_{target}_permutation.csv', index = False)
    print(idx)

len(permutation_sample[target])
len(ppi_net.PPI_dict.keys()) # 15065692
np.save(r"Data\Code_04_Outcome\PPI_dict.npy",ppi_net.PPI_dict)
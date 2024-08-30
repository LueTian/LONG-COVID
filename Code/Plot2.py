'''
生成 comorbidity network 所需要的数据
'''
import pandas as pd 
import numpy as np 
from sklearn.preprocessing import StandardScaler
import os 
net = pd.read_csv(r"Data\Code_01_Outcome\network_pairwise(2).csv")
net['source'] = net['source'].str.split(' ').str.get(1)
net['target'] = net['target'].str.split(' ').str.get(1)

len(set(net['source'].tolist() + net['target'].tolist()))
net['Pair'] = net['source'].str.cat(net['target'],sep = '_')
net_node = pd.read_csv(r"Data\Code_01_Outcome\node_detail(2).csv")
color_ls = ['#FFB09B','#CBCD84','#FFBF74','#FFADC9','#FFBAF6','#A6D8FF','#EAD55E','#55F5FF','#59F3D1','#EDBBCE','#85E49F','#64CACC','#B2E779','#A4D7C2']

nodes_ls = list(sorted(set(net_node['Detail_extend'])))
color_dict = {i:j for i,j in zip(nodes_ls,color_ls)}
net_node['color'] = [color_dict[i] for i in net_node['Detail_extend']]

net_node.columns = ['ID', 'Detail', 'Detail_extend', 'color']
net_node['ID'] = net_node['ID'].str.split(' ').str.get(1)
net_node.to_csv(r"Data\Network Data\comorbidity_net_node.csv",index=False)
net = net[['source', 'target', 'Diff_RR', 'Diff_Corr', 'Pair']]
net.to_csv(r"Data\Network Data\comorbidity_net_edge.csv",index=False)
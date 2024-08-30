import pandas as pd 
import numpy as np 
import networkx as nx 
import itertools
import os 
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
np.random.seed(0)

#################
### Data Load 
PPI = pd.read_csv(r'Data\Code_Original\Protein_Protein_EntrezID.csv')
COVID = pd.read_csv(r'Data\Code_Original\COVID_Protein_EntrezID.csv')
ICD = pd.read_csv(r'Data\Code_00_Outcome\ICD_gene_map_df.csv')

#################
### PPI network
Network = PPI[['proteinA_entrezid','proteinB_entrezid']]
G = nx.Graph()
G.add_edges_from(np.array(Network[['proteinA_entrezid','proteinB_entrezid']]))
c = max(nx.connected_components(G), key=len)
G = G.subgraph(c)
G_node = list(G.nodes)

###########
### PPI LCC COVID Target Protein
COVID_LCC = COVID[COVID['EntrezID'].isin(G_node)].reset_index(drop = True)
###########
### PPI LCC Interaction   327924 ---> 327868 
PPI_LCC = PPI[(PPI['proteinA_entrezid'].isin(G_node)) & (PPI['proteinB_entrezid'].isin(G_node))]
###########
### PPI LCC ICD Target Protein
def Protein_Overlap(input,ref):
    '''
    keep the overlapping part between input and ref
    '''
    input = eval(input)
    outcome = sorted(set(input).intersection(set(ref)))
    return str(outcome)
ICD['geneId_LCC'] = ICD['geneID'].apply(Protein_Overlap,ref = G_node)
ICD['Num_geneId_LCC'] = [len(eval(i)) for i in ICD['geneId_LCC']]
ICD_LCC = ICD[ICD['Num_geneId_LCC']!= 0 ].reset_index(drop = True)

def Node_Fliter(input):
    '''
    drop all code including 'V' & 'E'
    '''
    Keep_list =['E' + str(i + 950) for i in range(10)]
    index_list = []
    for i in input.index:
        source = input.loc[i,'ICDID']
        if  ('V' in source)|(('E' in source)&(source not in Keep_list)):
            index_list.append(i)
    input = input.drop(index = index_list).reset_index(drop = True)
    return input
ICD_LCC = Node_Fliter(ICD_LCC) # 537 

COVID_LCC.to_csv(r'Data\Code_02_Outcome\COVID_LCC.csv',index = False)
PPI_LCC.to_csv(r'Data\Code_02_Outcome\PPI_LCC.csv',index = False)
ICD_LCC.to_csv(r'Data\Code_02_Outcome\ICD_LCC.csv',index = False)

###########################################################################
# Negative Sample
ICD_LCC = pd.read_csv(r'Data\Code_02_Outcome\ICD_LCC.csv')
def Int2Str(input_str1):
    input_str1 = str(input_str1)
    if len(input_str1) < 3:
        input_str1 = '0'*(3-len(input_str1)) + input_str1
    return input_str1
ICD_LCC['ICDID'] = ICD_LCC['ICDID'].apply(Int2Str)
ICD_ls = ICD_LCC['ICDID'].tolist()
len(ICD_ls) # 537 

# 145
comornet = pd.read_csv(r'Data\Code_01_Outcome\network_pairwise(2).csv')
comornet = comornet[(comornet['source'].isin(ICD_ls))&(comornet['target'].isin(ICD_ls))].reset_index(drop = True) # 145
comornet.value_counts('target').shape

len(set(comornet['source'].tolist() + comornet['target'].tolist()))
comornet.head()
comornet.columns

ICD_pair = [[i,j] for i in ICD_ls for j in ICD_ls]
ICD_pair = pd.DataFrame(ICD_pair,columns=['source','target'])
ICD_pair = ICD_pair[ICD_pair['source']!=ICD_pair['target']].reset_index(drop = True) # 287832
ICD_pair = ICD_pair[ICD_pair['source'].isin(comornet['source'])].reset_index(drop = True) # 39664
ICD_pair['edge'] = ICD_pair[['source', 'target']].apply(lambda x: str([x[0],x[1]]), axis = 1)

ICD_pair_T = ICD_pair[ICD_pair['edge'].isin(comornet['edge'])]# 145
ICD_pair_F = ICD_pair.drop(index = ICD_pair_T.index) # 39519

def ICD_str(input):
    outcome = 'ICD9 ' + input
    return outcome

ICD_LCC = pd.read_csv(r'Data\Code_02_Outcome\ICD_LCC.csv')
# ICD_LCC['ICDID'] = ICD_LCC['ICDID'].apply(Int2Str)
# ICD_LCC['ICDID'] = ICD_LCC['ICDID'].apply(ICD_str)
# ICD_LCC.to_csv(r'Data\Code_02_Outcome\ICD_LCC.csv',index = False)
ICD_pair_T.to_csv(r'Data\Code_02_Outcome\ICD_pair_T(2).csv',index = False)
ICD_pair_F.to_csv(r'Data\Code_02_Outcome\ICD_pair_F(2).csv',index = False)






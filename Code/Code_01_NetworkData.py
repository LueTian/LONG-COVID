### In this code, we need to filter the data 
import pandas as pd 
import numpy  as np 
from scipy import stats

Net1 = pd.read_csv(r"F:\7_3_Code\network_outcome\network_edges_ed1(2).csv")
Net1 = Net1[(Net1['RR (treat)']>Net1['RR (control)'])&(Net1['Correlation (treat)']>Net1['Correlation (control)']) & (Net1['edge num (treat)'] > 10) & (Net1['p_value'] < 0.05)]
Net1 = Net1.reset_index(drop = True)
Net1.shape #(209, 16)
Net1['source (int)'] = Net1['source'].str.split(' ').str.get(1)
Net1['target (int)'] = Net1['target'].str.split(' ').str.get(1) 

Net1[Net1['source (int)'] == '519']['node num (source) (control)']

def Node_Fliter(input):
    Keep_list =['E' + str(i + 950) for i in range(10)]
    index_list = []
    for i in input.index:
        source = input.loc[i,'source (int)']
        target = input.loc[i,'target (int)']
        if ('V' in target) | ('V' in source) | (('E' in target)&(target not in Keep_list)) | (('E' in source)&(source not in Keep_list)):
            index_list.append(i)
    input = input.drop(index = index_list).reset_index(drop = True)
    return input

Net1 = Node_Fliter(Net1) # 161
Net1['Diff_RR'] = Net1['RR (treat)'] / Net1['RR (control)']
Net1['Diff_Corr'] = Net1['Correlation (treat)'] - Net1['Correlation (control)']

Net1.to_csv(r'Data\Code_01_Outcome\network_pairwise(2).csv',index=False)


Net1 = pd.read_csv(r'Data\Code_01_Outcome\network_pairwise(2).csv')

from icd9cms.icd9 import search
Node_list = list(set(Net1['source (int)'].tolist()+ Net1['target (int)'].tolist())) # 96
len(Node_list) # 96
Node_list = ['0'*(3-len(str(i))) + str(i) if len(str(i)) < 3  else str(i) for i in Node_list]
Node_DF = pd.DataFrame(Node_list,columns=['Conditions'])
Node_DF['Detail'] = [search(str(i)) for i in Node_list]
Node_DF['Detail'] = Node_DF['Detail'].astype('str')
Node_DF['Detail_0'] = Node_DF['Detail'].str.split(':').str.get(0)
Node_DF['Detail_1'] = Node_DF['Detail'].str.split(':').str.get(1)
Node_DF['Detail_2'] = Node_DF['Detail'].str.split(':').str.get(2)
Node_DF = Node_DF[Node_DF['Detail_0'] != 'None'].reset_index(drop = True)
Node_DF = Node_DF[['Conditions','Detail_1']]
Node_DF.columns = ['Conditions','Detail']
Node_DF['Detail_extend'] = [str(search(str(i)).parent).rsplit(':',1)[0] for i in Node_DF['Conditions']]
Node_DF =  Node_DF.sort_values('Conditions',ascending=True)

def find_top_parent(node):
    node = search(node)
    while node and node.parent:
        node = node.parent 
    node = str(node).rsplit(':',1)[0]
    return node

Node_DF['Detail_extend'] = Node_DF['Conditions'].astype(str).apply(find_top_parent)
Node_DF.value_counts('Detail_extend').shape

Node_DF['Conditions'] = ['ICD9 ' + i for i in Node_DF['Conditions'].astype(str)]
# Node_DF.to_csv(r'Data\Code_01_Outcome\node_detail(2).csv',index=False)


### compute the frequency of each diseases
### count num 【frequency】
DF_count_S = Net1.value_counts('source').reset_index(drop = False)
DF_count_S.columns = ['Conditions','Count Num']
DF_count_T = Net1.value_counts('target').reset_index(drop = False)
DF_count_T.columns = ['Conditions','Count Num']
### count num【according to correlation】
DF_count_S_ = Net1.groupby('source')['Diff_Corr'].sum().reset_index(drop = False)
DF_count_S_.columns = ['Conditions','Count Num(Weighted Corr)']
DF_count_T_ = Net1.groupby('target')['Diff_Corr'].sum().reset_index(drop = False)
DF_count_T_.columns = ['Conditions','Count Num(Weighted Corr)']

DF_count_S = pd.merge(DF_count_S,DF_count_S_,on = 'Conditions',how = 'left')
DF_count_T = pd.merge(DF_count_T,DF_count_T_,on = 'Conditions',how = 'left')
### count num【according to Relative Risk】 
DF_count_S_ = Net1.groupby('source')['Diff_RR'].sum().reset_index(drop = False)
DF_count_S_.columns = ['Conditions','Count Num(Weighted RR)']
DF_count_T_ = Net1.groupby('target')['Diff_RR'].sum().reset_index(drop = False)
DF_count_T_.columns = ['Conditions','Count Num(Weighted RR)']

DF_count_S = pd.merge(DF_count_S,DF_count_S_,on = 'Conditions',how = 'left')
DF_count_T = pd.merge(DF_count_T,DF_count_T_,on = 'Conditions',how = 'left')

DF_count_T = DF_count_T.sort_values('Conditions',ascending=True)
DF_count_S = DF_count_S.sort_values('Conditions',ascending=True)

# DF_count_T.to_csv(r'Data\Code_01_Outcome\node_count_T(2).csv',index=False)
# DF_count_S.to_csv(r'Data\Code_01_Outcome\node_count_S(2).csv',index=False)

DF_count_T.columns = ['Conditions', 'Count Num Target', 'Count Num(Weighted Corr) Target', 'Count Num(Weighted RR) Target']
DF_count_S.columns = ['Conditions', 'Count Num Source', 'Count Num(Weighted Corr) Source', 'Count Num(Weighted RR) Source']

Condition_ls = DF_count_T['Conditions'].tolist() + DF_count_S['Conditions'].tolist()
Condition_ls = pd.DataFrame(Condition_ls,columns= ['Conditions'])
Condition_ls = pd.merge(Condition_ls,DF_count_S,on ='Conditions',how ='left')
Condition_ls = pd.merge(Condition_ls,DF_count_T,on ='Conditions',how ='left')
Condition_ls = Condition_ls.fillna(0)
Condition_ls = Condition_ls.drop_duplicates().reset_index(drop = True)
Condition_ls = pd.merge(Condition_ls,Node_DF, on = 'Conditions',how = 'left')
# Condition_ls.to_csv(r'Data\Code_01_Outcome\node_frequency(2).csv',index=False)

### correlation test 【Net1】 0.551598413042954
correlation,p_value = stats.pearsonr(Net1['Diff_Corr'],Net1['Diff_RR'])
# (0.36704038444975384, 1.6748341398024643e-06)
correlation,p_value = stats.spearmanr(Net1['Diff_Corr'],Net1['Diff_RR'])
# (0.5310674028065333, 4.2571352930472584e-13)
correlation,p_value = stats.kendalltau(Net1['Diff_Corr'],Net1['Diff_RR'])
# (0.3666149068322981, 5.104369805535729e-12)



### construct Group - Group Disease 
Node_DF_Dict = {i:j for i,j in zip(Node_DF['Conditions'],Node_DF['Detail_extend'])}
Node_DF_Dict.keys()
Net1['source G'] = [Node_DF_Dict[i] for i in Net1['source']]
Net1['target G'] = [Node_DF_Dict[i] for i in Net1['target']]
Net1['G_G_Type'] = Net1['source G'].str.cat(Net1['target G'],sep = '_')
DF_1_count1 = Net1.value_counts('G_G_Type').reset_index(drop = False)
DF_1_count2 = Net1.groupby('G_G_Type')['Diff_Corr'].sum().reset_index(drop = False)
DF_1_count3 = Net1.groupby('G_G_Type')['Diff_RR'].sum().reset_index(drop = False)
DF_1_count1 = pd.merge(DF_1_count1,DF_1_count2,on = 'G_G_Type')
DF_1_count1 = pd.merge(DF_1_count1,DF_1_count3,on = 'G_G_Type')
DF_1_count1.columns = ['G_G_Type', 'Frequency', 'Diff_Corr', 'Diff_RR']
DF_1_count1['source'] = DF_1_count1['G_G_Type'].str.split('_').str.get(0)
DF_1_count1['target'] = DF_1_count1['G_G_Type'].str.split('_').str.get(1)
DF_1_count1 = DF_1_count1.sort_values(['source','target'])
# DF_1_count1.to_csv(r'Data\Code_01_Outcome\group_count(2).csv',index=False)

### Group Correlation  (0.8707024318964619, 2.1391962095335627e-21)
correlation,p_value = stats.pearsonr(DF_1_count1['Diff_Corr'],DF_1_count1['Diff_RR'])
# (0.9019309604757212, 1.7069115789679767e-26)
correlation,p_value = stats.spearmanr(DF_1_count1['Diff_Corr'],DF_1_count1['Diff_RR'])
# (0.8280815326743067, 9.376293103852769e-19)
correlation,p_value = stats.kendalltau(DF_1_count1['Diff_Corr'],DF_1_count1['Diff_RR'])
# (0.6530020703933748, 1.2973056033218199e-15)
DF_1_count1.value_counts('source').shape
DF_1_count1.value_counts('target').shape
##################################################################################

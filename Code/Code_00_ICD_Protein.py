'''
Disease-Protein Mapping 
'''
import pandas as pd
import numpy as np 
import pdfplumber
import re

#######################################################
# Function
#######################################################
def Item_list(Group,str): 
    '''
    output list
    '''
    return sorted(set(Group[str].tolist()))
def DF_Dict(DF):
    DF.columns = ['key','value']
    outcome = {i:j for i,j in zip(DF['key'],DF['value'])}
    return outcome 
def disease_ICD_Dict(input):
    '''将范围拓展成单个 ICD code'''
    v1 = input[0]
    v2 = input[1]
    value = input[2]
    outcome = []
    for i in range(int(v1),int(v2)+1):
        outcome.append([i,value])
    return str(outcome)
def Dict_Merge(Dict1,Dict2):
    result = {}
    for key_ in Dict1.keys():
        value_ = Dict1[key_]
        outcome = []
        for key__ in value_:
            if key__ in Dict2.keys():
                outcome += Dict2[key__]
        result[key_] = sorted(set(outcome))
    return result
def Int2Str(input_str1):
    input_str1 = str(input_str1)
    if len(input_str1) < 3:
        input_str1 = '0'*(3-len(input_str1)) + input_str1
    return input_str1
def Phenotype_Split(input):
    phenotype = input[0]
    phenotypeMimNumber = input[1]
# Long phenotype
    matcher = re.match(r'^(.*),\s(\d{6})\s\((\d)\)(|, (.*))$', phenotype)
    if matcher:

        # Get the fields
        phenotype = matcher.group(1)
        phenotypeMimNumber = matcher.group(2)
        phenotypeMappingKey = matcher.group(3)
        inheritances = matcher.group(5)

        # Get the inheritances, may or may not be there
        if inheritances:
            for inheritance in inheritances.split(','):
                inheritance = inheritance.strip()

    # Short phenotype
    else:
        # Short phenotype
        matcher = re.match(r'^(.*)\((\d)\)(|, (.*))$', phenotype)
        if matcher:
            # Get the fields
            phenotype = matcher.group(1)
            phenotypeMappingKey = matcher.group(2)
            inheritances = matcher.group(3)

            # Get the inheritances, may or may not be there
            if inheritances:
                for inheritance in inheritances.split(','):
                    inheritance = inheritance.strip()
    return phenotypeMimNumber

#######################################################
# Function
#######################################################
### 使用 DisgeNet 数据库
### geneId --- diseaseId map
gene_disease_map = pd.read_csv(r"Data\Code_Original\curated_gene_disease_associations.tsv",sep = '\t')
gene_disease_map = gene_disease_map[['geneId', 'geneSymbol', 'diseaseId', 'diseaseName','diseaseType']]
disease_gene_map = gene_disease_map.groupby('diseaseId').apply(Item_list,str = 'geneId').reset_index(drop = False)
disease_gene_map = DF_Dict(disease_gene_map)

### diseaseId --- ICD9 code map
disease_map  = pd.read_csv(r"Data\Code_Original\disease_mappings.tsv",sep = '\t')
disease_map.value_counts('vocabulary')
disease_OMIM_map = disease_map[disease_map['vocabulary'] == 'OMIM'].reset_index(drop = True)
disease_OMIM_map = disease_OMIM_map[['code','diseaseId']]
disease_OMIM_map = disease_OMIM_map.groupby('code').apply(Item_list,str = 'diseaseId').reset_index(drop = False)
OMIM_disease_map = DF_Dict(disease_OMIM_map)
OMIM_gene_map = Dict_Merge(OMIM_disease_map,disease_gene_map)
'''
[OMIM,geneId]: OMIM_gene_map 
[diseaseId,geneId]: disease_gene_map
'''
disease_ICD_map = disease_map[disease_map['vocabulary']  == 'ICD9CM'].reset_index(drop = True)
disease_ICD_map['Range'] = [1 if '-' in i else 0  for i in disease_ICD_map['code']]
disease_ICD_map_T = disease_ICD_map[disease_ICD_map['Range'] == 1].reset_index(drop = True) # 41
disease_ICD_map_F = disease_ICD_map[disease_ICD_map['Range'] != 1].reset_index(drop = True) # 4030

disease_ICD_map_T['code1'] = disease_ICD_map_T['code'].str.split('-').str.get(0) # 41
disease_ICD_map_T['code2'] = disease_ICD_map_T['code'].str.split('-').str.get(1)
disease_ICD_map_T['code2'] = disease_ICD_map_T['code2'].str.split('.').str.get(0)

disease_ICD_map_T1 = disease_ICD_map_T[disease_ICD_map_T['code1'] == disease_ICD_map_T['code2']].reset_index(drop = True)
disease_ICD_map_T2 = disease_ICD_map_T[disease_ICD_map_T['code1'] != disease_ICD_map_T['code2']].reset_index(drop = True)

disease_ICD_map_T2.loc[:35,'dict_ls'] = disease_ICD_map_T2.loc[:35,['code1','code2','diseaseId']].apply(disease_ICD_Dict,axis= 1)
disease_ICD_map_T2.loc[36,'dict_ls'] = str([['E' + str(i),'C0000921'] for i in range(880,889)])
disease_ICD_ls = [j for i in disease_ICD_map_T2['dict_ls'] for j in eval(i)]
disease_ICD_map_T2 = pd.DataFrame(disease_ICD_ls,columns=['code','diseaseId'])

disease_ICD_map_T1 = disease_ICD_map_T1[['code1','diseaseId']]
disease_ICD_map_T1.columns = ['code','diseaseId']

disease_ICD_map_T = pd.concat([disease_ICD_map_T1,disease_ICD_map_T2],axis=0).reset_index(drop = True)

disease_ICD_map_F['code'] = disease_ICD_map_F['code'].str.split('.').str.get(0)
disease_ICD_map_F = disease_ICD_map_F[['code','diseaseId']]

disease_ICD_map = pd.concat([disease_ICD_map_T,disease_ICD_map_F],axis = 0).reset_index(drop =True)
disease_ICD_map['code'] = disease_ICD_map['code'].apply(Int2Str)
disease_ICD_map = disease_ICD_map.drop_duplicates().reset_index(drop = True)
disease_ICD_map = disease_ICD_map.groupby('code').apply(Item_list,str = 'diseaseId').reset_index(drop = False)
ICD_disease_map = DF_Dict(disease_ICD_map)
ICD_gene_map = Dict_Merge(ICD_disease_map,disease_gene_map)
'''
[ICD,geneId]: ICD_gene_map
'''
disease_gene_map_df = {i:str(j) for i,j in disease_gene_map.items()}
OMIM_gene_map_df = {i:str(j) for i,j in OMIM_gene_map.items()}
ICD_gene_map_df = {i:str(j) for i,j in ICD_gene_map.items()}

disease_gene_map_df = pd.DataFrame.from_dict(disease_gene_map_df,orient='index').reset_index(drop = False)
OMIM_gene_map_df = pd.DataFrame.from_dict(OMIM_gene_map_df,orient='index').reset_index(drop = False)
ICD_gene_map_df = pd.DataFrame.from_dict(ICD_gene_map_df,orient='index').reset_index(drop = False)
disease_gene_map_df.columns = ['diseaseID','geneID']
OMIM_gene_map_df.columns = ['OMIMID','geneID']
ICD_gene_map_df.columns = ['ICDID','geneID']

disease_gene_map_df.to_csv(r'Data\Code_00_Outcome\disease_gene_map_df.csv',index = False)
OMIM_gene_map_df.to_csv(r'Data\Code_00_Outcome\OMIM_gene_map_df.csv',index = False)
ICD_gene_map_df.to_csv(r'Data\Code_00_Outcome\ICD_gene_map_df_p1.csv',index = False)


###################################################################################################
#### Additional Map Data     [DisGeNET] + [Manual Map Dataset] 
#### Manual Map Dataset
###################################################################################################
Manual_Map = pd.read_excel(r"Data\Code_Original\Map.xlsx")
Manual_Map['ICD-9-CM code'] =  Manual_Map['ICD-9-CM code'].str.strip('[]')
Manual_Map['ICD-9-CM code'] = Manual_Map['ICD-9-CM code'].str.split('.').str.get(0)
Manual_Map = Manual_Map[['Disease Code for OMIM from Goh et al. (2007)', 'ICD-9-CM code']]
Manual_Map.columns = ['Disease ID', 'code']

Manual_Map_G = Manual_Map.groupby('code').apply(Item_list,str = 'Disease ID').reset_index(drop = False)
ICD_Manual_Map = DF_Dict(Manual_Map_G)
'''
[ICD,Manual]: ICD_Manual_Map
'''
df = pd.DataFrame()
pdf =  pdfplumber.open(r'Data\Code_Original\01361table1.pdf') 
for i in range(62):
    page = pdf.pages[i]
    table = page.extract_table()
    if i == 0:
        column = table[1]
        add_df = table[2:]
    else: 
        add_df = table
    add_df = pd.DataFrame(add_df,columns=column)
    df = df.append(add_df,ignore_index=True)

df['Disease OMIM ID'] = df[['Disorder name','OMIM ID']].apply(Phenotype_Split,axis = 1)
df_map = df[['Disease ID','Disease OMIM ID']].drop_duplicates().reset_index(drop = True)
df_map['Disease ID'] = df_map['Disease ID'].astype('int')
df_map.shape # (2393, 2)
df_map_G = df_map.groupby('Disease ID').apply(Item_list,str = 'Disease OMIM ID').reset_index(drop = False)
Manual_OMIM_Map = DF_Dict(df_map_G)
ICD_OMIM_Map = Dict_Merge(ICD_Manual_Map,Manual_OMIM_Map)
ICD_gene_Map_2 = Dict_Merge(ICD_OMIM_Map,OMIM_gene_map)
'''
[Manual,OMIM ID]: Manual_OMIM_Map
[ICD,OMIM ID]: ICD_OMIM_Map
[ICD,geneId]: ICD_gene_Map_2
'''
ICD_gene_Map_2_df = {i:str(j) for i,j in ICD_gene_Map_2.items()}
ICD_gene_Map_2_df = pd.DataFrame.from_dict(ICD_gene_Map_2_df,orient='index').reset_index(drop = False)
ICD_gene_Map_2_df.columns = ['ICDID','geneID']
ICD_gene_Map_2_df.to_csv(r'Data\Code_00_Outcome\ICD_gene_map_df_p2.csv',index = False)


#####################################################################################################
ICD_gene_Map_df = pd.read_csv(r'Data\Code_00_Outcome\ICD_gene_map_df_p1.csv')
ICD_gene_Map_2_df = pd.read_csv(r'Data\Code_00_Outcome\ICD_gene_map_df_p2.csv')
ICD_gene_Map_df = pd.concat([ICD_gene_Map_df,ICD_gene_Map_2_df],axis = 0).reset_index(drop = True)
ICD_gene_Map_df = ICD_gene_Map_df.groupby('ICDID').apply(Item_list,str = 'geneID').reset_index(drop = False)
ICD_gene_Map_df.columns = ICD_gene_Map_2_df.columns
ICD_gene_Map_df['geneID'] = [sorted(list(set([k for j in i for k in eval(j)]))) for i in ICD_gene_Map_df['geneID']]
ICD_gene_Map_df['geneID_num'] = [len(i) for i in ICD_gene_Map_df['geneID']]
ICD_gene_Map_df = ICD_gene_Map_df[ICD_gene_Map_df['geneID_num'] != 0].reset_index(drop = True) # 539
ICD_gene_Map_df_add = pd.read_csv(r'Data\Code_Original\ICD_Protein_EntrezID_PPI_ed1.csv')

DF1_Dict = {i:j for i,j in zip(ICD_gene_Map_df['ICDID'],ICD_gene_Map_df['geneID'])}
DF2_Dict = {i:eval(j) for i,j in zip(ICD_gene_Map_df_add['ICD'],ICD_gene_Map_df_add['Protein'])}
ICD_set = set(ICD_gene_Map_df['ICDID'].tolist()).intersection(set(ICD_gene_Map_df_add['ICD'].tolist()))
print(len(set(ICD_gene_Map_df['ICDID'].tolist()) - set(ICD_gene_Map_df_add['ICD'].tolist())))
outcome = {}
for i in ICD_set:
    new = set(DF1_Dict[i]) - set(DF2_Dict[i])
    outcome[i] = new

ICD_gene_Map_df['geneID_new'] = 0
ICD_gene_Map_df['geneID_new'] = [outcome[i] if i in outcome.keys() else set() for i in ICD_gene_Map_df['ICDID']]
ICD_gene_Map_df['geneID_new_num'] = [len(i) for i in ICD_gene_Map_df['geneID_new']]
ICD_gene_Map_df[ICD_gene_Map_df['geneID_new_num'] > 0]
ICD_gene_Map_df['geneID_new_num'].mean()

ICD_gene_Map_df[['ICDID','geneID']]
ICD_gene_Map_df[['ICD','Protein']].dtypes

ICD_gene_Map_df['geneID'].tolist()[0]
ICD_gene_Map_df_add['Protein'].tolist()[0]
ICD_gene_Map_df.to_csv(r'Data\Code_00_Outcome\ICD_gene_map_df.csv',index = False)



##### Acknowledgement 
# OMIM DATASETS RELATED FUNCTION





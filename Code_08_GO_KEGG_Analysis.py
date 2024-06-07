'''
In this part, we compute the importance to evaluate the KEGG and GO Functions
'''
import pandas as pd 
import numpy as np 
import os

path1 = "Data\\Code_08_Outcome\\Source Protein"
path2 = "Data\\Code_08_Outcome\\Source Protein Add"
path3 = "Data\\Code_08_Outcome\\Target Protein"


def Folder_Analysis(Path):
    KEGG_path = Path + "\\KEGG"
    GO_path =  Path + "\\GO"

    #### Folder KEGG
    KEGG_List_DF = pd.DataFrame()
    Edge_KEGG_DF = []
    KEGG_file = os.listdir(KEGG_path)
    for KEGG_name in KEGG_file:
        KEGG_df = pd.read_csv(os.path.join(KEGG_path,KEGG_name))
        if KEGG_df.shape[0] == 0:
            continue
        KEGG_df = KEGG_df[KEGG_df['p.adjust']<0.05].reset_index(drop = True)
        KEGG_df['bg.gene'] = KEGG_df['BgRatio'].str.split('/').str.get(0).astype('int')
        #KEGG_df = KEGG_df[(KEGG_df['bg.gene']>20)&(KEGG_df['bg.gene']<200)].reset_index(drop = True)
        Edge_KEGG_DF += [[KEGG_name,i] for i in KEGG_df['ID']]
        ### KEGG_detail_dataframe
        KEGG_List_DF = pd.concat([KEGG_List_DF,KEGG_df[['category','subcategory','ID','Description']]],axis = 0).drop_duplicates().reset_index(drop=True)

    #### Folder KEGG
    GO_List_DF = pd.DataFrame()
    Edge_GO_DF = []
    GO_file = os.listdir(GO_path)
    for GO_name in GO_file:
        GO_df = pd.read_csv(os.path.join(GO_path,GO_name))
        if GO_df.shape[0] == 0:
            continue
        GO_df = GO_df[GO_df['p.adjust']<0.05].reset_index(drop = True)
        GO_df['bg.gene'] = GO_df['BgRatio'].str.split('/').str.get(0).astype('int')
        # GO_df = GO_df[(GO_df['bg.gene']>20)&(GO_df['bg.gene']<200)].reset_index(drop = True)
        Edge_GO_DF += [[GO_name,i] for i in GO_df['ID']]
        ### GO_detail_dataframe
        GO_List_DF = pd.concat([GO_List_DF,GO_df[['ONTOLOGY','ID','Description']]],axis = 0).drop_duplicates().reset_index(drop=True)

    Edge_GO_DF = pd.DataFrame(Edge_GO_DF,columns= ['Edge','GO'])
    Edge_KEGG_DF = pd.DataFrame(Edge_KEGG_DF,columns= ['Edge','KEGG'])

    return Edge_GO_DF,Edge_KEGG_DF,GO_List_DF,KEGG_List_DF

Edge_GO_DF_S,Edge_KEGG_DF_S,GO_List_DF_S,KEGG_List_DF_S = Folder_Analysis(path1)
Edge_GO_DF_S_add,Edge_KEGG_DF_S_add,GO_List_DF_S_add,KEGG_List_DF_S_add = Folder_Analysis(path2)
Edge_GO_DF_T,Edge_KEGG_DF_T,GO_List_DF_T,KEGG_List_DF_T = Folder_Analysis(path3)

GO_List = pd.concat([GO_List_DF_S,GO_List_DF_S_add,GO_List_DF_T],axis = 0).drop_duplicates().reset_index(drop = True)
KEGG_List = pd.concat([KEGG_List_DF_S,KEGG_List_DF_S_add,KEGG_List_DF_T],axis = 0).drop_duplicates().reset_index(drop =True)

def Edge_Columen(DF):
    DF['Edge1'] = DF['Edge'].str.split('_').str.get(1)
    DF['Edge2'] = DF['Edge'].str.split('_').str.get(2)
    DF['Edge'] = DF['Edge1'].str.cat(DF['Edge2'],sep = '_')
    DF = DF.drop(columns = ['Edge1','Edge2'])
    DF.columns = ['Edge','ID']
    DF['Ref'] = DF['Edge'].str.cat(DF['ID'],sep = '_')
    return DF

Edge_GO_DF_S = Edge_Columen(Edge_GO_DF_S)
Edge_KEGG_DF_S = Edge_Columen(Edge_KEGG_DF_S)
Edge_GO_DF_S_add = Edge_Columen(Edge_GO_DF_S_add)
Edge_KEGG_DF_S_add = Edge_Columen(Edge_KEGG_DF_S_add)
Edge_GO_DF_T = Edge_Columen(Edge_GO_DF_T)
Edge_KEGG_DF_T = Edge_Columen(Edge_KEGG_DF_T)

Edge_GO_DF_S_add = Edge_GO_DF_S_add.drop(index = Edge_GO_DF_S_add[Edge_GO_DF_S_add['Ref'].isin(Edge_GO_DF_S['Ref'])].index).reset_index(drop = True)
GO_overlap = Edge_GO_DF_S_add[Edge_GO_DF_S_add['Ref'].isin(Edge_GO_DF_T['Ref'])].reset_index(drop = True)
GO_overlap

Edge_KEGG_DF_S_add = Edge_KEGG_DF_S_add.drop(index = Edge_KEGG_DF_S_add[Edge_KEGG_DF_S_add['Ref'].isin(Edge_KEGG_DF_S['Ref'])].index).reset_index(drop = True)
KEGG_overlap = Edge_KEGG_DF_S_add[Edge_KEGG_DF_S_add['Ref'].isin(Edge_KEGG_DF_T['Ref'])].reset_index(drop = True)

KEGG_overlap = pd.merge(KEGG_overlap,KEGG_List,on = 'ID',how = 'left')
GO_overlap = pd.merge(GO_overlap,GO_List,on = 'ID',how = 'left')


comornet = pd.read_csv('Data\\Code_03_Outcome\\network_pairwise.csv')
corr_dict = {i:j for i,j in zip(comornet['edge'],comornet['Diff_Corr'])}
rr_dict = {i:j for i,j in zip(comornet['edge'],comornet['Diff_RR'])}
KEGG_overlap['Diff_Corr'] = [corr_dict[i] for i in KEGG_overlap['Edge']]
KEGG_overlap['Diff_RR'] = [rr_dict[i] for i in KEGG_overlap['Edge']]
GO_overlap['Diff_Corr'] = [corr_dict[i] for i in GO_overlap['Edge']]
GO_overlap['Diff_RR'] = [rr_dict[i] for i in GO_overlap['Edge']]


KEGG_overlap_RR = KEGG_overlap.groupby('ID')['Diff_RR'].sum().reset_index(drop = False)
KEGG_overlap_Corr = KEGG_overlap.groupby('ID')['Diff_Corr'].sum().reset_index(drop = False)
KEGG_overlap_Fre = KEGG_overlap.value_counts('ID').reset_index(drop = False)
KEGG_overlap_Fre.columns =['ID','Frequency']
KEGG_overlap_RR = pd.merge(KEGG_overlap_RR,KEGG_overlap_Corr,on = 'ID',how ='left')
KEGG_overlap_RR = pd.merge(KEGG_overlap_RR,KEGG_overlap_Fre,on = 'ID',how= 'left')

GO_overlap_RR = GO_overlap.groupby('ID')['Diff_RR'].sum().reset_index(drop = False)
GO_overlap_Corr = GO_overlap.groupby('ID')['Diff_Corr'].sum().reset_index(drop = False)
GO_overlap_Fre = GO_overlap.value_counts('ID').reset_index(drop = False)
GO_overlap_Fre.columns =['ID','Frequency']
GO_overlap_RR = pd.merge(GO_overlap_RR,GO_overlap_Corr,on = 'ID',how ='left')
GO_overlap_RR = pd.merge(GO_overlap_RR,GO_overlap_Fre,on = 'ID',how= 'left')

def Element_Index(input,reverse_ = True):
    ref_ls = sorted(input,reverse = reverse_) 
    outcome = [ref_ls.index(i) for i in input]
    return outcome

GO_overlap_RR['RR_idx'] = Element_Index(GO_overlap_RR['Diff_RR'])
GO_overlap_RR['Corr_idx'] = Element_Index(GO_overlap_RR['Diff_Corr'])
GO_overlap_RR['mean_idx'] = GO_overlap_RR[['RR_idx','Corr_idx']].mean(axis= 1)
GO_overlap_RR['idx'] = Element_Index(GO_overlap_RR['mean_idx'],reverse_=False)
GO_overlap_RR = GO_overlap_RR.sort_values('idx',ascending=True).reset_index(drop = True)
GO_overlap_RR = pd.merge(GO_overlap_RR,GO_List,on = 'ID',how = 'left')

KEGG_overlap_RR['RR_idx'] = Element_Index(KEGG_overlap_RR['Diff_RR'])
KEGG_overlap_RR['Corr_idx'] = Element_Index(KEGG_overlap_RR['Diff_Corr'])
KEGG_overlap_RR['mean_idx'] = KEGG_overlap_RR[['RR_idx','Corr_idx']].mean(axis= 1)
KEGG_overlap_RR['idx'] = Element_Index(KEGG_overlap_RR['mean_idx'],reverse_=False)
KEGG_overlap_RR = KEGG_overlap_RR.sort_values('idx',ascending=True).reset_index(drop = True)
KEGG_overlap_RR = pd.merge(KEGG_overlap_RR,KEGG_List,on = 'ID',how = 'left')

KEGG_overlap_RR.to_csv("Data\\Code_08_Outcome\\KEGG_overlap_v2.csv",index = False)
GO_overlap_RR.to_csv("Data\\Code_08_Outcome\\GO_overlap_v2.csv",index = False)
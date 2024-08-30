'''
In this part, we compute the importance to evaluate the KEGG and GO Functions
'''
import pandas as pd 
import numpy as np 
import os

path1 = "Source Protein"
path2 = "Source Protein Add"
path3 = "Target Protein"

def Folder_Analysis(Path):
    KEGG_path = "Data\\Code_06_Outcome\\KEGG\\" + Path
    GO_path = "Data\\Code_06_Outcome\\GO\\" + Path

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
        # KEGG_df = KEGG_df[(KEGG_df['bg.gene']>20)&(KEGG_df['bg.gene']<200)].reset_index(drop = True)
        Edge_KEGG_DF += [[KEGG_name.split('_')[1],i] for i in KEGG_df['ID']]
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
        Edge_GO_DF += [[GO_name.split('_')[1],i] for i in GO_df['ID']]
        ### GO_detail_dataframe
        GO_List_DF = pd.concat([GO_List_DF,GO_df[['ONTOLOGY','ID','Description']]],axis = 0).drop_duplicates().reset_index(drop=True)
    Edge_GO_DF = pd.DataFrame(Edge_GO_DF,columns= ['Edge','GO'])
    Edge_KEGG_DF = pd.DataFrame(Edge_KEGG_DF,columns= ['Edge','KEGG'])
    return Edge_GO_DF,Edge_KEGG_DF,GO_List_DF,KEGG_List_DF

Edge_GO_DF_S,Edge_KEGG_DF_S,GO_List_DF_S,KEGG_List_DF_S = Folder_Analysis(path1)
Edge_GO_DF_S_add,Edge_KEGG_DF_S_add,GO_List_DF_S_add,KEGG_List_DF_S_add = Folder_Analysis(path2)
Edge_GO_DF_T,Edge_KEGG_DF_T,GO_List_DF_T,KEGG_List_DF_T = Folder_Analysis(path3)

### GO / KEGG Information
GO_List = pd.concat([GO_List_DF_S,GO_List_DF_S_add,GO_List_DF_T],axis = 0).drop_duplicates().reset_index(drop = True)
KEGG_List = pd.concat([KEGG_List_DF_S,KEGG_List_DF_S_add,KEGG_List_DF_T],axis = 0).drop_duplicates().reset_index(drop =True)

Edge_GO_DF_S['Ref'] = Edge_GO_DF_S[['Edge', 'GO']].apply(lambda x: x[0] + '_' + x[1], axis = 1)
Edge_KEGG_DF_S['Ref'] = Edge_KEGG_DF_S[['Edge', 'KEGG']].apply(lambda x: x[0] + '_' + x[1], axis = 1)
Edge_GO_DF_S_add['Ref'] = Edge_GO_DF_S_add[['Edge', 'GO']].apply(lambda x: x[0] + '_' + x[1], axis = 1)
Edge_KEGG_DF_S_add['Ref'] = Edge_KEGG_DF_S_add[['Edge', 'KEGG']].apply(lambda x: x[0] + '_' + x[1], axis = 1)
Edge_GO_DF_T['Ref'] = Edge_GO_DF_T[['Edge', 'GO']].apply(lambda x: x[0] + '_' + x[1], axis = 1)
Edge_KEGG_DF_T['Ref'] = Edge_KEGG_DF_T[['Edge', 'KEGG']].apply(lambda x: x[0] + '_' + x[1], axis = 1)

Edge_GO_DF_S_add = Edge_GO_DF_S_add.drop(index = Edge_GO_DF_S_add[Edge_GO_DF_S_add['Ref'].isin(Edge_GO_DF_S['Ref'])].index).reset_index(drop = True)
GO_overlap = Edge_GO_DF_S_add[Edge_GO_DF_S_add['Ref'].isin(Edge_GO_DF_T['Ref'])].reset_index(drop = True)
GO_overlap

Edge_KEGG_DF_S_add = Edge_KEGG_DF_S_add.drop(index = Edge_KEGG_DF_S_add[Edge_KEGG_DF_S_add['Ref'].isin(Edge_KEGG_DF_S['Ref'])].index).reset_index(drop = True)
KEGG_overlap = Edge_KEGG_DF_S_add[Edge_KEGG_DF_S_add['Ref'].isin(Edge_KEGG_DF_T['Ref'])].reset_index(drop = True)

KEGG_overlap = pd.merge(KEGG_overlap,KEGG_List,on = 'ID',how = 'left') # 0 
GO_overlap = pd.merge(GO_overlap,GO_List,left_on = 'GO', right_on = 'ID',how = 'left') # 323
GO_overlap.value_counts('Edge')


comornet = pd.read_csv(r'Data\Code_01_Outcome\network_pairwise(2).csv')
comornet.rename(columns = {"edge":"Edge"}, inplace = True)
GO_overlap = pd.merge(GO_overlap, comornet[['Edge','Diff_RR', 'Diff_Corr']], on = 'Edge', how = 'left')
GO_overlap.value_counts('GO') # 145

class TF_IDF():
    def __init__(self,DF) -> None:
        self.DF = DF

    def TF(self):
        DF = self.DF
        DF_ = DF.value_counts('Edge').reset_index(drop = False)
        DF_.columns = ['Edge','Count']
        DF = pd.merge(DF,DF_,on = 'Edge',how = 'left')
        DF['TF(Corr)'] = DF['Diff_Corr']/DF['Count']
        DF['TF(RR)'] = DF['Diff_RR']/DF['Count']
        DF['TF(Count)'] = 1/DF['Count']
        self.TF = DF

    def IDF(self):
        DF = self.DF

        DF_Corr = DF.groupby('ID')['Diff_Corr'].sum().reset_index(drop = False)
        DF_RR = DF.groupby('ID')['Diff_RR'].sum().reset_index(drop = False)
        DF_Count = DF.value_counts('ID').reset_index(drop = False)    

        DF_ = pd.merge(DF_Corr,DF_RR,on = 'ID',how ='left')
        DF_ = pd.merge(DF_,DF_Count,on = 'ID',how = 'left')
        DF_.columns = ['ID','Diff_Corr','Diff_RR','Frequency']
        # smoothing
        # DF__ = DF[['Edge', 'Diff_Corr', 'Diff_RR']].drop_duplicates().reset_index(drop = True)
        # DF_['IDF(Corr)'] = DF_[['Diff_Corr', 'Frequency']].apply(lambda x: np.log(DF__['Diff_Corr'].sum()/(x[0] * (1 + 1/x[1]))),axis=1) 
        # DF_['IDF(RR)'] = DF_[['Diff_RR','Frequency']].apply(lambda x: np.log(DF__['Diff_RR'].sum()/(x[0] * (1 + 1/x[1]))),axis=1)
        # DF_['IDF(Count)'] = DF_['Frequency'].apply(lambda x: np.log(DF__.shape[0]/(x + 1)))
        DF__ = DF[['Edge', 'Diff_Corr', 'Diff_RR']].drop_duplicates().reset_index(drop = True)
        DF_['IDF(Corr)'] = DF_['Diff_Corr'].apply(lambda x: np.log(DF__['Diff_Corr'].sum()/x)) 
        DF_['IDF(RR)'] = DF_['Diff_RR'].apply(lambda x: np.log(DF__['Diff_RR'].sum()/x))
        DF_['IDF(Count)'] = DF_['Frequency'].apply(lambda x: np.log(DF__.shape[0]/x))
        self.IDF = DF_

    def TFIDF_DF(self):
        TF = self.TF
        IDF = self.IDF

        DF = pd.merge(TF, IDF[['ID', 'IDF(Corr)', 'IDF(RR)','IDF(Count)']], on='ID', how='left')
        DF['TFIDF(Corr)'] = DF[['TF(Corr)','IDF(Corr)']].apply(lambda x: x[0] * x[1], axis =1)
        DF['TFIDF(RR)'] = DF[['TF(RR)','IDF(RR)']].apply(lambda x: x[0] * x[1], axis =1)
        DF['TFIDF(Count)'] = DF[['TF(Count)','IDF(Count)']].apply(lambda x: x[0] * x[1], axis =1)

        self.TFIDF = DF

GO_TFIDF = TF_IDF(GO_overlap)
GO_TFIDF.TF()
GO_TFIDF.IDF()
GO_TFIDF.TFIDF_DF()
GO_tfidf = GO_TFIDF.TFIDF
GO_tfidf.to_csv("Data\\Code_07_Outcome\\GO_TFIDF.csv",index = False)


GO_Count_G = GO_tfidf.groupby('ID')['TFIDF(Count)'].sum().reset_index(drop = False)
GO_RR_G = GO_tfidf.groupby('ID')['TFIDF(RR)'].sum().reset_index(drop = False)
GO_Corr_G = GO_tfidf.groupby('ID')['TFIDF(Corr)'].sum().reset_index(drop = False)
GO_G = pd.merge(GO_Count_G,GO_RR_G,on = 'ID',how = 'left')
GO_G = pd.merge(GO_G,GO_Corr_G,on = 'ID',how = 'left')

def Element_Index(input,reverse_ = True):
    ref_ls = sorted(input,reverse = reverse_) 
    outcome = [ref_ls.index(i) for i in input]
    return outcome

GO_G['RR_idx'] = Element_Index(GO_G['TFIDF(RR)'])
GO_G['Corr_idx'] = Element_Index(GO_G['TFIDF(Corr)'])
GO_G['Count_idx'] = Element_Index(GO_G['TFIDF(Count)'])
GO_G['mean_idx'] = GO_G[['RR_idx','Corr_idx','Count_idx']].mean(axis= 1)
GO_G['idx'] = Element_Index(GO_G['mean_idx'],reverse_=False)
GO_G = GO_G.sort_values('idx',ascending=True).reset_index(drop = True)
GO_G.rename(columns = {'ID':'GO Term', 'mean_idx':'average of idx (value)'},inplace = True)
GO_G = pd.merge(GO_G,GO_List,left_on = 'GO Term', right_on='ID',how = 'left')
GO_G.to_csv("Data\\Code_07_Outcome\\GO_Total_TFIDF.csv",index = False)
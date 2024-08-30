import pandas as pd
import numpy as np 
import os 
comornet = pd.read_csv(r'Data\Code_01_Outcome\network_pairwise(2).csv')
disease_pairs = pd.read_csv('Data/Code_03_Outcome/ICD_pair_T_add(2).csv')
disease_pairs.drop(columns = ['outcome1','outcome2','outcome3','outcome1_add','outcome2_add','outcome3_add'], inplace = True)
disease_pairs = pd.merge(disease_pairs, comornet[['edge', 'Diff_RR', 'Diff_Corr']], on = 'edge', how = 'left')
disease_pairs.rename(columns = {'source_portein':'source_protein' ,'source_portein_add':'source_protein_add' ,'target_portein':'target_protein'}, inplace = True)

disease_pairs['overlap'] = disease_pairs[['source_protein_add','target_protein']].apply(lambda x: str(list(set(eval(x[0])).intersection(set(eval(x[1]))))),axis = 1)
disease_pairs = disease_pairs[disease_pairs['overlap'] !='[]'].reset_index(drop=True) # 47 disease pair 

class TF_IDF():
    def __init__(self,DF) -> None:
        self.DF = DF

    def TF(self):
        TF_dict = {}
        TF_sum_dict = {}
        DF = self.DF
        for idx in DF.index:
            protein_dict = {}
            overlap_protein = eval(DF.loc[idx,'overlap'])
            Diff_RR = DF.loc[idx,'Diff_RR']
            Diff_Corr = DF.loc[idx,'Diff_Corr']
            for protein in overlap_protein:
                protein_dict[protein] = {'TF(count)': 1/len(overlap_protein), 'TF(RR)': Diff_RR/len(overlap_protein), 'TF(Corr)':Diff_Corr/len(overlap_protein)}
                if protein in TF_sum_dict.keys():
                    TF_sum_dict[protein]['Sum (count)'] += 1
                    TF_sum_dict[protein]['Sum (RR)'] += Diff_RR
                    TF_sum_dict[protein]['Sum (Corr)'] += Diff_Corr
                else:
                    TF_sum_dict[protein] = {'Sum (count)': 1, 'Sum (RR)': Diff_RR, 'Sum (Corr)': Diff_Corr}
            TF_dict[DF.loc[idx,'edge']] = protein_dict
        self.TF_dict = TF_dict
        self.TF_sum_dict = TF_sum_dict

    def IDF(self):
        DF = self.DF
        TF_sum_dict = self.TF_sum_dict
        total_rr = DF['Diff_RR'].sum()
        total_corr = DF['Diff_Corr'].sum()
        total_num = DF.shape[0]
        IDF_dict = {}
        for protein in TF_sum_dict.keys():
            ### smoothing
            # IDF_dict[protein] = {'IDF(count)': np.log(total_num/(TF_sum_dict[protein]['Sum (count)']+1)),
            #                      'IDF(RR)': np.log(total_rr/(TF_sum_dict[protein]['Sum (RR)'] * (1 + 1/TF_sum_dict[protein]['Sum (count)']))),
            #                      'IDF(Corr)': np.log(total_corr/(TF_sum_dict[protein]['Sum (Corr)'] * (1 + 1/TF_sum_dict[protein]['Sum (count)'])))}
            IDF_dict[protein] = {'IDF(count)': np.log(total_num/(TF_sum_dict[protein]['Sum (count)'])),
                                 'IDF(RR)': np.log(total_rr/TF_sum_dict[protein]['Sum (RR)']),
                                 'IDF(Corr)': np.log(total_corr/TF_sum_dict[protein]['Sum (Corr)'])}            
        self.IDF_dict = IDF_dict

    def TFIDF(self):
        TF_dict = self.TF_dict
        IDF_dict = self.IDF_dict
        TFIDF_dict = {}
        for edge in TF_dict.keys():
            TFIDF_dict[edge] = {}
            for protein in TF_dict[edge].keys():
                TFIDF_dict[edge][protein] = {'TFIDF(count)': TF_dict[edge][protein]['TF(count)'] * IDF_dict[protein]['IDF(count)'],
                                             'TFIDF(RR)': TF_dict[edge][protein]['TF(RR)'] * IDF_dict[protein]['IDF(RR)'],
                                             'TFIDF(Corr)': TF_dict[edge][protein]['TF(Corr)'] * IDF_dict[protein]['IDF(Corr)']}
        self.TFIDF_dict = TFIDF_dict

    def TFIDF_DF(self):
        TFIDF_dict = self.TFIDF_dict
        TFIDF_df = []
        for edge in TFIDF_dict.keys():
            for protein in TFIDF_dict[edge].keys():
                TFIDF_df.append([edge, protein, TFIDF_dict[edge][protein]['TFIDF(count)'],
                                 TFIDF_dict[edge][protein]['TFIDF(RR)'], TFIDF_dict[edge][protein]['TFIDF(Corr)']])
        TFIDF_df = pd.DataFrame(TFIDF_df, columns=['edge', 'protein', 'TFIDF(count)', 'TFIDF(RR)','TFIDF(Corr)'])
        self.TFIDF_df = TFIDF_df
import mygene
def entrez_id_to_gene_symbol(entrez_id):
    mg = mygene.MyGeneInfo()
    result = mg.getgene(entrez_id, fields='symbol')
    gene_symbol = result.get('symbol')
    return gene_symbol

tf_idf = TF_IDF(disease_pairs)
tf_idf.TF()
tf_idf.IDF()
tf_idf.TFIDF() 
tf_idf.TFIDF_DF()
tdidf_df = tf_idf.TFIDF_df
tdidf_df = tdidf_df.sort_values(['edge','TFIDF(count)'],ascending=False).reset_index(drop = True)
tdidf_df['symbol'] = [entrez_id_to_gene_symbol(str(i)) for i in tdidf_df['protein']]
tdidf_df.to_csv("Data/Code_06_Outcome/Protein_TFIDF.csv",index = False)

Protein_Count_G = tdidf_df.groupby('protein')['TFIDF(count)'].sum().reset_index(drop = False)
Protein_RR_G = tdidf_df.groupby('protein')['TFIDF(RR)'].sum().reset_index(drop = False)
Protein_Corr_G = tdidf_df.groupby('protein')['TFIDF(Corr)'].sum().reset_index(drop = False)

Protein_G = pd.merge(Protein_Count_G,Protein_RR_G,on = 'protein',how = 'left')
Protein_G = pd.merge(Protein_G,Protein_Corr_G,on = 'protein',how = 'left')

Protein_G = pd.merge(Protein_G, tdidf_df[['protein', 'symbol']].drop_duplicates(), on = 'protein', how = 'left')
Protein_G.to_csv("Data/Code_06_Outcome/Protein_Frequency_TFIDF.csv",index = False)
########################################################################################################
# Protein Index 
def Element_Index(input,reverse_ = True):
    ref_ls = sorted(input,reverse = reverse_)
    outcome = [ref_ls.index(i) for i in input]
    return outcome
Protein_G['RR_idx'] = Element_Index(Protein_G['TFIDF(RR)'])
Protein_G['Corr_idx'] = Element_Index(Protein_G['TFIDF(Corr)'])
Protein_G['Count_idx'] = Element_Index(Protein_G['TFIDF(count)'])
Protein_G['mean_idx'] = Protein_G[['RR_idx','Corr_idx','Count_idx']].mean(axis= 1)
Protein_G['idx'] = Element_Index(Protein_G['mean_idx'],reverse_=False)
Protein_G =Protein_G.sort_values('idx',ascending=True).reset_index(drop = True)
Protein_G.rename(columns = {'protein':'entrez_id', 'mean_idx': 'average value (idx)'},inplace = True)
Protein_G.to_csv("Data/Code_06_Outcome/Protein_Total_TFIDF.csv",index = False)
########################################################################################################
# Based on Protein_Corr_G into Source Protein, Source Protein Add, Target Protein 
path1 = "Data\\Code_06_Outcome\\Source protein"
path2 = "Data\\Code_06_Outcome\\Source protein Add"
path3 = "Data\\Code_06_Outcome\\Target Protein"

def List2DF(list,col):
    DF = pd.DataFrame(list,columns=col)
    return DF

def DF_Split(DF,path1,path2,path3):
    for idx in DF.index:
        edge = DF.loc[idx,'edge']
        source_protein = [[edge,i] for i in eval(DF.loc[idx,'source_protein'])]
        source_protein_add =[[edge,i] for i in eval( DF.loc[idx,'source_protein_add'])] 
        target_protein = [[edge,i] for i in eval(DF.loc[idx,'target_protein'])]
        source_protein_df = List2DF(source_protein,['edge','protein'])
        source_protein_add_df = List2DF(source_protein_add,['edge','protein'])       
        target_protein_df = List2DF(target_protein,['edge','protein'])
        source_protein_df.to_csv(os.path.join(path1,f"{edge}_source_protein.csv"),index = False)
        source_protein_add_df.to_csv(os.path.join(path2,f"{edge}_source_protein_add.csv"),index = False)
        target_protein_df.to_csv(os.path.join(path3,f"{edge}_target_protein.csv"),index = False)

DF_Split(disease_pairs,path1,path2,path3)
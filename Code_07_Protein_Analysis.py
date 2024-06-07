import pandas as pd
import numpy as np 
import os 

Pair = pd.read_csv("Data\\Code_06_Outcome\\Pair_T_significant_Diff.csv")
COVID_LCC = pd.read_csv("Data\Code_02_Outcome\COVID_LCC.csv")
comornet = pd.read_csv(r'Data\Code_01_Outcome\network_pairwise.csv')
comornet['edge'] = comornet['edge'].apply(lambda x: str(eval(eval(x)[0])[0])).str.cat(comornet['edge'].apply(lambda x: str(eval(x)[1])), sep = '_')
Pair = pd.merge(Pair,comornet[['edge','Correlation_treat','RR_treat','Diff_RR','Diff_Corr']],on = 'edge',how = 'left')
Pair.rename(columns = {'source_portein_add':'source_protein_add', 'target_portein':'target_protein','source_portein':'source protein'}, inplace = True)
def Overlap_Protein(input):
    s = eval(input[0])
    t = eval(input[1])
    overlap = sorted(list(set(s).intersection(set(t))))
    return str(overlap)

Pair['overlap'] = Pair[['source_protein_add','target_protein']].apply(Overlap_Protein,axis = 1)
Pair= Pair[Pair['overlap'] !='[]'].reset_index(drop=True) # 75 disease pair 

Protein_RR = [[k,j] for i,j in zip(Pair['overlap'],Pair['Diff_RR']) for k in eval(i)]
Protein_RR = pd.DataFrame(Protein_RR,columns = ['Protein','RR'])
Protein_RR_G = Protein_RR.groupby('Protein')['RR'].sum().reset_index(drop = False)

Protein_Corr = [[k,j] for i,j in zip(Pair['overlap'],Pair['Diff_Corr']) for k in eval(i)]
Protein_Corr = pd.DataFrame(Protein_Corr,columns = ['Protein','Corr'])
Protein_Corr_G = Protein_Corr.groupby('Protein')['Corr'].sum().reset_index(drop = False)

Protein_Corr_G = Protein_Corr_G.sort_values('Corr',ascending=False).reset_index(drop =True)
Protein_RR_G = Protein_RR_G.sort_values('RR',ascending=False).reset_index(drop =True)
Protein_Corr_G = pd.merge(Protein_Corr_G,Protein_RR_G,on = 'Protein',how = 'left')

import mygene
def entrez_id_to_gene_symbol(entrez_id):
    mg = mygene.MyGeneInfo()
    result = mg.getgene(entrez_id, fields='symbol')
    gene_symbol = result.get('symbol')
    return gene_symbol
Protein_Corr_G['symbol'] = [entrez_id_to_gene_symbol(str(i)) for i in Protein_Corr_G['Protein']]
Protein_Corr_G.to_csv("Data\\Code_07_Outcome\\Protein_Frequency.csv",index = False)

########################################################################################################
# Protein Index 
def Element_Index(input,reverse_ = True):
    ref_ls = sorted(input,reverse = reverse_) 
    outcome = [ref_ls.index(i) for i in input]
    return outcome

Protein_Corr_G = pd.read_csv("Data\\Code_07_Outcome\\Protein_Frequency.csv")
Protein_Corr_G['RR_idx'] = Element_Index(Protein_Corr_G['RR'])
Protein_Corr_G['Corr_idx'] = Element_Index(Protein_Corr_G['Corr'])
Protein_Corr_G['mean_idx'] = Protein_Corr_G[['RR_idx','Corr_idx']].mean(axis= 1)
Protein_Corr_G['idx'] = Element_Index(Protein_Corr_G['mean_idx'],reverse_=False)
Protein_Corr_G =Protein_Corr_G.sort_values('idx',ascending=True).reset_index(drop = True)
Protein_Corr_G.to_csv("Data\\Code_07_Outcome\\Protein_Frequency.csv",index = False)

########################################################################################################
# Based on Protein_Corr_G into Source Protein, Source Protein Add, Target Protein 

Pair = pd.read_csv(".\\Data\\Code_06_Outcome\\Pair_T_significant.csv")
path1 = "Data\\Code_07_Outcome\\Source Protein"
path2 = "Data\\Code_07_Outcome\\Source Protein Add"
path3 = "Data\\Code_07_Outcome\\Target Protein"


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

DF_Split(Pair,path1,path2,path3)



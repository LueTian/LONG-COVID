'''
2023-12-3
In this part, we analysis Positive PPI distance; Negative PPI distance; Repeatation PPI distance for positive sample
'''

import pandas as pd 
import numpy as np 
import scipy.stats as stats
import statsmodels.stats.weightstats as sw 
import matplotlib.pyplot as plt
import seaborn as sns
import os 

##############################################
# data load
Pair_T= pd.read_csv("Data\\Code_03_Outcome\\ICD_pair_T_add.csv")

len(set(Pair_T['source'].tolist() + Pair_T['target'].tolist()))
Pair_T.value_counts('source').shape
Pair_T.value_counts('target').shape

Pair_F= pd.read_csv("Data\\Code_03_Outcome\\ICD_pair_F_add.csv")
comornet = pd.read_csv('Data\\Code_01_Outcome\\network_pairwise.csv')

Pair_T['Diff_outcome1'] = Pair_T['outcome1'] - Pair_T['outcome1_add']
Pair_F['Diff_outcome1'] = Pair_F['outcome1'] - Pair_F['outcome1_add']

comornet['edge'] = comornet['edge'].apply(lambda x: str(eval(eval(x)[0])[0])).str.cat(comornet['edge'].apply(lambda x: str(eval(x)[1])), sep = '_')
Pair_T = pd.merge(Pair_T,comornet[['edge','Correlation_treat','RR_treat','Diff_RR','Diff_Corr']],on = 'edge',how = 'left')

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
### [Part One]overall statistic test  基于 negative / positive disease pairs 的 PPI distance 进行假设检验 
Pair_F['Diff_outcome1'].mean()
Pair_T['Diff_outcome1'].mean()
stats.mannwhitneyu(Pair_F['Diff_outcome1'].tolist(),Pair_T['Diff_outcome1'].tolist(),alternative='less') 
# MannwhitneyuResult(statistic=4537864.0, pvalue=0.04920759367425814)

boxplot_data = [Pair_F['Diff_outcome1'].tolist(),Pair_T['Diff_outcome1'].tolist()]
boxplot_labels = ['Negative Sample','Positive Sample']
sns.boxplot(data = boxplot_data)
plt.xlabel('Sample')
plt.xticks([0, 1], ['Negative Sample', 'Positive Sample'])

plt.ylabel('PPI Distance Reduction')
plt.title('Boxplot about PPI Distance Reduction in Negative and Positive Sample')
plt.legend(boxplot_labels)
plt.show()
plt.savefig('Data\\Code_06_Outcome\\boxplot.svg')

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
### [Part Three] overall statistic test  基于  positive disease pairs 的 PPI distance 和 重复检验的结果 
file_path = 'Data\\Code_05_Outcome'
file_ls = os.listdir('Data\\Code_05_Outcome')
ref_list = []
for file_name in file_ls:
    df = pd.read_csv(os.path.join(file_path,file_name))    
    edge = df.loc[0,'edge']
    df['Diff_outcome1'] = df['outcome1'] - df['outcome1_add']
    mean_= df['Diff_outcome1'].mean()
    std_  = df['Diff_outcome1'].std(ddof=1)
    value_ = df['Diff_outcome1'].tolist()
    ref_list.append([edge,mean_,std_,str(value_)])
ref_list = pd.DataFrame(ref_list,columns=['edge','bg_mean','bg_std','Diff_repeat'])

Pair_T= pd.read_csv("Data\\Code_03_Outcome\\ICD_pair_T_add.csv")
Pair_T = pd.merge(Pair_T,ref_list,on = 'edge',how = 'left')
Pair_T['Diff_outcome1'] = Pair_T['outcome1'] - Pair_T['outcome1_add']
def Z_Score(input):
    real = input[0]
    m = input[1]
    std = input[2]
    repeat = eval(input[3])
    p = np.mean([1 if i > real else 0 for i in repeat])
    z = (real - m)/std
    return z, p

Pair_T[['z_score','p_value']] = Pair_T[['Diff_outcome1','bg_mean','bg_std','Diff_repeat']].apply(Z_Score,axis = 1,result_type='expand')
Pair_T[Pair_T['p_value'] < 0.05]
Pair_T.to_csv("Data\\Code_06_Outcome\\Pair_T_significant_Diff.csv",index = False)


file_path = 'Data\\Code_05_Outcome'
file_ls = os.listdir('Data\\Code_05_Outcome')
ref_list = []
for file_name in file_ls:
    df = pd.read_csv(os.path.join(file_path,file_name))    
    edge = df.loc[0,'edge']
    df['Diff_outcome1'] = df['outcome1'] - df['outcome1_add']
    mean_= df['outcome1_add'].mean()
    std_  = df['outcome1_add'].std(ddof=1)
    value_ = df['outcome1_add'].tolist()
    ref_list.append([edge,mean_,std_,str(value_)])
ref_list = pd.DataFrame(ref_list,columns=['edge','bg_mean','bg_std','Diff_repeat'])
Pair_T= pd.read_csv("Data\\Code_03_Outcome\\ICD_pair_T_add.csv")
Pair_T = pd.merge(Pair_T,ref_list,on = 'edge',how = 'left')
Pair_T['Diff_outcome1'] = Pair_T['outcome1'] - Pair_T['outcome1_add']
def Z_Score(input):
    real = input[0]
    m = input[1]
    std = input[2]
    repeat = eval(input[3])
    p = np.mean([1 if i > real else 0 for i in repeat])
    z = (real - m)/std # 如果 Z 小于 0 说明 加入 COVID Protein 之后 显著变小
    return z, p

Pair_T[['z_score','p_value']] = Pair_T[['outcome1_add','bg_mean','bg_std','Diff_repeat']].apply(Z_Score,axis = 1,result_type='expand')
Pair_T['p_value'] = 1 - Pair_T['p_value']
Pair_T[Pair_T['p_value'] < 0.05]
### 包括P-value 和 Z-score
Pair_T.to_csv("Data\\Code_06_Outcome\\Pair_T_significant_add.csv",index = False)


'''
2023-12-3
In this part, we analysis Positive PPI distance; Negative PPI distance; Repeatation PPI distance for positive sample
'''

import pandas as pd 
import numpy as np 
import scipy.stats as stats
import statsmodels.stats.weightstats as sw 
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt
import seaborn as sns
import openchord as ocd
import os 

##############################################
# data load
Pair_T= pd.read_csv("Data\\Code_03_Outcome\\ICD_pair_T_add(2).csv")
Pair_F= pd.read_csv("Data\\Code_03_Outcome\\ICD_pair_F_add(2).csv")

Pair_T = Pair_T[['source', 'target', 'edge', 'source_portein', 'target_portein','outcome1', 'source_portein_add','outcome1_add']]
Pair_F = Pair_F[['source', 'target', 'edge', 'source_portein', 'target_portein','outcome1', 'source_portein_add','outcome1_add']]
Pair_T.rename(columns = {'outcome1':'distance (ST)', 'outcome1_add': 'distance (SCT)'}, inplace = True)
Pair_F.rename(columns = {'outcome1':'distance (ST)', 'outcome1_add': 'distance (SCT)'}, inplace = True)

# Distance Reduction
Pair_T['distance (DIFF)'] = Pair_T['distance (ST)'] - Pair_T['distance (SCT)']
Pair_F['distance (DIFF)'] = Pair_F['distance (ST)'] - Pair_F['distance (SCT)']
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
### For each disease pair, compare it with negative disease pairs with the same source disease
def Compare_Single(treat,control,str1,str2):
    '''
    str1: the column for comparision
    str2: the parameter for one-tail test ( 'less': the values in control group is smaller than it in treatment group;
                                            'greater': the values in control group is karger than it in treatment group)
    compare each comorbidity patterns in the treatment group with comorbidity patterns which have the same pre-existing disease in control group 
    '''
    outcome = [] 
    control_G = control.groupby('source')
    for idx in treat.index:
        edge = treat.loc[idx, 'edge']
        item = treat.loc[idx, 'source']
        control_p = control_G.get_group(item)
        value = treat.loc[idx, str1]
        statistic,pvalue = wilcoxon(control_p[str1].tolist(),[value] * control_p.shape[0],alternative= str2)
        outcome.append([edge, pvalue,'W',value,control_p[str1].mean()])
    outcome = pd.DataFrame(outcome, columns = ['edge', 'pvalue', 'method', 'treat value', 'control value'])
    outcome.sort_values('pvalue', ascending = True, inplace  = True)
    outcome['source'] = outcome['edge'].apply(lambda x: eval(x)[0])
    outcome['target'] = outcome['edge'].apply(lambda x: eval(x)[1])
    return outcome

pvalues = Compare_Single(Pair_T,Pair_F,'distance (DIFF)','less')
pvalues[pvalues['pvalue'] < 0.05].shape # (54, 5)
len(pvalues[pvalues['pvalue'] < 0.05]['source'].unique()) # 42 
pvalues.to_csv('Data/Code_05_Outcome/diff_distance_test.csv', index = False) 

### Plot
### violin plot
def BoxPlot_Data(Treat, Control, str1):
    Treat = Treat[['source',str1]]
    Control = Control[['source',str1]]
    outcome = pd.DataFrame()    
    target = Treat['source'].unique().tolist()
    Treat_G = Treat.groupby('source')
    Control_G = Control.groupby('source')
    for disease in target:
        Treat_p = Treat_G.get_group(disease)
        Control_p = Control_G.get_group(disease)
        Treat_p = Treat_p[['source',str1]]
        Control_p = Control_p[['source',str1]]
        Treat_p['type'] = 'positive'
        Control_p['type'] = 'negative'
        res = pd.concat([Treat_p, Control_p], axis = 0)
        outcome = pd.concat([outcome,res], axis = 0)  
    Treat['type'] = 'positive'
    Control['type'] = 'negative'
    Treat['source'] = 'Total'
    Control['source'] = 'Total'
    outcome = pd.concat([Treat,Control,outcome], axis = 0).reset_index(drop = True)
    outcome = outcome.reset_index(drop = True)
    return outcome
Pair_T_significant = Pair_T[Pair_T['edge'].isin(pvalues[pvalues['pvalue'] < 0.05]['edge'])]
Pair_T_significant = Pair_T_significant[Pair_F.columns]
boxplot_data = BoxPlot_Data(Pair_T_significant, Pair_F, 'distance (DIFF)')  
boxplot_data.to_csv('Data/Code_05_Outcome/box_plot_data(significant).csv', index = False) # 62042

boxplot_data_total = BoxPlot_Data(Pair_T, Pair_F, 'distance (DIFF)') 
boxplot_data_total.to_csv('Data/Code_05_Outcome/box_plot_data(total).csv', index = False) # 79328
def Boxplot(boxplot_data, path):
    boxplot_data.sort_values(by = ['source', 'type'], ascending = True, inplace = True)
    rgb_values = [(129, 184, 223),(254, 129, 125)]
    # Convert RGB values to normalized values between 0 and 1
    colors = [(r / 255, g / 255, b / 255) for r, g, b in rgb_values]
    sns.set(font = 'Arial')
    sns.set_style('white')
    plt.figure(figsize = (20,10))
    sns.boxplot(data = boxplot_data, x = 'source', y = 'distance (DIFF)', hue = 'type', palette = colors)
    plt.xticks(rotation = 90, fontsize = 15)
    plt.yticks(fontsize = 15)
    # Add labels and title
    plt.xlabel('Disease', fontsize = 20)
    plt.ylabel('Distance Difference', fontsize = 20)
    # plt.title('Distance Difference in Positive and Negative Sample')
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    # plt.show()

Boxplot(boxplot_data,'Data/Code_05_Outcome/box_plot_data(significant).svg')
Boxplot(boxplot_data_total,'Data/Code_05_Outcome/box_plot_data(total).svg')
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
### [Part Three] Z_score Computation
def Z_Score(input):
    real = input[0]
    m = input[1]
    std = input[2]
    repeat = eval(input[3])
    p = np.mean([0 if i > real else 1 for i in repeat])
    z = (real - m)/std
    return z, p
files = os.listdir('Data/Code_04_Outcome')
files = [i for i in files if 'permutation' in i]
ref_dict = {}
for file in files:
    file_content = pd.read_csv(f'Data\Code_04_Outcome\{file}')
    if file_content.shape[0] != 1000:
        print(file)
    edge = file_content.loc[0, 'edge']
    file_content['distance (DIFF)'] = file_content['outcome1'] - file_content['outcome1_add']

    mean_sct = file_content['outcome1_add'].mean()
    std_sct = file_content['outcome1_add'].std()
    values_sct = file_content['outcome1_add'].tolist() 

    mean = file_content['distance (DIFF)'].mean()
    std = file_content['distance (DIFF)'].std()
    values = file_content['distance (DIFF)'].tolist() 

    ref_dict[edge] = {'mean distance (SCT)':mean_sct, 'std distance (SCT)':std_sct, 'values distance (SCT)': str(values_sct), 
                      'mean distance (DIFF)':mean, 'std distance (DIFF)':std, 'values distance (DIFF)': str(values)}

Ref_DF = pd.DataFrame.from_dict(ref_dict, orient = 'index').reset_index(drop = False)
Ref_DF.rename(columns = {'index':'edge'}, inplace = True)
Pair_T= pd.read_csv("Data\\Code_03_Outcome\\ICD_pair_T_add(2).csv")
Pair_T = pd.merge(Pair_T, Ref_DF, how = 'left', on = 'edge')

Pair_T[['z_score (distance (DIFF))','p_value (distance (DIFF))']] = Pair_T[['outcome1_add','mean distance (DIFF)','std distance (DIFF)','values distance (DIFF)']].apply(Z_Score,axis = 1,result_type='expand')
Pair_T['p_value (distribution) (distance (DIFF))'] = Pair_T['z_score (distance (DIFF))'].apply(lambda x: 1-stats.norm.cdf(x))
Pair_T['p_value (distance (DIFF))'] = 1 - Pair_T['p_value (distance (DIFF))']
Pair_T[Pair_T['p_value (distribution) (distance (DIFF))'] < 0.05] # 90
### 包括P-value 和 Z-score
Pair_T.to_csv("Data\\Code_05_Outcome\\z_score.csv",index = False)
len(Pair_T['source'].unique()) # 74
len(Pair_T['target'].unique()) # 37
len(set(Pair_T['source'].unique().tolist()).union(set(Pair_T['target'].unique().tolist()))) # 84

def AM(DF):
    items = sorted(list(set(DF['source'].unique().tolist() + DF['target'].unique().tolist())))
    items_dict = {j:i for i,j in enumerate(items)}
    adjacency_matrix = np.zeros((len(items), len(items)))
    for idx in DF.index:
        source = DF.loc[idx, 'source']
        target = DF.loc[idx, 'target']
        z_score = DF.loc[idx, 'z_score (distance (DIFF))']
        adjacency_matrix[items_dict[source], items_dict[target]] = z_score
        adjacency_matrix[items_dict[target], items_dict[source]] = z_score/4
    return adjacency_matrix,items

Z_socres_significant = Pair_T[Pair_T['p_value (distribution) (distance (DIFF))'] < 0.05][['source', 'target', 'edge','z_score (distance (DIFF))','p_value (distance (DIFF))','p_value (distribution) (distance (DIFF))']] 
Z_socres_total = Pair_T[Z_socres_significant.columns]
Z_scores_total_positive = Z_socres_total[Z_socres_total['z_score (distance (DIFF))'] > 0]
Z_scores_total_negative = Z_socres_total[Z_socres_total['z_score (distance (DIFF))'] <= 0]
Z_scores_total_negative['z_score (distance (DIFF))'] = Z_scores_total_negative['z_score (distance (DIFF))'].apply(lambda x: np.abs(x))

[Z_socres_significant.shape, Z_socres_total.shape] # [(90, 6), (145, 6)]
[Z_scores_total_positive.shape, Z_scores_total_negative.shape] # [(117, 6), (28, 6)]

Z_socres_significant.to_csv("Data\\Code_05_Outcome\\chord_significant.csv",index = False)
Z_socres_total.to_csv("Data\\Code_05_Outcome\\chord_total.csv",index = False)
Z_scores_total_positive.to_csv("Data\\Code_05_Outcome\\chord_total(positive).csv",index = False)
Z_scores_total_negative.to_csv("Data\\Code_05_Outcome\\chord_total(negative).csv",index = False)

net_node = pd.read_csv(r"Data\Code_01_Outcome\node_detail(2).csv")
group_ls = list(sorted(net_node['Detail_extend'].unique().tolist()))
color_ls = ['#FFB09B','#CBCD84','#FFBF74','#FFADC9','#FFBAF6','#A6D8FF','#EAD55E','#55F5FF','#59F3D1','#EDBBCE','#85E49F','#64CACC','#B2E779','#A4D7C2']
color_dict = {i:j for i,j in zip(group_ls,color_ls)}
net_node['color'] = [color_dict[i] for i in net_node['Detail_extend']]
node_color_dict = {i:j for i,j in zip(net_node['Conditions'].tolist(), net_node['color'].tolist())}

def Chord_Plot(data,color_dict,path):
    adjacency_matrix, labels = AM(data)
    labels_color = [color_dict[i] for i in labels]
    fig = ocd.Chord(adjacency_matrix, labels)
    fig.radius = 600
    fig.padding = 600
    fig.font_size = 20
    fig.gap_size = 0.03
    fig.font_family = "Arial"
    fig.bg_transparancy = 0.05
    fig.colormap = labels_color
    # fig.show()
    fig.save_svg(path)

Chord_Plot(Z_socres_total, node_color_dict, 'Data/Code_05_Outcome/chord_plot_data(total).svg')
Chord_Plot(Z_socres_significant, node_color_dict, 'Data/Code_05_Outcome/chord_plot_data(significant).svg')
Chord_Plot(Z_scores_total_positive, node_color_dict, 'Data/Code_05_Outcome/chord_plot_data(total)(positive).svg')
Chord_Plot(Z_scores_total_negative, node_color_dict, 'Data/Code_05_Outcome/chord_plot_data(total)(negative).svg')

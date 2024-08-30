'''
绘制论文中涉及到的图像
绘制 sankey plot
'''
'''
Sankey Diagrams: source condition ICD category ---> target condition ICD9 catrgory 
'''
import plotly.io as pio
import plotly.graph_objects as go
import plotly.express as pex
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns

net_node = pd.read_csv(r"Data\Code_01_Outcome\group_count(2).csv")
net_node.head()
# 设置标签
net_node['Source Group Label'] = ['S ' + i.split(':')[0] for i in net_node['source'].tolist()]
net_node['Target Group Label'] = ['T ' + i.split(':')[0] for i in net_node['target'].tolist()]

set_nodes =  sorted(list(set(net_node['Source Group Label'].tolist() + net_node['Target Group Label'].tolist())))
color_ls = ['#FFB09B','#CBCD84','#FFBF74','#FFADC9','#FFBAF6','#A6D8FF','#EAD55E','#55F5FF','#59F3D1','#EDBBCE','#85E49F','#64CACC','#B2E779','#A4D7C2']
alpha = 1

def hex_to_rgba(hex_code, alpha):
    # 去除 # 号并将十六进制代码转换为RGB元组
    hex_code = hex_code.lstrip('#')
    rgb = tuple(int(hex_code[i:i+2], 16) for i in (0, 2, 4))
    # 添加透明度分量并返回RGBA表示
    return rgb + (alpha,)

color_ls = [hex_to_rgba(item_.strip('#'), alpha) for item_ in color_ls]
color_ls = ['rgba'+str(item_) for item_ in color_ls]
nodes_ls = sorted(set([i.split(' ')[1] for i in set_nodes]))
len(nodes_ls)
### 为每一类疾病设置配色 
color_dict = {i:j for i,j in zip(nodes_ls,color_ls)}
color_df = pd.DataFrame.from_dict(color_dict,orient='index').reset_index(drop = False)
color_df.columns = ['Group Label','color']
color_df_1 = pd.DataFrame()
color_df_2 = pd.DataFrame()
color_df_1['Group Label'] = ['S ' + item_ for item_ in color_df['Group Label']]
color_df_2['Group Label'] = ['T ' + item_ for item_ in color_df['Group Label']]
color_df_1['color'] = color_df['color']
color_df_2['color'] = color_df['color']
color_df = pd.concat([color_df_1,color_df_2],axis=0).reset_index(drop  = True)
color_dict = {i:j for i,j in zip(color_df['Group Label'].tolist(),color_df['color'].tolist())}

################################################# 
# 图表美化
## 对于每一个source选择top3的关系 设置大的alpha
## 对于其他的 设置小的alpha
#################################################
net_node['edge color'] = [color_dict[item_] for item_ in net_node['Source Group Label']]
value_str = 'Diff_Corr'
net_node = net_node.sort_values(value_str,ascending=False).reset_index(drop = True)
net_node['Group Group Ref'] = net_node['Source Group Label'].str.cat(net_node['Target Group Label'],sep='_')
important_ls = net_node['Group Group Ref'].tolist()[:10]

def Judge_Important_Edge(input_,important_ls):
    gg_ = input_[0]
    color_ = input_[1] 
    if gg_ in important_ls:
        color_ = color_.replace(', 1',', 0.5')
    else:
        color_ = color_.replace(', 1',', 0.1')
    return color_

net_node['edge color'] = net_node[['Group Group Ref','edge color']].apply(Judge_Important_Edge,important_ls = important_ls,axis = 1)

len(color_dict.keys())
## All chart nodes
all_nodes = sorted(list(set(net_node['Source Group Label'].tolist() + net_node['Target Group Label'].tolist())))
len(all_nodes)
x_nodes = {i:j for i,j in zip(all_nodes,[0.1]*14+[0.9]*14)}
y_nodes = {i:j for i,j in zip(all_nodes,np.linspace(0.1,0.9,14).tolist()+ np.linspace(0.1,0.9,14).tolist())}
sorted(list(set(net_node['Source Group Label'].str.split(' ').str.get(1).tolist() + net_node['Target Group Label'].str.split(' ').str.get(1).tolist())))
## Indices of sources and destinations
source_indices = [all_nodes.index(item_) for item_ in net_node['Source Group Label'].tolist()] ## Retrieve source nodes indexes as per all nodes list.
target_indices = [all_nodes.index(item_) for item_ in net_node['Target Group Label'].tolist()] ## Retrieve destination nodes indexes as per all nodes list.


fig = go.Figure(data=[go.Sankey(
                        # Define nodes
    arrangement = "snap",
    node = dict(
        # pad = 12.5,
        thickness = 20, 
        label = all_nodes,
        color =  [color_dict[item_] for item_ in all_nodes],
        x = [x_nodes[item_] for item_ in all_nodes],
        y = [y_nodes[item_] for item_ in all_nodes],
    ),
    # Add links
    link = dict(
        source = source_indices,
        target = target_indices,
        value =  net_node[value_str].tolist(),
        color =  net_node['edge color']
    )
    )])

fig.update_layout(title={
        'x': 0.5,  # 设置标题的水平位置为中央
        'xanchor': 'center',  # 锚点设置为中央
        'y': 0.95}, # 设置标题的垂直位置，0为底部，1为顶部    
        font = dict(family = 'Arial', size = 20)

    )
fig.update_layout(height=800, width=1200) 
fig.show()
fig.write_image(r"Data\SankeyPlot\Diff_Sum.svg")



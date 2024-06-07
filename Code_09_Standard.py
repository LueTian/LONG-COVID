import pandas as pd 
import numpy  as np
import os 
path1 = r"C:\Users\31435\Desktop\Appendix\pdf"
def Stand(path_):
    for df_name in os.listdir(path_):
        df = pd.read_csv(os.path.join(path_,df_name))
        float_columns = df.select_dtypes(include=['float']).columns
        df[float_columns] = df[float_columns].round(3)
        df_name = df_name.split('.')[0]
        df_name = df_name + '.xlsx'
        df.to_excel(os.path.join(path_,df_name),index = False)
Stand(path1)

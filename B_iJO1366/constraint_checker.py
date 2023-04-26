import os
import sys

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 
 
from Ob_Met_Net_solmethods import MILP_sol_OP
from support_functions import brute_check_cons
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'C_Templates')))
from myDF import dataframe

import pandas as pd
from MN_ijo1366 import MN_ijo1366

mn = MN_ijo1366

df = dataframe

df_analysis = df[['C_Che','IC_Che','D_Che','Strain','K','Method','Tgt']][(df.IC_Che == 2000) & (df.Method == 'M') & (df.Strain == mn.Name) & (df.D_Che > 1e-2) & (df.Tgt >= 0.1) & (df.Tgt < 1)]

file_name = f"../Results/Cts_Check{mn.Name[:3]}.txt"

with open(file_name,"w+") as fb:
    fb.write(f"=========== {mn.Name} ===========")
    fb.write(f"\n")
    fb.write(df_analysis.to_markdown())
    fb.write(f"\n") 

if len(df_analysis) == 0:
    with open(file_name,"w+") as fb:

        fb.write(f'=========== No Analysys ===========')
        fb.write(f"\n")
else: 
    for i in range(len(df_analysis.head())):
        target = df_analysis.iloc[i].Tgt
        k = df_analysis.iloc[i].K
        mn.target = target
        print(mn.FVA[mn.chemical])
        c = MILP_sol_OP(network=mn,k=k,log=False)
        brute_check_cons(network=mn,solution=c,file=file_name)

import os
import sys

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 
 
from Ob_Met_Net_solmethods import MILP_sol_OP
from support_functions import brute_check_cons
import pandas as pd
from MN_ijr904 import MN_ijr904
from tqdm import tqdm
mn = MN_ijr904

df = pd.read_csv(f"../Results/Envelopes/MC_2.csv")

df_analysis = df[['C_Che','IC_Che','D_Che','Strain','K','Method','Tgt']][(df.IC_Che == 2000)&(df.Method == 'M') & (df.Strain == mn.Name) & (df.D_Che > 1e-2) & (df.Tgt >= 0.1) & (df.Tgt < 1)]


file_name = f"../Results/Cts_Check_{mn.Name[:3]}_3.txt"


values = []

for i in range(len(df_analysis)):
    tgt_df = float(df_analysis.iloc[i].Tgt)
    k_df = int(df_analysis.iloc[i].K)
    values.append((tgt_df,k_df))

with open(file_name,"w+") as fb:
    fb.write(f"=========== {mn.Name} ===========\n")
    fb.write(f"=========== FBA =================\n")
    fb.write(f"v[b] = {mn.FBA[mn.biomass]} v[c] = {mn.FBA[mn.chemical]}")
    fb.write(f"\n")
    fb.write(df_analysis.to_markdown())
    fb.write(f"\n") 

if len(df_analysis) == 0:
    with open(file_name,"a") as fb:
        fb.write(f'=========== No Analysys ===========')
        fb.write(f"\n")
else: 
    for pairs in tqdm(values,desc='Progress'):
        mn.target , ks = pairs
        c = MILP_sol_OP(network=mn,k=ks,log=False)
        brute_check_cons(network=mn,solution=c,file=file_name)


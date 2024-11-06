import pandas as pd

from itertools import product


# def mydataframes(strain:str=None,k:int=None):
#     df = pd.read_csv(f"FD_OP_k{k}_{strain}_w.csv")

#     new_r = {'Biomass':df.MB.tolist() + df.CB.tolist() + df.PB.tolist(),
#       'Chemical':df.MC.tolist() + df.CC.tolist() + df.PC.tolist(),
#       'Method': ['M' for _ in range(len(df.MB))] + ['C' for _ in range(len(df.CB))] + ['P' for _ in range(len(df.PB))],
#       'Time': df.Mt.tolist() + df.Ct.tolist() + df.Pt.tolist(),
#       'Tgt': df.tgt.tolist() + df.tgt.tolist() + df.tgt.tolist(),
#       'Strain':[strain for _ in range(len(df)*3)],
#       'K': [k for _ in range(len(df)*3)],
#       'Strategy':df.Mys.tolist() + df.Cys.tolist() + df.Pys.tolist()}
    
#     return pd.DataFrame.from_dict(new_r)


def mymasterframe(strain:str=None,k:int=None):
    df = pd.read_csv(f"FD_OP_k{k}_{strain}_w.csv")

    new_r = {'C_Bio':df.MB.tolist() + df.CB.tolist() + df.PB.tolist(),
      'C_Chel':df.MC.tolist() + df.CC.tolist() + df.PC.tolist(),
      'Method': ['M' for _ in range(len(df.MB))] + ['O' for _ in range(len(df.CB))] + ['P' for _ in range(len(df.PB))],
      'Time': df.Mt.tolist() + df.Ct.tolist() + df.Pt.tolist(),
      'Tgt': df.tgt.tolist() + df.tgt.tolist() + df.tgt.tolist(),
      'Strain':[strain for _ in range(len(df)*3)],
      'K': [k for _ in range(len(df)*3)],
      'Strategy':df.Mys.tolist() + df.Cys.tolist() + df.Pys.tolist(),
      'E_Bio':df.Bio.tolist() + df.Bio.tolist() + df.Bio.tolist(),
      'E_Che':df.Che.tolist() + df.Che.tolist() + df.Che.tolist(),
      'IC_Bio':df.ICB_mb.tolist() + df.ICB_cb.tolist() + df.ICB_pb.tolist(),
      'IC_Che':df.ICB_mc.tolist() + df.ICB_cc.tolist() + df.ICB_pc.tolist() }
    
    return pd.DataFrame.from_dict(new_r)


strains = ['iAF1260','iJO1366','iJR904']
ks = [1,2,3]
dfs = []
for k,strain in product(ks,strains):
    print(f"Compiling k= {k} & strain {strain}")
    di = mymasterframe(strain=strain,k=k) 
    dfs.append(di)
    print(f"appended...")
    

f = pd.concat(dfs,sort=False)

file_excel = f"Master_complete_2.xlsx"
file_csv = f"FD_OP_complete_2.csv" 

f.to_csv(file_csv)
with pd.ExcelWriter(file_excel, mode='w') as writer:
    f.to_excel(writer,sheet_name='Summary')




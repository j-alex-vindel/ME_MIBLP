import pandas as pd
from itertools import product


def compiledfs(K:int=None,strain:str=None):

    df = pd.read_csv(f"C:/Users/alexa/OneDrive - University of Strathclyde/ME_PHD/ME_MIBLP/Results/Envelopes/Revised/Revised_{strain[:3]}_k{K}.csv")

    new_r = {'C_Bio':df.MB.tolist() + df.OB.tolist() + df.PB.tolist(),
        'C_Che':df.MC.tolist() + df.OC.tolist() + df.PC.tolist(),
        'Method': ['M' for _ in range(len(df.MB))] + ['O' for _ in range(len(df.OB))] + ['P' for _ in range(len(df.PB))],
        'Time': df.Mt.tolist() + df.Ot.tolist() + df.Pt.tolist(),
        'Tgt': df.tgt.tolist() + df.tgt.tolist() + df.tgt.tolist(),
        'Strain':[strain for _ in range(len(df)*3)],
        'K': [K for _ in range(len(df)*3)],
        'Ys':df.Mys.tolist() + df.Oys.tolist() + df.Pys.tolist(),
        'E_Bio':df.Bio.tolist() + df.Bio.tolist() + df.Bio.tolist(),
        'E_Che':df.Che.tolist() + df.Che.tolist() + df.Che.tolist(),
        'IB_Bio':df.ICB_mb.tolist() + df.ICB_ob.tolist() + df.ICB_pb.tolist(),
        'IC_Bio':df.ICC_mb.tolist() + df.ICC_ob.tolist() + df.ICC_pb.tolist(),
        'IB_Che':df.ICB_mc.tolist() + df.ICB_oc.tolist() + df.ICB_pc.tolist(),
        'IC_Che': df.ICC_mc.tolist() + df.ICC_oc.tolist() + df.ICC_pc.tolist(),
        'GKO': df.MGKO.tolist() + df.OGKO.tolist() + df.PGKO.tolist()}

    return pd.DataFrame.from_dict(new_r)

def unifieddf(dfs,file:bool=None):
    f = pd.concat(dfs,sort=False)
    txt = 'Rev'
    for i in f.Strain.unique():
        txt += f"_{i[:3]}"
    f_xlx = f"{txt}_m.xlsx"
    f_csv = f"{txt}_m.csv" 
    f.to_csv(f_csv)
    
    if file:
        with pd.ExcelWriter(f_xlx, mode='w') as writer:
            f.to_excel(writer,sheet_name='Summary')
    return f


strains = ['iAF1260','iJO1366','iJR904']

ks = [1,2,3]
dfs = []

for K,strain in product(ks,strains):
    print(strain,'->',K)
    temp_df = compiledfs(K=K,strain=strain)
    dfs.append(temp_df)

mdf = unifieddf(dfs=dfs,file=False)

print(mdf.Strain.unique())
    



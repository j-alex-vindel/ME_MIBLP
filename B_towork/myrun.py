import os
import sys
import time
import pandas as pd
import datetime
import matplotlib.pyplot as plt
from itertools import product

def main(arg1=None,arg2=None,name:str=None,sleep:int=None,file:str=None):
    
    # os.system(f"python ../A_modules/{sys.argv[1]}.py")
    print(f">> Folder {name}")
    print(f">> Running -> ../B_{name}/{file} {arg1} {arg2}.py")
    if arg2 == None:
        os.system(f"python ../B_{name}/{file}.py {arg1}")
    else:   
        os.system(f"python ../B_{name}/{file}.py {arg1} {arg2}")

    print(f'>> Resting for {sleep}')
    time.sleep(sleep)
    print(f">> Ready to start over")
    
    print('Run Completed')


if __name__ == '__main__':
    strains = ['ijr']
    checks = ['A','B','C']
    sleep = 10
    file = f"full_run"
    ks = [1,2,3]
    
    strains = {"iJR904":[366,376],'iJO1366':[340,422,91,87],'iAF1260':[845,772,707,764]}
    for k in ks:
        for strain,indeces in strains.items():
            for index in indeces:
                main(arg1=k,arg2=index,name=strain,file=file,sleep=sleep)
        

# =========================================================================================
# =========================================================================================
# =========================================================================================


# def latex_tables(bacteria:str=None):
#     fsol1 = f"../Results/EQ/Methods_K1_{bacteria}.csv"
#     fsol2 = f"../Results/EQ/Methods_K2_{bacteria}.csv"
#     fsol3 = f"../Results/EQ/Methods_K3_{bacteria}.csv"

#     fgeq1 = f"../Results/EQ/IC_K1_{bacteria}.csv"
#     fgeq2 = f"../Results/EQ/IC_K2_{bacteria}.csv"
#     fgeq3= f"../Results/EQ/IC_K3_{bacteria}.csv"
    
    
#     s1,s2,s3 = pd.read_csv(fsol1),pd.read_csv(fsol2),pd.read_csv(fsol3)

#     st = pd.concat([s1,s2,s3],sort=False).round(decimals=4)

#     i1,i2,i3 = pd.read_csv(fgeq1),pd.read_csv(fgeq2),pd.read_csv(fgeq3)

#     it = pd.concat([i1,i2,i3],sort=False)
#     table = pd.pivot_table(it,index=['From','FBA_Objct','K'],columns=['Criteria'],values=['V_biomass','V_chemical'])
    
#     # output
#     file_dir = f"../Results/Latex_tables/EQ/{bacteria}.txt"

#     with open(file_dir,'a') as txt:
#         txt.write(f'>> Table Solution {bacteria} <<\n')
#         txt.write(st.style.to_latex())
#         txt.write(f'\n')
#         txt.write(f'\n')
#         txt.write(f'>> Inner Checks EQ {bacteria} <<\n')
#         txt.write(table.style.to_latex())
    
#     print('Tables Writen!!!')

# def latex_tables_pes(bacteria:str=None):
#     fsol1 = f"../Results/Methods_OP2_{bacteria}.csv"
#     # fsol2 = f"../Results/Pessimistic/Methods_Kn_{bacteria}.csv"
#     # fsol3 = f"../Results/Pessimistic/Methods_Kn_{bacteria}.csv"

#     st =pd.read_csv(fsol1)
#     # s1,s2,s3 = pd.read_csv(fsol1),pd.read_csv(fsol2),pd.read_csv(fsol3)
    
#     # st = pd.concat([s1,s2,s3],sort=False).round(decimals=4)
    
#     # output
#     file_dir = f"../Results/Latex_tables/OP2_{bacteria}.txt"

#     with open(file_dir,'a') as txt:
#         txt.write(f'>> Table Solution {bacteria} <<\n')
#         txt.write(st.style.to_latex())
#         txt.write(f'\n')
#         txt.write(f'\n')   
#     print('Tables Writen!!!')

# def pessimistic_run(bacteria:str=None,target:float=None,k:int=None,log:bool=False,version:int=None):

#     dw = {0:'mo',1:'tu',2:'we',3:'th',4:'fr',5:'sa',5:'su'}

#     print(f"\n>> Running-> ../B_towork/pessimistic_{bacteria}.py")
#     print(f">> Target -> {target} Biomass")
#     print(f">> K -> {k}")
#     print(f">> Log file -> {log}")
#     now = datetime.datetime.now()
#     txt_file = f"../Results/LOG_P/log_{bacteria}_{k}{int(100*target)}_{dw[now.weekday()]}{version}.txt"

#     if log:
#         os.system(f"python ../B_towork/pessimistic_{bacteria}.py {target} {k} > {txt_file}")
#     else:
#         os.system(f"python ../B_towork/pessimistic_{bacteria}.py {target} {k}")
    
#     print(f"File Executed!! \n")

# def tables_md_tex(entry=None):
#     df,name,ko = entry
#     df[f'biomass pct'] = [f"{int(i*100)}%" for i in df['tgt']]
#     df_b = df[['biomass pct','Bio','Mys','MB','ICB_mb','ICC_mb','Oys','OB','ICB_ob','ICC_ob','Pys','PB','ICB_pb','ICC_pb']]
#     df_c = df[['biomass pct','Che','Mys','MC','ICB_mc','ICC_mc','Oys','OC','ICB_oc','ICC_oc','Pys','PC','ICB_pc','ICC_pc']]

#     file_b = f"../Results/Latex_tables/Full Detail/Revised/Biomass_{name[:3]}_{ko}.txt"
#     isfile_b = os.path.exists(file_b)
#     if not isfile_b:
#         with open(file_b,'w+') as fb:
#             fb.write(f"======= {name} ============")
#             fb.write(f"======= Biomass ===========")
#             fb.write(" ")
#             fb.write(f"\n")
#             fb.write(df_b.to_markdown())
#             fb.write(" ")
#             fb.write(f"\n")
#             fb.write(df_b.to_latex())
    
#     file_c = f"../Results/Latex_tables/Full Detail/Revised/Chemical_{name[:3]}_{ko}.txt"
#     isfile_c = os.path.exists(file_c)
#     if not isfile_c:
#         with open(file_c,"w+") as fc:
#             fc.write(f"======= {name} ============")
#             fc.write(f"======= Chemical ===========")
#             fc.write(" ")
#             fc.write(f"\n")
#             fc.write(df_c.to_markdown())
#             fc.write(" ")
#             fc.write(f"\n")
#             fc.write(df_c.to_latex())
            
#     fx = f"../Results/Excel_files/{name[:3]}_summary_{ko}.xlsx"
    
#     with pd.ExcelWriter(fx,mode='w') as writer:

#         df.to_excel(writer,sheet_name="Summary")
#         df_b.to_excel(writer,sheet_name='Biomass')
#         df_c.to_excel(writer,sheet_name='Chemical')

# def txttables(entry,file):
#     df,name,ko = entry
#     dfp = df[["biomass pct",'PB','PC','ICB_pb','ICB_pc','ICC_pb','ICC_pc','Pys','Pt']]
#     dfc = df[["biomass pct",'CB','CC','ICB_cb','ICB_cc','ICC_cb','ICC_cc','Cys','Ct']]
#     dfm = df[["biomass pct",'MB','MC','ICB_mb','ICB_mc','ICC_mb','ICC_mc','Mys','Mt']]
#     with open(file,mode='a') as f:
#         f.write(f"\n>>> Bacteria {name}\n")
#         f.write(f">>> Callbacks - Pessimistic")
#         f.write(f"  \n")
#         f.write(dfp.to_markdown())
#         f.write(f" \n")

#         f.write(f">>> Bacteria {name}\n")
#         f.write(f">>> Callbacks - Optimistic")
#         f.write(f"  \n")
#         f.write(dfc.to_markdown())
#         f.write(f" \n")

#         f.write(f">>> Bacteria {name}\n")
#         f.write(f">>> MILP - Optimistic")
#         f.write(f"  \n")
#         f.write(dfm.to_markdown())
#         f.write(f" \n")
#         fx = f"../Results/Excel_files/{name}_summary_{ko}.xlsx"
#         with pd.ExcelWriter(fx, mode='a', engine="openpyxl") as writer:
#             dfp.to_excel(writer,sheet_name = 'Pessimistic-CB')
#             dfc.to_excel(writer,sheet_name = 'Optimistic-CB')
#             dfm.to_excel(writer,sheet_name = 'MILP')

# def normgraphs(entry):
#     df,name,ko = entry
#     nb = [i/max(df.Bio) for i in df.Bio]
#     nc = [i/max(df.Che) for i in df.Che]

#     nmb = [i/max(df.MB) for i in df.MB]
#     nmc = [i/max(df.MC) for i in df.MC]

#     ncb = [i/max(df.CB) for i in df.CB]
#     ncc = [i/max(df.CC) for i in df.CC]
    
#     npb = [i/max(df.PB) for i in df.PB]
#     npc = [i/max(df.PC) for i in df.PC]

#     fig,ax = plt.subplots()
#     plt.plot(nb,nc,linestyle='dotted',c='grey')

#     ax.set(xlim=(0,1.01),ylim=(0,1.01))
#     ax.set_title(f"Strain {name} ")
#     plt.xlabel(r""" % Biomass production $mmol/g(Dw)h$""")
#     plt.ylabel(r""" % Chemical production $mmol/g(Dw)h$""")
#     plt.scatter(nmb,nmc,marker='x',label="MILP",s=65,c='black')
#     plt.scatter(ncb,ncc,marker='d',label="CB_O",s=60,c='teal')
#     plt.scatter(npb,npc,marker='o',label="CB_P",s=60,c='red')
#     ax.legend()
#     plt.savefig(f"../Results/Graphs/NN_FE_{name}_k{ko}.png")
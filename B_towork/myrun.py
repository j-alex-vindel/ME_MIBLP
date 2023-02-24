import os
import sys
import time
import pandas as pd

def main(ks:list=None,name:str=None,sleep:int=None):
    for k in ks:
        # os.system(f"python ../A_modules/{sys.argv[1]}.py")
        print(f">> Folder {name}")
        print(f">> Running -> ../{name}/run_experiment_k{k}")

        os.system(f"python ../B_{name}/run_experiment_k{k}.py")
        print(f'>> Resting for {sleep}')
        time.sleep(sleep)
        print(f">> Ready to start over")
    
    print('Run Completed')


def latex_tables(bacteria:str=None):
    fsol1 = f"../Results/EQ/Methods_K1_{bacteria}.csv"
    fsol2 = f"../Results/EQ/Methods_K2_{bacteria}.csv"
    fsol3 = f"../Results/EQ/Methods_K3_{bacteria}.csv"

    fgeq1 = f"../Results/EQ/IC_K1_{bacteria}.csv"
    fgeq2 = f"../Results/EQ/IC_K2_{bacteria}.csv"
    fgeq3= f"../Results/EQ/IC_K3_{bacteria}.csv"
    
    
    s1,s2,s3 = pd.read_csv(fsol1),pd.read_csv(fsol2),pd.read_csv(fsol3)

    st = pd.concat([s1,s2,s3],sort=False).round(decimals=4)

    i1,i2,i3 = pd.read_csv(fgeq1),pd.read_csv(fgeq2),pd.read_csv(fgeq3)

    it = pd.concat([i1,i2,i3],sort=False)
    table = pd.pivot_table(it,index=['From','FBA_Objct','K'],columns=['Criteria'],values=['V_biomass','V_chemical'])
    
    # output
    file_dir = f"../Results/Latex_tables/EQ/{bacteria}.txt"

    with open(file_dir,'a') as txt:
        txt.write(f'>> Table Solution {bacteria} <<\n')
        txt.write(st.style.to_latex())
        txt.write(f'\n')
        txt.write(f'\n')
        txt.write(f'>> Inner Checks EQ {bacteria} <<\n')
        txt.write(table.style.to_latex())
    
    print('Tables Writen!!!')

def latex_tables_pes(bacteria:str=None):
    fsol1 = f"../Results/Methods_OP2_{bacteria}.csv"
    # fsol2 = f"../Results/Pessimistic/Methods_Kn_{bacteria}.csv"
    # fsol3 = f"../Results/Pessimistic/Methods_Kn_{bacteria}.csv"

    st =pd.read_csv(fsol1)
    # s1,s2,s3 = pd.read_csv(fsol1),pd.read_csv(fsol2),pd.read_csv(fsol3)
    
    # st = pd.concat([s1,s2,s3],sort=False).round(decimals=4)
    
    # output
    file_dir = f"../Results/Latex_tables/OP2_{bacteria}.txt"

    with open(file_dir,'a') as txt:
        txt.write(f'>> Table Solution {bacteria} <<\n')
        txt.write(st.style.to_latex())
        txt.write(f'\n')
        txt.write(f'\n')   
    print('Tables Writen!!!')





if __name__ == '__main__':
    ks = ['n']
    names = ['iJR904']
    sleep = 5
    
    for name in names:
        main(ks=ks,name=name,sleep=sleep)
        latex_tables(bacteria=name)
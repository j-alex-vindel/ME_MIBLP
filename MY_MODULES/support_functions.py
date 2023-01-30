import gurobipy as gp
from gurobipy import GRB
from typing import List,NewType
import copy
import pandas as pd


FBA = NewType('FBA_vector',List[float])

def set_constructor(list:List[str])->List[int]:
    if list != None:
        return [i for i in range(len(list))]

def wildtype_FBA(obj)->FBA:
    LB_wt = copy.deepcopy(obj.LB)
    UB_wt = copy.deepcopy(obj.UB)

    wt = gp.Model()
    v_wt = wt.addVars(obj.M, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='v_wt')
    
    wt.setObjective(1*v_wt[obj.biomass],GRB.MAXIMIZE)
    wt.addMConstr(obj.S,v_wt,'=',obj.b,name='Stoi')
    wt.addConstrs((LB_wt[j] <= v_wt[j] for j in obj.M), name='LBwt')
    wt.addConstrs((UB_wt[j] >= v_wt[j] for j in obj.M), name='UBwt')
    
    wt.Params.OptimalityTol = obj.infeas
    wt.Params.IntFeasTol = obj.infeas
    wt.Params.FeasibilityTol = obj.infeas
    wt.Params.OutputFlag = 0
    wt.optimize()
    if wt.status == GRB.OPTIMAL:
        wt_vs =  [wt.getVarByName('v_wt[%s]'%a).x for a in obj.M] 

    return wt_vs

def maketables(name:str=None):
    if name is None:
        return f"Enter a name iJO1366 - iJR904 - iAF1260"

    r1 = pd.read_csv(f'Results_NOP_k1_{name}.csv',index_col=False)
    r2 = pd.read_csv(f'Results_NOP_k2_{name}.csv',index_col=False)
    r3 = pd.read_csv(f'Results_NOP_k3_{name}.csv',index_col=False)

    r1 = r1[['Method','Biomass','Chemical','Time','Strategy']].round(decimals=5)
    r2 = r2[['Method','Biomass','Chemical','Time','Strategy']].round(decimals=5)
    r3 = r3[['Method','Biomass','Chemical','Time','Strategy']].round(decimals=5)

    ri1 = pd.read_csv(f'Results_NOP_IC_k1_{name}.csv')
    ri2 = pd.read_csv(f'Results_NOP_IC_k2_{name}.csv')
    ri3 = pd.read_csv(f'Results_NOP_IC_k3_{name}.csv')
    table1 = pd.pivot_table(ri1.round(decimals=5),index=['From','FBA_Objct'],columns=['Criteria'],values=['V_biomass','V_chemical'])
    table2 = pd.pivot_table(ri2.round(decimals=5),index=['From','FBA_Objct'],columns=['Criteria'],values=['V_biomass','V_chemical'])
    table3 = pd.pivot_table(ri3.round(decimals=5),index=['From','FBA_Objct'],columns=['Criteria'],values=['V_biomass','V_chemical'])
    
    
    print(f"======={name} K=1 ===========")
    print(' ')
    print(r1.to_markdown())
    print(' ')
    print('>> Inner Check')
    print(table1.to_string())
    print('===============================')
    print(' ')
    print(f"======={name} K=2 ===========")
    print(' ')
    print(r2.to_markdown())
    print(' ')
    print('>> Inner Check')
    print(table2.to_string())
    print('===============================')
    print(' ')
    print(f"======={name} K=3 ===========")
    print(' ')
    print(r3.to_markdown())
    print(' ')
    print('>> Inner Check')
    print(table3.to_string())
    print('===============================')
    
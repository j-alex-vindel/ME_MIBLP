import gurobipy as gp
from gurobipy import GRB
from typing import List,NewType,Type
import copy
import pandas as pd
from dict2xml import dict2xml
from collections import namedtuple
import matplotlib.pyplot as plt
import seaborn as sns
from bs4 import BeautifulSoup
import os

FBA = NewType('FBA_vector',List[float])
Result = namedtuple('Result',['MetNet','Strategy','Ys','Vs','Time','Soltype','Method'])
Result_cb = namedtuple('Result_cb',['MetNet','Strategy','Ys','Vs','Vij','Time','Soltype','Method'])
R_xml = namedtuple('Result_xml',['MetNet','Strategy','Ys','Vs','Time','Method','K'])
RM = Type[Result]
RC = Type[Result_cb]


def result_parser(file):
    k = int(file[-5])
    with open(file,'r') as f:
        data = f.read()

    soup = BeautifulSoup(data,"xml")

    vs = [float(flows.text) for flows in soup.find_all('Vs')]
    ys = [float(Bi.text) for Bi in soup.find_all('Ys')]
    time = float(soup.find('Time').text)
    
    if k ==1:
        strategy = [soup.find('Strategy')]
    else:
        strategy = [strat.text for strat in soup.find_all("Strategy")]

    method = soup.find('Method').text
    microbe = soup.find('MetNet').text

    return R_xml(microbe,strategy,ys,vs,time,method,k)


def set_constructor(list:List[str])->List[int]:
    if list != None:
        return [i for i in range(len(list))]

def wildtype_FBA(obj, wildtype:bool=True,mutant:bool=False)->FBA:
    LB_wt = copy.deepcopy(obj.LB)
    UB_wt = copy.deepcopy(obj.UB)

    if wildtype and not mutant:
        objective = obj.biomass
        FVA = False
    elif mutant and not wildtype:
        objective = obj.chemical
        FVA = True
    else:
        raise Exception("Both values of wildtype and mutant can't be both True or False at the same Time, pick one or the other")
        

    wt = gp.Model()
    v_wt = wt.addVars(obj.M, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='v_wt')
    
    wt.setObjective(1*v_wt[objective],GRB.MAXIMIZE)
    wt.addMConstr(obj.S,v_wt,'=',obj.b,name='Stoi')
    wt.addConstrs((LB_wt[j] <= v_wt[j] for j in obj.M), name='LBwt')
    wt.addConstrs((UB_wt[j] >= v_wt[j] for j in obj.M), name='UBwt')
    if FVA:
        wt.addConstr((v_wt[obj.biomass] >= obj.minprod), name='minprod')

    wt.Params.OptimalityTol = obj.infeas
    wt.Params.IntFeasTol = obj.infeas
    wt.Params.FeasibilityTol = obj.infeas
    wt.Params.OutputFlag = 0
    wt.Params.LogToConsole = 0
    wt.optimize()
    if wt.status == GRB.OPTIMAL:
        wt_vs =  [wt.getVarByName('v_wt[%s]'%a).x for a in obj.M] 
    else:
        wt_vs = [-200 for _ in obj.M]

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
    

def save_results(c:RC=None,m:RM=None):
    if c is None and m is None:
        return f"Nothing to save"
    elif c != None and m is None:
        file = f"../Results/XML/{c.Method}_{c.MetNet}_k{len(c.Strategy)}.xml"
        xml = dict2xml(dict(c._asdict()),wrap='root',indent="   ")

    elif m != None and c is None:
        file = f"../Results/XML/{m.Method}_{m.MetNet}_k{len(m.Strategy)}.xml"
        xml = dict2xml(dict(m._asdict()),wrap='root',indent="   ")
    else:
        return f"cannot handle two sol at the same time"

    is_file = os.path.exists(file)

    if not is_file:
        with open(file,'w+') as f:
            f.write(xml)

    return f" Solution saved!"

def do_them_graphs(obj,df,aproach):
    timeplot = f"../Results/Graphs/{aproach}_Time_bar_{obj.Name}.png"
    bioplot = f"../Results/Graphs/{aproach}_Biomass_{obj.Name}.png"
    chemplot = f"../Results/Graphs/{aproach}_Chemical_{obj.Name}.png"
    metplot = f"../Results/Graphs/{aproach}_Time_bymethod_{obj.Name}.png"
    timeline = f"../Results/Graphs/{aproach}_Time_line_{obj.Name}.png"


    is_time_linegraph = os.path.exists(timeline)

    is_time_bargraph = os.path.exists(timeplot)

    is_bioplot = os.path.exists(bioplot)

    is_chemplot = os.path.exists(chemplot)

    is_metplot = os.path.exists(metplot)

    if not is_time_linegraph:

        sns.catplot(x='Method',y='Time',hue='K',kind='point',data=df,palette='flare',linestyles=['dashed','dotted','dashdot'],markers=['^','o','s'])
        plt.savefig(timeline)

    if not is_time_bargraph:
        sns.catplot(data=df,x='Method',y='Time',hue='K',kind='bar',palette='flare')
        plt.savefig(timeplot)

    if not is_bioplot:
        bp = sns.catplot(data=df, x="Biomass", y="Method", hue='K',kind='swarm',palette='flare')
        plt.savefig(bioplot)

    if not is_chemplot:
        cp = sns.catplot(data=df, x="Chemical", y="Method", hue='K',kind='swarm',palette='flare')
        plt.savefig(chemplot)

    if not is_metplot:  
            sns.catplot(data=df,x='K',y='Time',hue='Method',kind='bar',palette='flare')
            plt.savefig(metplot)

    return f" Graphs completed!"
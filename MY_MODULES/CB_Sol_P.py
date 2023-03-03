import gurobipy as gp
from gurobipy import GRB 
from typing import List, Type, NewType
from collections import namedtuple
import copy
from Ob_Met_Net import Metabolic_Network
from itertools import combinations
import math 

M_Network = Type[Metabolic_Network]
Ks = NewType('K Strategies',int)
Vector = List[int]
Inner_obj = NewType('Inner Objective',str)
Result_cb = namedtuple('Result_cb',['MetNet','Strategy','Ys','Vs','Vij','Time','Soltype','Method'])


def CB_P(network:M_Network=None, k:Ks=None,log:bool=None,speed:bool=False,threads:bool=False) -> Result_cb:
    '''
        cb = CB_solve_2_NOP(network=network,k=k,log=True,speed=False,threads=False)
           
        cb.MetNet    = Metabolic Network's name
        cb.Stragtegy = List of rxn to knockout
        cb.Ys        = Binary solution as a vector
        cb.Vs        = Optimal bilevel flows
        cb.Vij       = Flows in the inner problem
        cb.Time      = Solving time 
        cb.Soltype   = Type of solution [optimal, timelimit , infeasible]
        cb.Method    = Solving Method - set to CBP (Callbacks-Pessimistic)

        '''
    print(f'**** Solving Callbacks k={k} ****')
    print(f'# Variables (reactions in the network): {len(network.M)}')
    print('Current Infeasibility:',network.infeas,sep=' -> ')
    print('KO set: ',len(network.KO), ' reactions')
    print(f"MN: {network.Name}")
    print(f"-- Pessimistic Approach --")

    lb = copy.deepcopy(network.LB)
    ub = copy.deepcopy(network.UB)
    minprod = copy.deepcopy(network.minprod)
    lb[network.biomass] = minprod

    def inner(imodel, yoj:Vector):
        global vij

        imodel.setAttr('LB',imodel.getVars(),[lb[j]*yoj[j] for j in network.M])
        imodel.setAttr('UB',imodel.getVars(),[ub[j]*yoj[j] for j in network.M])
        imodel.Params.OptimalityTol = network.infeas
        imodel.Params.IntFeasTol = network.infeas
        imodel.Params.FeasibilityTol = network.infeas
        imodel.Params.Presolve = 0
        imodel.optimize()
        status = imodel.status
        if status == GRB.OPTIMAL:
            vij = [imodel.getVarByName('fv[%s]'%a).x for a in network.M] # rounded inner values
        elif status in (GRB.INFEASIBLE, GRB.UNBOUNDED, GRB.INF_OR_UNBD):
            vij = [2000 if i == network.biomass else yoj[i] for i in network.M]
        return vij,status

    def lazycall(model,where):
        if where == GRB.Callback.MIPSOL:
            model._voj = model.cbGetSolution(model._vars)
            model._yoj = model.cbGetSolution(model._varsy)
            knockset =  [i for i,y in enumerate(model._yoj) if model._yoj[i] < 1e-6]

            if len(knockset) != k:
                return
            cur_obj = round(model.cbGet(GRB.Callback.MIPSOL_OBJBST),6)
            cur_bd = round(model.cbGet(GRB.Callback.MIPSOL_OBJBND),6)
            print(f"Vbio: {model._voj[network.biomass]}")
            print(f"Vche: {model._voj[network.chemical]}")
            model._vi, inner_status = inner(model._inner, model._yoj)

# ============================ Checking Inner Optimality Status ===================================
            # print('MIPSOL')
            # print(f"Objective = {cur_obj}")
            # print(f"Best Bound = {cur_bd}")
            # print(f"Curnt Pbnd = {model._pbnd}")
            print(f"inner value{model._vi[network.chemical]}")
            if inner_status != GRB.OPTIMAL:
                # print(f"optimality cuts inner not optinal MIPNODE")
                model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
            else:
                vi_biom_val = round(model._vi[network.biomass],6)
                vi_chem_val = round(model._vi[network.chemical],6)
                knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.KO)
                ki = (i for i in combinations(knockset_inner,k))

                # bestknownchem = cur_obj
                
                if model._pbnd - model._vi[network.chemical] >= -1e-6: # try with vij instead of model._vij to access the inner values, changed the value to -1e-6
                    # print(f'optimality cuts pbnd - vi[chemical')
                    model.cbLazy(sum(model._varsy[j] for j in knockset) >= 1)
                
                elif model._vi[network.chemical] - model._voj[network.chemical] >= 1e-6:
                    model.cbLazy(vi_chem_val <= model._vars[network.chemical] + network.BM(sum(model._varsy[i] for i in knockset)))
                    model._pbnd = model._vi[network.chemical] 

                elif abs((model._vi[network.biomass] - model._voj[network.biomass])) > 1e-6:
                    # print(f"Big M cut")
                    for comb in ki:
                        model.cbLazy(vi_biom_val <= model._vars[network.biomass] +
                            (math.ceil(vij[network.biomass]*10)/10) *(sum(model._varsy[f] for f in comb)))

                else:
                    # print(f"pbnd = cur_obj")
                    model._pbnd = cur_obj
# =============================================================================================================================                        

        elif where == GRB.Callback.MIPNODE:
            status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
            if status == GRB.OPTIMAL:
                # print('MIPNODE')
                mipobj = model.cbGet(GRB.Callback.MIPNODE_OBJBST)
                mipbnd = model.cbGet(GRB.Callback.MIPNODE_OBJBND)
                # print(f"Obj= {mipobj}")
                # print(f"Bnd= {mipbnd}")
                # print(f"Curnt set solution v: {model._sv[network.chemical]:.6f}")
                # print(f"Curnt set solution len y: {len(model._sy)}")
                # print(f"Curnt PBND: {model._pbnd}")
                model._cbcnt += 1
                if model._cbcnt % 10 == 0:
                    return
                # print(f"CB counts : {model._cbcnt}")
                model._voj = model.cbGetNodeRel(model._vars)
                model._ryo = model.cbGetNodeRel(model._varsy) #  
                for i,y in enumerate(model._ryo):
                    if model._ryo[y] >= 0.8:
                        model._ryo[y] = 1.0
                    elif model._ryo[y] <= 0.2:
                        model._ryo[y] = 0.0
                    else:
                        model._ryo[y] = 1.0
                knock = [i for i,y in enumerate(model._ryo) if model._ryo[i] < 1e-6]
                if sum(model._ryo.values()) != len(model._ryo)-k:
                    return
                else:
                    model._vi, inner_status = inner(model._inner,model._ryo)
                    # print(f"Optimality Code - {inner_status}")

                    if inner_status != GRB.OPTIMAL:
                        # print('Optimality cuts - inner not optimal')
                        model.cbCut(sum([model._varsy[f] for f in knock]) >= 1)
                    
                    elif (model._vi[network.chemical] > model._pbnd) and (model._vi[network.chemical] > model._sv[network.chemical]): # added condition to set solution always better
                        model._sv = model._vi
                        model._sy = model._ryo

                        # print(f"Set Solution ")
                        # print(f"biomas {model._vi[network.biomass]:.6f}")
                        # print(f"chemical {model._vi[network.chemical]:.6f}")
                        model.cbSetSolution(model._vars, model._sv)
                        model.cbSetSolution(model._varsy, model._sy)
                    
                    else:
                        # print(f"Optimality cuts - vi not larger or equal than pnbd")
                        model.cbCut(sum([model._varsy[f] for f in knock]) >= 1)
    
    m = gp.Model()
    cbv = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='cbv')
    cby = m.addVars(network.M,vtype=GRB.BINARY,name='cby')

    m.setObjective(1*cbv[network.chemical],GRB.MAXIMIZE)
    
    m.addMConstr(network.S,cbv,'=',network.b,name='Stoi')
    # m.addConstrs((gp.quicksum(network.S[i,j]*cbv[j] for j in network.M) == 0 for i in network.N),name='Stoichiometry')

    m.addConstr(cbv[network.biomass] >= minprod, name='target')

    m.addConstrs((lb[j]*cby[j] <= cbv[j] for j in network.M),name='LB')
    m.addConstrs((cbv[j] <= ub[j]*cby[j] for j in network.M),name='UB')
    
    m.addConstr(sum(1-cby[j] for j in network.KO) == k, name='knapsack')
    m.addConstrs((cby[j] == 1 for j in network.M if j not in network.KO))


    imodel = gp.Model()
    fv = imodel.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='fv')
    imodel.params.LogToConsole = 0
    
    # Inner model Objective
    imodel.setObjective(2000*fv[network.biomass] - fv[network.chemical], GRB.MAXIMIZE)
    
    imodel.addMConstr(network.S,fv,'=',network.b,name='Stoi')
    # imodel.addConstrs((gp.quicksum(network.S[i,j]*fv[j] for j in network.M) == 0 for i in network.N),name='S2')
    
    imodel.addConstr(fv[network.biomass] >= minprod, name='target2')

    imodel.update()
    
    m._inner = imodel.copy()
    # m._innerv = fv
    m._vars = cbv
    m._varsy = cby
    m.Params.lazyConstraints = 1
    m._pbnd = -1000
    m._cbcnt = 0
    m._sv = [0 for i in network.M]
    m._sy = []

    m.Params.OptimalityTol = network.infeas
    m.Params.IntFeasTol = network.infeas
    m.Params.FeasibilityTol = network.infeas
    m.Params.Presolve = 0
    if not log: m.Params.OutputFlag = 0
    if speed: m.Params.NodefileStart = 0.5
    if threads: m.Params.Threads = 6

    m.Params.TimeLimit = network.time_limit

    m.optimize(lazycall)
    
    cb_time = m.Runtime

    if m.status == GRB.OPTIMAL:
        ys = [m.getVarByName('cby[%d]'%j).x for j in network.M]
        vs = [m.getVarByName('cbv[%d]'%j).x for j in network.M]
        del_strat_cb = [network.Rxn[i] for i in network.M if ys[i] < .5]
        soltype = 'Optimal'
    elif m.status == GRB.TIME_LIMIT:
        ys = [m.getVarByName('cby[%d]'%j).x for j in network.M]
        vs = [m.getVarByName('cbv[%d]'%j).x for j in network.M]
        del_strat_cb = [network.Rxn[i] for i in network.M if ys[i] < .5]
        soltype = 'Time_Limit'

    elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
        ys = ['all' for i in network.M]
        vs = ['~' for i in network.M]
        del_strat_cb = ['all']
        soltype = 'Infeasible'

    return Result_cb(network.Name,del_strat_cb,ys,vs,vij,cb_time,soltype,'CBP')

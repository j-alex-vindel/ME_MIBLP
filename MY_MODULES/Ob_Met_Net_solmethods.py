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
Result_cb = namedtuple('Result_cb',['MetNet','Strategy','Ys','Vs','Vij','Time','Soltype'])
Result = namedtuple('Result',['MetNet','Strategy','Ys','Vs','Time','Soltype'])
Result_inner = namedtuple('Result_IC',['Biomass','Chemical','Soltype'])
Pareto_point = namedtuple('P_point',['Biomass','Chemical'])



# ================================================ NO Presolve ==================================================================================

def CB_sol_2_NOP(network:M_Network=None,k:Ks=None,log:bool=True,speed:bool=False,threads:bool=False) -> Result_cb:
    '''
    cb = CB_solve_2_NOP(network=network,k=k,log=True,speed=False,threads=False)
        returns the bilevel optimistic solution as a named tupple 
    cb.MetNet    = Metabolic Network's name
    cb.Stragtegy = List of rxn to knockout
    cb.Ys        = Binary solution as a vector
    cb.Vs        = Optimal bilevel flows
    cb.Vij       = Flows in the inner problem
    cb.Time      = Solving time 
    cb.Soltype   = Type of solution [optimal, timelimit , infeasible]

    '''
    print(f'**** Solving Callbacks k={k} ****')
    print(f'# Variables (reactions in the network): {len(network.M)}')
    print('Current Infeasibility:',network.infeas,sep=' -> ')
    print('KO set: ',len(network.KO), ' reactions')
    print(f"MN: {network.Name}")

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

            model._vi, inner_status = inner(model._inner, model._yoj)

# ============================ Checking Inner Optimality Status ===================================
            # print('MIPSOL')
            # print(f"Objective = {cur_obj}")
            # print(f"Best Bound = {cur_bd}")
            # print(f"Curnt Pbnd = {model._pbnd}")

            if inner_status != GRB.OPTIMAL:
                # print(f"optimality cuts inner not optinal MIPNODE")
                model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
            else:
                vinner_biomas_value = round(model._vi[network.biomass],6)
                knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.KO)
                ki = (i for i in combinations(knockset_inner,k))

                # bestknownchem = cur_obj
                
                if model._pbnd - model._vi[network.chemical] >= -1e-6: # try with vij instead of model._vij to access the inner values, changed the value to -1e-6
                    # print(f'optimality cuts pbnd - vi[chemical')
                    model.cbLazy(sum(model._varsy[j] for j in knockset) >= 1)

                elif (abs(vij[network.biomass] - model._voj[network.biomass]) > 1e-6):
                    # print(f"Big M cut")
                    for comb in ki:
                        model.cbLazy(vinner_biomas_value <= model._vars[network.biomass] +
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
                model._ryoj = model.cbGetNodeRel(model._varsy) #  
                for i,y in enumerate(model._ryoj):
                    if model._ryoj[y] >= 0.8:
                        model._ryoj[y] = 1.0
                    elif model._ryoj[y] <= 0.2:
                        model._ryoj[y] = 0.0
                    else:
                        model._ryoj[y] = 1.0
                knock = [i for i,y in enumerate(model._ryoj) if model._yoj[i] < 1e-6]
                if sum(model._ryoj.values()) != len(model._ryoj)-k:
                    return
                else:
                    model._vi, inner_status = inner(model._inner,model._ryoj)
                    # print(f"Optimality Code - {inner_status}")

                    if inner_status != GRB.OPTIMAL:
                        # print('Optimality cuts - inner not optimal')
                        model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
                    
                    elif (model._vi[network.chemical] > model._pbnd) and (model._vi[network.chemical] > model._sv[network.chemical]): # added condition to set solution always better
                        model._sv = model._vi
                        model._sy = model._ryoj

                        # print(f"Set Solution ")
                        # print(f"biomas {model._vi[network.biomass]:.6f}")
                        # print(f"chemical {model._vi[network.chemical]:.6f}")
                        model.cbSetSolution(model._vars, model._sv)
                        model.cbSetSolution(model._varsy, model._sy)
                    
                    else:
                        # print(f"Optimality cuts - vi not larger or equal than pnbd")
                        model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
    
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
    imodel.setObjective(2000*fv[network.biomass] + fv[network.chemical], GRB.MAXIMIZE)
    
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

    return Result_cb(network.Name,del_strat_cb,ys,vs,vij,cb_time,soltype)

def MILP_solve_NOP(network:M_Network=None,k:Ks=None,log:bool=True,speed:bool=False,threads:bool=False) -> Result:
    '''
    MILP_solve(network=network,k=k)
        return Result[MetNet,Strategy,Flows,Time,Soltype]
    '''
    print(f'**** Solving ReacKnock k={k} ****')
    print(f'# Variables (reactions in the network): {len(network.M)}')
    print('Current Infeasibility:',network.infeas,sep=' -> ')
    print('KO set: ',len(network.KO), ' reactions')
    print(f"MN: {network.Name}")

    lb = copy.deepcopy(network.LB)
    lb[network.biomass] = network.minprod

    m = gp.Model()

    # Variables
    v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')
    y = m.addVars(network.M,vtype=GRB.BINARY,name='y')

    # Dual Variables
    l = m.addVars(network.N,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='lambda')
    a1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha1')
    b1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta1')
    a2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha1')
    b2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta1')
    a = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha')
    b = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta')
    
    # Objective
    m.setObjective((1*v[network.chemical]),GRB.MAXIMIZE)

    # Knapsack Constrs
    m.addConstrs((y[j] == 1 for j in network.M if j not in network.KO), name='y_essentials')

    m.addConstr(sum(1-y[j] for j in network.KO) == k, name='knapsack')

    # Stoichimetric Constrs
    m.addMConstr(network.S,v,'=',network.b,name='Stoi')
    # m.addConstrs((gp.quicksum(network.S[i,j] * v[j] for j in network.M) == network.b[i] for i in network.N),name='Stoichiometry')
    
    # Dual Objective
    m.addConstr((v[network.biomass] >= (sum(a1[j]*network.UB[j] - b1[j]*lb[j] for j in network.M)
     + sum(a2[j]*network.UB[j] - b2[j]*lb[j] for j in network.M))),name='dual-objective')
    
    # Dual Constraints
    m.addConstrs((gp.quicksum(network.S.transpose()[i,j]*l[j] for j in network.N)
              - b[i]
              + a[i] - b2[i] + a2[i]
               == network.c[i] for i in network.M)
             ,name='S_dual')

    # m.addConstr((gp.quicksum(network.S.transpose()[network.biomas,j]*l[j] for j in network.N)
    #         - b[network.biomas]
    #         + a[network.biomas]
    #         - b2[network.biomas] + a2[network.biomas] == 1), name='Sdual_t')
    
    # Linearization
    m.addConstrs((a1[j] <= network.BM*y[j] for j in network.M),name='l1_a1')

    m.addConstrs((a1[j] >= - network.BM*y[j] for j in network.M),name='l2_a1')

    m.addConstrs((a1[j] <= a[j] + network.BM*(1-y[j]) for j in network.M),name='l3_a1')

    m.addConstrs((a1[j] >= a[j] - network.BM*(1-y[j]) for j in network.M),name='l4_a1')

    m.addConstrs((b1[j] <= network.BM*y[j] for j in network.M),name='l1_b1')

    m.addConstrs((b1[j] >= -network.BM*y[j] for j in network.M),name='l2_b1')

    m.addConstrs((b1[j] <= b[j] + network.BM*(1-y[j]) for j in network.M),name='l3_b1')

    m.addConstrs((b1[j] >= b[j] - network.BM*(1-y[j]) for j in network.M),name='l4_b1')

    # Bounds
    m.addConstrs((lb[j]*y[j] <= v[j] for j in network.M), name='LB')
    m.addConstrs((v[j] <= network.UB[j]*y[j] for j in network.M), name='UB')

    m.addConstrs((lb[j] <= v[j] for j in network.M),name='lb')
    m.addConstrs((v[j] <= network.UB[j] for j in network.M),name='ub')

    m.Params.OptimalityTol = network.infeas
    m.Params.IntFeasTol = network.infeas
    m.Params.FeasibilityTol = network.infeas
    # m.Params.NodefileStart = 0.5
    m.Params.Presolve = 0
    if not log: m.Params.OutputFlag = 0
    if speed: m.Params.NodefileStart = 0.5
    if threads: m.Threads = 2
    m.optimize()

    s = m.Runtime
    if m.status == GRB.OPTIMAL:
        chem = m.getObjective().getValue()
        ys = [m.getVarByName('y[%d]'%j).x for j in network.M]
        vs = [m.getVarByName('v[%d]'%j).x for j in network.M]
        soltype = 'Optimal'
        del_strat = [network.Rxn[i] for i in network.M if ys[i] <.5]

    elif m.status == GRB.TIME_LIMIT:
        ys = [m.getVarByName('my[%d]'%j).x for j in network.M]
        vs = [m.getVarByName('mv[%d]'%j).x for j in network.M]
        del_strat = [network.Rxn[i] for i in network.M if ys[i] <.5]
        soltype = 'Time_limit'
    
    if m.status in (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED):
        # print('Model status: *** INFEASIBLE or UNBOUNDED ***')
        ys = ['$' for i in network.M]
        vs = ['~' for i in network.M]
        # print('Chemical:',vs[network.chemical],sep=' ^ ')
        # print('Biomass:',vs[network.biomass],sep=' ^ ')
        del_strat = 'all'

    print('*'*4,' FINISHED!!! ','*'*4)

    return  Result(network.Name,del_strat,ys, vs, s,soltype)

def Inner_check_vs_ys_NOP(network:M_Network=None,result_cb:Result_cb=None,result_milp:Result=None,criteria:str='both',milp:bool=False,cb:bool=False,objective:Inner_obj=None,log:bool=True) -> Result_inner:

    lb = copy.deepcopy(network.LB)
    ub = copy.deepcopy(network.UB)
    lb[network.biomass] = network.minprod

    if (milp and cb) or (not milp and not cb):
        return f"Choose milp=True,cb=False or milp=False,cb=True"
    elif milp and not cb:
        print('>> Checking on milp')
        vs = copy.deepcopy(result_milp.Vs)
        ys = copy.deepcopy(result_milp.Ys)
    elif cb and not milp:
        print('>> Checking on CB')
        vs = copy.deepcopy(result_cb.Vs)
        ys = copy.deepcopy(result_cb.Ys)

    if objective == 'biomass':
        objct = network.biomass
    elif objective == 'chemical':
        objct = network.chemical
    else:
        return f"No Inner Objective ('biomass' or 'chemical')"
    


    if criteria == 'both':

        m = gp.Model()

        v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')

        m.setObjective((1*v[objct]),GRB.MAXIMIZE)
        m.addMConstr(network.S,v,'=',network.b,name='Stoi')

        m.addConstrs((lb[j]*ys[j]<= v[j] for j in network.M),name='lb')

        m.addConstrs((ub[j]*ys[j] >= v[j] for j in network.M),name='ub')

        m.addConstrs((v[j] == vs[j] for j in network.M),name='hards_vs')

        m.Params.OptimalityTol = network.infeas
        m.Params.IntFeasTol = network.infeas
        m.Params.FeasibilityTol = network.infeas
        m.Params.Presolve = 0
        if not log: m.Params.OutputFlag = 0
        m.optimize()

        if m.status == GRB.OPTIMAL:
            vs =  [m.getVarByName('v[%s]'%a).x for a in network.M]
            soltype = 'optimal'
        elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
            vs = [1 if i == objct else 2000 for i in network.M]
            soltype = 'Infeasible'
    
    elif criteria == 'vs':
        m = gp.Model()

        v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')

        m.setObjective((1*v[objct]),GRB.MAXIMIZE)
        m.addMConstr(network.S,v,'=',network.b,name='Stoi')

        m.addConstrs((lb[j] <= v[j] for j in network.M),name='lb')

        m.addConstrs((ub[j] >= v[j] for j in network.M),name='ub')

        m.addConstrs((v[j] == vs[j] for j in network.M),name='hards_vs')

        m.Params.OptimalityTol = network.infeas
        m.Params.IntFeasTol = network.infeas
        m.Params.FeasibilityTol = network.infeas
        m.Params.Presolve = 0
        if not log: m.Params.OutputFlag = 0
        m.optimize()

        if m.status == GRB.OPTIMAL:
            vs =  [m.getVarByName('v[%s]'%a).x for a in network.M]
            soltype = 'optimal'
        elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
            vs = [1 if i == objct else 2000 for i in network.M]
            soltype = 'Infeasible'

    elif criteria == 'ys':
        m = gp.Model()

        v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')

        m.setObjective((1*v[objct]),GRB.MAXIMIZE)
        m.addMConstr(network.S,v,'=',network.b,name='Stoi')

        m.addConstrs((lb[j] * ys[j] <= v[j] for j in network.M),name='lb')

        m.addConstrs((ub[j] * ys[j] >= v[j] for j in network.M),name='ub')

        m.Params.OptimalityTol = network.infeas
        m.Params.IntFeasTol = network.infeas
        m.Params.FeasibilityTol = network.infeas
        m.Params.Presolve = 0
        if not log: m.Params.OutputFlag = 0
        m.optimize()

        if m.status == GRB.OPTIMAL:
            vs =  [m.getVarByName('v[%s]'%a).x for a in network.M]
            soltype = 'optimal'
        elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
            vs = [1 if i == objct else 2000 for i in network.M]
            soltype = 'Infeasible'
    else:
        return f"Choose the type if inner check [vs,ys,both]"

    return Result_inner(vs[network.biomass],vs[network.chemical],soltype) 

def Pareto_Frontier(network:M_Network=None,coef:int=None) -> Pareto_point:
    lb = copy.deepcopy(network.LB)
    ub = copy.deepcopy(network.UB)

    m = gp.Model()

    v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')

    m.setObjective((coef*v[network.chemical] + (1-coef)*v[network.biomass]),GRB.MAXIMIZE)

    m.addMConstr(network.S,v,'=',network.b,name='Stoi')
    m.addConstrs((lb[j] <= v[j] for j in network.M), name='LBwt')
    m.addConstrs((ub[j] >= v[j] for j in network.M), name='UBwt')
    
    m.Params.OptimalityTol = network.infeas
    m.Params.IntFeasTol = network.infeas
    m.Params.FeasibilityTol = network.infeas
    m.Params.Outputmag = 0
    m.optimize()
    if m.status == GRB.OPTIMAL:
        vs =  [m.getVarByName('v[%s]'%a).x for a in network.M] 

    bio_obj = vs[network.biomass]
    che_obj = vs[network.chemical]
    
    return Pareto_point(bio_obj,che_obj)

def CB_sol_2_P(network:M_Network=None, k:Ks=None,log:bool=None,speed:bool=False,threads:bool=False) -> Result_cb:
    '''
        cb = CB_solve_2_NOP(network=network,k=k,log=True,speed=False,threads=False)
            returns the bilevel pesimistic solution as a named tupple 
        cb.MetNet    = Metabolic Network's name
        cb.Stragtegy = List of rxn to knockout
        cb.Ys        = Binary solution as a vector
        cb.Vs        = Optimal bilevel flows
        cb.Vij       = Flows in the inner problem
        cb.Time      = Solving time 
        cb.Soltype   = Type of solution [optimal, timelimit , infeasible]

        '''
    print(f'**** Solving Callbacks k={k} ****')
    print(f'# Variables (reactions in the network): {len(network.M)}')
    print('Current Infeasibility:',network.infeas,sep=' -> ')
    print('KO set: ',len(network.KO), ' reactions')
    print(f"MN: {network.Name}")

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
                vinner_biomas_value = round(model._vi[network.biomass],6)
                knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.KO)
                ki = (i for i in combinations(knockset_inner,k))

                # bestknownchem = cur_obj
                
                if model._pbnd - model._vi[network.chemical] >= -1e-6: # try with vij instead of model._vij to access the inner values, changed the value to -1e-6
                    # print(f'optimality cuts pbnd - vi[chemical')
                    model.cbLazy(sum(model._varsy[j] for j in knockset) >= 1)

                elif (abs(vij[network.biomass] - model._voj[network.biomass]) > 1e-6):
                    # print(f"Big M cut")
                    for comb in ki:
                        model.cbLazy(vinner_biomas_value <= model._vars[network.biomass] +
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
                model._ryoj = model.cbGetNodeRel(model._varsy) #  
                for i,y in enumerate(model._ryoj):
                    if model._ryoj[y] >= 0.8:
                        model._ryoj[y] = 1.0
                    elif model._ryoj[y] <= 0.2:
                        model._ryoj[y] = 0.0
                    else:
                        model._ryoj[y] = 1.0
                knock = [i for i,y in enumerate(model._ryoj) if model._yoj[i] < 1e-6]
                if sum(model._ryoj.values()) != len(model._ryoj)-k:
                    return
                else:
                    model._vi, inner_status = inner(model._inner,model._ryoj)
                    # print(f"Optimality Code - {inner_status}")

                    if inner_status != GRB.OPTIMAL:
                        # print('Optimality cuts - inner not optimal')
                        model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
                    
                    elif (model._vi[network.chemical] > model._pbnd) and (model._vi[network.chemical] > model._sv[network.chemical]): # added condition to set solution always better
                        model._sv = model._vi
                        model._sy = model._ryoj

                        # print(f"Set Solution ")
                        # print(f"biomas {model._vi[network.biomass]:.6f}")
                        # print(f"chemical {model._vi[network.chemical]:.6f}")
                        model.cbSetSolution(model._vars, model._sv)
                        model.cbSetSolution(model._varsy, model._sy)
                    
                    else:
                        # print(f"Optimality cuts - vi not larger or equal than pnbd")
                        model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
    
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

    return Result_cb(network.Name,del_strat_cb,ys,vs,vij,cb_time,soltype)


# =============================== Old Funcs ====================================================================================================


# def MILP_solve(network:M_Network=None,k:Ks=None) -> Result:
#     '''
#     MILP_solve(network=network,k=k)
#         return Result[MetNet,Strategy,Flows,Time,Soltype]
#     '''
#     print(f'**** Solving ReacKnock k={k} ****')
#     print(f'# Variables (reactions in the network): {len(network.M)}')
#     print('Current Infeasibility:',network.infeas,sep=' -> ')
#     print('KO set: ',len(network.KO), ' reactions')
#     print(f"MN: {network.Name}")

#     lb = copy.deepcopy(network.LB)
#     lb[network.biomass] = network.minprod

#     m = gp.Model()

#     # Variables
#     v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')
#     y = m.addVars(network.M,vtype=GRB.BINARY,name='y')

#     # Dual Variables
#     l = m.addVars(network.N,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='lambda')
#     a1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha1')
#     b1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta1')
#     a2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha1')
#     b2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta1')
#     a = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha')
#     b = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta')
    
#     # Objective
#     m.setObjective((1*v[network.chemical]),GRB.MAXIMIZE)

#     # Knapsack Constrs
#     m.addConstrs((y[j] == 1 for j in network.M if j not in network.KO), name='y_essentials')

#     m.addConstr(sum(1-y[j] for j in network.KO) == k, name='knapsack')

#     # Stoichimetric Constrs
#     m.addConstrs((gp.quicksum(network.S[i,j] * v[j] for j in network.M) == network.b[i] for i in network.N),name='Stoichiometry')
    
#     # Dual Objective
#     m.addConstr((v[network.biomass] >= (sum(a1[j]*network.UB[j] - b1[j]*lb[j] for j in network.M)
#      + sum(a2[j]*network.UB[j] - b2[j]*lb[j] for j in network.M))),name='dual-objective')
    
#     # Dual Constraints
#     m.addConstrs((gp.quicksum(network.S.transpose()[i,j]*l[j] for j in network.N)
#               - b[i]
#               + a[i] - b2[i] + a2[i]
#                == network.c[i] for i in network.M)
#              ,name='S_dual')

#     # m.addConstr((gp.quicksum(network.S.transpose()[network.biomas,j]*l[j] for j in network.N)
#     #         - b[network.biomas]
#     #         + a[network.biomas]
#     #         - b2[network.biomas] + a2[network.biomas] == 1), name='Sdual_t')
    
#     # Linearization
#     m.addConstrs((a1[j] <= network.BM*y[j] for j in network.M),name='l1_a1')

#     m.addConstrs((a1[j] >= - network.BM*y[j] for j in network.M),name='l2_a1')

#     m.addConstrs((a1[j] <= a[j] + network.BM*(1-y[j]) for j in network.M),name='l3_a1')

#     m.addConstrs((a1[j] >= a[j] - network.BM*(1-y[j]) for j in network.M),name='l4_a1')

#     m.addConstrs((b1[j] <= network.BM*y[j] for j in network.M),name='l1_b1')

#     m.addConstrs((b1[j] >= -network.BM*y[j] for j in network.M),name='l2_b1')

#     m.addConstrs((b1[j] <= b[j] + network.BM*(1-y[j]) for j in network.M),name='l3_b1')

#     m.addConstrs((b1[j] >= b[j] - network.BM*(1-y[j]) for j in network.M),name='l4_b1')

#     # Bounds
#     m.addConstrs((lb[j]*y[j] <= v[j] for j in network.M), name='LB')
#     m.addConstrs((v[j] <= network.UB[j]*y[j] for j in network.M), name='UB')

#     m.addConstrs((lb[j] <= v[j] for j in network.M),name='lb')
#     m.addConstrs((v[j] <= network.UB[j] for j in network.M),name='ub')

#     m.Params.OptimalityTol = network.infeas
#     m.Params.IntFeasTol = network.infeas
#     m.Params.FeasibilityTol = network.infeas
#     m.Params.NodefileStart = 0.5
#     m.optimize()

#     s = m.Runtime
#     if m.status == GRB.OPTIMAL:
#         chem = m.getObjective().getValue()
#         ys = [m.getVarByName('y[%d]'%j).x for j in network.M]
#         vs = [m.getVarByName('v[%d]'%j).x for j in network.M]
#         soltype = 'Optimal'
#         del_strat = [network.Rxn[i] for i in network.M if ys[i] <.5]

#     elif m.status == GRB.TIME_LIMIT:
#         ys = [m.getVarByName('my[%d]'%j).x for j in network.M]
#         vs = [m.getVarByName('mv[%d]'%j).x for j in network.M]
#         del_strat = [network.Rxn[i] for i in network.M if ys[i] <.5]
#         soltype = 'Time_limit'
    
#     if m.status in (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED):
#         print('Model status: *** INFEASIBLE or UNBOUNDED ***')
#         ys = ['$' for i in network.M]
#         vs = ['~' for i in network.M]
#         print('Chemical:',vs[network.chemical],sep=' ^ ')
#         print('Biomass:',vs[network.biomass],sep=' ^ ')
#         del_strat = 'all'

#     print('*'*4,' FINISHED!!! ','*'*4)

#     return  Result(network.Name,del_strat,ys, vs, s,soltype)

# def CB_sol(network:M_Network=None,k:Ks=None) -> Result_cb:
#     '''
#     CB_solve(network=network,k=k)
#         return Result_cb['MetNet','Strategy','Ys','Voj','Vij','Time','Soltype']
#     '''
#     print(f'**** Solving Callbacks k={k} ****')
#     print(f'# Variables (reactions in the network): {len(network.M)}')
#     print('Current Infeasibility:',network.infeas,sep=' -> ')
#     print('KO set: ',len(network.KO), ' reactions')
#     print(f"MN: {network.Name}")

#     lb = copy.deepcopy(network.LB)
#     ub = copy.deepcopy(network.UB)
#     minprod = copy.deepcopy(network.minprod)

#     def inner(imodel, yoj:Vector):
#         global vij

#         imodel.setAttr('LB',imodel.getVars(),[lb[j]*yoj[j] for j in network.M])
#         imodel.setAttr('UB',imodel.getVars(),[ub[j]*yoj[j] for j in network.M])
#         imodel.Params.OptimalityTol = network.infeas
#         imodel.Params.IntFeasTol = network.infeas
#         imodel.Params.FeasibilityTol = network.infeas
#         imodel.optimize()
#         status = imodel.status
#         if status == GRB.OPTIMAL:
#             vij = [round(imodel.getVarByName('fv[%s]'%a).x,6) for a in network.M] # rounded inner values
#         elif status in (GRB.INFEASIBLE, GRB.UNBOUNDED, GRB.INF_OR_UNBD):
#             vij = [2000 if i == network.biomass else yoj[i] for i in network.M]
#         return vij,status

#     def lazycall(model,where):
#         if where == GRB.Callback.MIPSOL:
#             model._voj = model.cbGetSolution(model._vars)
#             model._yoj = model.cbGetSolution(model._varsy)
#             knockset =  [i for i,y in enumerate(model._yoj) if model._yoj[i] < 1e-6]

#             if len(knockset) != k:
#                 return
#             cur_obj = round(model.cbGet(GRB.Callback.MIPSOL_OBJBST),6)
#             cur_bd = round(model.cbGet(GRB.Callback.MIPSOL_OBJBND),6)
#             print('MIPSOL')
#             print(f"Obj = {cur_obj}")
#             print(f"Pnbd = {model._pbnd}")
#             print(f"Cur Bd = {cur_bd}")

#             model._vij, inner_status = inner(model._inner, model._yoj)
# # ============================ Checking Inner Optimality Status ===================================
#             if inner_status != GRB.OPTIMAL: # if the inner model is infeasible 
#                 model.cbLazy(sum(model._varsy[j] for j in knockset) >=1) # both ys can't be zero
#             else:
#                 vinner_biomas_value = model._vij[network.biomass]
#                 knockset_inner = (i for i,y in enumerate(model._vij) if abs(model._vij[i]) < 1e-6 and i in network.KO)
#                 ki = (i for i in combinations(knockset_inner,k))

#                 # bestknownchem = cur_obj
                
#                 if model._pbnd - model._vij[network.chemical] >= -1e-6: 
#                     model.cbLazy(sum(model._varsy[j] for j in knockset) >= 1) #  both ys can't be zero

#                 elif (abs(model._vij[network.biomass] - model._voj[network.biomass]) > 1e-6):
#                     for comb in ki:
#                         model.cbLazy(vinner_biomas_value <= model._vars[network.biomass] +
#                             (math.ceil(model._vij[network.biomass]*10)/10) *(sum(model._varsy[f] for f in comb)))

#                 else:
#                     model._pbnd = cur_obj
# # =============================================================================================================================                        

#         elif where == GRB.Callback.MIPNODE:
#             status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
#             if status == GRB.OPTIMAL:
#                 print('MIPNODE')
#                 mipobj = model.cbGet(GRB.Callback.MIPNODE_OBJBST)
#                 mipbnd = model.cbGet(GRB.Callback.MIPNODE_OBJBND)
#                 print(f"Obj= {mipobj}")
#                 print(f"Bnd= {mipbnd}")
#                 print(f"Pbnd={model._pbnd}")
#                 model._cbcnt += 1
#                 if model._cbcnt % 10 == 0:
#                     return
#                 model._voj = model.cbGetNodeRel(model._vars)
#                 model._ryoj = model.cbGetNodeRel(model._varsy) #  
#                 for i,y in enumerate(model._ryoj):
#                     if model._ryoj[y] >= 0.8:
#                         model._ryoj[y] = 1.0
#                     elif model._ryoj[y] <= 0.2:
#                         model._ryoj[y] = 0.0
#                     else:
#                         model._ryoj[y] = 1.0
#                 knock = [i for i,y in enumerate(model._ryoj) if model._yoj[i] < 1e-6]
#                 if sum(model._ryoj.values()) != len(model._ryoj)-k:
#                     return
#                 else:
#                     model._vij, inner_status = inner(model._inner,model._ryoj)
#                     print(f"Optimality Code - {inner_status}")

#                     if inner_status != GRB.OPTIMAL:
#                         model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
                    
#                     elif model._vij[network.chemical] > model._pbnd:

#                         print(f"Set Solution ")
#                         print(f"biomas {model._vij[network.biomass]}")
#                         print(f"chemical {model._vij[network.chemical]}")
#                         model.cbSetSolution(model._vars, model._vij)
#                         model.cbSetSolution(model._varsy, model._ryoj)
                    
#                     else:
#                         model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
    
#     m = gp.Model()
#     cbv = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='cbv')
#     cby = m.addVars(network.M,vtype=GRB.BINARY,name='cby')

#     m.setObjective(1*cbv[network.chemical],GRB.MAXIMIZE)
    
#     m.addMConstr(network.S,cbv,'=',network.b,name='Stoi')
#     # m.addConstrs((gp.quicksum(network.S[i,j]*cbv[j] for j in network.M) == 0 for i in network.N),name='Stoichiometry')

#     m.addConstr(cbv[network.biomass] >= minprod, name='target')

#     m.addConstrs((lb[j]*cby[j] <= cbv[j] for j in network.M),name='LB')
#     m.addConstrs((cbv[j] <= ub[j]*cby[j] for j in network.M),name='UB')
    
#     m.addConstr(sum(1-cby[j] for j in network.KO) == k, name='knapsack')
#     m.addConstrs((cby[j] == 1 for j in network.M if j not in network.KO))

#     m._vars = cbv
#     m._varsy = cby
#     m.Params.lazyConstraints = 1
#     m.Params.NodefileStart = 0.5

#     imodel = gp.Model()
#     fv = imodel.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='fv')
#     imodel.params.LogToConsole = 0
#     imodel.setObjective(2000*fv[network.biomass] + fv[network.chemical], GRB.MAXIMIZE)
    
#     imodel.addMConstr(network.S,fv,'=',network.b,name='Stoi')
#     # imodel.addConstrs((gp.quicksum(network.S[i,j]*fv[j] for j in network.M) == 0 for i in network.N),name='S2')
    
#     imodel.addConstr(fv[network.biomass] >= minprod, name='target2')

#     imodel.update()
    
#     m._inner = imodel.copy()
#     # m._innerv = fv
#     m._pbnd = -1000
#     m._cbcnt = 0

#     m.Params.OptimalityTol = network.infeas
#     m.Params.IntFeasTol = network.infeas
#     m.Params.FeasibilityTol = network.infeas
#     # m.Params.OutputFlag = 0
#     m.Params.TimeLimit = network.time_limit
#     # m.Params.NodefileStart = 0
#     # m.Params.Threads = 4
    
#     m.optimize(lazycall)
    
#     # m.setParam(GRB.Param.PoolSolutions, 10)
#     # m.setParam(GRB.Param.PoolSearchMode, 2)
#     # m.setParam(GRB.Param.PoolGap, 0.01)
#     cb_time = m.Runtime
#     # nsolutions = m.SolCount

#     if m.status == GRB.OPTIMAL:
#         ys = [m.getVarByName('cby[%d]'%j).x for j in network.M]
#         vs = [m.getVarByName('cbv[%d]'%j).x for j in network.M]
#         del_strat_cb = [network.Rxn[i] for i in network.M if ys[i] < .5]
#         soltype = 'Optimal'
#     elif m.status == GRB.TIME_LIMIT:
#         ys = [m.getVarByName('cby[%d]'%j).x for j in network.M]
#         vs = [m.getVarByName('cbv[%d]'%j).x for j in network.M]
#         del_strat_cb = [network.Rxn[i] for i in network.M if ys[i] < .5]
#         soltype = 'Time_Limit'

#     elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
#         ys = ['all' for i in network.M]
#         vs = ['~' for i in network.M]
#         del_strat_cb = ['all']
#         soltype = 'Infeasible'

#     return Result_cb(network.Name,del_strat_cb,ys,vs,vij,cb_time,soltype)


# def CB_sol_2(network:M_Network=None,k:Ks=None) -> Result_cb:
#     '''
#     CB_solve_2(same as CB_solve but with the difference in the inner model)
#     '''
#     print(f'**** Solving Callbacks k={k} ****')
#     print(f'# Variables (reactions in the network): {len(network.M)}')
#     print('Current Infeasibility:',network.infeas,sep=' -> ')
#     print('KO set: ',len(network.KO), ' reactions')
#     print(f"MN: {network.Name}")

#     lb = copy.deepcopy(network.LB)
#     ub = copy.deepcopy(network.UB)
#     minprod = copy.deepcopy(network.minprod)
#     lb[network.biomass] = minprod

#     def inner(imodel, yoj:Vector):
#         global vij

#         imodel.setAttr('LB',imodel.getVars(),[lb[j]*yoj[j] for j in network.M])
#         imodel.setAttr('UB',imodel.getVars(),[ub[j]*yoj[j] for j in network.M])
#         imodel.Params.OptimalityTol = network.infeas
#         imodel.Params.IntFeasTol = network.infeas
#         imodel.Params.FeasibilityTol = network.infeas
#         imodel.optimize()
#         status = imodel.status
#         if status == GRB.OPTIMAL:
#             vij = [imodel.getVarByName('fv[%s]'%a).x for a in network.M] # rounded inner values
#         elif status in (GRB.INFEASIBLE, GRB.UNBOUNDED, GRB.INF_OR_UNBD):
#             vij = [2000 if i == network.biomass else yoj[i] for i in network.M]
#         return vij,status

#     def lazycall(model,where):
#         if where == GRB.Callback.MIPSOL:
#             model._voj = model.cbGetSolution(model._vars)
#             model._yoj = model.cbGetSolution(model._varsy)
#             knockset =  [i for i,y in enumerate(model._yoj) if model._yoj[i] < 1e-6]

#             if len(knockset) != k:
#                 return
#             cur_obj = round(model.cbGet(GRB.Callback.MIPSOL_OBJBST),6)
#             cur_bd = round(model.cbGet(GRB.Callback.MIPSOL_OBJBND),6)

#             model._vi, inner_status = inner(model._inner, model._yoj)

# # ============================ Checking Inner Optimality Status ===================================
#             print('MIPSOL')
#             print(f"Objective = {cur_obj}")
#             print(f"Best Bound = {cur_bd}")
#             print(f"Curnt Pbnd = {model._pbnd}")

#             if inner_status != GRB.OPTIMAL:
#                 print(f"optimality cuts inner not optinal MIPNODE")
#                 model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
#             else:
#                 vinner_biomas_value = round(model._vi[network.biomass],6)
#                 knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.KO)
#                 ki = (i for i in combinations(knockset_inner,k))

#                 # bestknownchem = cur_obj
                
#                 if model._pbnd - model._vi[network.chemical] >= -1e-6: # try with vij instead of model._vij to access the inner values, changed the value to -1e-6
#                     print(f'optimality cuts pbnd - vi[chemical')
#                     model.cbLazy(sum(model._varsy[j] for j in knockset) >= 1)

#                 elif (abs(vij[network.biomass] - model._voj[network.biomass]) > 1e-6):
#                     print(f"Big M cut")
#                     for comb in ki:
#                         model.cbLazy(vinner_biomas_value <= model._vars[network.biomass] +
#                             (math.ceil(vij[network.biomass]*10)/10) *(sum(model._varsy[f] for f in comb)))

#                 else:
#                     print(f"pbnd = cur_obj")
#                     model._pbnd = cur_obj
# # =============================================================================================================================                        

#         elif where == GRB.Callback.MIPNODE:
#             status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
#             if status == GRB.OPTIMAL:
#                 print('MIPNODE')
#                 mipobj = model.cbGet(GRB.Callback.MIPNODE_OBJBST)
#                 mipbnd = model.cbGet(GRB.Callback.MIPNODE_OBJBND)
#                 print(f"Obj= {mipobj}")
#                 print(f"Bnd= {mipbnd}")
#                 print(f"Curnt set solution v: {model._sv[network.chemical]:.6f}")
#                 print(f"Curnt set solution len y: {len(model._sy)}")
#                 print(f"Curnt PBND: {model._pbnd}")
#                 model._cbcnt += 1
#                 if model._cbcnt % 10 == 0:
#                     return
#                 print(f"CB counts : {model._cbcnt}")
#                 model._voj = model.cbGetNodeRel(model._vars)
#                 model._ryoj = model.cbGetNodeRel(model._varsy) #  
#                 for i,y in enumerate(model._ryoj):
#                     if model._ryoj[y] >= 0.8:
#                         model._ryoj[y] = 1.0
#                     elif model._ryoj[y] <= 0.2:
#                         model._ryoj[y] = 0.0
#                     else:
#                         model._ryoj[y] = 1.0
#                 knock = [i for i,y in enumerate(model._ryoj) if model._yoj[i] < 1e-6]
#                 if sum(model._ryoj.values()) != len(model._ryoj)-k:
#                     return
#                 else:
#                     model._vi, inner_status = inner(model._inner,model._ryoj)
#                     print(f"Optimality Code - {inner_status}")

#                     if inner_status != GRB.OPTIMAL:
#                         print('Optimality cuts - inner not optimal')
#                         model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
                    
#                     elif (model._vi[network.chemical] > model._pbnd) and (model._vi[network.chemical] > model._sv[network.chemical]): # added condition to set solution always better
#                         model._sv = model._vi
#                         model._sy = model._ryoj

#                         print(f"Set Solution ")
#                         print(f"biomas {model._vi[network.biomass]:.6f}")
#                         print(f"chemical {model._vi[network.chemical]:.6f}")
#                         model.cbSetSolution(model._vars, model._sv)
#                         model.cbSetSolution(model._varsy, model._sy)
                    
#                     else:
#                         print(f"Optimality cuts - vi not larger or equal than pnbd")
#                         model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
    
#     m = gp.Model()
#     cbv = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='cbv')
#     cby = m.addVars(network.M,vtype=GRB.BINARY,name='cby')

#     m.setObjective(1*cbv[network.chemical],GRB.MAXIMIZE)
    
#     m.addMConstr(network.S,cbv,'=',network.b,name='Stoi')
#     # m.addConstrs((gp.quicksum(network.S[i,j]*cbv[j] for j in network.M) == 0 for i in network.N),name='Stoichiometry')

#     m.addConstr(cbv[network.biomass] >= minprod, name='target')

#     m.addConstrs((lb[j]*cby[j] <= cbv[j] for j in network.M),name='LB')
#     m.addConstrs((cbv[j] <= ub[j]*cby[j] for j in network.M),name='UB')
    
#     m.addConstr(sum(1-cby[j] for j in network.KO) == k, name='knapsack')
#     m.addConstrs((cby[j] == 1 for j in network.M if j not in network.KO))

#     m._vars = cbv
#     m._varsy = cby
#     m.Params.lazyConstraints = 1
#     m.Params.NodefileStart = 0.5

#     imodel = gp.Model()
#     fv = imodel.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='fv')
#     imodel.params.LogToConsole = 0
#     imodel.setObjective(2000*fv[network.biomass] + fv[network.chemical], GRB.MAXIMIZE)
    
#     imodel.addMConstr(network.S,fv,'=',network.b,name='Stoi')
#     # imodel.addConstrs((gp.quicksum(network.S[i,j]*fv[j] for j in network.M) == 0 for i in network.N),name='S2')
    
#     imodel.addConstr(fv[network.biomass] >= minprod, name='target2')

#     imodel.update()
    
#     m._inner = imodel.copy()
#     # m._innerv = fv
#     m._pbnd = -1000
#     m._cbcnt = 0
#     m._sv = [0 for i in network.M]
#     m._sy = []

#     m.Params.OptimalityTol = network.infeas
#     m.Params.IntFeasTol = network.infeas
#     m.Params.FeasibilityTol = network.infeas
#     # m.Params.OutputFlag = 0
#     m.Params.TimeLimit = network.time_limit
#     # m.Params.NodefileStart = 0
#     # m.Params.Threads = 4
    
#     m.optimize(lazycall)
    
#     # m.setParam(GRB.Param.PoolSolutions, 10)
#     # m.setParam(GRB.Param.PoolSearchMode, 2)
#     # m.setParam(GRB.Param.PoolGap, 0.01)
#     cb_time = m.Runtime
#     # nsolutions = m.SolCount

#     if m.status == GRB.OPTIMAL:
#         ys = [m.getVarByName('cby[%d]'%j).x for j in network.M]
#         vs = [m.getVarByName('cbv[%d]'%j).x for j in network.M]
#         del_strat_cb = [network.Rxn[i] for i in network.M if ys[i] < .5]
#         soltype = 'Optimal'
#     elif m.status == GRB.TIME_LIMIT:
#         ys = [m.getVarByName('cby[%d]'%j).x for j in network.M]
#         vs = [m.getVarByName('cbv[%d]'%j).x for j in network.M]
#         del_strat_cb = [network.Rxn[i] for i in network.M if ys[i] < .5]
#         soltype = 'Time_Limit'

#     elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
#         ys = ['all' for i in network.M]
#         vs = ['~' for i in network.M]
#         del_strat_cb = ['all']
#         soltype = 'Infeasible'

#     return Result_cb(network.Name,del_strat_cb,ys,vs,vij,cb_time,soltype)


# def Inner_check_ys(network:M_Network=None,result_ys:Result=None, objective:Inner_obj=None)->Result_inner:
#     '''
#     Inner_check_sol(network:Metabolic Network,ys:Ys Vector)
#         return soltype()
#     '''
#     lb = copy.deepcopy(network.LB)
#     ub = copy.deepcopy(network.UB)
#     lb[network.biomass] = network.minprod

#      # Objective
#     if objective == 'biomass':
#         objct = network.biomass
#     elif objective == 'chemical':
#         objct = network.chemical
#     else:
#         return f"No Inner Objective ('biomass' or 'chemical')"
    
#     m = gp.Model()

#     v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')

#     m.setObjective((1*v[objct]),GRB.MAXIMIZE)
#     m.addMConstr(network.S,v,'=',network.b,name='Stoi')

#     m.addConstrs((lb[j]*result_ys.Ys[j] <= v[j] for j in network.M),name='lb')

#     m.addConstrs((ub[j]*result_ys.Ys[j] >= v[j] for j in network.M),name='ub')

#     # m.addConstrs((cbv[j] <= ub[j]*cby[j] for j in network.M),name='UB')

#     m.Params.OptimalityTol = network.infeas
#     m.Params.IntFeasTol = network.infeas
#     m.Params.FeasibilityTol = network.infeas
#     # m.Params.OutputFlag = 0
#     m.optimize()
#     if m.status == GRB.OPTIMAL:
#         vs =  [m.getVarByName('v[%s]'%a).x for a in network.M]
#         soltype = 'optimal'
#     elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
#         vs = [1 if i == objct else 2000 for i in network.M]
#         soltype = 'Infeasible'
    
#     return Result_inner(vs[network.biomass],vs[network.chemical],soltype)

# def Inner_check_vs(network:M_Network=None,result_cb:Result_cb=None,result_milp:Result=None,milp:bool=None,cb:bool=None,objective:Inner_obj=None) -> Result_inner:

#     lb = copy.deepcopy(network.LB)
#     ub = copy.deepcopy(network.UB)
#     lb[network.biomass] = network.minprod

#     if (milp and cb) or (not milp and not cb):
#         return f"Choose milp=True,cb=False or milp=False,cb=True"
#     elif milp and not cb:
#         vs = copy.deepcopy(result_milp.Vs)
#     elif cb and not milp:
#         vs = copy.deepcopy(result_cb.Voj)
    
#     if objective == 'biomass':
#         objct = network.biomass
#     elif objective == 'chemical':
#         objct = network.chemical
#     else:
#         return f"No Inner Objective ('biomass' or 'chemical')"
#     m = gp.Model()

#     v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')

#     m.setObjective((1*v[objct]),GRB.MAXIMIZE)
#     m.addMConstr(network.S,v,'=',network.b,name='Stoi')

#     m.addConstrs((lb[j]<= v[j] for j in network.M),name='lb')

#     m.addConstrs((ub[j] >= v[j] for j in network.M),name='ub')

#     m.addConstrs((v[j] == vs[j] for j in network.M),name='hards_vs')

#     m.Params.OptimalityTol = network.infeas
#     m.Params.IntFeasTol = network.infeas
#     m.Params.FeasibilityTol = network.infeas
#     # m.Params.OutputFlag = 0
#     m.optimize()
    
#     if m.status == GRB.OPTIMAL:
#         vs =  [m.getVarByName('v[%s]'%a).x for a in network.M]
#         soltype = 'optimal'
#     elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
#         vs = [1 if i == objct else 2000 for i in network.M]
#         soltype = 'Infeasible'
    
#     return Result_inner(vs[network.biomass],vs[network.chemical],soltype) 

# def Inner_check_vs_ys(network:M_Network=None,result_cb:Result_cb=None,result_milp:Result=None,milp:bool=None,cb:bool=None,objective:Inner_obj=None) -> Result_inner:

#     lb = copy.deepcopy(network.LB)
#     ub = copy.deepcopy(network.UB)
#     lb[network.biomass] = network.minprod

#     if (milp and cb) or (not milp and not cb):
#         return f"Choose milp=True,cb=False or milp=False,cb=True"
#     elif milp and not cb:
#         vs = copy.deepcopy(result_milp.Vs)
#         ys = copy.deepcopy(result_milp.Ys)
#     elif cb and not milp:
#         vs = copy.deepcopy(result_cb.Voj)
#         ys = copy.deepcopy(result_cb.Ys)

#     if objective == 'biomass':
#         objct = network.biomass
#     elif objective == 'chemical':
#         objct = network.chemical
#     else:
#         return f"No Inner Objective ('biomass' or 'chemical')"
#     m = gp.Model()

#     v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')

#     m.setObjective((1*v[objct]),GRB.MAXIMIZE)
#     m.addMConstr(network.S,v,'=',network.b,name='Stoi')

#     m.addConstrs((lb[j]*ys[j]<= v[j] for j in network.M),name='lb')

#     m.addConstrs((ub[j]*ys[j] >= v[j] for j in network.M),name='ub')

#     m.addConstrs((v[j] == vs[j] for j in network.M),name='hards_vs')

#     m.Params.OptimalityTol = network.infeas
#     m.Params.IntFeasTol = network.infeas
#     m.Params.FeasibilityTol = network.infeas
#     # m.Params.OutputFlag = 0
#     m.optimize()

#     if m.status == GRB.OPTIMAL:
#         vs =  [m.getVarByName('v[%s]'%a).x for a in network.M]
#         soltype = 'optimal'
#     elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
#         vs = [1 if i == objct else 2000 for i in network.M]
#         soltype = 'Infeasible'

#     return Result_inner(vs[network.biomass],vs[network.chemical],soltype) 

    
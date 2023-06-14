import sys
import gurobipy as gp
from gurobipy import GRB
from typing import List, Type
from collections import namedtuple
import copy
from Ob_Met_Net import Met_Net
from itertools import combinations
import math 

MethodResult = namedtuple('Result',['MetNet','Strategy','Vs','Ys','KO_index','Time','Soltype','Method'])
Check = namedtuple('Check',['Biomass','Chemical','Soltype','of','Criteria'])
M_Network = Type[Met_Net]
Vector = List[int]

class BilevelMethods:
    '''
    A class for the bi-level methods used in the metabolic engineering context for optimistic - pesimistic optknock
    
    '''
    r_m = None
    r_o = None
    r_p = None
    r_c = None
    log = None

    def __init__(self,log:bool=None):
        self.log = log
    
    def info(self):

        text = ''' Bilevel methods to solve Gene Knockouts  '''

        return text

    def MILP(self,network:M_Network=None,K:int=None) -> MethodResult:
        mub = copy.deepcopy(network.UB)
        mlb = copy.deepcopy(network.LB)
        mlb[network.biomass] = network.minprod

        m = gp.Model()

        # Variables
        v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')
        vs = [v[i] for i in network.M]
    # ============================= Commented Section ===================
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
    # ============================= Commented Section ===================
        if network.KO is not None:
        # Knapsack Constrs
            m.addConstr((sum(1-y[j] for j in network.KO) == K), name='nonEssen') 
            m.addConstrs((y[j] == 1 for j in network.M if j not in network.KO),name='Essen')
            KOSum = 0
            # for i in network.M:
            #     if i not in network.KO and yss[i] < 1 -1e-6: 
            #         print("notko", i)
            #     if i not in network.KO:
            #         m.addConstr(y[i] >= 1)
            #     if  i in network.KO: 
            #         KOSum += yss[i]
            # print ("KOsum", KOSum)
            alpha = 1-1e-6
            # m.addConstr(sum(y[j] for j in network.M if j in network.KO) >= len(network.KO)-K - alpha, name='knapsack1')
            # m.addConstr(sum(y[j] for j in network.M if j in network.KO) <= len(network.KO)-K + alpha, name='knapsack2')

    # ============================= Commented Section ===================
        # Stoichimetric Constrs
        m.addMConstr(network.S,vs,'=',network.b,name='Stoi')
        
        # Dual Objective
        m.addConstr((v[network.biomass] >= (sum(a1[j]*mub[j] - b1[j]*mlb[j] for j in network.M)
        + sum(a2[j]*mub[j] - b2[j]*mlb[j] for j in network.M))),name='dual-objective')
        
        # Dual Constraints
        m.addConstrs((gp.quicksum(network.S.transpose()[i,j]*l[j] for j in network.N)
                - b[i]
                + a[i] - b2[i] + a2[i]
                == network.c[i] for i in network.M)
                ,name='S_dual')
        
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
        m.addConstrs((mlb[j]*y[j] <= v[j] for j in network.M), name='LB')
        m.addConstrs((v[j] <= mub[j]*y[j] for j in network.M), name='UB')

        m.addConstrs((mlb[j] <= v[j] for j in network.M),name='lb')
        m.addConstrs((v[j] <= mub[j] for j in network.M),name='ub')

        m.Params.OptimalityTol = network.infeas
        m.Params.IntFeasTol = network.infeas
        m.Params.FeasibilityTol = network.infeas
        m.Params.Presolve = 0
        if not self.log: m.Params.OutputFlag = 0
        
        m.optimize()
        
        m_time = m.Runtime
        if m.status == GRB.OPTIMAL:
            mys = [m.getVarByName('y[%d]'%j).x for j in network.M]
            vss = [m.getVarByName('v[%d]'%j).x for j in network.M]
            soltype = 'Optimal'
            del_strat = [network.Rxn[i] for i in network.M if mys[i] <.5]

        elif m.status == GRB.TIME_LIMIT:
            mys = [m.getVarByName('my[%d]'%j).x for j in network.M]
            vss = [m.getVarByName('mv[%d]'%j).x for j in network.M]
            del_strat = [network.Rxn[i] for i in network.M if mys[i] <.5]
            soltype = 'Time_limit'
        
        if m.status in (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED):
         
            mys = [0 for i in network.M]
            vss = [2000 if i in [network.biomass,network.chemical] else 0 for i in network.M]
            soltype = 'Infeasible or Unbounded'
            del_strat = [i for i in network.M]
        rys = self.__ysheuristic(mys)
        
        ko_index = [network.Rxn.index(i) for i in del_strat]

        self.r_m = MethodResult(network.Name,del_strat,vss,mys,ko_index,m_time,soltype,'M')
        return self.r_m

    def CB_O(self,network:M_Network=None,K:int=None) -> MethodResult:
        
        olb = copy.deepcopy(network.LB)
        oub = copy.deepcopy(network.UB)
        ominprod = copy.deepcopy(network.minprod)
        olb[network.biomass] = ominprod
        def inner(imodel, yoj:Vector):

            imodel.setAttr('LB',imodel.getVars(),[olb[j]*yoj[j] for j in network.M])
            imodel.setAttr('UB',imodel.getVars(),[oub[j]*yoj[j] for j in network.M])
            imodel.Params.OptimalityTol = network.infeas
            imodel.Params.IntFeasTol = network.infeas
            imodel.Params.FeasibilityTol = network.infeas
            imodel.Params.Presolve = 0
            imodel.optimize()
            status = imodel.status
            if status == GRB.OPTIMAL:
                vij = [imodel.getVarByName('fv[%s]'%a).x for a in network.M] 
            elif status in (GRB.INFEASIBLE, GRB.UNBOUNDED, GRB.INF_OR_UNBD):
                vij = [2000 if i in [ network.biomass,network.chemical] else yoj[i] for i in network.M]
            return vij,status

        def lazycall(model,where):
            if where == GRB.Callback.MIPSOL:
                # print(f"\n MIPSOL")
                model._voj = model.cbGetSolution(model._vars)
                model._yoj = model.cbGetSolution(model._varsy)
                knockset =  [i for i,y in enumerate(model._yoj) if model._yoj[i] < 1e-6]

                if len(knockset) != K:
                    return
                cur_obj = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
                cur_bd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
                model._vi, inner_status = inner(model._inner, model._yoj)

                if inner_status != GRB.OPTIMAL:

                    model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
                    return
                
                else:
                    vi_biom_val = model._vi[network.biomass]
                    vi_chem_val = model._vi[network.chemical]
                    knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.KO)
                    ki = (i for i in combinations(knockset_inner,K))
                    
                    if model._pbnd - model._vi[network.chemical] >= -1e-6: 

                        model.cbLazy(sum(model._varsy[j] for j in knockset) >= 1)
                        return
                    
                    elif (abs(model._vi[network.biomass] - model._voj[network.biomass]) > 1e-6):
                        model.cbLazy(vi_biom_val <= model._vars[network.biomass] + (math.ceil(model._vi[network.biomass]*10)/10) *(sum(model._varsy[f] for f in knockset)))
                        for comb in ki:
                            # print(f"{' '*4}{vi_biom_val} <= v[b] +{math.ceil(model._vi[network.biomass]*10)/10} *(sum{['y[%d]'%g for g in comb]})")
                            model.cbLazy(vi_biom_val <= model._vars[network.biomass] +
                                (math.ceil(model._vi[network.biomass]*10)/10) *(sum(model._varsy[f] for f in comb)))
                        return
                    else:

                        model._pbnd = cur_obj
                        return
    # =============================================================================================================================                        
        m = gp.Model("Optimistic")

        m.Params.OptimalityTol = network.infeas
        m.Params.IntFeasTol = network.infeas
        m.Params.FeasibilityTol = network.infeas
        m.Params.Presolve = 0
        m.Params.PreCrush = 1

        cbv = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='cbv')
        cby = m.addVars(network.M,vtype=GRB.BINARY,name='cby')
        cbvs = [cbv[i] for i in network.M]

        m.setObjective(1*cbv[network.chemical],GRB.MAXIMIZE)
        
        m.addMConstr(network.S,cbvs,'=',network.b,name='Stoi')

        m.addConstr(cbv[network.biomass] >= ominprod, name='target')

        m.addConstrs((olb[j]*cby[j] <= cbv[j] for j in network.M),name='LB')
        m.addConstrs((cbv[j] <= oub[j]*cby[j] for j in network.M),name='UB')
        
        if network.KO is not None:
            m.addConstr(sum(1-cby[j] for j in network.KO) == K, name='knapsack')
            m.addConstrs((cby[j] == 1 for j in network.M if j not in network.KO),name='Essen')
            # return to y[i] =1 for i essentials i not in KO
        elif network.KO is None:
            m.addConstr(sum(1-cby[j] for j in network.M) == K, name='knapsack')


        imodel = gp.Model()
        fv = imodel.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='fv')
        fvs = [fv[i] for i in network.M]
        imodel.params.LogToConsole = 0
        
        imodel.setObjective(2000*fv[network.biomass] + fv[network.chemical], GRB.MAXIMIZE)
        imodel.addMConstr(network.S,fvs,'=',network.b,name='Stoi')
        imodel.addConstr(fv[network.biomass] >= ominprod, name='target2')
        imodel.update()
        
        m._inner = imodel.copy()
        m._vars = cbv
        m._varsy = cby
        m.Params.lazyConstraints = 1
        m._pbnd = -1000
        m._cbcnt = 0
        m._sv = [0 for i in network.M]
        m._sy = []
        m._vi = None 
    
        m.Params.TimeLimit = network.time_limit
        if not self.log: m.Params.OutputFlag = 0
        m.optimize(lazycall)
        cb_time = m.Runtime

        if m.status == GRB.OPTIMAL:
            yss = [m.getVarByName('cby[%d]'%j).x for j in network.M]
            vss = [m.getVarByName('cbv[%d]'%j).x for j in network.M]
            del_strat_cb = [network.Rxn[i] for i in network.M if yss[i] < .5]
            soltype = 'Optimal'
        elif m.status == GRB.TIME_LIMIT:
            yss = [m.getVarByName('cby[%d]'%j).x for j in network.M]
            vss = [m.getVarByName('cbv[%d]'%j).x for j in network.M]
            del_strat_cb = [network.Rxn[i] for i in network.M if yss[i] < .5]
            soltype = 'Time_Limit'

        elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
            yss = [0 for i in network.M]
            vss = [2000 if i in [network.biomass,network.chemical] else 0 for i in network.M]
            del_strat_cb = [_ for _ in network.M]
            soltype = 'Infeasible'
        
        rys = self.__ysheuristic(yss)
        ko_index = [network.Rxn.index(i) for i in del_strat_cb]

        self.r_o = MethodResult(network.Name,del_strat_cb,vss,yss,ko_index,cb_time,soltype,'O')
        return self.r_o

    def CB_P(self,network:M_Network=None,K:int=None) -> MethodResult:
        
        plb = copy.deepcopy(network.LB)
        pub = copy.deepcopy(network.UB)
        pminprod = copy.deepcopy(network.minprod)
        plb[network.biomass] = pminprod

        def inner(imodel, yoj:Vector):

            imodel.setAttr('LB',imodel.getVars(),[plb[j]*yoj[j] for j in network.M])
            imodel.setAttr('UB',imodel.getVars(),[pub[j]*yoj[j] for j in network.M])
            imodel.Params.OptimalityTol = network.infeas
            imodel.Params.IntFeasTol = network.infeas
            imodel.Params.FeasibilityTol = network.infeas
            imodel.Params.Presolve = 0
            imodel.optimize()
            status = imodel.status
            if status == GRB.OPTIMAL:
                vij = [imodel.getVarByName('fv[%s]'%a).x for a in network.M] # rounded inner values
            elif status in (GRB.INFEASIBLE, GRB.UNBOUNDED, GRB.INF_OR_UNBD):
                vij = [2000 if i in [network.biomass,network.chemical] else yoj[i] for i in network.M]
            return vij,status

        def lazycall(model,where):
            
            if where == GRB.Callback.MIPSOL:
                model._voj = model.cbGetSolution(model._vars)
                model._yoj = model.cbGetSolution(model._varsy)
                knockset =  [i for i,y in enumerate(model._yoj) if model._yoj[i] < 1e-6]
                
                if len(knockset) != K:
                    return
                
                cur_obj = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
                cur_bd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
                if cur_bd > 100:
                    cur_bd = 100
                    
                model._vi, inner_status = inner(model._inner, model._yoj)

                if inner_status != GRB.OPTIMAL:
                    model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
                    return

                else:
                    vi_biom_val = model._vi[network.biomass]
                    vi_chem_val = model._vi[network.chemical]
                    
                    if network.KO is not None:
                        knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.KO)
                    else:
                        knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.M)
                    
                    ki = (i for i in combinations(knockset_inner,K))

                    flag = True
                    if abs(model._vi[network.biomass] - model._voj[network.biomass]) > 1e-6:
                            flag = False
                            for comb in ki:
                                model.cbLazy(vi_biom_val <= model._vars[network.biomass] +
                                                    (math.ceil(model._vi[network.biomass]*10)/10) *(sum(model._varsy[f] for f in comb)))
                    
                    if cur_bd + vi_chem_val + model._voj[network.chemical] == 0:
                        model.cbLazy(sum(model._varsy[j] for j in knockset)>=1)
                    
                    elif (cur_obj - vi_chem_val < 1e-6) and (flag):
 
                        model._pbnd = cur_obj
                        return
                    
                    elif model._pbnd < vi_chem_val: 
       
                        model.cbLazy(vi_chem_val >= model._vars[network.chemical] - math.ceil(cur_bd)*(sum(model._varsy[g] for g in knockset)))
                        if vi_chem_val > 1e-3:
                            model._sv = model._vi
                            model._sy = model._yoj

                            model.cbSetSolution(model._vars, model._sv)
                            model.cbSetSolution(model._varsy, model._sy)
                            model.cbUseSolution()
                        return

                    elif model._pbnd > vi_chem_val:

                        model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
                        return
    # =============================================================================================================================                        
        p = gp.Model('Pessimistic')

        p.Params.OptimalityTol = network.infeas
        p.Params.IntFeasTol = network.infeas
        p.Params.FeasibilityTol = network.infeas
        p.Params.Presolve = 0
        p.Params.PreCrush = 1

        cbv = p.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='cbv')
        cby = p.addVars(network.M,vtype=GRB.BINARY,name='cby')
        cbvs = [cbv[i] for i in network.M]
        
        p.setObjective(1*cbv[network.chemical],GRB.MAXIMIZE)
        
        p.addMConstr(network.S,cbvs,'=',network.b,name='Stoi')

        p.addConstr(cbv[network.biomass] >= pminprod, name='target')

        p.addConstrs((plb[j]*cby[j] <= cbv[j] for j in network.M),name='LB')
        p.addConstrs((cbv[j] <= pub[j]*cby[j] for j in network.M),name='UB')
        
        if network.KO is not None:
            p.addConstr(sum(1-cby[j] for j in network.KO) == K, name='knapsack')
            # p.addConstrs((cby[j] for j in network.M) == len(network.M)-K,name='Essen')
            p.addConstrs((cby[j]==1 for j in network.M if j not in network.KO),name='Essen')
        elif network.KO is None:
            p.addConstr(sum(1-cby[j] for j in network.M) == K, name='knapsack')


        imodel = gp.Model()
        fv = imodel.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='fv')
        fvs = [fv[i] for i in network.M]
        imodel.params.LogToConsole = 0
        
        # Inner model Objective
        imodel.setObjective(2000*fv[network.biomass] - fv[network.chemical], GRB.MAXIMIZE)
        
        imodel.addMConstr(network.S,fvs,'=',network.b,name='Stoi')
        
        imodel.addConstr(fv[network.biomass] >= pminprod, name='target2')

        imodel.update()
        
        p._inner = imodel.copy()
        p._vars = cbv
        p._varsy = cby
        p.Params.lazyConstraints = 1
        p._pbnd = -1000
        p._cbcnt = 0
        p._sv = [0 for i in network.M]
        p._sy = []
        p._vi = None 

        p.Params.TimeLimit = network.time_limit
        if not self.log: p.Params.OutputFlag = 0
        p.optimize(lazycall)
        
        pcb_time = p.Runtime
        
        if p.status == GRB.OPTIMAL:
            yss = [p.getVarByName('cby[%d]'%j).x for j in network.M]
            vss = [p.getVarByName('cbv[%d]'%j).x for j in network.M]
            del_strat_cb = [network.Rxn[i] for i in network.M if yss[i] < .5]
            soltype = 'Optimal'
        elif p.status == GRB.TIME_LIMIT:
            yss = [p.getVarByName('cby[%d]'%j).x for j in network.M]
            vss = [p.getVarByName('cbv[%d]'%j).x for j in network.M]
            del_strat_cb = [network.Rxn[i] for i in network.M if yss[i] < .5]
            soltype = 'Time_Limit'

        elif p.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
            yss = [0 for i in network.M]
            vss = [2000 if i in [network.biomass,network.chemical] else 0 for i in network.M]
            del_strat_cb = [_ for _ in network.M]
            soltype = 'Infeasible'

        rys = self.__ysheuristic(yss)

        ko_index = [network.Rxn.index(i) for i in del_strat_cb]
        self.r_p = MethodResult(network.Name,del_strat_cb,vss,yss,ko_index,pcb_time,soltype,'P')
        return self.r_p

    def FBA_check(self,network:M_Network,solution:MethodResult=None,obj_v:str=None,c_params:str=None) -> Check:
        clb = copy.deepcopy(network.LB)
        cub = copy.deepcopy(network.UB)
        
        clb[network.biomass] = network.minprod
        vs = copy.deepcopy(solution.Vs)
        cys = copy.deepcopy(solution.Ys)

        # print(f"{len(cys)} -> {len(clb)} -> {len(cub)}")
        if obj_v == 'biomass':
            f_objective = network.biomass
        elif obj_v == 'chemical':
            f_objective = network.chemical
        else:
            raise Exception("No Inner Objective ('biomass' or 'chemical')")
        
        m = gp.Model()

        v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')
        v_s = [v[i] for i in network.M]
        m.setObjective((1*v[f_objective]),GRB.MAXIMIZE)
       
        # print(len(m.getVars()))
        # print(len([clb[j] * cys[j] for j in network.M]))    
        # out_s = [0 if abs(i)<1e-6 else (i,j) for i,j in enumerate([network.S.dot(vs) - network.b])]
        # print(network.S.dot(vs) - network.b)
        # print(out_s)
        m.addMConstr(network.S,v_s,'=',network.b,name='Stoi')
        # m.addConstr(v[network.biomass] >= network.minprod,name='minprod')
        m.update()
        if c_params == 'both':
            m.setAttr('LB',m.getVars(),[clb[j]*cys[j] for j in network.M])
            m.setAttr('UB',m.getVars(),[cub[j]*cys[j] for j in network.M])

            # m.addConstrs((clb[j]*cys[j]<= v[j] for j in network.M),name='lb')

            # m.addConstrs((cub[j]*cys[j] >= v[j] for j in network.M),name='ub')

            m.addConstrs((v[j] == vs[j] for j in network.M),name='hards_vs')
        
        elif c_params == 'vs':
            m.setAttr('LB',m.getVars(),[clb[j] for j in network.M])
            m.setAttr('UB',m.getVars(),[cub[j] for j in network.M])

            m.addConstrs((v[j] == vs[j] for j in network.M),name='hards_vs')

        elif c_params == 'ys':
            m.setAttr('LB',m.getVars(),[clb[j]*cys[j] for j in network.M])
            m.setAttr('UB',m.getVars(),[cub[j]*cys[j] for j in network.M])
            # for j in network.M:
            #     if (lb[j]*ys[j] - vs[j] > 1e-6): print ("errorlb", j, lb[j], ys[j], vs[j])
            #     if (ub[j]*ys[j] - vs[j] < -1e-6): print ("errorub", j, ub[j], ys[j], vs[j])
            
        
        else:
            raise Exception("Choose the type if inner check [vs,ys,both]")
        
        m.Params.OptimalityTol = network.infeas
        m.Params.IntFeasTol = network.infeas
        m.Params.FeasibilityTol = network.infeas
        m.Params.Presolve = 0
        if not self.log: m.Params.OutputFlag = 0
        m.optimize()

        if m.status == GRB.OPTIMAL:
            vss =  [m.getVarByName('v[%s]'%a).x for a in network.M]
            soltype = 'optimal'
        elif m.status in (GRB.INFEASIBLE,GRB.UNBOUNDED, GRB.INF_OR_UNBD):
            vss = [network.BM if i in[network.biomass,network.chemical] else 0 for i in network.M]
            soltype = 'Infeasible'
        self.r_c = Check(vss[network.biomass],vss[network.chemical],soltype,obj_v,c_params)

        return self.r_c
    
    def __ysheuristic(self,ys:Vector=None) -> Vector:
        rounded_ys = [0]*len(ys)
        for i,y in enumerate(ys):
            if y < .5:
                rounded_ys[i] = 0
            else:
                rounded_ys[1] = 1
        return rounded_ys

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
Result_cb = namedtuple('Result_cb',['MetNet','Strategy','Ys','Vs','Vij','KO_index','Time','Soltype','Method'])


def CB_P(network:M_Network=None, k:Ks=None,log:bool=None,speed:bool=False,threads:bool=False,lp:bool=False) -> Result_cb:
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
        print(f'\n**** Solving Callbacks k={k} ****')
        print(f'# Variables (reactions in the network): {len(network.M)}')
        print('Current Infeasibility:',network.infeas,sep=' -> ')
        if network.KO is not None:
            print('KO set: ',len(network.KO), ' reactions')
        else:
            print(f"KO set: all reactions in KO -> {len(network.M)}")
        print(f"MN: {network.Name}")
        print(f"Chemical: {network.Rxn[network.chemical]} -> {network.chemical}")
        print(f"Growth: {network.Rxn[network.biomass]} -> {network.biomass}")
        print(f"Target: {network.target}")
        print(f"Minprod: {network.FVA[network.biomass]}")
        print(f"FBA [b]: {network.FBA[network.biomass]}")
        print(f"-- Pessimistic Approach -- \n")

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
                print(f'\nMIPSOL')
                model._voj = model.cbGetSolution(model._vars)
                model._yoj = model.cbGetSolution(model._varsy)
                knockset =  [i for i,y in enumerate(model._yoj) if model._yoj[i] < 1e-6]
                
                if len(knockset) != k:
                    return
                
                cur_obj = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
                cur_bd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)

                if cur_bd > 100:
                    cur_bd = 100
                print(f"KO set: {knockset} -> {[network.Rxn[i] for i in knockset]}")
                print(f"Vo[b]: {model._voj[network.biomass]}")
                print(f"Vo[c]: {model._voj[network.chemical]}")
                model._vi, inner_status = inner(model._inner, model._yoj)

    # ============================ Checking Inner Optimality Status ===================================
                
                print(f"MIPSOL Incbnt = {cur_obj}")
                print(f"MIPSOL Bound = {cur_bd}")
                print(f"Curnt Pbnd = {model._pbnd}")
                print(f"Vi[c]: {model._vi[network.chemical]}")
                print(f"Vi[b]: {model._vi[network.biomass]}")
                print(f"Sum: {cur_obj+model._voj[network.chemical]}")
                print(f">> Algorithm response: \n")

                if inner_status != GRB.OPTIMAL:
                    print(f"{' '*3}feasibility cuts inner not optinal MIPSOL \n")
                    print(f"{' '*4} (sum{['y[%d]'%g for g in knockset]}) >= 1")
                    model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
                    if lp:
                        ysum = [model._varsy[j] for j in knockset]
                        expr = sum(ysum)
                        sense = '>='
                        rhs = 1
                        lazycts.append((expr,sense,rhs))
                    return
         
                else:
                    vi_biom_val = model._vi[network.biomass]
                    vi_chem_val = model._vi[network.chemical]
                    
                    if network.KO is not None:
                        knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.KO)
                    else:
                        knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.M)
                    
                    
                    ki = (i for i in combinations(knockset_inner,k))

                 
                    flag = True
                    if abs(model._vi[network.biomass] - model._voj[network.biomass]) > 1e-6:
                            print(f"{' '*3}|vi[b] - vo[b]| > 1e-6 \n")
                            flag = False
                            for comb in ki:
                                model.cbLazy(vi_biom_val <= model._vars[network.biomass] +
                                                    (math.ceil(model._vi[network.biomass]*10)/10) *(sum(model._varsy[f] for f in comb)))
                                print(f"{' '*4}{vi_biom_val} <= v[b] +{math.ceil(model._vi[network.biomass]*10)/10} *(sum{['y[%d]'%g for g in comb]})")
                                if lp:
                                    sense = '>='
                                    rhs = vi_biom_val
                                    ysum = [model._varsy[j] for j in comb]
                                    expr = model._vars[network.biomass] + (math.ceil(model._vi[network.biomass]*10)/10) *sum(ysum)
                                    lazycts.append((expr,sense,rhs))
                    
                    print(f"This is a new print")  

# New Lines to avoid zeroes
                    if cur_obj + vi_chem_val + model._voj[network.chemical] == 0:
                        print(f"{' '*3}Cur Incubent = 0 & Vo = 0")
                        print(f"{' '*4} (sum{['y[%d]'%g for g in knockset]}) >= 1")
                        model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)

                    elif (cur_obj - vi_chem_val < 1e-6) and (flag):
                        print(f"{' '*3}curobj - vi[c] < 1e-6 & flag \n")
                        print(f"{' '*4}Update: pbdn = {vi_chem_val}")
                        model._pbnd = cur_obj
                        return
                    
                    elif model._pbnd < vi_chem_val: 
                        print(f"\n{' '*3}pbnd < vi[c] \n")
                        print(f"{' '*4}{vi_chem_val} >= v[c] - {math.ceil(cur_bd)}*(sum{['y[%d]'%g for g in knockset]})")
                       
                        
                        model.cbLazy(vi_chem_val >= model._vars[network.chemical] - math.ceil(cur_bd)*(sum(model._varsy[g] for g in knockset)))
                        # Updating the model variables
                        
                        if vi_chem_val > 1e-3: 
                            model._sv = model._vi
                            model._sy = model._yoj
                            # Setting solution
                            model.cbSetSolution(model._vars, model._sv)
                            model.cbSetSolution(model._varsy, model._sy)
                            model.cbUseSolution()
                            print(f"{' '*4}Set Solution")
                        if lp:
                            sense = '<='
                            rhs = vi_chem_val
                            ysum = [model._varsy[g] for g in knockset]
                            expr =  model._vars[network.chemical] - cur_obj*(sum(ysum))
                            lazycts.append((expr,sense,rhs))                       
                        return
                    
                    elif model._pbnd > vi_chem_val:
                        print(f"\n{' '*3}pbnd > vi_chem \n")
                        print(f"{' '*4}(sum{['y[%d]'%g for g in knockset]}) >= 1")
                        model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
                        if lp:
                            sense = '>='
                            rhs = 1
                            ysum = [model._varsy[j] for j in knockset]
                            expr = sum(ysum)
                            lazycts.append((expr,sense,rhs))
                        return
    # =============================================================================================================================                        
    # ============================================== MIPNODE =====================================       
            # elif where == GRB.Callback.MIP:
            #     # print(f'\nMIP')
            #     mip_objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
            #     mip_objbst = model.cbGet(GRB.Callback.MIP_OBJBST)

            #     # print(f"MIP_objbnd -> {mip_objbnd}")
            #     # print(f"MIP_objbst -> {mip_objbst}")


            # elif where == GRB.Callback.MIPNODE:
            elif 1 ==0:
                status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
                if status == GRB.OPTIMAL:
                    # print(f'\nMIPNODE')
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
                        # print(f"Viche: {model._vi[network.chemical]:.6}")
                        # print(f"Algorithm response: \n")
                        if inner_status != GRB.OPTIMAL:
                            # print(f"{' '*3}Feasibility cuts - inner not optimal")
                            model.cbCut(sum([model._varsy[f] for f in knock]) >= 1)
                            if lp:
                                sense = '>='
                                rhs = 1
                                ysum = [model._varsy[j] for j in knock]
                                expr = sum(ysum)
                                lazycts.append((expr,sense,rhs))
                        else:

                            if (model._vi[network.chemical] >= model._pbnd) and (model._vi[network.chemical] > model._sv[network.chemical]): # added condition to set solution always better
                                model._sv = model._vi
                                model._sy = model._ryo
                                model._pbnd = model._vi[network.chemical]
                                # print(f"{' '*3} vi[chemical] > pbnd & vi[chemical] > sv[chemical]")
                                # print(f"{' '*3} Set Solution ")
                                # print(f"{' '*3} pbnd = vi[chemical] ")
                                model.cbSetSolution(model._vars, model._sv)
                                model.cbSetSolution(model._varsy, model._sy)
                            
                            elif model._vi[network.chemical] < model._pbnd:
                                # print(f"{' '*3} Optimality cuts - vi[chemical] < pbnd")
                                model.cbCut(sum([model._varsy[f] for f in knock]) >= 1)
                                if lp:
                                    sense = '>='
                                    rhs = 1
                                    ysum = [model._varsy[f] for f in knock]
                                    expr = sum(ysum)
                                    lazycts.append((expr,sense,rhs))
                            else:
                                model.cbSetSolution(model._vars, model._sv)
                                model.cbSetSolution(model._varsy, model._sy)
        
       
        m = gp.Model("Pessimistic")

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
        # m.addConstrs((gp.quicksum(network.S[i,j]*cbv[j] for j in network.M) == 0 for i in network.N),name='Stoichiometry')

        m.addConstr(cbv[network.biomass] >= pminprod, name='target')

        m.addConstrs((plb[j]*cby[j] <= cbv[j] for j in network.M),name='LB')
        m.addConstrs((cbv[j] <= pub[j]*cby[j] for j in network.M),name='UB')
        
        if network.KO is not None:
            m.addConstr(sum(1-cby[j] for j in network.KO) == k, name='knapsack')
            # m.addConstr(sum(cby[j]for j in network.M) == len(network.M)-k,name='Essen')
            m.addConstrs((cby[j]==1 for j in network.M if j not in network.KO),name='essen')
        elif network.KO is None:
            m.addConstr(sum(1-cby[j] for j in network.M) == k, name='knapsack')


        imodel = gp.Model()
        fv = imodel.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='fv')
        fvs = [fv[i] for i in network.M]
        imodel.params.LogToConsole = 0
        
        # Inner model Objective
        imodel.setObjective(2000*fv[network.biomass] - fv[network.chemical], GRB.MAXIMIZE)
        
        imodel.addMConstr(network.S,fvs,'=',network.b,name='Stoi')
        # imodel.addConstrs((gp.quicksum(network.S[i,j]*fv[j] for j in network.M) == 0 for i in network.N),name='S2')
        
        imodel.addConstr(fv[network.biomass] >= pminprod, name='target2')

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
        lazycts = []
        m._vi = None 

        if not log: m.Params.OutputFlag = 0
        if speed: m.Params.NodefileStart = 0.5
        if threads: m.Params.Threads = 6

        m.Params.TimeLimit = network.time_limit

        m.optimize(lazycall)
        
        cb_time = m.Runtime
        vij = m._vi
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
            ys = [0 for i in network.M]
            vs = [2000 for i in network.M]
            del_strat_cb = ['all']
            soltype = 'Infeasible'
        ko_index = [network.Rxn.index(i) for i in del_strat_cb]
        if lp:
            for l in range(len(lazycts)):
                if lazycts[l][1] == '=':
                    c = m.addConstr(lazycts[l][0] == lazycts[l][2],'lzy%d'%l)
                    c.Lazy = 1
                elif lazycts[l][1] == '>=':
                    c = m.addConstr(lazycts[l][0] >= lazycts[l][2],'lzy%d'%l)
                    c.Lazy = 1
                elif lazycts[l][1] == '<=':
                    c = m.addConstr(lazycts[l][0] <= lazycts[l][2],'lzy%d'%l)
                    c.Lazy = 1

            filename = f"pess_log_k{k}_{network.Name}.lp"
            m.write(filename)

        return Result_cb(network.Name,del_strat_cb,ys,vs,vij,ko_index,cb_time,soltype,'P')


##### Pessimistic Callback Function with changes from Monday's meeting no flag added ####
# def CB_P_t(network:M_Network=None, k:Ks=None,log:bool=None,speed:bool=False,threads:bool=False,lp:bool=False) -> Result_cb:
#     '''
#         cb = CB_solve_2_NOP(network=network,k=k,log=True,speed=False,threads=False)
           
#         cb.MetNet    = Metabolic Network's name
#         cb.Stragtegy = List of rxn to knockout
#         cb.Ys        = Binary solution as a vector
#         cb.Vs        = Optimal bilevel flows
#         cb.Vij       = Flows in the inner problem
#         cb.Time      = Solving time 
#         cb.Soltype   = Type of solution [optimal, timelimit , infeasible]
#         cb.Method    = Solving Method - set to CBP (Callbacks-Pessimistic)

#         '''
#     print(f'**** Solving Callbacks k={k} ****')
#     print(f'# Variables (reactions in the network): {len(network.M)}')
#     print('Current Infeasibility:',network.infeas,sep=' -> ')
#     print('KO set: ',len(network.KO), ' reactions')
#     print(f"MN: {network.Name}")
#     print(f"-- Pessimistic Approach --")

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
#         imodel.Params.Presolve = 0
#         imodel.optimize()
#         status = imodel.status
#         if status == GRB.OPTIMAL:
#             vij = [imodel.getVarByName('fv[%s]'%a).x for a in network.M] # rounded inner values
#         elif status in (GRB.INFEASIBLE, GRB.UNBOUNDED, GRB.INF_OR_UNBD):
#             vij = [2000 if i == network.biomass else yoj[i] for i in network.M]
#         return vij,status

#     def lazycall(model,where):
#         # cur_obj1 = model.cbGet(GRB.Callback.MIP_OBJBST)
#         if where == GRB.Callback.MIPSOL:
#             print(f'\nMIPSOL')
#             model._voj = model.cbGetSolution(model._vars)
#             model._yoj = model.cbGetSolution(model._varsy)
#             knockset =  [i for i,y in enumerate(model._yoj) if model._yoj[i] < 1e-6]
#             ysum = [model._varsy[g] for g in knockset]

#             if len(knockset) != k:
#                 return
#             cur_obj = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
#             # cur_obj1 = model.cbGet(GRB.Callback.MIP_OBJBST)
#             cur_bd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
#             print(f"Vbio: {model._voj[network.biomass]}")
#             print(f"Vche: {model._voj[network.chemical]}")
#             model._vi, inner_status = inner(model._inner, model._yoj)

# # ============================ Checking Inner Optimality Status ===================================
            
#             print(f"Objective = {cur_obj}")
#             # print(f"Objective = {cur_obj1}")
#             print(f"Best Bound = {cur_bd}")
#             print(f"Curnt Pbnd = {model._pbnd}")
#             print(f"Viche: {model._vi[network.chemical]}")
#             print(f"KO Set{knockset}")
#             print(f"Algorithm response: \n")
#             if inner_status != GRB.OPTIMAL:
#                 print(f"{' '*3}feasibility cuts inner not optinal MIPSOL")
#                 model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
#                 ysum = [model._varsy[j] for j in knockset]
#                 expr = sum(ysum)
#                 sense = '>='
#                 rhs = 1
#                 lazycts.append((expr,sense,rhs))

#             else:
#                 vi_biom_val = model._vi[network.biomass]
#                 vi_chem_val = model._vi[network.chemical]
#                 knockset_inner = (i for i,y in enumerate(model._vi) if abs(model._vi[i]) < 1e-6 and i in network.KO)
#                 ki = (i for i in combinations(knockset_inner,k))

#                 # bestknownchem = cur_obj
                
#                 if model._pbnd < vi_chem_val: 
#                     print(f"{' '*3}Optimality cuts, pbnd < vi[chemical]")
#                     print(f"{vi_chem_val} >= vchem + {cur_obj}*(sum{['y[%d]'%g for g in knockset]})")
#                     print(f"Update: pbdn = {vi_chem_val}")
#                     model.cbLazy(vi_chem_val >= model._vars[network.chemical] - cur_obj*(sum(model._varsy[i] for i in knockset)))
#                     model._pbnd = vi_chem_val
#                     model._sv = model._vi
#                     model._sy = model._yoj

#                     model.cbSetSolution(model._vars, model._sv)
#                     model.cbSetSolution(model._varsy, model._sy)
#                     model.cbUseSolution()

#                     sense = '<='
#                     rhs = vi_chem_val
#                     expr =  model._vars[network.chemical] - cur_obj*(sum(ysum))
#                     lazycts.append((expr,sense,rhs))

#                     if abs(model._vi[network.biomass] - model._voj[network.biomass]) > 1e-6:
#                         for comb in ki:
#                             model.cbLazy(vi_biom_val <= model._vars[network.biomass] +
#                             (math.ceil(vij[network.biomass]*10)/10) *(sum(model._varsy[f] for f in comb)))
#                             print(f"{vi_biom_val} <= vbiom +{math.ceil(vij[network.biomass]*10)/10} *(sum{['y[%d]'%g for g in comb]})")
#                             sense = '>='
#                             rhs = vi_biom_val
#                             ysumc = [model._varsy[j] for j in comb]
#                             expr = model._vars[network.biomass] + (math.ceil(vij[network.biomass]*10)/10) *sum(ysumc)
#                             lazycts.append((expr,sense,rhs))

                
#                 elif model._pbnd >= vi_chem_val:
#                     print(f"{' '*3}PBND > vi_chem")
#                     model.cbLazy(sum(model._varsy[j] for j in knockset) >=1)
#                     expr = sum(ysum)
#                     sense = '>='
#                     rhs = 1
#                     lazycts.append((expr,sense,rhs))




# # =============================================================================================================================                        

#         # elif where == GRB.Callback.MIPNODE:
#         elif 1==0:
#             status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
#             if status == GRB.OPTIMAL:
#                 print(f'\nMIPNODE')
#                 mipobj = model.cbGet(GRB.Callback.MIPNODE_OBJBST)
#                 mipbnd = model.cbGet(GRB.Callback.MIPNODE_OBJBND)
#                 print(f"Obj= {mipobj}")
#                 print(f"Bnd= {mipbnd}")
#                 # print(f"Curnt set solution v: {model._sv[network.chemical]:.6f}")
#                 # print(f"Curnt set solution len y: {len(model._sy)}")
#                 print(f"Curnt PBND: {model._pbnd}")
#                 model._cbcnt += 1
#                 if model._cbcnt % 10 == 0:
#                     return
#                 # print(f"CB counts : {model._cbcnt}")
#                 model._voj = model.cbGetNodeRel(model._vars)
#                 model._ryo = model.cbGetNodeRel(model._varsy) #  
#                 for i,y in enumerate(model._ryo):
#                     if model._ryo[y] >= 0.8:
#                         model._ryo[y] = 1.0
#                     elif model._ryo[y] <= 0.2:
#                         model._ryo[y] = 0.0
#                     else:
#                         model._ryo[y] = 1.0
#                 knock = [i for i,y in enumerate(model._ryo) if model._ryo[i] < 1e-6]
#                 ysum = [model._varsy[g] for g in knock]
#                 if sum(model._ryo.values()) != len(model._ryo)-k:
#                     return
#                 else:
#                     model._vi, inner_status = inner(model._inner,model._ryo)
#                     # print(f"Optimality Code - {inner_status}")
#                     print(f"Viche: {model._vi[network.chemical]:.6}")
#                     print(f"Algorithm response: \n")
#                     if inner_status != GRB.OPTIMAL:
#                         print(f"{' '*3}Feasibility cuts - inner not optimal")
#                         model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
#                         expr = sum(ysum)
#                         sense = '>='
#                         rhs = 1
#                         lazycts.append((expr,sense,rhs))

#                     elif (model._vi[network.chemical] > model._pbnd) and (model._vi[network.chemical] > model._sv[network.chemical]): # added condition to set solution always better
#                         model._sv = model._vi
#                         model._sy = model._ryo
#                         model._pbnd = model._vi[network.chemical]
#                         print(f"{' '*3} Set Solution ")
#                         # print(f"biomas {model._vi[network.biomass]:.6f}")
#                         # print(f"chemical {model._vi[network.chemical]:.6f}")
#                         model.cbSetSolution(model._vars, model._sv)
#                         model.cbSetSolution(model._varsy, model._sy)
#                         # new pbnd = vichem
                    
#                     else:
#                         print(f"{' '*3}Optimality cuts - vi not larger or equal than pnbd")
#                         model.cbLazy(sum([model._varsy[f] for f in knock]) >= 1)
#                         expr = sum(ysum)
#                         sense = '>='
#                         rhs = 1
#                         lazycts.append((expr,sense,rhs))
    
#     m = gp.Model()
#     cbv = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='cbv')
#     cby = m.addVars(network.M,vtype=GRB.BINARY,name='cby')
#     cbvs = [cbv[i] for i in network.M]

#     m.setObjective(1*cbv[network.chemical],GRB.MAXIMIZE)
    
#     m.addMConstr(network.S,cbvs,'=',network.b,name='Stoi')
#     # m.addConstrs((gp.quicksum(network.S[i,j]*cbv[j] for j in network.M) == 0 for i in network.N),name='Stoichiometry')

#     m.addConstr(cbv[network.biomass] >= minprod, name='target')

#     m.addConstrs((lb[j]*cby[j] <= cbv[j] for j in network.M),name='LB')
#     m.addConstrs((cbv[j] <= ub[j]*cby[j] for j in network.M),name='UB')
    
#     m.addConstr(sum(1-cby[j] for j in network.KO) == k, name='knapsack')
#     m.addConstrs((cby[j] == 1 for j in network.M if j not in network.KO),name='Essen')


#     imodel = gp.Model()
#     fv = imodel.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='fv')
#     imodel.params.LogToConsole = 0
#     fvs = [fv[i] for i in network.M]
#     # Inner model Objective
#     imodel.setObjective(2000*fv[network.biomass] - fv[network.chemical], GRB.MAXIMIZE)
    
#     imodel.addMConstr(network.S,fvs,'=',network.b,name='Stoi')
#     # imodel.addConstrs((gp.quicksum(network.S[i,j]*fv[j] for j in network.M) == 0 for i in network.N),name='S2')
    
#     imodel.addConstr(fv[network.biomass] >= minprod, name='target2')

#     imodel.update()
    
#     m._inner = imodel.copy()
#     # m._innerv = fv
#     m._vars = cbv
#     m._varsy = cby
#     m.Params.lazyConstraints = 1
#     m._pbnd = -1000
#     m._cbcnt = 0
#     m._sv = [0 for i in network.M]
#     m._sy = []
#     lazycts = []
#     m.Params.OptimalityTol = network.infeas
#     m.Params.IntFeasTol = network.infeas
#     m.Params.FeasibilityTol = network.infeas
#     m.Params.Presolve = 0
#     if not log: m.Params.OutputFlag = 0
#     if speed: m.Params.NodefileStart = 0.5
#     if threads: m.Params.Threads = 6

#     m.Params.TimeLimit = network.time_limit

#     m.optimize(lazycall)
    
#     cb_time = m.Runtime

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
#         ys = [0 for i in network.M]
#         vs = [2000 for i in network.M]
#         del_strat_cb = ['all']
#         soltype = 'Infeasible'

#     if lp:
#         for l in range(len(lazycts)):
#             if lazycts[l][1] == '=':
#                 c = m.addConstr(lazycts[l][0] == lazycts[l][2],'lzy%d'%l)
#                 c.Lazy = 1
#             elif lazycts[l][1] == '>=':
#                 c = m.addConstr(lazycts[l][0] >= lazycts[l][2],'lzy%d'%l)
#                 c.Lazy = 1
#             elif lazycts[l][1] == '<=':
#                 c = m.addConstr(lazycts[l][0] <= lazycts[l][2],'lzy%d'%l)
#                 c.Lazy = 1

#         filename = f"pess_log_k{k}.lp"
#         m.write(filename)

#     return Result_cb(network.Name,del_strat_cb,ys,vs,vij,cb_time,soltype,'CBP')

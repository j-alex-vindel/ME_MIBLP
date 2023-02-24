import sys
import os

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from Ob_Met_Net_solmethods import MILP_solve_NOP,CB_sol_2_NOP, Inner_check_vs_ys_NOP
from Ob_Met_Net_ijr904 import MetNet_ijr904

import pandas as pd

metnet = MetNet_ijr904
metnet.target = .1
# =========================================== K=1 ===========================

m = MILP_solve_NOP(network=metnet,k=1,log=False)
cb = CB_sol_2_NOP(network=metnet,k=1,log=False)

r  = {
            'Method':['MILP','CB'],
            'Biomass':[m.Vs[metnet.biomass],cb.Vs[metnet.biomass]],
            'Chemical':[m.Vs[metnet.chemical],cb.Vs[metnet.chemical]],
            'Time':[m.Time,cb.Time],
            'Strategy':[m.Strategy,cb.Strategy]
}

file_name_r = f"../Results/EQ/Methods_K{sys.argv[0][-4]}_{metnet.Name}.csv"

isfile_r = os.path.exists(file_name_r)

if not isfile_r:
    rdf = pd.DataFrame.from_dict(r)
    rdf.to_csv(file_name_r)

# Check on the inner problem vs - MILP
m_vs_b = Inner_check_vs_ys_NOP(network= metnet, result_milp = m, criteria= 'vs', milp= True, cb = False, objective = 'biomass',log=False)
m_vs_c = Inner_check_vs_ys_NOP(network= metnet, result_milp = m, criteria= 'vs', milp= True, cb = False, objective = 'chemical',log=False)


# Check on the inner problem ys - MILP
m_ys_b = Inner_check_vs_ys_NOP(network= metnet, result_milp = m, criteria= 'ys', milp= True, cb = False, objective = 'biomass',log=False)
m_ys_c = Inner_check_vs_ys_NOP(network= metnet, result_milp = m, criteria= 'ys', milp= True, cb = False, objective = 'chemical',log=False)


# Check on the inner problem vs - CB

cb_vs_b = Inner_check_vs_ys_NOP(network= metnet, result_cb = cb,  criteria= 'vs', milp= False, cb = True, objective = 'biomass',log=False)
cb_vs_c = Inner_check_vs_ys_NOP(network= metnet, result_cb = cb,  criteria= 'vs', milp= False, cb = True, objective = 'chemical',log=False)

# Check on the inner problem ys - CB

cb_ys_b = Inner_check_vs_ys_NOP(network= metnet, result_cb = cb, criteria= 'ys', milp= False, cb = True, objective = 'biomass',log=False)
cb_ys_c = Inner_check_vs_ys_NOP(network= metnet, result_cb = cb, criteria= 'ys', milp= False, cb = True, objective = 'chemical',log=False)


r_inner = {'From':      ['CB'            ,'MILP'         ,'CB'            ,'MILP'         ,'CB'            ,'MILP'         ,'CB'            ,'MILP'],
           'Criteria':  ['vs'            ,'vs'           ,'ys'            ,'ys'           ,'vs'            ,'vs'           ,'ys'            ,'ys'],
           'FBA_Objct': ['B'             ,'B'            ,'B'             ,'B'            ,'C'             ,'C'            ,'C'             ,'C'],
           'V_biomass': [cb_vs_b.Biomass ,m_vs_b.Biomass ,cb_ys_b.Biomass ,m_ys_b.Biomass ,cb_vs_c.Biomass ,m_vs_c.Biomass ,cb_ys_c.Biomass ,m_ys_c.Biomass],
           'V_chemical':[cb_vs_b.Chemical,m_vs_b.Chemical,cb_ys_b.Chemical,m_ys_b.Chemical,cb_vs_c.Chemical,m_vs_c.Chemical,cb_ys_c.Chemical,m_ys_c.Chemical],
           'Soltype':   [cb_vs_b.Soltype ,m_vs_b.Soltype ,cb_ys_b.Soltype ,m_ys_b.Soltype ,cb_vs_c.Soltype ,m_vs_c.Soltype ,cb_ys_c.Soltype ,m_ys_c.Soltype],
           'Strategy':  [cb.Strategy     ,m.Strategy     ,cb.Strategy     ,m.Strategy     ,cb.Strategy     ,m.Strategy     ,cb.Strategy     ,m.Strategy    ],
           'K': [len(cb.Strategy)     ,len(m.Strategy)     ,len(cb.Strategy)     ,len(m.Strategy)     ,len(cb.Strategy)     ,len(m.Strategy)     ,len(cb.Strategy)     ,len(m.Strategy)    ]
    }


file_name_ri = f"../Results/EQ/IC_K{sys.argv[0][-4]}_{metnet.Name}.csv"
isfile_ri = os.path.exists(file_name_ri)

if not isfile_ri:
    ridf = pd.DataFrame.from_dict(r_inner)
    ridf.to_csv(file_name_ri)
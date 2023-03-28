import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from CB_Sol_P import CB_P,CB_P_t
from Ob_Met_Net_solmethods import Inner_check_vs_ys_NOP

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'B_iJO1366')))
from MN_ijo1366 import MN_ijo1366
import pandas as pd

metnet = MN_ijo1366




# cp = CB_P(network=metnet,k=2,log=True)


# ccc = Inner_check_vs_ys_NOP(network=metnet,result_cb=cp,criteria='ys',objective='chemical')

# ccb = Inner_check_vs_ys_NOP(network=metnet,result_cb=cp,criteria='ys',objective='biomass')

# r = {'V_index':[metnet.biomass,     metnet.biomass,       metnet.chemical,         metnet.chemical],
#      'Name':[metnet.Rxn[metnet.biomass],metnet.Rxn[metnet.biomass],metnet.Rxn[metnet.chemical],metnet.Rxn[metnet.chemical]],
#     'CB_P':[cp.Vs[metnet.biomass],  cp.Vs[metnet.biomass], cp.Vs[metnet.chemical], cp.Vs[metnet.chemical]   ],
#      'IOF':[ccb.OF,                 ccc.OF,                ccb.OF                , ccc.OF],
#      'IC':[ccb.Biomass,             ccc.Biomass,           ccb.Chemical,            ccc.Chemical]}



c1p = CB_P_t(network=metnet,k=1,log=True)

cc1c = Inner_check_vs_ys_NOP(network=metnet,result_cb=c1p,criteria='ys',objective='chemical')

cc1b = Inner_check_vs_ys_NOP(network=metnet,result_cb=c1p,criteria='ys',objective='biomass')

r1 = {'V_index':[metnet.biomass,     metnet.biomass,       metnet.chemical,         metnet.chemical],
      'Name':[metnet.Rxn[metnet.biomass],metnet.Rxn[metnet.biomass],metnet.Rxn[metnet.chemical],metnet.Rxn[metnet.chemical]],
    'CB_P':[c1p.Vs[metnet.biomass],  c1p.Vs[metnet.biomass], c1p.Vs[metnet.chemical], c1p.Vs[metnet.chemical]   ],
     'IOF':[cc1b.OF,                 cc1c.OF,                cc1b.OF                , cc1c.OF],
     'IC':[cc1b.Biomass,             cc1c.Biomass,           cc1b.Chemical,            cc1c.Chemical]}

# df = pd.DataFrame.from_dict(r)
# df.round(decimals=5)
# print(f" ")
# print(f">>Strategy {cp.Strategy}\n")
# print(df.to_markdown())
# print(f" ")


df1 = pd.DataFrame.from_dict(r1)
df1.round(decimals=5)
print(f" ")
print(f">>Strategy {c1p.Strategy}\n")
print(df1.to_markdown())
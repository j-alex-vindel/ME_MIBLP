import sys
import os
from Ob_Met_Net import Met_Net



def strainsele(strain:str=None) -> Met_Net:

    if strain == 'ijo':
        sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'B_iJO1366'))) 
        from MN_ijo1366 import MN_ijo1366
        met = MN_ijo1366
    elif strain == 'ijr':
        sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'B_iJR904'))) 
        from MN_ijr904 import MN_ijr904
        met = MN_ijr904
    elif strain == 'iaf':
        sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'B_iAF1260'))) 
        from MN_iaf1260 import MN_iaf1260
        met = MN_iaf1260
    return met


def strain_id(name:str=None)->int:
    if name == 'iJO':
        idstrain = 1000
    elif name == 'iJR':
        idstrain = 2000
    elif name == 'iAF':
        idstrain = 3000
    return idstrain

def method_id(method:str=None)->int:
    if method == 'O':
        meid = 100
    elif method == 'P':
        meid = 200
    elif method == 'M':
        meid = 300 
    return meid


if __name__=="__main__":
    metnet = strainsele('iaf')
    print(metnet.Name)
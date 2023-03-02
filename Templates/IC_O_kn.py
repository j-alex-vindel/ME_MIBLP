import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from support_functions import result_parser
from Ob_Met_Net_solmethods import Inner_check_vs_ys_NOP


file = f"../Results/XML/OPCB_iAF1260_k1.xml"

k1 = result_parser(file)

print(f"Met Net : {k1.MetNet}")

import os
import sys
import time
import pandas as pd
import datetime
import matplotlib.pyplot as plt
from itertools import product

def main(arg1=None,arg2=None,name:str=None,sleep:int=None):
    
    # os.system(f"python ../A_modules/{sys.argv[1]}.py")
    print(f">> Folder {name}")
    file = name[:3].lower()


    print(f">> Running -> ../B_{name}/{file} {arg1} {arg2}.py")
    if arg2 == None:
        os.system(f"python ../B_{name}/{file}.py {arg1}")
    else:   
        os.system(f"python ../B_{name}/{file}.py {arg1} {arg2}")

    print(f'>> Resting for {sleep}')
    time.sleep(sleep)
    print(f">> Ready to start over")
    
    print('Run Completed')


if __name__ =="__main__":
    strains = ['iAF1260']
    args1 = [1,2,3]
    args2 = [10,20,30,40,50,60,70,80,90]

    for strain in strains:
        for arg1 in args1:
            for arg2 in args2:
                main(arg1=arg1,arg2=arg2,name=strain,sleep=5)

  


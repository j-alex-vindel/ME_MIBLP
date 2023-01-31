import sys
import os
import time
import shutil


def copyfile(names:list=None,t:int=None):
    path = os.getcwd()
    print(f"Current Directory <{path}>")

    parent_directory = os.path.abspath(os.path.join(path,os.pardir))
    print(f"Parent Directory <{parent_directory}>")

    for name in names:
        folder = f"B_{name}"
        destination = os.path.join(parent_directory,folder)
        sizedir = len(os.listdir(destination))
        for k in [1,2,3]:
            file = f"../Metabolic_Engineering_Bilevel/Templates/run_experiment_k{k}.py"
            
            source = os.path.abspath(os.path.join(parent_directory,file))

            print(f">> Copying from <{source}>")
            print(f">> Pasting to <{destination}>")
            time.sleep(t)
         
            shutil.copy(source,destination)
            
            if sizedir < len(os.listdir(destination)):
                
                print('Successful')
                print(os.listdir(destination))
                time.sleep(t+1)
    print(f'''
        TERMINATED!!!!
    
    ''')
if __name__ == '__main__':
    copyfile(names=sys.argv[1:len(sys.argv)-1],t=int(sys.argv[-1]))
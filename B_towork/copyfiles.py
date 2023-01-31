import sys
import os
import time
import shutil


def copyfile_name(names:list=None,t:int=None):
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


def copy_file_chname(cfile:str=None,bacteria:list=None,nfiles:int=None,t:int=None):

    path = os.getcwd()
    parent_directory = os.path.abspath(os.path.join(path,os.pardir))
    print(f"Parent Directory <{parent_directory}>")
    sfile = f"Templates/{cfile}"
    source = os.path.abspath(os.path.join(parent_directory,sfile))

    folder = f"B_{bacteria}"
    dest_folder = os.path.join(parent_directory,folder)
    cursizedir = len(os.listdir(dest_folder))
    print(dest_folder)

    ks = [i+1 for i in range(nfiles)]
    
    print(f'From {source} ->')

    for k in ks:
        file = f"run_experiment_p_k{k}.py"
        dest_file = os.path.join(dest_folder,file)
        print(dest_file)
        shutil.copyfile(source,dest_file)

    finsizedir = len(os.listdir(dest_folder))

    print(f"# {finsizedir-cursizedir} Files Added!!")


if __name__ == '__main__':
    # copyfile_name(names=sys.argv[1:len(sys.argv)-1],t=int(sys.argv[-1]))
    file_to_copy = sys.argv[1]
    to_bacteria = sys.argv[2]
    n_copies = sys.argv[3]
    copy_file_chname(cfile=file_to_copy,bacteria=to_bacteria,nfiles=int(n_copies))


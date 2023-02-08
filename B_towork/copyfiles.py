import sys
import os
import shutil

def file_copier(folder:str=None,o_file:str=None,d_folder:str=None,c_file:str=None):
    path = os.getcwd()
    parent_directory = os.path.abspath(os.path.join(path,os.pardir))
    
    print(f">> Information about Source <<")
    folder = input(">> Enter name of the folder: \n -> ")
    original_folder = os.path.join(parent_directory,folder)
    print(f"Files in the source folder : \n {os.listdir(original_folder)}")
    
    
    
    o_file = input(">> Enter the name of the file: \n -> ")
    source = f"{folder}/{o_file}"
    source_dir = os.path.abspath(os.path.join(parent_directory,source))



    print(f">> Information about Destination <<")
    d_folder = input(">> Enter name of the folder: \n -> ") 
    d_dir = os.path.join(parent_directory,d_folder)

    print(f"Files in the destinarion folder : \n {os.listdir(d_dir)}")

    c_file = input(">> Enter new name of file: \n -> ")
    newfile = os.path.join(d_dir,c_file)

    shutil.copyfile(source_dir,newfile)

    print("New file copied!!")



if __name__ == '__main__':
    file_copier()









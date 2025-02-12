import os
import sys

def list_all_files(folder_path):
    all_files = []

    for root, dirs, files in os.walk(folder_path):
        # Include hidden files and directories
        files = [f for f in files if not f.startswith('.')]

        for file in files:
            file_path = os.path.join(root, file)
            all_files.append(file_path)

    return all_files
def get_ftype(filename):
    idx = filename.rfind(".")
    file_type = filename[idx+1:]
    return file_type


folder_path = None #full path to the zip folder including /zip/
all_files = list_all_files(folder_path)
new_paths = []
for file in all_files:
    old_f = file
    new_fname = file.split("/zip/")[1].replace("/", "-")
    if ".zip" not in new_fname:
        ftype = get_ftype(new_fname)
        new_folder_destination = folder_path[:-3]+f"{ftype}"
        if not os.path.exists(new_folder_destination): 
            os.makedirs(new_folder_destination) 
        file = f"{folder_path}/{new_fname}"
        file = file.replace("/zip/", f"/{ftype}/")
        new_paths.append((old_f, file))
for i in range(len(new_paths)):
    n = new_paths[i]
    os.rename(n[0], n[1])
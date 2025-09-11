'''
Move files froom asym_program's output to correct t_binned folders

uses an input csv as a guide for distiguishing t_bins 
'''
from asym_extraction_chamber.py_extractor import asym_extractor_py

import ROOT
import numpy as np

import pexpect as p
import os
import csv
import json


def has_header(file_path, header_list):
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return False
    
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        first_row = next(reader, None)
        return first_row == header_list

def move_files():
    t_bins = []
    path_in_progress = True
    is_nominal = False
    var_set = ""
    var = ""
    var_shift = ""
    
    while path_in_progress:
        var_set = str(input(f"Which variable set was processed?\n{os.listdir('../systematic_study/')} or \'exit\' to go quit: ")) # CTime, HMS, RF
        
        if var_set == "exit":
            print("Exiting...")
            path_in_progress = False
            return
        elif var_set not in os.listdir('../systematic_study/'):
            print(f"Variable set {var_set} does not exist in the systematic study folder. Please check your input or")
            
            create_new_var_folder = str(input("Would you like to create a new variable set folderr? (y/n): ")).lower()
            if create_new_var_folder == "y" and var_set != "exit":
                os.makedirs(f"../systematic_study/{var_set}", exist_ok=True)
                print(f"Created new variable set folder: {var_set}")
            else:
                continue
        elif var_set == "nominal":
            is_nominal = True
            print(f"Variable set {var_set} has been chosen.")
            print("Established a direct path.")
            path_in_progress = False
            continue

        
        print(f"Variable set {var_set} has been chosen.")
        print("Select the variable in the set. Options: ")
        
        var = str(input(f"{os.listdir(f'../systematic_study/{var_set}/')} or \'reset\' to change variable set: "))

        if var not in os.listdir(f'../systematic_study/{var_set}/') and var != "reset":
            print(f"Variable {var} does not exist in the variable set folder. Please check your input.")
            
            create_new_var_folder = str(input("Would you like to create a new variable folder under selected name? (y/n): ")).lower()
            if create_new_var_folder == "y" and var != "reset":
                os.makedirs(f'../systematic_study/{var_set}/{var}', exist_ok=True)
                print(f"Created new variable folder: {var}")
            else:
                continue
        elif var == "reset":
            print("Reseting variable set selection...")
            continue
        
        print(f"Variable {var} has been chosen.")
        print("Select the variable shift. Variable should have increased or decreased by user set amount. Options: ")
        
        var_shift = str(input(f"{os.listdir(f'../systematic_study/{var_set}/{var}')} or \'reset\' to change variable and variable set: ")).lower()
        
        if var_shift not in os.listdir(f'../systematic_study/{var_set}/{var}') and var_shift != "reset":
            print(f"Variable shift {var_shift} does not exist in the variable folder. Please check your input.")
            
            create_new_var_folder = str(input("Would you like to create a new variable shift folder under selected name? (y/n): ")).lower()
            if create_new_var_folder == "y" and var_shift != "reset":
                os.makedirs(f"../systematic_study/{var_set}/{var}/{var_shift}", exist_ok=True)
                print(f"Created new variable shift folder: {var_shift}")
            else:
                continue
        elif var_shift == "reset":
            print("Reseting variable selection...")
            continue
        
        print(f"Variable shift {var_shift} has been chosen.")
        path_in_progress = False
        
    print("Variables chosen for path construction")  
    
    path_to_store = f"../systematic_study/{var_set}/{var}/{var_shift}/"
    if is_nominal:
        path_to_store = f'../systematic_study/{var_set}/'
    
    # Read the input file to know the particle and t_bins in lists
    with open("input/rcl_wa_inp.csv", newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            t_bins.append(row[1])
        
    for t_bin in t_bins:
        graphs_name = f"output/*{t_bin}.root"
        yield_name = f"output/*{t_bin}.csv"
        kin_name = f"output/kinav_*"
        
        dest_dir = t_bin[1:]
        
        p.run(f"mkdir -p {path_to_store}/{dest_dir}")
        print(f"mkdir -p {path_to_store}/{dest_dir}")
        
        p.run(f"mkdir -p {path_to_store}/kinematics")
        print(f"mkdir -p {path_to_store}/kinematics")
        
        os.system(f"mv -f {graphs_name} {yield_name} {path_to_store}/{dest_dir}/")
        print(f"mv -f {graphs_name} {yield_name} {path_to_store}/{dest_dir}/")
        
        os.system(f"mv -f {kin_name} {path_to_store}/kinematics/")
        print(f"mv -f {kin_name} {path_to_store}/kinematics/")
    
    to_extractor = str(input("\nWould you like to copy WA plots and kinematics to the asym extractor? (y/n): ")).lower()
    
    if to_extractor == "y":
        wa = "WA_*"
        kin_name = "*_plots.root"
        
        for t_bin in t_bins:
            dest_dir = t_bin[1:]
            
            os.system(f"cp -f {path_to_store}/{dest_dir}/{wa} asym_extraction_chamber/")
            print(f"cp -f {path_to_store}/{dest_dir}/{wa} asym_extraction_chamber/")
            
            os.system(f"cp -f {path_to_store}/kinematics/{kin_name} asym_extraction_chamber/")
            print(f"cp -f {path_to_store}/kinematics/{kin_name} asym_extraction_chamber/")

def extract_LTp():
    ptcl = ""
    t_bins =[]
    t_bin_count = 0
    with open("input/kin_wa_inp.csv", newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            ptcl = row[0]
            if float(row[3]) in t_bins:
                break
            if t_bin_count < 1:
                t_bins.append(float(row[2]))
            t_bins.append(float(row[3]))
            t_bin_count += 1
    
    os.chdir("asym_extraction_chamber/")
    asym_extractor_py(ptcl, t_bin_count, t_bins)
    os.chdir("../")
    
    to_move = str(input("\nWould like to copy the LTp/asym vs t to other directory? \nEnter the path to that directory or \'n\' to skip: "))
    if to_move.lower() != "n":
        os.system(f"cp -f asym_extraction_chamber/LTp_vs_t_{ptcl}.root asym_extraction_chamber/asym_vs_vars_{ptcl}.root {to_move}")
        print(f"Copied LTp vs t to {to_move}")
    else:
        print("Skipping copying LTp vs t. Invalid input or \'n\' entered.")

def csv_LTp():
    print("\n Creating LTp csv for this shift.")
    
    t_bins = []
    
    nt = 0
    with open("input/kin_wa_inp.csv", newline='') as inpfile:
        reader = csv.reader(inpfile)
        next(reader)
        for row in reader:
            particle = row[0]
            if float(row[3]) in t_bins:
                break
            if nt < 1:
                t_bins.append(float(row[2]))
            t_bins.append(float(row[3]))
            nt += 1
    
    headers = ["t min", "t max", "var name", "shift", "LTp", "dLTp"]
    
    
    
    # Select the correct variables and folders to extract the LTP into a CSV.
    
    path_in_progress = True
    is_nominal = False
    var_set = ""
    var = ""
    var_shift = ""
    
    while path_in_progress:
        var_set = str(input(f"Which variable set was processed?\n{os.listdir('../systematic_study/')} or \'exit\' to go back: ")) # CTime, HMS, RF
        
        if var_set == "exit":
            print("Exiting...")
            path_in_progress = False
            return
        elif var_set not in os.listdir('../systematic_study/'):
            print(f"Variable set {var_set} does not exist in the systematic study folder. Please check your input or")
            
            create_new_var_folder = str(input("Would you like to create a new variable set folderr? (y/n): ")).lower()
            if create_new_var_folder == "y" and var_set != "exit":
                os.makedirs(f"../systematic_study/{var_set}", exist_ok=True)
                print(f"Created new variable set folder: {var_set}")
            else:
                continue
        elif var_set == "nominal":
            is_nominal = True
            print(f"Variable set {var_set} has been chosen.")
            print("Established a direct path.")
            path_in_progress = False
            continue
        
        print(f"Variable set {var_set} has been chosen.")
        print("Select the variable in the set. Options: ")
        
        var = str(input(f"{os.listdir(f'../systematic_study/{var_set}/')} or \'reset\' to change variable set: "))

        if var not in os.listdir(f'../systematic_study/{var_set}/') and var != "reset":
            print(f"Variable {var} does not exist in the variable set folder. Please check your input.")
            
            create_new_var_folder = str(input("Would you like to create a new variable folder under selected name? (y/n): ")).lower()
            if create_new_var_folder == "y" and var != "reset":
                os.makedirs(f'../systematic_study/{var_set}/{var}', exist_ok=True)
                print(f"Created new variable folder: {var}")
            else:
                continue
        elif var == "reset":
            print("Reseting variable set selection...")
            continue
        
        print(f"Variable {var} has been chosen.")
        print("Select the variable shift. Variable should have increased or decreased by user set amount. Options: ")
        
        var_shift = str(input(f"{os.listdir(f'../systematic_study/{var_set}/{var}')} or \'reset\' to change variable and variable set: ")).lower()
        
        if var_shift not in os.listdir(f'../systematic_study/{var_set}/{var}') and var_shift != "reset":
            print(f"Variable shift {var_shift} does not exist in the variable folder. Please check your input.")
            
            create_new_var_folder = str(input("Would you like to create a new variable shift folder under selected name? (y/n): ")).lower()
            if create_new_var_folder == "y" and var_shift != "reset":
                os.makedirs(f"../systematic_study/{var_set}/{var}/{var_shift}", exist_ok=True)
                print(f"Created new variable shift folder: {var_shift}")
            else:
                continue
        elif var_shift == "reset":
            print("Reseting variable selection...")
            continue
    
        print(f"Variable shift {var_shift} has been chosen.")
        path_in_progress = False
    
    print("Variables chosen for path construction")
    
    path_to_csv = f"../systematic_study/comp_Alu.csv"
    path_to_root = f"../systematic_study/{var_set}/{var}/{var_shift}/asym-ltp/asym_vs_vars_{particle}.root"
    if is_nominal:
        path_to_root = f"../systematic_study/{var_set}/asym-ltp/asym_vs_vars_{particle}.root"
    path_to_nom_cfg = f"../systematic_study/nominal/cut_config.json"
    path_to_sft_cfg = f"../systematic_study/{var_set}/{var}/{var_shift}/cut_config.json"
    
    nom_file = open(path_to_nom_cfg, "r")
    sft_file = open(path_to_sft_cfg, "r")
    nom_cfg = json.load(nom_file)
    sft_cfg = json.load(sft_file)
    
    shifted_var = "none"
    var_shift_value = 0.0
    
    for i in nom_cfg:
        if sft_cfg[i] - nom_cfg[i] != 0:
            shifted_var = i
            var_shift_value = '{0:.5f}'.format(sft_cfg[i] - nom_cfg[i])
            break
    
    if not has_header(path_to_csv, headers):
        mode = "w"
    else:
        mode = "a"
    
    with open(path_to_csv, mode, newline='') as outfile:
        writer = csv.writer(outfile)
        
        if mode == "w":
            writer.writerow(headers)
        
        LTp_file = ROOT.TFile(path_to_root, "READ")

        LTp_graph = LTp_file.Get("asym_vs_vars").FindObject("Graph")
        
        for i in range(nt):
            row = []

            LTp_pt = LTp_graph.GetY()[i]
            LTp_pt_err = LTp_graph.GetErrorY(i)
            
            row.append(t_bins[i])
            row.append(t_bins[i+1])
            row.append(shifted_var)
            row.append(var_shift_value)
            row.append(LTp_pt)
            row.append(LTp_pt_err)
            
            writer.writerow(row)
        
        LTp_file.Close()
    

ans = True
print("Welcome to systematic file system!")
while ans:
    print("""
    Actions:
    1. Move graphs and yields to respective folders.
    2. Calculate LTp vs t.
    3. Extract csv of LTp per tbin.
    4. Exit
    """)
    
    
    ans = int(input("What would you like to do?: "))
    if ans == 1:
        move_files()
    elif ans == 2:
        extract_LTp()
    elif ans == 3:
        csv_LTp()
    elif ans == 4:
        print("Exiting...")
        ans = None
        print("Exited. C-Ya")
    else:
        print("Invalid Entry! Try again.")

    
# dir_serch = p.spawn(f"find ../systematic_study/ -path \"*{variable_shift}/{dest_dir}\"")
# dir_found = dir_serch.expect([
#     f"{variable_shift}/{dest_dir}",
#     p.EOF, 
#     p.TIMEOUT
# ])
# dir_serch.close()

# if dir_found == 0:
#     x = "send file to dest_dir of the right t_bin"
# elif dir_found == 1:
#     x = "directory might not exists, create new directory for a t bin and send the right files to it"
# elif dir_found == 2:
#     x = "timeout for the expect(). unexpected output/error in find"
        

# mypassword = ""

# p = pexpect.spawn ('scp output/kinav_pion_center.csv ivanz@loon.phys.uregina.ca:/home/ivanz/bsa2025/ivan/asym_program/output')
# p.expect ('password: ')
# p.sendline (mypassword)
# p.expect(pexpect.EOF)


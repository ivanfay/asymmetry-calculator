import subprocess
import csv

action = input("What program do you want to run (Asym, WA, Kin)?: ")

if action == "asym" or action == "Asym":
    cpp_exe = "./asym_calc"

    with open("input/iterator_inp.csv", newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            cmd = [cpp_exe] + row 
            print("Running:" + " ".join(cmd))
            result = subprocess.run(cmd)
            
            
elif action == "WA" or action == "wa":
    
    with open("input/rcl_wa_inp.csv", newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            ptcl, t_bin = row[0], row[1]
            cmd = [
                "root",
                "-l",
                "-q",
                "-b",
                f'error_weighted_average.C("{ptcl}", "{t_bin}")'
            ]
            print("Running:" + " ".join(cmd))
            result = subprocess.run(cmd)
            
elif action == "Kin" or action == "kin":
    tcount = 0
    tlist = []
    ptcl = ""
    
    cpp_exe = "./kin_average"
    
    with open("input/kin_wa_inp.csv", newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        
        for row in reader:
            ptcl, shms_pos, t_min, t_max = row[0], row[1], row[2], row[3]
            
            if t_min not in tlist:
                tlist.append(t_min)
                tcount += 1
            
            cmd = [cpp_exe, ptcl, shms_pos, t_min, t_max]
            
            print("Running:" + " ".join(cmd))
            result = subprocess.run(cmd)

        cmd = [
                "root",
                "-l",
                "-q",
                "-b",
                f'kinav_reader.C("{ptcl}", {tcount})'
            ]
        print("Running:" + " ".join(cmd))
        result = subprocess.run(cmd)

else:
    print("Wrong inout. Try different action.")
        
    
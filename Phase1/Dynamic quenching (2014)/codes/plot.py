import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import os
import subprocess

data_path = "/mnt/Gayathri's Work drive/MSc_thesis_9th_sem/Vishniac-flux-9th-sem/Fortran codes 2014/run_files"
fig_path = "/mnt/Gayathri's Work drive/MSc_thesis_9th_sem/Vishniac-flux-9th-sem/Fortran codes 2014/figures"
data_save_path = "/mnt/Gayathri's Work drive/MSc_thesis_9th_sem/Vishniac-flux-9th-sem/Fortran codes 2014/data_files"

passed_args = argparse.ArgumentParser()

passed_args.add_argument(
    '-t',
    '--trial',
    type=int,
    help='trial number',
)

args = passed_args.parse_args()

if args.trial is None or type(args.trial) != int:
    save_path = Path(data_path)
    exisiting_trials = [d for d in save_path.iterdir() if d.is_dir()]
    trial_numbers = [int(d.name.split('trial_')[-1]) for d in exisiting_trials]
    last_trial = max(trial_numbers)
    trial_num = last_trial + 1
else:
    trial_num = args.trial
    
# check if the trial directory exists
if os.path.exists(f'{data_save_path}/trial_{trial_num}'):
    print(f"trial_{trial_num} directory already exists")
    permission = input("Do you want to overwrite the existing directory? (y/n): ")
    if permission == 'n' or permission == 'N':
        print("Exiting the program")
        exit()
    elif permission == 'y' or permission == 'Y':
        print("Overwriting the existing directory")
        

fig_path = f'{fig_path}/trial_{trial_num}'
data_save_path = f'{data_save_path}/trial_{trial_num}'

os.makedirs(data_save_path, exist_ok=True)
os.makedirs(fig_path, exist_ok=True)

#import txt file 
filename = 'eta_fz_values.txt'
file_path = f"{data_path}/{filename}"

with open(file_path,'r') as f:
    lines=f.readlines()
# print(lines)

eta=np.array(lines,dtype=float)
#import txt file

filename = 'z_values.txt'
file_path = f"{data_path}/{filename}"

with open(file_path,'r') as f:
    lines=f.readlines()
# print(lines)

z=np.array(lines,dtype=float)
np.savetxt(f'{data_save_path}/z_values.txt', z)
plt.plot(z,eta)
plt.xlabel('z')
plt.ylabel('eta')
plt.title('eta vs z')
plt.xlim(-15,15)
plt.savefig(f'{fig_path}/eta_vs_z.png')
plt.close()
#import txt file
filename = 'alpha_values.txt'

file_path = f"{data_path}/{filename}"

try:
    with open(file_path,'r') as f:
        lines=f.readlines()
    # print(lines)

    alpha=np.array(lines,dtype=float)
    np.savetxt(f'{data_save_path}/alpha_values.txt', alpha)

    plt.plot(z,alpha)
    plt.xlabel('z')
    plt.ylabel('alpha')
    plt.title('alpha vs z')
    plt.xlim(-15,15)
    plt.savefig(f'{fig_path}/alpha_vs_z.png')
    plt.close()
except:
    print('alpha_values.txt file not found')

#import txt file
filename1 = 'Br_ini.txt'
filename2 = 'B_phi_ini.txt'
file_path1 = f"{data_path}/{filename1}"
file_path2 = f"{data_path}/{filename2}"

try:
    with open(file_path1,'r') as f:
        lines1=f.readlines()
    # print(lines1)
    with open(file_path2,'r') as f:
        lines2=f.readlines()
    # print(lines2)
    Br=np.array(lines1,dtype=float)
    Bphi=np.array(lines2,dtype=float)

    plt.plot(z,Br)
    plt.plot(z,Bphi)
    plt.xlabel('z')
    plt.ylabel('Br, Bphi')
    plt.title('Br, Bphi vs z - initial')
    plt.savefig(f'{fig_path}/Br_Bphi_vs_z_initial.png')
    plt.close()
except:
    print('Br_ini.txt or B_phi_ini.txt file not found')

#import txt file
filename1 = 'Br_final.txt'
filename2 = 'B_phi_final.txt'
filename3 = 'time.txt'
file_path1 = f"{data_path}/{filename1}"
file_path2 = f"{data_path}/{filename2}"
file_path3 = f"{data_path}/{filename3}"


try:
    with open(file_path1,'r') as f:
        lines1=f.readlines()
    # print(lines1)
    with open(file_path2,'r') as f:
        lines2=f.readlines()
    # print(lines2)
    with open(file_path3,'r') as f:
        lines3=f.readlines()
    # print(lines3)

    Br_list = []
    for line in lines1:
        line = line.strip()
        line = line.split()
        curr = np.array(line, dtype=float)
        Br_list.append(curr)
    Br_list = np.array(Br_list)

    Bphi_list = []
    for line in lines2:
        line = line.strip()
        line = line.split()
        curr = np.array(line, dtype=float)
        Bphi_list.append(curr)
    Bphi_list = np.array(Bphi_list)


    # Br_list=np.array(lines1,dtype=float)
    # Bphi_list=np.array(lines2,dtype=float)
    time_list=np.array(lines3,dtype=float)
    
    np.savetxt(f'{data_save_path}/time.txt', time_list)
    np.savetxt(f'{data_save_path}/Br_final.txt', Br_list)
    np.savetxt(f'{data_save_path}/B_phi_final.txt', Bphi_list)

    plt.plot(z, Br_list[-1], label='Br')
    plt.plot(z, Bphi_list[-1], label='Bphi')
    plt.xlabel('z')
    plt.ylabel('Br, Bphi')
    plt.title(f'Br, Bphi vs z at t={time_list[-1]}')
    plt.legend()
    plt.savefig(f'{fig_path}/Br_Bphi_vs_z_final.png')
    plt.close()
except:
    print('Br_final.txt or B_phi_final.txt or time.txt file not found')


#plot of Br and Bphi at against time
plt.plot(time_list, Br_list[:,25], label='Br')
plt.plot(time_list, Bphi_list[:,25], label='Bphi')
plt.xlabel('time')
plt.ylabel('Br, Bphi')
plt.title('Br, Bphi vs time')
plt.legend()
plt.savefig(f'{fig_path}/Br_Bphi_vs_time.png')
plt.close()


B_strength = np.sqrt(Br_list**2 + Bphi_list**2)
np.savetxt(f'{data_save_path}/B_strength.txt', B_strength)

plt.plot(time_list, B_strength[:,-1])
plt.xlabel('time')
plt.ylabel('B_strength')
plt.title('B_strength vs time')
plt.yscale('log')
plt.savefig(f'{fig_path}/B_strength_vs_time.png')
plt.close()
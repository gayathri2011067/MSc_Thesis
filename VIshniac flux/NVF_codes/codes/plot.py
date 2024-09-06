import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import os
import subprocess

data_path = "/home/gayathri/MSc_thesis/VIshniac flux/NVF_codes/run_files"
fig_path = "/home/gayathri/MSc_thesis/VIshniac flux/NVF_codes/figures"
data_save_path = "/home/gayathri/MSc_thesis/VIshniac flux/NVF_codes/data_files"

passed_args = argparse.ArgumentParser()

passed_args.add_argument(
    '-t',
    '--trial',
    type=int,
    default=None,
    help='trial number',
)

args = passed_args.parse_args()

if args.trial is None or type(args.trial) != int:
    save_path = Path(data_save_path)
    exisiting_trials = [d for d in save_path.iterdir() if d.is_dir()]
    trial_numbers = []
    for d in exisiting_trials:
        if d.name.startswith('trial_'):
            trial_numbers.append(int(d.name.split('trial_')[-1]))
    # print(f"Existing trials: {trial_numbers}")
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
else:
    print(f"trial_{trial_num} directory does not exist. Creating new directory")

fig_path = f'{fig_path}/trial_{trial_num}'
data_save_path = f'{data_save_path}/trial_{trial_num}'

os.makedirs(data_save_path, exist_ok=True)
os.makedirs(fig_path, exist_ok=True)

#import txt file 
# filename = 'eta_fz_values.txt'
# file_path = f"{data_path}/{filename}"

# with open(file_path,'r') as f:
#     lines=f.readlines()
# # print(lines)

# eta=np.array(lines,dtype=float)
# #import txt file

filename = 'z_values.txt'
file_path = f"{data_path}/{filename}"

with open(file_path,'r') as f:
    lines=f.readlines()
# print(lines)

z=np.array(lines,dtype=float)
np.savetxt(f'{data_save_path}/z_values.txt', z)
# plt.plot(z,eta)
# plt.xlabel('z')
# plt.ylabel('eta')
# plt.title('eta vs z')
# plt.xlim(-15,15)
# plt.savefig(f'{fig_path}/eta_vs_z.png')
# plt.close()
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
    # print(Br_list.shape)
    plt.plot(z, Br_list[-1], label='Br')
    plt.plot(z, Bphi_list[-1], label='Bphi')
    # print(Br_list[1])
    # print(Bphi_list[1])
    # plt.plot(z, Br_list[1], label='Br')
    # plt.plot(z, Bphi_list[1], label='Bphi')
    # plt.plot(z, Br_list[2], label='Br')
    # plt.plot(z, Bphi_list[2], label='Bphi')
    plt.xlabel('z')
    plt.ylabel('Br, Bphi')
    plt.title(f'Br, Bphi vs z at t={time_list[0]}')
    # plt.xlim(-0.25,0.25)
    plt.axvline(x=1, color='r', linestyle='--')
    plt.axvline(x=-1, color='r', linestyle='--')
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

# plt.plot(time_list,alpha)
# plt.xlabel('time')
# plt.ylabel('alpha')
# plt.title('alpha vs time')
# plt.xlim(-15,15)
# plt.savefig(f'{fig_path}/alpha_vs_time.png')
# plt.close()
B_strength = np.sqrt(Br_list**2 + Bphi_list**2)
np.savetxt(f'{data_save_path}/B_strength.txt', B_strength)

plt.plot(time_list, B_strength[:,-1])
plt.xlabel('time')
plt.ylabel('B_strength')
plt.title('B_strength vs time')
plt.yscale('log')

plt.savefig(f'{fig_path}/B_strength_vs_time.png')
plt.close()

# plot of alpha_m final vs z
filename = 'alpham_final.txt'
file_path = f"{data_path}/{filename}"

try:
    with open(file_path,'r') as f:
        lines=f.readlines()

    alpha_list = []
    for line in lines:
        line = line.strip()
        line = line.split()
        curr = np.array(line, dtype=float)
        alpha_list.append(curr)
    alpha_m=np.array(alpha_list,dtype=float)
    np.savetxt(f'{data_save_path}/alpha_m.txt', alpha_m)

    n1 = 500
    Nt = 50000
    total_t = 50
    n2 = total_t*Nt/n1
    t_print = [10, 20, 30, 40, 50]
    t_print = [int(t_val*n1/total_t) for t_val in t_print]

    for t_val in t_print:
        plt.plot(z, alpha_m[t_val - 1], label=f't={t_val*total_t/n1:.2f}') # t_val - 1 because the index starts from 0
    plt.xlabel('z')
    plt.ylabel('alpha_m')
    plt.title('alpha_m vs z')
 
    plt.legend()
    plt.savefig(f'{fig_path}/alpha_m_vs_z.png')
    plt.close()
except:
    print('alpha_m.txt file not found')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os, sys
import subprocess

data_path = "/home/gayathri/MSc_thesis/VIshniac flux/Sharanyas_notes_codes/data_files"
fig_path = "/home/gayathri/MSc_thesis/VIshniac flux/Sharanyas_notes_codes/figures"

trial_numbers = []
labels = []

t_num = 1
while True:
    try:
        trial_numbers.append(int(sys.argv[t_num]))
        labels.append(sys.argv[t_num+1])
        t_num += 2
    except:
        break

folder_name = f'trials_({"_".join([str(i) for i in trial_numbers])})'
    
# check if the trial directory exists
if os.path.exists(f'{fig_path}/{folder_name}'):
    print(f"{folder_name} directory already exists")
    permission = input("Do you want to overwrite the existing directory? (y/n): ")
    if permission == 'n' or permission == 'N':
        print("Exiting the program")
        exit()
    elif permission == 'y' or permission == 'Y':
        print("Overwriting the existing directory")

fig_path = f'{fig_path}/{folder_name}'

os.makedirs(fig_path, exist_ok=True)

z_lists = []
filename = 'z_values.txt'

for trial_num in trial_numbers:
    file_path = f"{data_path}/trial_{trial_num}/{filename}"

    z = np.loadtxt(file_path)
   
    z_lists.append(np.copy(z))

#import txt file
filename1 = 'Br_final.txt'
filename2 = 'B_phi_final.txt'
filename3 = 'time.txt'


Brs = []
Bphis = []
B_strengths = []
times = []


for i in range(len(trial_numbers)):
    
    trial_num = trial_numbers[i]

    file_path1 = f"{data_path}/trial_{trial_num}/{filename1}"
    file_path2 = f"{data_path}/trial_{trial_num}/{filename2}"
    file_path3 = f"{data_path}/trial_{trial_num}/{filename3}"
    

    Br_list = np.loadtxt(file_path1)
    Bphi_list = np.loadtxt(file_path2)
    time_list = np.loadtxt(file_path3)
   
    
    # Br_list=np.array(lines1,dtype=float)
    # Bphi_list=np.array(lines2,dtype=float)
    B_strength = np.sqrt(Br_list**2 + Bphi_list**2)
    
    B_strengths.append(np.copy(B_strength))
    
    times.append(np.copy(time_list))
    
    Brs.append(np.copy(Br_list))
    Bphis.append(np.copy(Bphi_list))
   

    plt.plot(z_lists[i], Br_list[-1], label=f'Br: {labels[i]}')
    plt.plot(z_lists[i], Bphi_list[-1], label=f'Bphi: {labels[i]}')
    plt.xlim(-1,1)
    
plt.xlabel('z')
plt.ylabel('Br, Bphi')
plt.title(f"Br, Bphi vs z at t=({','.join([f'trial{trial_numbers[i]}:{times[i][-1]:.2e}' for j in range(len(trial_numbers))])})")

plt.legend()
plt.savefig(f'{fig_path}/Br_Bphi_vs_z_final.png')
plt.close()
        
#plot of Br and Bphi at against time
print(Br_list.shape)
space_indices = []
for i in range(len(trial_numbers)):

    lab = labels[i]
    
    space_idx = (Br_list.shape[1]-1)//2
    space_indices.append(space_idx)
    
    plt.plot(times[i], Brs[i][:,space_idx], label=f'Br: {lab}')
    plt.plot(times[i], Bphis[i][:,space_idx], label=f'Bphi: {lab}')
print(space_indices)   

plt.xlabel('time')
plt.ylabel('Br, Bphi')

plt.title(f"Br, Bphi vs time at z=({','.join([f'trial{trial_numbers[i]}:{int(z_lists[i][space_indices[i]])}' for i in range(len(trial_numbers))])})")
plt.legend()
plt.savefig(f'{fig_path}/Br_Bphi_vs_time.png')
plt.close()

for i in range(len(trial_numbers)):

    lab = labels[i]
    space_idx = 107//2

    plt.plot(times[i], B_strengths[i][:,-1], label=f'{lab}')
plt.xlim(0,7)

plt.yscale('log')
plt.ylim(10**-4,1)
# plt.axhline(y=1, color='b', linestyle='--')
plt.xlabel('time(in Gyr)')
plt.ylabel(r'$B_{\text{strength}} / B_0$')

plt.title(f"Dynamical quenching")
plt.text(0.90, 1.05, r'$B_0=0.82\mu G$', transform=plt.gca().transAxes, verticalalignment='top')

plt.legend()
plt.savefig(f'{fig_path}/B_strength_vs_time.png')
plt.close()

#plot alpha_m vs z similarly
filename = 'alpha_m.txt'

alpha_m_values = []
for i in range(len(trial_numbers)):
    
    trial_num = trial_numbers[i]

    file_path = f"{data_path}/trial_{trial_num}/{filename}"
    
    alpha_m_list = np.loadtxt(file_path)
    
    alpha_m_values.append(np.copy(alpha_m_list))
    
    plt.plot(z_lists[i], alpha_m_list[-1], label=f'alpha_m: {labels[i]}')
    plt.xlim(-1,1)
    
plt.xlabel('z')
plt.ylabel('alpha_m')
plt.title(f"alpha_m vs z at t=({','.join([f'trial{trial_numbers[i]}:{times[i][-1]:.2e}' for i in range(len(trial_numbers))])})")
plt.legend()
plt.savefig(f'{fig_path}/alpha_m_vs_z_final.png')
plt.close()





import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import os
import subprocess


data_path = "/home/gayathri/MSc_thesis/VIshniac flux/FOSA/run_files"
fig_path = "/home/gayathri/MSc_thesis/VIshniac flux/FOSA/figures"
data_save_path = "/home/gayathri/MSc_thesis/VIshniac flux/FOSA/data_files"

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
    # print(f"Existing trials: {trial_numbers}")♥
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

filename = 'turb_vel.txt'
file_path = f"{data_path}/{filename}"

with open(file_path,'r') as f:
    lines=f.readlines()
# print(lines)

turb_vel=np.array(lines,dtype=float)
np.savetxt(f'{data_save_path}/turb_vel.txt', turb_vel)
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
    # plt.xlim(-15,15)
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
    # print(time_list)
    np.savetxt(f'{data_save_path}/time.txt', time_list)
    np.savetxt(f'{data_save_path}/Br_final.txt', Br_list)
    np.savetxt(f'{data_save_path}/B_phi_final.txt', Bphi_list)
    # print(Br_list.shape)
    plt.plot(z, Br_list[-1], label='Br')
    plt.plot(z, Bphi_list[-1], label='Bphi')

    # plt.plot(z, Br_list[-4], label='Br')
    # plt.plot(z, Bphi_list[-4], label='Bphi')
    print(Br_list.shape)
    # print(Bphi_list[1])
    # plt.plot(z, Br_list[1], label='Br')
    # plt.plot(z, Bphi_list[1], label='Bphi')
    # plt.plot(z, Br_list[2], label='Br')
    # plt.plot(z, Bphi_list[2], label='Bphi')
    plt.xlabel('z')
    plt.ylabel('Br, Bphi')
    plt.title(f'Br, Bphi vs z at t={time_list[-1]}')
    # plt.xlim(-0.25,0.25)
    plt.axvline(x=1, color='r', linestyle='--')
    plt.axvline(x=-1, color='r', linestyle='--')
    plt.axhline(y=0, color='g', linestyle='--')
    plt.legend()
    plt.savefig(f'{fig_path}/Br_Bphi_vs_z_final.png')
    plt.close()
except:
    print('Br_final.txt or B_phi_final.txt or time.txt file not found')


#plot of Br and Bphi at against time
plt.plot(time_list, Br_list[:,-1], label='Br')
# print(time_list)
plt.plot(time_list, Bphi_list[:,-1], label='Bphi')
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

# Plot B_strength at a specific z index
plt.plot(time_list, B_strength[:,51])
plt.xlabel('time')
# plt.xlim(0,10)

# plt.yscale('log')
plt.ylabel('B_strength/B0')
plt.title('B_strength vs time')
plt.yscale('log')


# plt.axhline(y=0.001, color='g', linestyle='--')
# plt.axvline(x=2.5, color='r', linestyle='--')
plt.savefig(f'{fig_path}/B_strength_vs_time.png')
plt.close()

# Plot field strength averaged over z as a function of time
B_strength_avg = np.sqrt(np.mean(B_strength**2, axis=1))
np.savetxt(f'{data_save_path}/B_strength_avg.txt', B_strength_avg)
plt.plot(time_list, B_strength_avg, color='red',linewidth=3)
plt.xlabel('time(Gyr)',fontsize=17)
plt.yscale('log')
# plt.xlim(4,5)
# plt.ylim(1,10)
plt.ylim(10**-2,1)
plt.ylabel(r'$B_{strength}/B_0$',fontsize=17)
# plt.title(r'$B_{strength}$ vs time')
plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
plt.savefig(f'{fig_path}/B_strength_avg_vs_time.png')
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
    n1 = 5000
    total_t = 10
    t_values = [2.5, 5,10]  # List of t values to plot

    # Initialize an empty list to store the t indices
    t_indices = []

    # Calculate t indices and append to the list
    for t_val in t_values:
        t_index = int(t_val * n1 / total_t)
        t_indices.append(t_index - 1)  # Append the adjusted index

    # Now you can use the t_indices list for plotting
    for t_index, t_val in zip(t_indices, t_values):
        plt.plot(z, alpha_m[t_index], label=f't={t_val:.2f}')  # Use the t_index directly from the list

    plt.axhline(y=0, color='grey', linestyle='--', linewidth=1, alpha=0.5)
    plt.axvline(x=0, color='grey', linestyle='--', linewidth=1, alpha=0.5)
    plt.xlim(-1, 1)
    plt.ylim(-100, 100)
    plt.xlabel('z')
    plt.ylabel(r'$\alpha_m$')
    plt.title(r'$\alpha_m$ vs z at different times')
    plt.legend()
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)

    # Plot turb velocity
    plt.plot(z, turb_vel, linestyle='--', color='grey', alpha=0.8, label='Turbulence Velocity')

    plt.legend()
    plt.savefig(f'{fig_path}/alpha_m_vs_z_multiple_times.png')
    plt.close()
except:
    print('alpha_m.txt file not found')





#plot bstrength squared on y and -alpha m on x
# print(np.shape(B_strength))

plt.plot( -alpha_m[-1],B_strength[-1]**2)
plt.ylabel('B_strength_avg')
plt.xlabel('-alpha_m')
plt.title('B_strength_avg vs -alpha_m')
plt.savefig(f'{fig_path}/B_strength_avg_vs_alpha_m.png')
plt.close()






filename = 'alpham_final.txt'
file_path = f"{data_path}/{filename}"
try:
    with open(file_path,'r') as f:
        lines=f.readlines()

    alpha_mm_list = []
    for line in lines:
        line = line.strip()
        line = line.split()  
        curr = np.array(line, dtype=float)
        alpha_mm_list.append(curr)
    alpha_mm=np.array(alpha_mm_list,dtype=float)
    np.savetxt(f'{data_save_path}/aallppmm.txt', alpha_mm)
    # print(alpha_tot)
    plt.plot(z, alpha_mm[-1],label='alpha_m',color='red',linewidth=3)
    plt.xlabel(r'z($\times$ 0.5 kpc)',fontsize=17)
    plt.ylabel(r'$\alpha_m$ (km s$^{-1}$)',fontsize=17)
    # plt.title('aallppmm vs z')
    plt.savefig(f'{fig_path}/aallppmm_vs_z.png')
    plt.close()
except:
    print('alpham_final.txt file not found')

#read and save alpha_tot similarly
filename = 'alp_tot.txt'
file_path = f"{data_path}/{filename}"

try:
    with open(file_path,'r') as f:
        lines=f.readlines()

    alpha_tot_list = []
    for line in lines:
        line = line.strip()
        line = line.split()  
        curr = np.array(line, dtype=float)
        alpha_tot_list.append(curr)
    alpha_tot=np.array(alpha_tot_list,dtype=float)
    np.savetxt(f'{data_save_path}/alpha_tot.txt', alpha_tot)
    # print(alpha_tot)
    plt.plot(z, alpha_tot[-1])
    plt.xlabel('z')
    plt.ylabel('alpha_tot')
    plt.title('alpha_tot vs z')
    plt.savefig(f'{fig_path}/alpha_tot_vs_z.png')
    plt.close()
except:
    print('alpha_tot.txt file not found')
    

  
# #  plot alpha and alpha total in same plot against time
plt.plot(time_list,alpha_mm[:,51],label='alpha_m')
plt.plot(time_list,alpha_tot[:,51],label='alpha_total')
plt.xlabel('time')
plt.ylabel('alpha, alpha_tot')
plt.title('alpha, alpha_tot vs time')
plt.legend()
plt.savefig(f'{fig_path}/alpha_alpha_tot_vs_time.png')
plt.close()
   
    

#write a code to find and plot pitch angles |p| = arctan |Br/Bphi|

# Calculate the pitch angle and limit the result to the range [-pi/2, pi/2]
p = np.degrees(np.arctan((Br_list / Bphi_list)))  # Convert the angle to degrees
np.savetxt(f'{data_save_path}/pitch_angle.txt', p)

# Plot the pitch angle vs z
plt.plot(z, p[-1])
plt.xlabel('z')
plt.ylabel('Pitch Angle (degrees)')
plt.title('Pitch Angle vs z')
plt.savefig(f'{fig_path}/pitch_angle_vs_z.png')
plt.close()




















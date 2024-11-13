import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os, sys
import subprocess

# Paths to data and figures
data_path = "/home/gayathri/MSc_thesis/VIshniac flux/FOSA/data_files"
fig_path = "/home/gayathri/MSc_thesis/VIshniac flux/FOSA/figures"

trial_numbers = []
labels = []

# Read trial numbers and labels from command line arguments
t_num = 1
while True:
    try:
        trial_numbers.append(int(sys.argv[t_num]))
        labels.append(sys.argv[t_num + 1])
        t_num += 2
    except:
        break

# Create folder name for trials
folder_name = f'trials_({"_".join([str(i) for i in trial_numbers])})'

# Check if directory exists
if os.path.exists(f'{fig_path}/{folder_name}'):
    print(f"{folder_name} directory already exists")
    permission = input("Do you want to overwrite the existing directory? (y/n): ")
    if permission.lower() == 'n':
        print("Exiting the program")
        exit()
    elif permission.lower() == 'y':
        print("Overwriting the existing directory")

# Update figure path with folder name
fig_path = f'{fig_path}/{folder_name}'
os.makedirs(fig_path, exist_ok=True)

# Initialize lists to store data
z_lists = []
filename = 'z_values.txt'

# Load z values from the files
for trial_num in trial_numbers:
    file_path = f"{data_path}/trial_{trial_num}/{filename}"
    z = np.loadtxt(file_path)
    z_lists.append(np.copy(z))

# Initialize lists for storing magnetic field components and times
Brs = []
Bphis = []
B_strengths = []
B_average = []
times = []

# Load magnetic field and time data
for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]
    file_path1 = f"{data_path}/trial_{trial_num}/Br_final.txt"
    file_path2 = f"{data_path}/trial_{trial_num}/B_phi_final.txt"
    file_path3 = f"{data_path}/trial_{trial_num}/time.txt"

    Br_list = np.loadtxt(file_path1)
    Bphi_list = np.loadtxt(file_path2)
    time_list = np.loadtxt(file_path3)

    # Calculate magnetic field strength
    B_strength = np.sqrt(Br_list**2 + Bphi_list**2)
    B_average.append(np.mean(B_strength))

    # Append data to lists
    B_strengths.append(np.copy(B_strength))
    B_average.append(np.mean(B_strength))
    times.append(np.copy(time_list))
    Brs.append(np.copy(Br_list))
    Bphis.append(np.copy(Bphi_list))

    # Plot Br and Bphi versus z for the last time step
    plt.figure(figsize=(10, 5))  # Set figure size to rectangular
    plt.plot(z_lists[i], Br_list[-1], label=f'Br: {labels[i]}')
    plt.plot(z_lists[i], Bphi_list[-1], label=f'Bphi: {labels[i]}')
    plt.xlim(-1, 1)

# Label and title the plot, ensuring proper escaping for LaTeX syntax
plt.xlabel(r'$z$')
plt.ylabel(r'$B_r, B_\phi$')
plt.title(r'$B_r, B_\phi \ \mathrm{vs} \ z \ \mathrm{at} \ t=(' +
          ', '.join([f'trial{trial_numbers[j]}:{times[j][-1]:.2e}' for j in range(len(trial_numbers))]) + r')$')
# plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(f'{fig_path}/Br_Bphi_vs_z_final.png', bbox_inches='tight')
plt.close()

# Plot Br and Bphi versus time at the midpoint of z
space_indices = []
plt.figure(figsize=(10, 5))  # Set figure size to rectangular
for i in range(len(trial_numbers)):
    space_idx = (Br_list.shape[1] - 1) // 2
    space_indices.append(space_idx)
    plt.plot(times[i], Brs[i][:, space_idx], label=f'Br: {labels[i]}')
    plt.plot(times[i], Bphis[i][:, space_idx], label=f'Bphi: {labels[i]}')

plt.xlabel(r'$t$')
plt.ylabel(r'$B_r, B_\phi$')
plt.title(r'$B_r, B_\phi \ \mathrm{vs} \ t \ \mathrm{at} \ z=(' +
          ', '.join([f'trial{trial_numbers[i]}:{int(z_lists[i][space_indices[i]])}' for i in range(len(trial_numbers))]) + r')$')
# plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(f'{fig_path}/Br_Bphi_vs_time.png', bbox_inches='tight')
plt.close()

# Plot magnetic field strength versus time at z midpoint
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 5))  # Set figure size to rectangular
B_strength_avg = np.mean(B_strengths, axis=2)
line_styles = ['-', '--', '-.', ':']  # Define different line styles
colors = ['#f39237', '#9c3848', '#06a77d', '#3f8efc']  # Darker color scheme

for i in range(len(trial_numbers)):
    lab = labels[i]
    plt.plot(times[i], B_strength_avg[i], label=f'{lab}', linestyle='-', linewidth=2, color=colors[i % len(colors)])  # Increase line thickness

plt.xlim(0, 7)
plt.yscale('log')
plt.ylim(10**-2, 2)
plt.xlabel('time(in Gyr)')
plt.ylabel(r'${B}(0) / B_0$')
plt.title(f"Dynamical quenching")
plt.text(0.90, 1.05, r'$B_0=0.82\mu G$', transform=plt.gca().transAxes, verticalalignment='top')
plt.legend(loc='lower right', labels=[f'${label}$' for label in labels])
plt.savefig(f'{fig_path}/B_strength_vs_time.png')
plt.close()

plt.figure(figsize=(10, 5))  # Set figure size to rectangular
for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]
    
    plt.plot(times[i], B_strength_avg[i], label=f'{labels[i]}')
plt.xlim(0, 14)
plt.ylim(0, 0.6)
plt.xlabel(r'$t$')
plt.ylabel(r'$\mathrm{B}_{\mathrm{strength, avg}}$')
plt.title(r'$\mathrm{B}_{\mathrm{strength, avg}} \ \mathrm{vs} \ t$')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Save or display the plot
plt.savefig(f'{fig_path}/B_strength_avg_vs_time.png', bbox_inches='tight')
# Alternatively, you can use plt.show() to display it in real-time


# Close the plot
plt.close()

# Plot alpha_m versus z for the last time step
plt.figure(figsize=(10, 5))  # Set figure size to rectangular
alpha_m_values = []
filename = 'alpha_m.txt'

for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]
    file_path = f"{data_path}/trial_{trial_num}/{filename}"
    alpha_m_list = np.loadtxt(file_path)
    alpha_m_values.append(np.copy(alpha_m_list))

    plt.plot(z_lists[i], alpha_m_list[-1], label=f'alpha_m: {labels[i]}')
    plt.xlim(-0.5, 0.5)
plt.axhline(y=0, color='grey', linestyle='--')
plt.axvline(x=0, color='grey', linestyle='--')
plt.xlabel(r'$z$')
plt.ylabel(r'$\alpha_m$')
plt.title(r'$\alpha_m \ \mathrm{vs} \ z$')
# plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(f'{fig_path}/alpha_m_vs_z_final.png', bbox_inches='tight')
plt.close()


# Bstr_values = []
# filename = 'B_strength_avg.txt'

# for i in range(len(trial_numbers)):
#     trial_num = trial_numbers[i]
#     file_path = f"{data_path}/trial_{trial_num}/{filename}"
#     BSTR_list = np.loadtxt(file_path)
#     Bstr_values.append(np.copy(BSTR_list))

#     plt.plot(times[i], BSTR_list[i], label=f'B_strength: {labels[i]}')

# plt.xlabel(r'$t$')
# plt.ylabel(r'$B_avg$')
# plt.title(r'$B_avg \ \mathrm{vs} \ t$')
# # plt.grid(True)
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.savefig(f'{fig_path}/Bavg.png', bbox_inches='tight')
# plt.close()

fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid of subplots
filename = 'alpha_m.txt'
alpha_m_values = []

# Flatten the 2D array of axes for easier iteration
axs = axs.flatten()

for i in range(4):  # Assuming only 4 trials
    trial_num = trial_numbers[i]
    file_path = f"{data_path}/trial_{trial_num}/{filename}"
    alpha_m_list = np.loadtxt(file_path)
    alpha_m_values.append(np.copy(alpha_m_list))
    
    # Plot on the corresponding subplot
    axs[i].plot(z_lists[i], alpha_m_list[-1], label=f'alpha_m: {labels[i]}')
    axs[i].set_xlim(-1,1)
    axs[i].axhline(y=0, color='grey', linestyle='--')
    axs[i].axvline(x=0, color='grey', linestyle='--')
    axs[i].set_xlabel(r'$z$')
    axs[i].set_ylabel(r'$\alpha_m$')
    axs[i].set_title(f'Trial {trial_num}')
    axs[i].legend(loc='upper right')
    
# Adjust the layout so subplots do not overlap
plt.tight_layout()

# Save the entire figure as a single image
plt.savefig(f'{fig_path}/alpha_m_vs_z_subplots.png', bbox_inches='tight')
plt.close()


# First set of subplots: B_strength vs t at specific z for each trial
fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid of subplots
axs = axs.flatten()

for i in range(4):  # Assuming only 4 trials
    trial_num = trial_numbers[i]
    
    # Plot B_strength vs time at the specific z location for each trial
    axs[i].plot(times[i], B_strengths[i][:, space_indices[i]], label=f'{labels[i]}')
    axs[i].set_xlim(0, 10)
    axs[i].set_yscale('log')
    axs[i].set_xlabel(r'$t$')
    axs[i].set_ylabel(r'$\mathrm{B}_{\mathrm{strength}}$')
    axs[i].set_title(f'Trial {trial_num} at z={int(z_lists[i][space_indices[i]])}')
    axs[i].legend(loc='upper right')

plt.tight_layout()
plt.savefig(f'{fig_path}/B_strength_vs_time_subplots.png', bbox_inches='tight')
plt.close()

# Second set of subplots: B_strength_avg vs t for each trial
fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid of subplots
axs = axs.flatten()

B_strength_avg = np.mean(B_strengths, axis=2)

for i in range(4):  # Assuming only 4 trials
    trial_num = trial_numbers[i]
    
    # Plot B_strength_avg vs time for each trial
    axs[i].plot(times[i], B_strength_avg[i], label=f'{labels[i]}', color='red')
    axs[i].set_xlim(0, 7)
    axs[i].set_ylim(0.01,1)
    axs[i].set_yscale('log')
    axs[i].set_xlabel(r'$t \ (\mathrm{computational \ units})$')
    axs[i].set_ylabel(r'$\log \ \langle B \rangle / B_0$')
    axs[i].set_title(f'Trial {trial_num}')
    axs[i].legend(loc='lower right')
    
    # Add ticklines on all 4 sides of the grid
    axs[i].tick_params(axis='both', which='both', direction='in', top=True, right=True)

plt.tight_layout()
plt.savefig(f'{fig_path}/B_strength_avg_vs_time_subplots.png', bbox_inches='tight')
plt.close()


fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid of subplots
alpha_m_filename = 'alpha_m.txt'
turb_vel_filename = 'turb_vel.txt'
alpha_m_values = []
turb_vel_values = []
n1 = 5000
total_t = 10
t_values = 2.5  #NOTE: Change this value to plot at different time steps
t_index = int(t_values * n1 / total_t)
# Flatten the 2D array of axes for easier iteration
axs = axs.flatten()

for i in range(4):  # Assuming only 4 trials
    trial_num = trial_numbers[i]
    
    # Load alpha_m data
    alpha_m_file_path = f"{data_path}/trial_{trial_num}/{alpha_m_filename}"
    alpha_m_list = np.loadtxt(alpha_m_file_path)
    alpha_m_values.append(np.copy(alpha_m_list))
    
    # Load turbulent velocity data
    turb_vel_file_path = f"{data_path}/trial_{trial_num}/{turb_vel_filename}"
    with open(turb_vel_file_path, 'r') as f:
        lines = f.readlines()
    turb_vel = np.array(lines, dtype=float)
    turb_vel_values.append(turb_vel)

    # Plot alpha_m on the corresponding subplot
    axs[i].plot(z_lists[i], alpha_m_list[t_index], label=f'alpha_m: {labels[i]}')
    axs[i].axhline(y=0, color='grey', linestyle='--')
    axs[i].axvline(x=0, color='grey', linestyle='--')

    # Plot turbulent velocity on the same subplot
    axs[i].plot(z_lists[i], turb_vel, linestyle='--', color='black', alpha=0.8, label='RMS turb velocity')

    # Set y-limits based on the subplot index
    if i == 0:
        axs[i].set_ylim(-100, 100)
        axs[i].set_xlim(-1, 1)
    else:
        axs[i].set_ylim(-40, 40)
        axs[i].set_xlim(-1, 1)

    # Add labels, title, and legend
    axs[i].set_xlabel(r'$z$')
    axs[i].set_ylabel(r'$\alpha_m$ / RMS turb velocity')
    axs[i].set_title(f'Trial {trial_num}')
    axs[i].legend(loc='upper left')

# Adjust the layout so subplots do not overlap
plt.tight_layout()

# Save the entire figure as a single image
plt.savefig(f'{fig_path}/alpha_m_and_turb_vel_vs_z_subplots.png', bbox_inches='tight')
plt.close()

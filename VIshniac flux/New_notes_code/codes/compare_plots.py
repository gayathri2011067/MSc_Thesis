import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os, sys
import subprocess

# Paths to data and figures
data_path = "/home/gayathri/MSc_thesis/VIshniac flux/New_notes_code/data_files"
fig_path = "/home/gayathri/MSc_thesis/VIshniac flux/New_notes_code/figures"

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
plt.figure(figsize=(10, 5))  # Set figure size to rectangular
for i in range(len(trial_numbers)):
    plt.plot(times[i], B_strengths[i][:, space_indices[i]], label=f'{labels[i]}')

plt.xlim(0, 14)
plt.ylim(0, 0.6)
plt.xlabel(r'$t$')
plt.ylabel(r'$\mathrm{B}_{\mathrm{strength}}$')  # Replacing \text with \mathrm
plt.title(r'$\mathrm{B}_{\mathrm{strength}} \ \mathrm{vs} \ t \ \mathrm{at} \ z=(' +
          ', '.join([f'trial{trial_numbers[i]}:{int(z_lists[i][space_indices[i]])}' for i in range(len(trial_numbers))]) + r')$')
# plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(f'{fig_path}/B_strength_vs_time.png', bbox_inches='tight')
plt.close()
B_strength_avg = np.mean(B_strengths, axis=2)

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
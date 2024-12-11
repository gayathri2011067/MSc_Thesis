import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import os
import subprocess
from plotly.subplots import make_subplots

data_path = "/home/gayathri/MSc_thesis/VIshniac flux/New_notes_code/run_files"
fig_path = "/home/gayathri/MSc_thesis/VIshniac flux/New_notes_code/figures"
data_save_path = "/home/gayathri/MSc_thesis/VIshniac flux/New_notes_code/data_files"

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
    # print(Br_list[1])
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
plt.plot(time_list, Br_list[:,25], label='Br')
# print(time_list)
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

# Plot B_strength at a specific z index
plt.plot(time_list, B_strength[:,51])
plt.xlabel('time')

# plt.ylim(0,0.6)
plt.ylabel('B_strength/B0')
plt.title('B_strength vs time')
plt.yscale('log')
# plt.xlim(0,7)

# plt.axhline(y=0.001, color='g', linestyle='--')
# plt.axvline(x=2.5, color='r', linestyle='--')
plt.savefig(f'{fig_path}/B_strength_vs_time.png')
plt.close()

# Plot field strength averaged over z as a function of time
B_strength_avg = np.sqrt(np.mean(B_strength**2, axis=1))
np.savetxt(f'{data_save_path}/B_strength_avg.txt', B_strength_avg)
plt.plot(time_list, B_strength_avg)
plt.yscale('log')
plt.xlabel('time')
plt.xlim(0,10)
plt.ylabel('Average B_strength/B0')
plt.title('Average B_strength vs time')
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
    n1= 5000
    total_t= 25
    t_val = 8.8
    t_index = int(t_val * n1 / total_t)

    plt.plot(z, alpha_m[t_index - 1], label=f't={t_val:.2f}')  # t_index - 1 because the index starts from 0
    plt.xlabel('z')
    plt.ylabel('alpha_m')
    plt.title(f'alpha_m vs z at t={t_val:.2f}')
    plt.legend()
    plt.savefig(f'{fig_path}/alpha_m_vs_z_t_{t_val:.2f}.png')
    plt.close()
except:
    print('alpha_m.txt file not found')
    #i have br_final in (t,z) array. want to plot a 2d heat map with z as y axis and time as x axis.
    # Plot 2D heatmap of Br_final with time on x-axis and z on y-axis
import plotly.io as pio
pio.kaleido.scope.plotlyjs = "https://cdn.plot.ly/plotly-latest.min.js"

import plotly.graph_objects as go



# Function to create contour plots
# def create_contour_plot(z_data, time_list, z, title, file_name):
#     fig = go.Figure()

#     # Add contour lines
#     fig.add_trace(go.Contour(
#         z=z_data.T,
#         x=time_list,
#         y=z,
#         colorscale='inferno',
#         showscale=True,
#         contours=dict(
#             start=z_data.min(),
#             end=z_data.max(),
#             size=(z_data.max() - z_data.min()) / 50,
#             # Uncomment the next line for only contour lines without fills
#             coloring='lines',
#             line=dict(color='black')  # Uncomment this to set contour lines to black
#         )
#     ))

#     # Update layout
#     fig.update_layout(
#         title=f'<b>{title}</b>',
#         title_font=dict(size=20),
#         xaxis_title='<b>Time</b>',
#         yaxis_title='<b>z</b>',
#         xaxis=dict(range=[9, 10]),
#         yaxis=dict(range=[-0.5, 0.5]),
#         width=800,  # Adjust width as needed
#         height=400,  # Adjust height to maintain rectangular aspect ratio
#         margin=dict(l=50, r=50, t=50, b=50),
#         font=dict(size=14),
#         template='plotly_white'  # Choose a clean template
#     )

#     # Save the figure
#     fig.write_image(f'{fig_path}/{file_name}')

# # Create contour plot for Br
# create_contour_plot(Br_list, time_list, z, '2D Heatmap of $B_r$ with Contour Lines', 'Br_final_heatmap_with_contours.png')

# # Create contour plot for Bphi
# create_contour_plot(Bphi_list, time_list, z, '2D Heatmap of $B_\\phi$ with Contour Lines', 'Bphi_final_heatmap_with_contours.png')

# # Create contour plot for B strength
# create_contour_plot(B_strength, time_list, z, '2D Heatmap of $B$ with Contour Lines', 'B_heatmap_with_contours.png')
# Function to create heatmap plots
# def create_heatmap_plot(z_data, time_list, z, title, file_name):
#     fig = go.Figure()

#     # Add heatmap
#     fig.add_trace(go.Heatmap(
#         z=z_data.T,
#         x=time_list,
#         y=z,
#         colorscale='inferno',
#         showscale=True
#     ))

#     # Add contour lines
#     fig.add_trace(go.Contour(
#         z=z_data.T,
#         x=time_list,
#         y=z,
#         showscale=False,
#         contours=dict(
#             coloring='lines',
#             # showlabels=True,  # Show labels on contour lines
#             start=z_data.min(),
#             end=z_data.max(),
#             size=(z_data.max() - z_data.min()) / 100  # More contours
#         ),
#         line=dict(color='black')
#     ))

#     # Update layout
#     fig.update_layout(
#         title=f'<b>{title}</b>',
#         title_font=dict(size=20),
#         xaxis_title='<b>Time</b>',
#         yaxis_title='<b>z</b>',
#         xaxis=dict(range=[9, 10]),
#         yaxis=dict(range=[-0.5, 0.5]),
#         width=800,  # Adjust width as needed
#         height=400,  # Adjust height to maintain rectangular aspect ratio
#         margin=dict(l=50, r=50, t=50, b=50),
#         font=dict(size=14),
#         template='plotly_white'  # Choose a clean template
#     )

#     # Save the figure
#     fig.write_image(f'{fig_path}/{file_name}')

# # Create heatmap plot for Br
# create_heatmap_plot(Br_list, time_list, z, '2D Heatmap of $B_r$', 'Br_final_heatmap.png')

# # Create heatmap plot for Bphi
# create_heatmap_plot(Bphi_list, time_list, z, '2D Heatmap of $B_\\phi$', 'Bphi_final_heatmap.png')

# # Create heatmap plot for B strength
# create_heatmap_plot(B_strength, time_list, z, '2D Heatmap of $B$', 'B_heatmap.png')

# # Create a combined plot with Br, Bphi, and B strength
# fig_combined = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.1,
#                              subplot_titles=('2D Heatmap of $B_r$', '2D Heatmap of $B_\\phi$', '2D Heatmap of $B$'))

# # Add Br heatmap
# fig_combined.add_trace(go.Heatmap(
#     z=Br_list.T,
#     x=time_list,
#     y=z,
#     colorscale='electric',
#     showscale=False
# ), row=1, col=1)

# # Add Bphi heatmap
# fig_combined.add_trace(go.Heatmap(
#     z=Bphi_list.T,
#     x=time_list,
#     y=z,
#     colorscale='electric',
#     showscale=False
# ), row=2, col=1)

# # Add B strength heatmap
# fig_combined.add_trace(go.Heatmap(
#     z=B_strength.T,
#     x=time_list,
#     y=z,
#     colorscale='electric',
#     showscale=True
# ), row=3, col=1)

# # Update layout
# fig_combined.update_layout(
#     height=1200,  # Adjust height as needed
#     width=800,  # Adjust width as needed
#     title_text='<b>Combined 2D Heatmaps of $B_r$, $B_\\phi$, and $B$</b>',
#     title_font=dict(size=20),
#     xaxis3_title='<b>Time</b>',
#     yaxis1_title='<b>z</b>',
#     yaxis2_title='<b>z</b>',
#     yaxis3_title='<b>z</b>',
#     font=dict(size=14),
#     template='plotly_white',  # Choose a clean template
#     xaxis=dict(range=[9, 10]),  # Set x-axis limit
#     yaxis=dict(range=[-0.5, 0.5]),  # Set y-axis limit for the first subplot
#     yaxis2=dict(range=[-0.5, 0.5]),  # Set y-axis limit for the second subplot
#     yaxis3=dict(range=[-0.5, 0.5])  # Set y-axis limit for the third subplot
# )

# # Add contour lines with sharp contrast
# for i, z_data in enumerate([Br_list, Bphi_list, B_strength], start=1):
#     fig_combined.add_trace(go.Contour(
#         z=z_data.T,
#         x=time_list,
#         y=z,
#         showscale=False,
#         contours=dict(
#             coloring='lines',
#             start=z_data.min(),
#             end=z_data.max(),
#             size=(z_data.max() - z_data.min()) / 20  # Fewer contours for sharper contrast
#         ),
#         line=dict(color='black')
#     ), row=i, col=1)

# # Save the combined figure
# fig_combined.write_image(f'{fig_path}/combined_heatmap.png')
#check if alpha_m is nan anywhere,print the index
if np.isnan(alpha_m).any():
    print('alpha_m has nan values')
    print(np.argwhere(np.isnan(alpha_m)))
# Plot alpha_m vs time
try:
    # Ensure alpha_m is loaded
    if 'alpha_m' not in locals():
        filename = 'alpham_final.txt'
        file_path = f"{data_path}/{filename}"
        with open(file_path, 'r') as f:
            lines = f.readlines()
        alpha_list = []
        for line in lines:
            line = line.strip()
            line = line.split()
            curr = np.array(line, dtype=float)
            alpha_list.append(curr)
        alpha_m = np.array(alpha_list, dtype=float)
        np.savetxt(f'{data_save_path}/alpha_m.txt', alpha_m)

    # Plot alpha_m vs time
    plt.plot(time_list, alpha_m[:, 25])  # Assuming alpha_m has the same shape as Br_list
    plt.xlabel('time')
    plt.ylabel('alpha_m')
    plt.title('alpha_m vs time')
    plt.savefig(f'{fig_path}/alpha_m_vs_time.png')
    plt.close()
except Exception as e:
    print(f'Error plotting alpha_m vs time: {e}')

#NOTE: - Rk finite case
#      - 2D heatmaps of Br, Bphi, and B strength with contour lines
#      - Peak discrepancies in Br, Bphi, and B strength
#      - averaging over z
#      - alpha_m vs z at a specific time
#      - dimensions
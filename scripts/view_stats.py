import numpy as np
import matplotlib.pyplot as plt
import sys

# define metrics
labels = ['mcc', 'sens', 'spec', 'accu', 'prec', 'FDR', 'F1']
num_vars = len(labels)

# read stats from file
def read_stats(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    # Extract the values for RM2, EG, and EDTA
    rm2_values = list(map(float, lines[1].strip().split(',')))
    eg_values = list(map(float, lines[4].strip().split(',')))
    edta_values = list(map(float, lines[7].strip().split(',')))
    
    return rm2_values, eg_values, edta_values

# read the stats data from the stats.txt file
stats_file = sys.argv[1]
rm2_values, eg_values, edta_values = read_stats(stats_file)

ideal_values = [1, 1, 1, 1, 0, 1]

# compute angle for each axis of the hexagon
angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

# close the loop (the last point is the same as the first)
rm2_values += rm2_values[:1]
eg_values += eg_values[:1]
edta_values += edta_values[:1]
ideal_values += ideal_values[:1]
angles += angles[:1]

# set up the figure and axis
fig, ax = plt.subplots(figsize=(6, 6), dpi=100, subplot_kw=dict(polar=True))

# plot the data
ax.plot(angles, rm2_values, color='b', linewidth=2, linestyle='solid', label='RM2')
ax.plot(angles, eg_values, color='g', linewidth=2, linestyle='solid', label='Earl Grey')
ax.plot(angles, edta_values, color='r', linewidth=2, linestyle='solid', label='EDTA')
#ax.plot(angles, ideal_values, color='r', linewidth=2, linestyle='solid')

# set the labels on the hexagon
ax.set_xticks(angles[:-1])
ax.set_xticklabels(labels, fontsize=12)

ax.set_ylim(0, 1)
ax.set_yticks([0.0, 0.3, 0.5, 0.8, 1.0])
ax.set_yticklabels(['0.0', '0.3', '0.5', '0.8', '1.0'], fontsize=10)

# set the range (0 to 1 for each metric)
ax.set_ylim(0, 1)
model=sys.argv[3]
coverage=sys.argv[4]
ax.set_title(f"{coverage} Metric Scores for {model}")
ax.legend(loc='center', bbox_to_anchor=(.2, .2), fontsize=12)

save_file = sys.argv[2]
plt.savefig(save_file)

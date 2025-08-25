import sys
import os
import re
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np


#define peak IO
# DB read 210GB, peak disk IO 2.6GBps
# 210/2.6 = 81 sec
peak_io_time = 81

time_to_files_dict = {
    1: "NP_670887.fa",
    2: "YP_244862.fa",
    3: "NP_308413.fa",
    4: "NP_931067.fa",
    5: "YP_182226.fa",
    6: "YP_054627.fa",
    7: "YP_209371.fa",
    8: "YP_141079.fa",
    9: "NP_436070.fa",
    10: "YP_000338.fa",
    11: "NP_349135.fa",
    12: "NP_566269.fa",
    13: "NP_894055.fa",
    14: "YP_137100.fa",
    15: "YP_062020.fa",
    16: "NP_864591.fa",
    17: "NP_214548.fa",
    18: "YP_069013.fa",
    19: "NP_872878.fa",
    20: "NP_345015.fa",
    21: "NP_085504.fa",
    22: "NP_174344.fa",
    23: "NP_611552.fa",
    24: "XP_610726.fa",
    25: "YP_189837.fa",
    26: "NP_245243.fa",
    27: "NP_818972.fa",
    28: "NP_633288.fa",
    29: "NP_781168.fa",
    30: "NP_784621.fa",
    31: "NP_952492.fa",
    32: "YP_121326.fa",
    33: "NP_965272.fa",
    34: "NP_706440.fa",
    35: "XP_499151.fa",
    36: "XP_237451.fa",
    37: "NP_733286.fa",
    38: "YP_129633.fa",
    39: "YP_081219.fa",
}

# convert time in "real 40m10.280s" format to seconds
def time_to_seconds(time_str):
    match = re.match(r'real\s+(\d+)m(\d+\.\d+)s', time_str)
    if match:
        minutes = int(match.group(1))
        seconds = float(match.group(2))
        return minutes * 60 + seconds
    return 0

dirname = sys.argv[1]

if not dirname.endswith("/"):
    dirname += "/"

#Ensure that dir is valid
if not(os.path.isdir(dirname)):
    print("perf dir is not valid.")

file_list = os.listdir(dirname)

files_starting_with_number = [f for f in file_list if re.match(r'^\d', f)]

#Throw an error if there are any missing files.
for i in range(1, 40):
    if (f"{i}_time_blaze" not in files_starting_with_number or
        f"{i}_time_st" not in files_starting_with_number or
        f"{i}_time_mt" not in files_starting_with_number):
        print(f"[ERROR] Missing time file for query {i}")
        exit(1)

time_dict = {}

for i in range(1, 40):
    exec_time = {}
    #Get query name from i
    query_name = time_to_files_dict.get(i)
    q_name = query_name.split('.')[0]

    for record in SeqIO.parse(dirname+query_name, "fasta"):
        exec_time['seq_length'] = len(record.seq)

    with open(dirname + f"{i}_time_blaze", 'r') as file:
        exec_time_str = next(line for line in file if line.startswith('real'))
        exec_time['blaze'] = (time_to_seconds(exec_time_str))
    
    with open(dirname + f"{i}_time_st", 'r') as file:
        exec_time_str = next(line for line in file if line.startswith('real'))
        exec_time['st'] = (time_to_seconds(exec_time_str))

    with open(dirname + f"{i}_time_mt", 'r') as file:
        exec_time_str = next(line for line in file if line.startswith('real'))
        exec_time['mt'] = (time_to_seconds(exec_time_str))

    exec_time['peak'] = peak_io_time

    time_dict[q_name] = exec_time

time_dict = dict(sorted(time_dict.items(), key=lambda item: item[1]['seq_length']))

queries = list(time_dict.keys())
st_time = [time_dict[query]['st'] for query in queries]
mt_time = [time_dict[query]['mt'] for query in queries]
peak_time = [time_dict[query]['peak'] for query in queries]
blaze_time = [time_dict[query]['blaze'] for query in queries]

peak_speedup = [cpu / peak for cpu, peak in zip(st_time, peak_time)]
multi_speedup = [cpu / multi for cpu, multi in zip(st_time, mt_time)]
blaze_speedup = [cpu / blaze for cpu, blaze in zip(st_time, blaze_time)]


#Plot graph of total time taken by each query in a plt
fig, axs = plt.subplots(1, 1, figsize=(12, 4))
# fig, axs = plt.subplots(1, 1)

axs.set_axisbelow(True)

#Plot the bars next to each other
barWidth = 0.3

# Set position of bar on X axis
r1 = np.arange(len(queries))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

# Plot the data
axs.bar(r1, peak_speedup, width = barWidth, label='Peak Speedup (Disk I/O)', color='orange')
axs.bar(r2, multi_speedup, width = barWidth, label='Multithreaded speedup', color='blue')
axs.bar(r3, blaze_speedup, width = barWidth, label='BLAZE', color='green')

axs.bar(len(queries), 0, width = barWidth, label='', color='orange')
axs.bar(len(queries)+barWidth, 0, width = barWidth, label='', color='blue')
axs.bar(len(queries)+2*barWidth, 0, width = barWidth, label='', color='green')

axs.bar(len(queries)+3, np.mean(peak_speedup), width = barWidth, color='orange')
axs.bar(len(queries)+3+barWidth, np.mean(multi_speedup), width = barWidth, color='blue')
axs.bar(len(queries)+3+2*barWidth, np.mean(blaze_speedup), width = barWidth, color='green')

axs.set_xticks([r + barWidth for r in range(len(queries))] + [barWidth+len(queries)+3])

group_labels = [str(i+1) for i in range(len(queries))] + ['Mean']

axs.set_xticklabels(group_labels)

axs.yaxis.grid(True, alpha=0.7)

axs.legend(prop={'size': 10})

axs.set_ylabel('Speedup', fontsize=12)
axs.set_xlabel('Benchmark', fontsize=12)
axs.set_title('Speedup of Peak possible, Multi-threaded and BLAZE over CPU baseline', fontsize=12)

plt.tight_layout(pad=0.1)
plt.savefig("main_blaze.png",bbox_inches='tight',dpi=300.0)

print("Generated graph and saved as main_blaze.png")
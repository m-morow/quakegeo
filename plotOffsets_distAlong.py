import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import numpy as np

"""
Plot offset measurements for Feb. 2023 Turkey earthquake
"""

project = 'optical_offset_measurements_b_UTM' #data name
df = pd.read_csv('/Users/mtan/Documents/Data/Turkey/08102023_testing/{}.csv'.format(project), index_col=None) #file path

#get rid of preferred =1000m, AKA piercing points that were identified but ultimately not measured
#get rid of partial =Yes, AKA piercing points that are not representative of full offset
df_adj = df[(df['preferred'] != 1000) & (df['partial'] == 'No')]

#identify Narli and EAF piercing points, assign colors
faults = list(df_adj['feature'])
color=[]
for i in range(0, len(faults)):
    if faults[i] == 'EAF':
        color.append('orange') 
    else:
        color.append('g') 

x = list(df_adj['dist_along']/1000) #convert to km
y_pref = list(df_adj['preferred'])

#calculate min/max values using errorbar module
y_lower = list(df_adj['min'])
y_upper = list(df_adj['max'])
y_l = (np.array(y_pref) - np.array(y_lower))
y_u = (np.array(y_upper) - np.array(y_pref))

#get rid of very large (~ +/-1000, indicative of unmeasured piercing points)
y_l_adj = [0 if (l > 100 or l < 0) else l for l in y_l] 
y_u_adj = [0 if (u > 100 or u < 0) else u for u in y_u]
errors = [y_l_adj, y_u_adj] #errorbar module format

#plot preferred values and corresponding fault color
for i in range(len(faults)):
    plt.scatter(x[i], y_pref[i], c=color[i], s=30, edgecolors='k', zorder=5)

#plot min/max values as bars
plt.errorbar(x, y_pref, yerr=errors, fmt="ob", ecolor='k', zorder=0)

#set up and plot legend
EAF = mlines.Line2D([], [], color='orange', marker='o', markersize=5, linestyle="None", label='East Anatolian Fault')
NF = mlines.Line2D([], [], color='green', marker='o', markersize=5, linestyle="None", label='Narli Fault')
plt.legend(handles=[EAF, NF], loc=2)

#plt.axhline(y=0)
plt.title("Offset measurements")
plt.xlabel("Distance along fault (km)")
plt.ylabel("Offset (m)")
plt.show()
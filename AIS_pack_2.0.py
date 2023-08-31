# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:30:00 2023

@author: yuhao
"""

#import plotly.graph_objects as go
import os
import numpy as np
import pandas as pd
import scipy.signal as signal
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

path=r"C:\Users\yuhao\Documents\Python Scripts"
os.chdir(path)

###to list all the files####
#for files in path:
    #if files endswith "csv":
        #f_lists=os.filelist(files)
###read the csv files for length measurment 
df = pd.read_csv("TRIM46_length.csv", usecols = ['gray_value'])
GV=df['gray_value']
###smooth the length profile by 4 ums
window_size = 80
GV_smo = np.convolve(GV, np.ones(window_size)/window_size, mode='same')
print("Max value",GV_smo.max())
print("mini value",GV_smo.min())
he=0.5*GV_smo.max()
###now detect peak from the length profile and extract the index and value 
indices = find_peaks(GV_smo, height=he, threshold=0.03, distance=1)[0]
peak_value = [GV_smo[j] for j in indices]
#print(indices)
#print(GV)
print("peak value is: ", peak_value)
#x=[i for i in range(len(GV))]
#y=[GV[j] for j in range(len(GV))]

for i in range(len(indices)):
    plt.annotate("", (indices[i], peak_value[i]), 
                 textcoords="offset points", xytext=(0,10), 
                 arrowprops=dict(facecolor='red', 
                                 edgecolor='none', shrink=0.01))

###find left edge and right edge of the length profile and store info in dic

for i, peak_idx in enumerate(indices):
    print(peak_idx)
    peak_value = GV_smo[peak_idx]
    thr_v = peak_value * 0.4
    print(thr_v)
    left_idx = peak_idx - 1
    print(left_idx)
    while left_idx >= 0 and abs(GV_smo[left_idx]) >= thr_v: 
        left_idx -= 1
    Lva=GV_smo[left_idx]
    
    right_idx = peak_idx + 1
    while right_idx < len(GV_smo) and abs(GV_smo[right_idx])>= thr_v:
        right_idx += 1
    Rva=GV_smo[right_idx]
    band_width=(right_idx-left_idx)
    peak_dict = {
                'peak_value': peak_value,
                'left_index': left_idx+1,
                'right_index': right_idx+1,
                'left_value': Lva,
                'right_value': Rva,
                'peak_width': band_width
            }
            
###left index would be the start of AIS and right index would be the end of AIS    
###the peak width would be the length of the AIS

#for li in peak_dict["left_index"]:
li=peak_dict["left_index"] 
ri=peak_dict["right_index"]
#print(GV_smo[li])
plt.plot(GV_smo)
plt.xlabel('Length (um)')
plt.ylabel('Gray Value')
plt.title('peak detection')
plt.annotate("", (li, GV_smo[li]), textcoords="offset points", 
             xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none', 
                                            shrink=0.01))
plt.annotate("", (ri, GV_smo[ri]), textcoords="offset points", 
             xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none', 
                                            shrink=0.01))
plt.show()
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 16:29:56 2023

@author: yuhao
"""

import os
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
#import plotly.graph_objects as go
import matplotlib.pyplot as plt

# Prompt the user for the input directory and cutoff value
input_path = ""      #input("Enter the input directory: ")
os.chdir(input_path)
cutoff =0.8        #float(input("Enter the cutoff value: "))   #cutoff for spacing calculation
X = 20              #float(input("value of your X axes"))  #pixel size
ave_ring_dis=[]
k_w="Value"         #value being extracted from csv file
# Prompt the user for the output directory and check if the directory exists
output_path = ""          #input("Enter the output directory: ")
os.makedirs(output_path, exist_ok=True)

# Loop through each CSV file in the directory
for files in os.listdir(input_path):
    # Check if the file is a CSV file
    if files.endswith(".csv"):
        # Read the data from the CSV file and extract gray value
        df = pd.read_csv(files, usecols=[k_w])
        GV = df[k_w]
        window_size = 5     #how many pixels for averaging
        GV_smo = np.convolve(GV, np.ones(window_size)/window_size, mode='same')
        
        he = 0.2 * GV_smo.max()
        
        # Find the peak indices and distances between peaks
        indices = find_peaks(GV_smo, height=he, threshold=0.03, distance=1)[0]
        peak_value = [GV_smo[j] for j in indices]
        peak_distances = []
        plt.plot(GV_smo)
        plt.xlabel('time (ms)')
        plt.ylabel('Gray Value')
        plt.title('peak detection')

        for i in range(len(indices)):
            plt.annotate("", (indices[i], peak_value[i]), 
            textcoords="offset points", xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none', shrink=0.01))
        #plt.savefig()
        plt.show()
        
        
        
        
        for i in range(len(indices) - 1):
            peak_dis = indices[i+1] - indices[i]
            peak_distances.append(peak_dis*X)   #convert pixel to microns
        avee=sum(peak_distances)/len(peak_distances)
        
        # Initialize a list to store the peak information
        peak_info = []
        left_edge=[]
        right_edge=[]

        # Iterate over the peak indices and calculate the peak values and left/right values
        for i, peak_idx in enumerate(indices):
            peak_value = GV_smo[peak_idx]
            thr_v = peak_value * cutoff
            left_idx = peak_idx - 1
            while left_idx >= 0 and abs(GV_smo[left_idx]- thr_v) >= abs(GV_smo[left_idx + 1] - thr_v):
                left_idx -= 1
            Lva=GV_smo[left_idx]
            left_edge.append(left_idx)
            
            right_idx = peak_idx + 1
            while right_idx < len(GV_smo) and abs(GV_smo[right_idx]- thr_v) >= abs(GV_smo[right_idx - 1] - thr_v):
                right_idx += 1
            Rva=GV_smo[right_idx]
            right_edge.append(right_idx)
            band_width=(right_idx-left_idx)*X      
            peak_dict = {
                'peak_value': peak_value,
                'peak_index': peak_idx+1,
                'left_index': left_idx+1,
                'right_index': right_idx+1,
                'left_value': Lva,
                'right_value': Rva,
                'peak_width': band_width
            }
            
            # Calculate the distance to the next peak (if there is one)
        
            #if i < len(indices) - 1:
            #    next_peak_idx = indices[i + 1]
            #    next_peak_value = GV_smo[next_peak_idx]
            #    next_left_idx = next_peak_idx - 1
            #    while next_left_idx >= 0 and GV_smo[next_left_idx] >= next_peak_value:
            #        next_left_idx -= 1
            #    next_left_idx += 1
            #    distance = (next_left_idx - right_idx)*X
            #    peak_dict['distance (nm)'] = distance

            # Append the dictionary to the peak_info list
            peak_info.append(peak_dict)
        # Save the peak information to an Excel file with the same name as the CSV file
        excel_file = os.path.splitext(files)[0] + '.xlsx'
        excel_path = os.path.join(output_path, excel_file)
        df = pd.DataFrame(data=peak_info)
        df.to_excel(excel_path, index=False)
        
        #ARD=[]
        #RD=0
        #for d in range(len(peak_info)-1):
        #    ARD.append(peak_info[d]["distance (nm)"])
        #    RD=sum(ARD)/len(ARD)
        #ring_distance.append(RD)
        #actin_ring_distance["dis"]=ring_distance
        #ARD=[]
        #rd=0
        #ave_ad=0
        #for indx in range(len(left_edge)):
    
        #    if indx < len(left_edge)-1:
        #rd=0
        #rd=right_edge(indx+1)-left_edge(indx)
        #ring_dis.append(rd)
        #print(indx)
        #print(indx, left_edge[indx+1], right_edge[indx])
        #        rd=(left_edge[indx+1]-right_edge[indx])*20
        #        ARD.append(rd)
        #ave_ad=sum(ARD)/len(ARD)
        #ring_distance.append(ave_ad)
        #ARD=[]
        #rd=0
        #ave_ad=0
    ave_ring_dis.append(avee)
fn="Average actin ring distance.xlsx"
#f_path=os.path.join(output_path, fn)
f_df=pd.DataFrame(ave_ring_dis)
f_df.to_excel(output_path+fn, index=False)
            
        

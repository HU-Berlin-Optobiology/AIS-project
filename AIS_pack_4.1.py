# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:59:35 2023
## Update: This is the latest version of AIS_pack. some naming issues are fixed
@author: yuhao
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import math
from itertools import zip_longest
import re

#go to the path where all profiles are stored
path = r"C:\Users\yuhao.HP640G8MIKH\Desktop\Final results revision\final results_AIS_plasticity\rapid\nonAcD\TTX"
save_exe=r"C:\Users\yuhao.HP640G8MIKH\Desktop\Final results revision\final results_AIS_plasticity\rapid\nonAcD"
os.chdir(path) #go to the directory
k_w="AnkG" #set the key word in your file name
window_size = 86 #set the distance for normalising the intensity profile
cut_off=0.4 #set the limit for length definition 
exp_name="TTX"


f_roots=[] #creat a list to store all roots with csv files
print("Checking files...")
nu_files=0 ## check  number of files
#loop over all files and roots in the directory
for roots, directory, files in os.walk(path):
    for f in files:        #focus on the files
        if f.endswith(".csv") and k_w in f: 
            #print(f)
            nu_files+=1
            
            if not roots in f_roots:
                f_roots.append(roots)  #add correct roots into the list
print("Number of csv files detected: ", str(nu_files))

AIS_ana={}
AIS_profile={}                

for rs in f_roots:
    print("Going to directory: ",rs)
    os.chdir(rs)
    file_list = os.listdir(rs)
    AIS_ana[rs]={}
    AIS_profile[rs]={}
    for fi in file_list:
        if fi.endswith(".csv") and k_w in fi:
            filename=str(fi)
            
            f_name=filename.split(".lif")[0]
            print("Processing file: ", filename)
# Read the CSV file for length measurement
            print("Loading csv file...")
            df = pd.read_csv(fi, usecols=['Value'])

# Smooth the length profile by 4 ums using convolution
            print("Signal normalisation...")
            smoothed_signal = np.convolve(df['Value'], np.ones(window_size) / window_size, mode='same')

# Normalize the smoothed signal values
            normalized_signal = smoothed_signal / smoothed_signal.max()

# Convert the normalized values to a numpy array
            GV_len = np.array(normalized_signal)
# append the normalised profile into a dictionary 
            AIS_profile[rs][filename]=GV_len
# extract the index of values in normalised AnkG profile
            indices=np.arange(len(GV_len))  
# extract the index of max value in normalised AnkG profile            
            max_index = np.argmax(GV_len)
# lock the max value and also set the threshold for detecting the border of AIS
            GV_max=GV_len.max()
            thr_v = GV_max * cut_off
            
            print("Index of max value:", max_index)
            print("threshold for length determination:", thr_v)
            
            print("Ploting intensity profile")
            plt.plot(GV_len)
            plt.xlabel('Length (pixel)')
            plt.ylabel('Gray Value (norm)')
            plt.title(filename)   
            plt.annotate("max", (max_index,GV_len[max_index]), textcoords="offset points", 
                         xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none',
                                                         shrink=0.01))
            #plt.show()        
# split the AnkG profile into two parts based on max value
            
            left_of_max = GV_len[:max_index][::-1]  # Reversed to start from max
            right_of_max = GV_len[max_index:]
# initialise the left and right border index and value
            left_index_below_threshold = 0
            right_index_below_threshold = 0
            Lva= None
            Rva= None
# loop over the velues on left and right side of max value
# and detect the border of AIS       
            for i, value in enumerate(left_of_max):

                if value < thr_v:
                    #Store the value
                    Lva=value
        # Store the index considering that left_of_max is reversed
                    left_index_below_threshold = max_index - i
                    break

# Find the index of the first value below 40% on the right
            for i, value in enumerate(right_of_max):
                if value < thr_v:
                    Rva=value
        # Store the index
                    right_index_below_threshold = max_index + i
                    break

# Return the results
            #print(f"The index of the first value below 40% to the left of max: {left_index_below_threshold}")
            #print(f"The index of the first value below 40% to the right of max: {right_index_below_threshold}")   
# using the right index minus left index to calculate the length of AIS
            
            # if left_index_below_threshold or right_index_below_threshold == None:
            #     print("Error!!", filename)
            print("left: ",left_index_below_threshold)
            print("right: ",right_index_below_threshold)
            #     print("#####")
            #     continue

            plt.annotate("left edge", (left_index_below_threshold,GV_len[left_index_below_threshold]), textcoords="offset points", 
                         xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none',
                                                         shrink=0.01))
            
            plt.annotate("right edge", (right_index_below_threshold,GV_len[right_index_below_threshold]), textcoords="offset points", 
                         xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none',
                                                         shrink=0.01))
            plt.show() 
            
            band_width= right_index_below_threshold - left_index_below_threshold
            print("AIS length: ", band_width)
# store results into dictionary  
            peak_dict = {
                
                'peak_value': GV_max,
                'AIS_start': left_index_below_threshold,
                'AIS_end': right_index_below_threshold,
                'left_value': Lva,
                'right_value': Rva,
                'AIS_length': band_width,
                'AIS_peak': max_index
                }
            
            AIS_ana[rs][filename]=peak_dict
            
            
#saving results as excel   
#writer=pd.ExcelWriter(save_exe+"\\"+exp_name+"_AIS_ana.xlsx", engine='openpyxl')
data_for_excel = []
for root, files in AIS_ana.items():
    for file, values in files.items():
        row = {'Root': root, 'File': file}
        row.update(values)
        data_for_excel.append(row)

# Convert the list of dictionaries to a DataFrame
df = pd.DataFrame(data_for_excel)

# Specify the file path for the Excel file
excel_file_path = os.path.join(save_exe, exp_name + "_AIS_ana.xlsx")

# Create an Excel writer object
with pd.ExcelWriter(excel_file_path, engine='openpyxl') as writer:
    # Write the DataFrame to an Excel sheet
    df.to_excel(writer, sheet_name='AIS_KCl', index=False)

print(f"Data saved to {excel_file_path}")


#plot all profiles to decide for cut-off
# for k in AIS_profile.keys():
#     print(k)
#     for p in AIS_profile[k].keys():
#         #print(p)
#         #print(AIS_ana[k][p]["AIS_start"])
#         plt.plot(AIS_profile[k][p])
#         li=AIS_ana[k][p]["AIS_start"]
#         ri=AIS_ana[k][p]["AIS_end"]
#         pi=AIS_ana[k][p]["AIS_peak"]
#         plt.annotate("", (li, AIS_profile[k][p][li]), textcoords="offset points", 
#                      xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none',
#                                                      shrink=0.01))
#         plt.annotate("", (ri, AIS_profile[k][p][ri]), textcoords="offset points", 
#                      xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none', 
#                                                      shrink=0.01))
#         plt.annotate("", (pi, AIS_profile[k][p][pi]), textcoords="offset points",
#                      xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none', 
#                                                      shrink=0.01))
# plt.xlabel('Length (pixel)')
# plt.ylabel('Gray Value')
# plt.title('peak detection')        
# plt.show()        



            # li=peak_dict["AIS_start"]
            # ri=peak_dict["AIS_end"]
            # pi=peak_dict["AIS_peak"]
            # plt.plot(GV_len)
            # plt.xlabel('Length (pixel)')
            # plt.ylabel('Gray Value')
            # plt.title('peak detection')
            # plt.annotate("", (li, GV_len[li]), textcoords="offset points", 
            #              xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none',
            #                                             shrink=0.01))
            # plt.annotate("", (ri, GV_len[ri]), textcoords="offset points", 
            #              xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none', 
            #                                             shrink=0.01))
            # plt.annotate("", (pi, GV_len[pi]), textcoords="offset points", 
            #              xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none', 
            #                                             shrink=0.01))




# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 16:01:37 2023

@author: yuhao
"""

import os
import glob
import re
import csv
import math
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
from itertools import zip_longest
from scipy.signal import find_peaks

### user input ###

file_path=r"E:\Yuhao_do not delet\images_tem\AIS_plasticity_Synaptopodin\20231018_AIS plasticity_synpo_DIV12_ctrl_02\AnkG_Synpo_profiles"   ## type in file path
os.chdir(file_path)  ## set directory to the input path
os.getcwd()  ## check ##



# ankg_pattern = re.compile(r'rt_AnkG_ROI_exp3_AnkG_MAX_AIS_chronic pla_ctrl_(\d+).lif.tif.csv') 
# synpo_pattern = re.compile(r'rt_Synpo_ROI_exp3_Synpo_MAX_AIS_chronic pla_ctrl_(\d+).lif.tif.csv')

ankg_pattern = re.compile('rt_AnkG_ROI_exp2_AnkG_MAX_AIS plasticity_synpo_DIV12_ctrl_(\d+).lif.tif.csv')
synpo_pattern = re.compile(r'rt_Synpo_ROI_exp2_Synpo_MAX_AIS plasticity_synpo_DIV12_ctrl_(\d+).lif.tif.csv')
# rt_Synpo_ROI_exp1_Synpo_MAX_AIS plasticity_synpo_DIV12_ctrl_(\d+).lif.tif.csv



## type in name of protein need to be detected ##
keyword_1="AnkG"
keyword_2="Synpo"
file_kw="rt"
#keyword_3="KCl"
exp_name="ctrl_DIV12"

save_exe=r"E:\Yuhao_do not delet\images_tem\AIS_plasticity_Synaptopodin\20231018_AIS plasticity_synpo_DIV12_ctrl_02"


## number of pixel for normalisation, threshold for detection and pixel size ##
window_size_synpo = 10
window_size_AnkG = 86
co_range=0.5
co_range_dy=600
cut_off=0.4
pixel_size=0.08

## define functions ##

## function_01: load coordinates and intensities of AnkG or Synpo into a dictionary 
## dict format: data_dict{AnkG_fileindex:{ coords: [[,], [,]] 
##                                        intensity: [, ,]}}

def read_csv_file(filename):
    data = pd.read_csv(filename)
    coords = data.iloc[:, :2].values.tolist()  # Convert DataFrame to a list of coordinates
    intensity = data.iloc[:, 2].tolist()  # Convert Series to a list of intensities

    result = {'coords': [], 'intensity': []}
    for coord, intens in zip(coords, intensity):
        result['coords'].append(coord)
        result['intensity'].append(intens)
    
    return result

## function_02: file name handling 
## check if the name of each files in the directory match the correct pattern. If not, correct it.

# def correct_file_names_in_folder(folder_path, keyword, keyword_2, pattern):
#     print("Protein name: ", keyword)
#     for filename in os.listdir(folder_path):
#         if filename.endswith(".csv") and keyword_2 in filename:
#             print("checking file: ", filename)
#             if keyword in filename and not pattern.match(filename):
#                 print("incorrect filename")
#                 print("correcting file name...")
#                 number = re.search(r'(\d+)', filename)
#                 if number:
#                     new_filename = f'rt_{keyword}_ROI_MAX_AIS_chronic pla_NaCl_{number.group(1)}.lif.csv'
#                     old_file_path = os.path.join(folder_path, filename)
#                     new_file_path = os.path.join(folder_path, new_filename)
#                     os.rename(old_file_path, new_file_path)
#                     print(f'Renamed: {filename} to {new_filename}')
#             if keyword in filename and pattern.match(filename):
#                 #print("filename: ", filename)
#                 print("filename correct")
#     print("check done!")
        
## step_01: check file names in the folder use function_02 ##
#correct_file_names_in_folder(file_path, keyword_1, file_kw, ankg_pattern)

#correct_file_names_in_folder(file_path, keyword_2, file_kw, synpo_pattern)


## step_02: list all the files in the directory and generate a list of files containing AnkG or Synpo

files=os.listdir(file_path)
fn_in=[]
for f in files:
    if f.endswith(".csv") and "rt" in f:
        #print(f)
        fn_in.append(f)
#print("fn_in: ", fn_in)

## step_03: use function_01 to extract coords and intensities and store into a dictionary with desired 
## file name

data_dict={}
for fn in fn_in:
    if "AnkG" in fn:
        print("file name: ", fn)
        file_ind=re.search(ankg_pattern, fn).group(1)
        #ankg_match=ankg_pattern.match(fn)
        #print("matching AnkG file: ", str(ankg_match))
        print("index: ", file_ind)
        data_dict[f'AnkG_{file_ind}'] = read_csv_file(fn)
        
        #synpo_file_name = f'rt_Synpo_ROI_exp3_Synpo_MAX_AIS_chronic pla_ctrl_{file_ind}.lif.tif.csv'
        synpo_file_name = f'rt_Synpo_ROI_exp2_Synpo_MAX_AIS plasticity_synpo_DIV12_ctrl_{file_ind}.lif.tif.csv'
        
        ## rt_Synpo_ROI_exp1_Synpo_MAX_AIS plasticity_synpo_DIV12_ctrl_(\d+).lif.tif.csv
        
        if synpo_file_name in fn_in:
            print("matching synpo file found")
            synpo_f=synpo_file_name
            print("synpo_file: ", synpo_f)
            print("extracting coords and intensity of the corresponding Synpo file...")
            data_dict[f'Synpo_{file_ind}'] = read_csv_file(synpo_f)
            print("done")
            print("moving to the next file")
            print()
            
        else:
            print("Error when searching for matching Synpo file")
            continue
        

## check how many files are loaded in the dic ##
print("Number of files loaded into the dict: ", len(data_dict))
print(data_dict.keys())


## initialise dictionary to store final results ##
final_results={}
error_file=[]

for k in data_dict.keys():       
    if keyword_2 in k:
        print("Detecting synaptopodin peaks...")
        print("name of profile: ",k)
        pro_ind=file_ind=re.search(r'\d+', k).group() ## extract file index
        cor_ankg="AnkG_"+str(pro_ind)  ## locate corresponding AnkG profile
        #print("corresponding AnkG file: ", cor_ankg)
        print("profile index: ", pro_ind)
        print("extracted intensity (first 10 values): ",data_dict[k]["intensity"][0:10])
        
        ## initialise dict with name: AIS + file index (e.g. AIS_01) ##
        final_results["AIS_"+str(pro_ind)]={} ## initialise dict with name: AIS + file index (e.g. AIS_01) ##
        
        ## smooth intensity using window size for synaptodin ##
        smoothed_signal = np.convolve(data_dict[k]["intensity"], np.ones(window_size_synpo) / window_size_synpo, mode='same')
        print("max_int of smoothed signal: ", smoothed_signal.max())
        normalized_signal = smoothed_signal / smoothed_signal.max()
        print(f"normalised intensity of {k}: ",normalized_signal[0:10])
        print("length of norm signal: ", len(normalized_signal))
        GV_len = np.array(normalized_signal)
        print("GV_len: ",GV_len[0:10])
        #print(len(GV_len))
        
        ## Find the peak indices and distances between peaks##
        he = 0.4 * GV_len.max()
        peak_indices = find_peaks(GV_len, height=he, threshold=0.001, distance=15)[0]
        peak_value = [GV_len[j] for j in peak_indices]
        print("peak indices: ", peak_indices, f" length of peak_indices: ", len(peak_indices))
        print("peak value: ", peak_value, f" length of peak values: ", len(peak_value))
        
        ## plot detected synpo clusters ##
        plt.plot(GV_len)
        plt.xlabel('distance (pixel)')
        plt.ylabel('Gray Value (normalised)')
        plt.title(k)
        for i in range(len(peak_indices)):
            plt.annotate("Synpo_" + str(peak_indices[i]), (peak_indices[i], peak_value[i]), 
            textcoords="offset points", xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none', shrink=0.01))
        #plt.savefig()
        #plt.show()
        
        ## extract synpo peak coords from synpo profile ##
        print("extracting synaptopodin peak coordinates...")
        synpo_coords = data_dict[k]["coords"]
        
        
        
        peak_coords = [synpo_coords[pc] for pc in peak_indices]
        print("Synpo peak coords: ", peak_coords)
        print("Number of synpo peak coords: ", len(peak_coords))
        if len(peak_coords) != len(peak_indices):
            raise ValueError("Error!! Number of detected peak does not match number of coordinates!")
            
        else:
            print("Number of detected cluster match the number of coordinates!")
            
        
        ## check corresponding AnkG file ##
        print()
        print("##############################################")
        print("checking corresponding AnkG file...")
        print("corresponding AnkG file: ", cor_ankg)
        ankg_coords = data_dict[cor_ankg]["coords"]
        print("coords from corresponding AnkG profile (first 10 coords): ", ankg_coords[0:10])
        print("length of corresponding AnkG profile coords: ", len(ankg_coords))
        
        
        
        # if not len(ankg_coords) > len(synpo_coords):   ## check if the length of AnkG line is shorter than Synpo line. if yes, skip this file
        #     print("length of AnkG coords in file ", k, " is shorter than Synpo coords")
        #     error_file.append(k)
        #     continue
        
        print()
        print("Searching Synpo cluster coords in AnkG file...")
        
        ## searching synpo peak coords in AnkG coordinates and store it in list det_synpo_coords ##
        det_synpo_coords=[]
        det_synpoc_ind=[]
        undet_synpo_coords=[]
        
        co_c= 0
        
        for co in peak_coords:
            #print("x", co[0], "y", co[1])  
            syn_x=co[0]
            syn_y=co[1]    ## to get x and y of synpo peak coords
            co_c+=1 ## count how many synpo cluster coords are processed
            print("calculating Synpo cluster ", co, " ...")
            print("cluster number: ", co_c)
            c_cal=0 
            print("finding matching AnkG coords based on dx...")
            label=0 ## initialise a label to indicate current state of the synpo coord
                    ## 0: unprocessed/unmatching ankg coord found;
                    ## 1: matching dx; 2: matching ankg coord found; 3: dy is not in range; 4: dx out of range
                    
            for a_co in ankg_coords:  ## loop through each ankg coords and assign x, y 
                ank_x = a_co[0]  ## get x and y of AnkG coords
                ank_y = a_co[1]  
                
                d_x= abs(syn_x - ank_x)  ## calculate dx of synpo coord to each ankg coords
                c_cal+=1 ## count how many ankg coords have been processed against synpo cluster coord
                
            ## search synpo peak coords in AnkG coords within range of 0.5 pixels in this case ##   
                #if ank_x-co_range <= syn_x <= ank_x+co_range or ank_y-co_range <= syn_y <= ank_y+co_range:
                    
     ## if d_x in range, calculate dy. If dy in range, append current ankg coord, then move to next Synpo coord and do the same. 
     ## if not, keep calculating dx and dy of synpo coord to ankg coords     
                if d_x <= co_range:                     
                    label=1 
                    d_y = abs(syn_y - ank_y)           
                    if d_y <= co_range_dy:             
                        det_synpo_coords.append(a_co)
                        det_synpoc_ind.append(int(c_cal)-1)
                        print("!!! matching coord in AnkG: ", a_co)
                        print("The index of matching AnkG coord: ", c_cal-1)
                        print("dx",d_x)
                        print("dy",d_y)
                        label=2
                        break
                    
    ## if dx is not in range. Add current synpo to the list if it has not been added. Otherwise continue to the next
                    else:
                        print("did not find matching coord in ankg profile for current synpo cluster: ", co)
                        print("dy "+ str(d_y) +" is not in range.")
                        print("coords calculated (synpo, ankg): ", co,  a_co)
                        print("adding synpo cluster to unidentified list and searching base on dy...")
                        label = 3
                        if not co in undet_synpo_coords:
                            undet_synpo_coords.append(co)
                            break
                        continue     
                    
            print("number of coords calculated: ", str(c_cal))
            
            #print()
            if label == 0 or label == 3:
                print("Did not find matching ankg coord for synpo cluster based on dx: ", co)
                print("searching again based on dy...")
                c_cal_dy=0
                rlabel = 0 ## set another label. same as label
                for ra_co in ankg_coords:  ## loop through each ankg coords and assign x, y 
                    rank_x = ra_co[0]  ## get x and y of AnkG coords
                    rank_y = ra_co[1]  
                    c_cal_dy+=1
                    rd_y = abs(syn_y - rank_y)  ## calculate dx of synpo coord to each ankg coords
                    if rd_y <= co_range:
                        rlabel=1
                        rd_x =  abs(syn_x - rank_x)
                        if rd_x <= co_range_dy:
                            if not ra_co in det_synpo_coords:
                                det_synpo_coords.append(ra_co)
                                det_synpoc_ind.append(int(c_cal_dy)-1)
                            print("!!! matching coord in AnkG based on dy: ", ra_co)
                            print("The index of matching AnkG coord: ", c_cal_dy-1)
                            print("dx",rd_x)
                            print("dy",rd_y)
                            rlabel=2
                            break
                        else:
                            print("did not find matching coord in ankg profile for current synpo cluster with dy: ", co)
                            print("coords calculated (synpo, ankg): ", co,  ra_co)
                            print("adding synpo cluster to unidentified list and searching base on dy...")
                            rlabel = 3
                            if not co in undet_synpo_coords:
                                undet_synpo_coords.append(co)
                                break
                            continue     
                            
                print()
            
            else:
                print("#######################################################")
                print()
            
            
            
        ## check if the length of detected synpo coords match the number of synpo clusters detected in synpo profile            
        if len(det_synpo_coords) == len(peak_indices):
            print("Synpo coords found in AnkG coords list")
            print(str(len(peak_indices))+" clusters detected in Synpo profile and " + str(len(det_synpo_coords))+
                  " coords detected in corresponding AnkG profile.")
            print("Coordinates of Synpo clusters detected in AnkG profile: ",det_synpo_coords)
            print("Calculating distance of cluster to axon start...")
        else:
            print("!!!#######################################")
            print("Error! Unmatching number of clusters detected in AnkG file")            
            print("Coordinates of Synpo clusters detected in AnkG profile: ", det_synpo_coords)
            print("Number of Synpo cluster coords detected in AnkG: ", len(det_synpo_coords))
            print("Number of Synpo clusters detected from Synpo profile: ", len(peak_coords))
            print("Coordinates of Synpo clusters detected in Synpo profile: ", peak_coords)
            print("length of Synpo coords: ", len(synpo_coords))
            print("length of AnkG coords: ", len(ankg_coords))
            print()
            print("!!!#######################################")
            raise ValueError("Error! Unmatching number of clusters detected in AnkG file")
        
        ## calculate the distance of synpo clusters in the AIS by using AnkG as referrence ##
        synpo_dis_to_axon=[]
        #synpo_cluster_dis=[]
        unit_dis = 0
        for a_co in range(len(ankg_coords)-1):
            xa, ya = ankg_coords[a_co]      #define x1,y1
            #print(xa,ya)                             #check
            xb, yb = ankg_coords[a_co+1] 
            dx, dy = xb - xa, yb - ya                 #calculate dx and dy
            unit_length = math.sqrt(dx**2 + dy**2)    
            unit_dis += unit_length                   # distance of each pixel to axon start
            if ankg_coords[a_co] in det_synpo_coords:
                print("synpo_coord", ankg_coords[a_co])
                print("distance to axon start (micron): ", unit_dis*pixel_size)
                synpo_dis_to_axon.append(unit_dis*pixel_size)
             
        print("Storing information of detected clusters...")
        synpo_dict= {"cluster number": len(peak_indices),
                     "cluster coords": peak_coords,
                     "cluster coords in AnkG profile": det_synpo_coords,
                     "distance to axon start (um)": synpo_dis_to_axon
                     }
        
        final_results["AIS_"+str(pro_ind)]["Synpo"]=synpo_dict
        print("##############################################")
        
        #############################################################################################################################
        ## measure AIS length using corresponding AnkG intensity ##
        print("Calculating corresponding AIS length and location...")
        ankg_int=data_dict[cor_ankg]["intensity"]
        smoothed_signal_ankg = np.convolve(ankg_int, np.ones(window_size_AnkG) / window_size_AnkG, mode='same')
        ## Normalize the smoothed signal values
        normalized_signal_ankg = smoothed_signal_ankg / smoothed_signal_ankg.max()
        ## Convert the normalized values to a numpy array
        GV_len_ankg = np.array(normalized_signal_ankg)
        
        ## plot AnkG profile ##
        plt.plot(GV_len_ankg )
        #plt.xlabel('distance (pixel)')
        #plt.ylabel('Gray Value (normalised)')
        #plt.title('AnkG length')
        
        ## extract the index of values in normalised AnkG profile
        indices_ankg=np.arange(len(GV_len_ankg))  
        ## extract the index of max value in normalised AnkG profile            
        max_index_ankg = np.argmax(GV_len_ankg)
        
        plt.annotate("max", (max_index_ankg, GV_len_ankg[max_index_ankg]), textcoords="offset points",
                     xytext=(0,10), arrowprops=dict(facecolor='pink', edgecolor='none', shrink=0.01))
        
        ## lock the max value and also set the threshold for detecting the border of AIS
        GV_max_ankg=GV_len_ankg.max()
        thr_v_ankg = GV_max_ankg * cut_off    
        #print("Index of max value:", max_index)
        #print("threshold for length determination:", thr_v)
        
        ## split the AnkG profile into two parts based on max value
        left_of_max = GV_len_ankg[:max_index_ankg][::-1]  # Reversed to start from max
        right_of_max = GV_len_ankg[max_index_ankg:]
        
        ## initialise the left and right border index and value
        left_index_below_threshold = 0
        right_index_below_threshold = 0
        Lva= None
        Rva= None
        
        ## loop over the velues on left and right side of max value
        ## and detect the border of AIS       
        
        ## Find the index of the first value below 40% on the left (AIS start)
        for i, value in enumerate(left_of_max):
            if value < thr_v_ankg:
                Lva=value
        # Store the index considering that left_of_max is reversed
                left_index_below_threshold = max_index_ankg - i
                break
        
        ## Find the index of the first value below 40% on the right (AIS end)
        for i, value in enumerate(right_of_max):
            if value < thr_v_ankg:
                Rva=value
                right_index_below_threshold = max_index_ankg + i
                break

        ## Return the results
        #print(f"The index of the first value below 40% to the left of max: {left_index_below_threshold}")
        #print(f"The index of the first value below 40% to the right of max: {right_index_below_threshold}")         
        print(f"AIS start: {left_index_below_threshold*pixel_size}")
        print(f"AIS end: {right_index_below_threshold*pixel_size}")  
        
        ## using the right index minus left index to calculate the length of AIS
        band_width= right_index_below_threshold - left_index_below_threshold
        print("length of AIS (AnkG, micron): ", band_width*pixel_size)
        
        
        
        plt.annotate("start", (left_index_below_threshold,GV_len_ankg[left_index_below_threshold]), textcoords="offset points", 
                     xytext=(0,10), arrowprops=dict(facecolor='pink', edgecolor='none',
                                                     shrink=0.01))
        plt.annotate("end", (right_index_below_threshold, GV_len_ankg[right_index_below_threshold]), textcoords="offset points", 
                     xytext=(0,10), arrowprops=dict(facecolor='pink', edgecolor='none', 
                                                     shrink=0.01))
        plt.annotate("Axon start", (0, GV_len_ankg[0]), textcoords="offset points", 
                     xytext=(0,10), arrowprops=dict(facecolor='black', edgecolor='none', 
                                                     shrink=0.01))
        
        for n in det_synpoc_ind:
            ind= n-1
            plt.annotate("", (ind,GV_len_ankg[ind-1]), textcoords="offset points", 
                         xytext=(0,10), arrowprops=dict(facecolor='blue', edgecolor='none',
                                                         shrink=0.01))
            
        
        plt.show()
        
        ## store information in a dictionary ##
        print("Storing information of corresponding AIS (AnkG)...")
        ankg_dict= {"AIS length (um)": band_width*pixel_size,
                    "AIS start (um)": left_index_below_threshold*pixel_size,
                    "AIS peak (um)": max_index_ankg*pixel_size,
                    "AIS end (um)": right_index_below_threshold*pixel_size
                     }
        
        final_results["AIS_"+str(pro_ind)]["AnkG"]=ankg_dict
        
        print("===================================== current profile finished! =====================================================")
        print("loading next profile...")
    print()
print("All files processed!")
print("analysis finished!")

## creating separated dict for Synpo and AnkG ##
print()
print("Creating individual dict for Synpo and AnkG...")
ankg_final={}
synpo_final={}
AIS_ind=0

s_dta=[]
s_cluster=[]
AIS_st=[]
AIS_pe=[]
AIS_en=[]
AIS_le=[]

for ks in final_results.keys():
    #print(ks)
    AIS_ind+=1
    #print("index", AIS_ind)
    for sk in final_results[ks].keys():
        #print(ks, sk)
        
        if sk == "Synpo":
            #print(ks, sk, final_results[ks][sk]['distance to axon start (um)'])
            dis_to_axon=final_results[ks][sk]['distance to axon start (um)']
            num_cluster=final_results[ks][sk]['cluster number']
            for values in dis_to_axon:
                s_dta.append(values)
                synpo_final["distance to axon (um)"]=s_dta
            
            s_cluster.append(num_cluster)
            synpo_final["number of cluster"]=s_cluster    
            
            
        if sk == "AnkG":
            #print(ks, sk, final_results[ks][sk].keys())
            A_S=final_results[ks][sk]['AIS start (um)']
            #print(A_S)
            A_P=final_results[ks][sk]['AIS peak (um)']
            A_E=final_results[ks][sk]['AIS end (um)']
            A_LE=final_results[ks][sk]['AIS length (um)']
            
            AIS_st.append(A_S)
            ankg_final["AIS_start"]=AIS_st
            
            AIS_pe.append(A_P)
            ankg_final["AIS_peak"]=AIS_pe
            
            AIS_en.append(A_E)
            ankg_final["AIS_end"]=AIS_en
            
            AIS_le.append(A_LE)
            ankg_final["AIS_length"]=AIS_le
ave_syn_clu=[]            
for i in range(len(synpo_final["number of cluster"])-1):
    #print(i, synpo_final["number of cluster"][i], ankg_final["AIS_length"][i])
    ave_synpo=(synpo_final["number of cluster"][i]/ankg_final["AIS_length"][i])*10
    ave_syn_clu.append(ave_synpo)
    ave_synpo=0
    ankg_final["num of Synpo clusters / 10 um"]=ave_syn_clu

print("Done! Saving results as excel...")           
print()

## saving results ##
#exp_name="D14_ctrl"

writer_synpo=pd.ExcelWriter(save_exe+"\\"+exp_name+"_synpo_AIS.xlsx", engine='openpyxl')
writer_ankg=pd.ExcelWriter(save_exe+"\\"+exp_name+"_AnkG_AIS.xlsx", engine='openpyxl')

df_synpo=pd.DataFrame.from_dict(synpo_final, orient='index').T
df_synpo.to_excel(writer_synpo, sheet_name=" ")

df_ankg=pd.DataFrame.from_dict(ankg_final, orient='index').T
df_ankg.to_excel(writer_ankg, sheet_name=" ")
#     #writer.book.create_sheet(sh_n[-4:-1])
# #writer.save()   

writer_synpo.close()
writer_ankg.close()

print("Results saved!!!")






# for sh_n, sh_data in final_results.items():    #d_n is name of directory, di is the sub_dictionary
#     #print(sh_n, sh_data)
#     #sh_n=sh_n.split("\\")[-2]
#     #print(sh_n)
    
#     #sheet=writer.book.add_worksheet(sh_n[-4:0])
#     #print(sh_n[0:-4])
#     df=pd.DataFrame.from_dict(sh_data, orient='index')
#     df.to_excel(writer, sheet_name=sh_n)
#     #writer.book.create_sheet(sh_n[-4:-1])
# #writer.save()   
# writer.close()
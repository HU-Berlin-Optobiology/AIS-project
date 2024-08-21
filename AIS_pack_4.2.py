import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import math
from itertools import zip_longest
import re

### user input ###
# set directories
path = r"Z:\Yuhao_do not delet\images_Marseille\Analysis\AIS plasticity"
save_exe = r"Z:\Yuhao_do not delet\images_Marseille\Analysis\AIS plasticity\final_output"

k_w = "ankg" # set the key word in the directory name that needs to be processed
exp_name = "nonacd" # set the name of file name for final results storage 

# set experimental condition in your experiment
exp_condition_1="ctrl" 
exp_condition_2="nmda"
exp_condition_3="none"

# set the distance for normalising the intensity profile
window_size = 54 
cut_off = 0.4 # set the limit for length definition 

# optional: generate averaged profile before analysing the length and shift of AIS protein
average_p = "no" # type yes if want to visualise averaged profile before set threshold for AIS analysis 
cut_off_ave= 0.2 # threshold after visualising average profile
align_max = "yes" # align at max for average profile
align_start = "no" # align at start of the AIS after analysis

######################## end of user input #############################
######################## start of actual processing #############################

### first sort csv files in the directory that have the correct keyword
#k_w = k_w.lower() 
os.chdir(path) #go to the directory
f_roots=[] #creat a list to store all roots with csv files
print("Checking files...")
nu_files=0 ## check  number of files
#loop over all files and roots in the directory
for roots, directory, files in os.walk(path):
    if k_w in roots:
        #print("roots: ", roots)
        f_roots.append(roots) # if keyword is in root, add this root into root list
        os.chdir(roots) # go to the directory
        file_list_r = os.listdir(roots) # list all files in the root
        for f in file_list_r:   # loop over all files in the root
            if f.endswith(".csv"): # if file ends with csv count the file
                #print(f)
                nu_files+=1
            
print("Number of csv files detected: ", str(nu_files))
print("")

### now we should have all the directories that contain AnkG profiles in csv format
### we then go into each directory to process these profiles
## before we start, we define useful functions to process profiles
## 1. define function to generate normalised intensity line profile
def normalise_profile (file_name, window_size):
    if file_name.endswith(".csv"):
        #print("loading csv file:", file_name)
        # Read the CSV file for length measurement
        
        df = pd.read_csv(file_name, usecols=['Value'])
        # Smooth the length profile by a specified window size using convolution
        #print("intensity normalisation...")
        smoothed_signal = np.convolve(df['Value'], np.ones(window_size) / window_size, mode='same')

        # Normalize the smoothed signal values
        normalized_signal = smoothed_signal / smoothed_signal.max()

        # Convert the normalized values to a numpy array
        GV_len = np.array(normalized_signal)
        
        return GV_len
    else:
        print("File does not end with .csv")
        return None
## define a function to plot all profiles, aligned profiles and averaged profiles
def plot_profiles (profiles_plt, experiment_name_condition):
    plt.figure(figsize=(10, 6))
    for p in profiles_plt:
        plt.plot(p)
    plt.xlabel('Length (pixel)')
    plt.ylabel('Gray Value (norm)')
    plt.title(experiment_name_condition+"_all profiles")   
    plt.show()

def plt_alignedprofiles (aligned_pro, experiment_name_condition):
    if not aligned_pro is None:
        plt.figure(figsize=(10, 6))
        for p in aligned_pro:
            plt.plot(p)
        plt.title(experiment_name_condition+"_aligned profiles")
        plt.xlabel("Position")
        plt.ylabel("Intensity")
        #plt.legend(loc=(0.7, 0.9))
        plt.show()
        
    if aligned_pro is None:
        print("no profiles to be plotted")

def plt_average_pro (average_profile_1, average_profile_2, average_profile_3, std_dev_profile_1, std_dev_profile_2,
                     std_dev_profile_3, label_1, label_2, label_3):
    
    plt.figure(figsize=(10, 6))
    
    average_profiles = [average_profile_1, average_profile_2, average_profile_3]
    std_li = [std_dev_profile_1, std_dev_profile_2, std_dev_profile_3]
    colors = ["black", "red","blue"]
    labels = [label_1, label_2, label_3]
    
    for index in range(len(average_profiles)):
        if average_profiles[index] is not None:
            plt.plot(average_profiles[index], color= colors[index], label= labels[index])
            plt.fill_between(range(len(std_li[index])),
                             average_profiles[index] - std_li[index],
                             average_profiles[index] + std_li[index],
                             color=colors[index], alpha=0.2, label="±SD")
        
    plt.title("Average Profile")
    plt.xlabel("length (pixels)")
    plt.ylabel("norm Intensity")
    plt.legend(loc=(0.8, 0.6))
    plt.show()   

def plt_aligned_average_pro (aligned_aveprofiles, std_dev_profiles, labels):
    
    plt.figure(figsize=(10, 6))    
    average_profiles = [item for item in aligned_aveprofiles if item is not None]
    stdv = [item for item in std_dev_profiles if item is not None]
    colors = ["black", "red","blue"]
    labels = labels
    
    for index in range(len(average_profiles)):
        if average_profiles[index] is not None:
            
            plt.plot(average_profiles[index], color= colors[index], label= labels[index])
            plt.fill_between(range(len(average_profiles[index])),
                             average_profiles[index] - stdv[index],
                             average_profiles[index] + stdv[index],
                             color=colors[index], alpha=0.2, label="±SD")
        
    plt.title("aligned Average Profile")
    plt.xlabel("length (pixels)")
    plt.ylabel("norm Intensity")
    plt.legend(loc=(0.85, 0.75))
    plt.show()  
    

## define function to generate aligned profiles from a list of profiles   
def align_profiles_max (profiles):
    maxindli = []
    aligned_profiles = []
    if len(profiles) > 0:
        for p in profiles:
            max_index = np.argmax(p)            
            maxindli.append(max_index)
    
        max_center = max(maxindli)
        print("")
        print("aligning profiles...")
        print("centre of alignment: ", str(max_center))
   
        for p in range(len(profiles)):
            profile = profiles[p]
            max_index_pro= maxindli[p]
            shift = max_center - max_index_pro
        
            if shift > 0:
                aligned_profile = np.pad(profile, (shift, 0), 'constant')
            elif shift < 0:
                aligned_profile = np.pad(profile, (0, -shift), 'constant')
            else:
                aligned_profile = profile
        
            aligned_profiles.append(aligned_profile)
    
        max_length = max(len(ap) for ap in aligned_profiles)
        aligned_profiles = [np.pad(ap, (0, max_length - len(ap)), 'constant') for ap in aligned_profiles]
    
        return aligned_profiles
    
    if len(profiles) == 0:
        print("")
        print("no profiles found")
    
## define function to generate average profile based on the aligned profiles
def average_profile(profiles):
    if profiles is not None and len(profiles) > 0:
        average_profile = np.mean(profiles, axis=0)
        std_dev_profile = np.std(profiles, axis=0)
        return average_profile, std_dev_profile
    else:
        return None, None
        
## define function to align averaged_profiles and stdv
def align_aveprofiles_max (profiles, stdvs):
    maxindli = []
    aligned_profiles = []
    aligned_stdv = []
    
    if len(profiles) > 0:
        for p in profiles:
            max_index = np.argmax(p)            
            maxindli.append(max_index)
    
        max_center = max(maxindli)
        print("")
        print("aligning averaged profiles...")
        print("centre of alignment: ", str(max_center))
   
        for p in range(len(profiles)):
            profile = profiles[p]
            max_index_pro= maxindli[p]
            shift = max_center - max_index_pro
        
            if shift > 0:
                aligned_profile = np.pad(profile, (shift, 0), 'constant')
            elif shift < 0:
                aligned_profile = np.pad(profile, (0, -shift), 'constant')
            else:
                aligned_profile = profile
        
            aligned_profiles.append(aligned_profile)
    
        max_length = max(len(ap) for ap in aligned_profiles)
        aligned_profiles = [np.pad(ap, (0, max_length - len(ap)), 'constant') for ap in aligned_profiles]
        
        
        ## align stdvs based on max_index of profile
        for sd in range(len(stdvs)):
            sd_profile = stdvs[sd]
            max_index_std = maxindli[sd]
            shift = max_center - max_index_std
        
            if shift > 0:
                aligned_sdprofile = np.pad(sd_profile, (shift, 0), 'constant')
            elif shift < 0:
                aligned_sdprofile = np.pad(sd_profile, (0, -shift), 'constant')
            else:
                aligned_sdprofile = sd_profile
        
            aligned_stdv.append(aligned_sdprofile)
    
        max_length_sd = max(len(apsd) for apsd in aligned_stdv)
        aligned_stdv = [np.pad(apsd, (0, max_length - len(apsd)), 'constant') for apsd in aligned_stdv]
    
    
        return aligned_profiles, aligned_stdv
    
    if len(profiles) == 0:
        print("")
        print("no profiles found")



## define a function to analyse AIS length and position
def AIS_analysis (p, cutoffs, name):
    
    indices=np.arange(len(p))  
# extract the index of max value in normalised AnkG profile            
    max_index = np.argmax(p)
# lock the max value and also set the threshold for detecting the border of AIS
    GV_max=p.max()
    thr_v = GV_max * cutoffs
            
    print("Index of max value:", max_index)
    print("threshold for length determination:", thr_v)        
    print("Ploting intensity profile")
    plt.plot(p)
    plt.xlabel('Length (pixel)')
    plt.ylabel('Gray Value (norm)')
    plt.title(name)   
    plt.annotate("max", (max_index,p[max_index]), textcoords="offset points", xytext=(0,10), 
                 arrowprops=dict(facecolor='red', edgecolor='none', shrink=0.01))   
# split the AnkG profile into two parts based on max value
          
    left_of_max = p[:max_index][::-1]  # Reversed to start from max
    right_of_max = p[max_index:]
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


    plt.annotate("left edge", (left_index_below_threshold,p[left_index_below_threshold]), textcoords="offset points", 
                 xytext=(0,10), arrowprops=dict(facecolor='red', edgecolor='none',
                                                 shrink=0.01))
            
    plt.annotate("right edge", (right_index_below_threshold,p[right_index_below_threshold]), textcoords="offset points", 
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
    print("")
    return peak_dict
            
        
        
    

############################# end of function definition #######################################

AIS_profile={} 
AIS_ana={}
print("normalising intensity profile to maximum intensity")
# loop over all roots in f_roots, if it contains keyword, go into it and process profiles
for rs in f_roots: 
    if k_w in rs:
        os.chdir(rs)
        file_list_f = os.listdir(rs) # list all files in the root
        AIS_profile[rs]={} #specify directories contain correct files
        AIS_ana[rs]={}
        
        for f in file_list_f: # loop over fll files in the root
            filename = str(f) # extract file name
            norm_pro = normalise_profile (f, window_size) # normalise profile as norm_pro
            AIS_profile[rs][filename] = norm_pro # store all normalised profile into a dictionary

print("profiles normalised! No errors occured :)")
print("") 


## bring user input to lower case
average_p =  average_p.lower()
align_max =  align_max.lower()
align_start = align_start.lower()    

## creat list to store profiles for alignment
profile_exp1=[]
profile_exp2=[]
profile_exp3=[]      


## if user decide to visualise averaged profiles 
if average_p == "yes":
    print("generating averaged profiles for", str(k_w))
    print("experiment name:", str(exp_name))
    
    for k in AIS_profile.keys():   # loop over all keys in AIS_profile
        if exp_name in k: # if experiment name is in the keys (e.g. nonacd is in key, do something)
            #print(k)
            for profiles in AIS_profile[k].keys(): # loop over profiles in the dataset 
                if exp_condition_1 in profiles.lower(): # go to profiles from experiment condition 1 
                    #print(profiles) # name of the profile
                    #print(AIS_profile[k][profiles][0:5]) # access the actual profile
                    profile_exp1.append(AIS_profile[k][profiles]) # append all profiles into the list
                    
                    
                if exp_condition_2 in profiles.lower():
                    #print(profiles)
                    #print(AIS_profile[k][profiles][0:5]) # access the actual profile
                    profile_exp2.append(AIS_profile[k][profiles])
                
                if exp_condition_3 in profiles.lower():
                    profile_exp3.append(AIS_profile[k][profiles])
   
    ## now all profiles are added into list according to experiment name (nonacd) and experimental condition
    ## we proceed by ploting profiles, align profiles and then visualise the aligned profiles
    # first we visualise all profiles according to conditions
    plot_profiles(profile_exp1, exp_name+"_"+exp_condition_1)
    plot_profiles(profile_exp2, exp_name+"_"+exp_condition_2)
    plot_profiles(profile_exp3, exp_name+"_"+exp_condition_3)
    
    # then we align profiles on maximum intensity and visualise them
    aligned_pro_exp1 = align_profiles_max(profile_exp1)
    aligned_pro_exp2 = align_profiles_max(profile_exp2)
    aligned_pro_exp3 = align_profiles_max(profile_exp3)
    
    plt_alignedprofiles(aligned_pro_exp1,  exp_name+"_"+exp_condition_1)
    plt_alignedprofiles(aligned_pro_exp2,  exp_name+"_"+exp_condition_2)
    plt_alignedprofiles(aligned_pro_exp3,  exp_name+"_"+exp_condition_3)
    
    # then generate averaged profile, calculate standard deviation and plot
    average_profile_exp1, sd_avepro_exp1 = average_profile(aligned_pro_exp1)
    average_profile_exp2, sd_avepro_exp2 = average_profile(aligned_pro_exp2)
    average_profile_exp3, sd_avepro_exp3 = average_profile(aligned_pro_exp3)
    
    #plt_average_pro(average_profile_exp1, average_profile_exp2, 
    #                average_profile_exp3, sd_avepro_exp1, sd_avepro_exp2, sd_avepro_exp3,
    #                exp_condition_1, exp_condition_2, exp_condition_3)
    
    ## align average profiles and corresponding sd to max intensity index and plot
    aveprofilelist = [average_profile_exp1, average_profile_exp2, average_profile_exp3]
    sd_aveprolist = [sd_avepro_exp1,  sd_avepro_exp2,  sd_avepro_exp3]
    aveprofilelist = [item for item in aveprofilelist if item is not None]
    sd_aveprolist = [item for item in sd_aveprolist if item is not None]

    aligned_aveprofiles, aligned_sdvs = align_aveprofiles_max(aveprofilelist, sd_aveprolist)

    labels = [exp_condition_1, exp_condition_2, exp_condition_3]
    plt_aligned_average_pro(aligned_aveprofiles, aligned_sdvs, labels)
    
    ################################################ end of making averaged profile ###################################
    ### now proceed with analysing AIS length and position
    
    cutoff_ave = 0.35
    #cutoff_ave = input("set threshold for AIS analysis!")
    #cutoff_ave = float(cutoff_ave)
    
    for k in AIS_profile.keys():   # loop over all keys in AIS_profile
        if exp_name in k: # if experiment name is in the keys (e.g. nonacd is in key, do something)
            #print(k)
            for profiles in AIS_profile[k].keys(): # loop over profiles in the dataset 
                if exp_condition_1 in profiles.lower(): # go to profiles from experiment condition 1 
                    # actual profile is AIS_profile[k][profiles]
                    print("processing", profiles)
                    AIS_results_1 = AIS_analysis(AIS_profile[k][profiles], cutoff_ave, profiles)
                    AIS_ana[k][profiles]=AIS_results_1
                    
                if exp_condition_2 in profiles.lower(): # go to profiles from experiment condition 1 
                    # actual profile is AIS_profile[k][profiles]
                    print("processing", profiles)
                    AIS_results_2 = AIS_analysis(AIS_profile[k][profiles], cutoff_ave, profiles)
                    AIS_ana[k][profiles]=AIS_results_2
                    
                if exp_condition_3 in profiles.lower(): # go to profiles from experiment condition 1 
                    # actual profile is AIS_profile[k][profiles]
                    print("processing", profiles)
                    AIS_results_3 = AIS_analysis(AIS_profile[k][profiles], cutoff_ave, profiles)
                    AIS_ana[k][profiles]=AIS_results_3
    
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
        df.to_excel(writer, sheet_name='AIS', index=False)

    print(f"Data saved to {excel_file_path}")
    
   
else:
    print("No averaged profile required...")
    print("Analysing AIS length and position...")
    print("Threshold for AIS length and position defination:", str(cut_off))
    
    cutoff_ave = cut_off
    for k in AIS_profile.keys():   # loop over all keys in AIS_profile
        if exp_name in k: # if experiment name is in the keys (e.g. nonacd is in key, do something)
            #print(k)
            for profiles in AIS_profile[k].keys(): # loop over profiles in the dataset 
                if exp_condition_1 in profiles.lower(): # go to profiles from experiment condition 1 
                    # actual profile is AIS_profile[k][profiles]
                    print("processing", profiles)
                    AIS_results_1 = AIS_analysis(AIS_profile[k][profiles], cutoff_ave, profiles)
                    AIS_ana[k][profiles]=AIS_results_1
                    
                if exp_condition_2 in profiles.lower(): # go to profiles from experiment condition 1 
                    # actual profile is AIS_profile[k][profiles]
                    print("processing", profiles)
                    AIS_results_2 = AIS_analysis(AIS_profile[k][profiles], cutoff_ave, profiles)
                    AIS_ana[k][profiles]=AIS_results_2
                    
                if exp_condition_3 in profiles.lower(): # go to profiles from experiment condition 1 
                    # actual profile is AIS_profile[k][profiles]
                    print("processing", profiles)
                    AIS_results_3 = AIS_analysis(AIS_profile[k][profiles], cutoff_ave, profiles)
                    AIS_ana[k][profiles]=AIS_results_3
    
    # save results
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
        df.to_excel(writer, sheet_name='AIS', index=False)

    print(f"Data saved to {excel_file_path}")
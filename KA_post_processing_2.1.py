# -*- coding: utf-8 -*-
"""
Created on Sun May  7 09:35:20 2023

@author: yuhao
"""

import os
import csv
import math
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
from itertools import zip_longest

### creat function to edit name of each track###################################
def name_mod(name):
    names=name.split(".")
    mo_name=names[0]
    return mo_name
#to calculate pause per minute
def calc_pause_permin(group):
    return (group['# of pause'] / group['TotalD (s)']) * 60
    
##############################################################################
label=0


file_dir=r"C:\Users\yuhao.HP640G8MIKH\Desktop\Final results revision\final results_LAMP1" # input your excel file directory

output_file = r"C:\Users\yuhao.HP640G8MIKH\Desktop\Final results revision\final results_LAMP1\KYMOA_PO.xlsx" # this is how you should put your final output



os.chdir(file_dir)
efile_list = os.listdir(file_dir)
#kymo_e="20211006kymo_ana.xlsx"
#f_re={}
for e_f in efile_list:
    #print(e_f)
    label+=1
    e_df=pd.read_excel(e_f, sheet_name=None, engine='openpyxl') #read excel file as df

    f_re={}

    for ks in e_df.keys(): # loop each sheet name in the file
        #print("current sheet: ", ks) #check
    #print(e_df[ks]) # check
        f_re[ks]={}
        sh_d=e_df[ks] #creat a new df equals data in the loaded excel file
    #sh_d.rename(columns={sh_d.columns[0]:"Name"}) #rename unnamed columns
        sh_d.columns.values[0]="Name"    #rename unnamed columns
    #print(sh_d.columns)
    #specify the columns that need to remove 0 values
        
        co_rep_zero=["aveAntD (s)","aveAntR (um)","aveAntV (um/s)","aveRetD (s)",
                    "aveRetR (um)","aveRetV (um/s)","avePauseD (s)","psv_dur (sec)","psv_dis (um)",
                    "psv_vel (um/s)","Reversals"]
        
        sh_d["Name"]=sh_d["Name"].apply(name_mod) #remove .tif in the name
            
        
        for c in sh_d.columns:
            if c in co_rep_zero:
                sh_d[c]=sh_d[c].replace(0,np.nan)
            else:
                sh_d[c]=sh_d[c]
        
        if len(sh_d[sh_d["Type"]=="Mobile"].groupby("Name"))>0:
            pau_permin=sh_d[sh_d["Type"]=="Mobile"].groupby("Name").apply(calc_pause_permin).reset_index().rename(columns={0:"pause_permin"})
            pau_permin=pau_permin.groupby("Name")["pause_permin"].mean()   
        else:
            pau_permin=sh_d.groupby("Name")["# of pause"].mean().reset_index().rename(columns={0:"pause_permin"})
            
        
        #sh_d=sh_d.replace(0,pd.np.nan) # replace all 0 values to n       

        co_to_sum=["# of anterograde run","# of retrograde run","# of pause","total # of run"]
        
        #co_to_mea=["aveAntD (s)","aveAntR (um)","aveAntV (um/s)","aveRetD (s)",
        #            "aveRetR (um)","aveRetV (um/s)","avePauseD (s)","psv_dur (sec)","psv_dis (um)",
        #            "psv_vel (um/s)","Total_Run_D (s)","Total_Pausing_D (s)","PercTimeAntRun (%)","PercTimeRetRun (%)",
        #            "Total travel len (um)","psv_dur (sec)","psv_dis (um)","psv_vel (um/s)","Displacement","TotalD (s)","PercTimePausing (%)","PercTimeRunning (%)",
        #            "PercTimepassivem (%)"] # creat a list to store column names in the df
        co_to_mea=[] # creat a list to store columns that need to be calculated
        
        for cs in sh_d.columns:  # loop through the columns from your excel sheet and extract those columns with numeric values
            #print(cs)
            if not cs in co_to_sum and cs!="Name" and pd.api.types.is_numeric_dtype(sh_d[cs]):
                co_to_mea.append(cs)
                
    #next step to store mobile track in mob_t and calculate mean for all columns excluding 0 and sum of the columns that were specified
        mob_t=sh_d[sh_d["Type"]=="Mobile"].groupby("Name")[co_to_mea].mean()
        mob_s=sh_d[sh_d["Type"]=="Mobile"].groupby("Name")[co_to_sum].sum()
    # now we merge 2 datasets    
        me_sum_ave=pd.merge(mob_t, mob_s, on='Name', how="left").merge(pau_permin,on="Name", how="left")
        
    ###count mobile tracks and stationary tracks
    #sort by mobile and stationary then count 
        co_mobile=sh_d[sh_d["Type"]=="Mobile"].groupby("Name")["Type"].count().reset_index().rename(columns={'Type': 'co_mobile'})
        co_stat=sh_d[sh_d["Type"]=="Stationary"].groupby("Name")["Type"].count().reset_index().rename(columns={'Type': 'co_stat'})
    #sort by ant and ret then count bcs only moving tracks have direction label, stat tracks dont
        co_ant=sh_d[sh_d["Direction"]=="Anterograde"].groupby("Name")["Direction"].count().reset_index().rename(columns={'Direction': 'co_ant'})
        co_ret=sh_d[sh_d["Direction"]=="Retrograde"].groupby("Name")["Direction"].count().reset_index().rename(columns={'Direction': 'co_ret'})
        co_reversal=sh_d[sh_d["Type"] == "Mobile"].groupby("Name")["Reversals"].apply(lambda x: (x > 0).sum()).reset_index().rename(columns={"Reversals": "count_reversals"})
        
    
        
    #merge previous df into one based on name of the tracks
        me_mobility=pd.merge(me_sum_ave, co_mobile, on='Name', how="left").merge(co_stat, on="Name", how="left") #merge mobilitiy count based on left edge
        me_final=pd.merge(me_mobility, co_ant, on="Name", how="left").merge(co_ret, on="Name", how="left") # merge directionality 
        me_final=pd.merge(me_final, co_reversal, on="Name", how="left")
    #now calculate percentage of mobility and directionality per cell
        me_final["total_track"]=me_final["co_mobile"].fillna(0)+me_final["co_stat"].fillna(0)
        me_final["per_mobile"]=(me_final["co_mobile"].fillna(0)/me_final["total_track"])*100
        me_final["per_stat"]=(me_final["co_stat"].fillna(0)/me_final["total_track"])*100
        me_final["per_ant"]=(me_final["co_ant"].fillna(0)/me_final["co_mobile"].fillna(0))*100
        me_final["per_ret"]=(me_final["co_ret"].fillna(0)/me_final["co_mobile"].fillna(0))*100
        me_final["per_reversals"]=(me_final["count_reversals"].fillna(0)/me_final["co_mobile"].fillna(0))*100
        
        f_re[ks]=me_final
  
    if label==1:
        with pd.ExcelWriter (output_file, engine='openpyxl') as writer: 
            #writer=pd.ExcelWriter(output_file, engine='openpyxl')
            for sh_name, sh_data in f_re.items():      #d_n is name of directory, di is the sub_dictionary
                #print(sh_name)
                df=pd.DataFrame.from_dict(sh_data)
                df.to_excel(writer, sheet_name=sh_name)        
        #writer.save()  
    
    #else:
    #    e_name=str(e_f.split(".")[0])
    #    writer_s=pd.ExcelWriter(r"Q:\Yuhao\Images\AIS project\AIS_trafficking\trafficking analysis\Rab3-GFP_kymo\results_cells"+"\\"+e_name+"_post.xlsx", engine='openpyxl')
    #    for sh_name, sh_data in f_re.items():      #d_n is name of directory, di is the sub_dictionary
        #print(sh_name)
    #        df=pd.DataFrame.from_dict(sh_data)
    #        df.to_excel(writer_s, sheet_name=sh_name)        
    #    writer_s.save()  
        

    # if label>1:
    #     excel_file=r".xlsx"
    #     o_df=pd.read_excel(excel_file, sheet_name=None)
    #     writer_a=pd.ExcelWriter(excel_file, engine='openpyxl', mode="w")
    #     for sh_name, sh_data in f_re.items():      
    # #        print(sh_name)
    #         n_df=pd.DataFrame.from_dict(f_re[sh_name])
    #         if sh_name in o_df:
    #             print(sh_name)
    #             o_df[sh_name]=o_df[sh_name].append(n_df, ignore_index=True)
    #         if not sh_name in o_df:
    #             print("new sheet",sh_name)
    #             o_df[sh_name]=n_df
    #         o_df[sh_name].to_excel(writer_a, sheet_name=sh_name, index=False)
    
    #     writer_a.save()
    #     writer_a.close() 
    
    
    f_re={}
    
    
    
    
    
    
    
    
    

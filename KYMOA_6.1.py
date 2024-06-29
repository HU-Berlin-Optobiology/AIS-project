#!/usr/bin/env python
# coding: utf-8

# In[3]:

########## User input for file directories
figure_save=r"" # type in the directory you want to save the final plot (e.g. r"c/s/c")
output_save=r"" # type in the directory you want to save the final output as excel sheet
file_directory=r"" # type in the directory where your tracks are
########## User input for imaging settings
pixel_size=0.13    #pixel size of the microscope (micro)
frame=0.2          #time for time laps  (s)
########## User input for defining runs and random movements
run_thr=5          # threshold to define run (number of pixels)
pause_thr=50         # upper threshold to define pause and stationary (number of frames)
pause_thr_l=3       # lower threshold to define pause or immobile
displacement_thr=0   #threshold to define the direction of the whole track

########## User input for file names
k_w="Resli"       # provide file name
exp_name="20230228"  # name of your experiment
Objective="100x"

# In[ ]:


import os
import csv
import math
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
from itertools import zip_longest

save_fig=figure_save           # type in where u wanna save the figures
save_exe=output_save   # type in where to save the final excel sheet
dir=file_directory  # directory to find your csv files
os.chdir(dir)
nu_files=0     #cunt how many files are processed
kymo_ana={}    #store all track information
#kymo={}        #store number of moving/stationary/total tracks
#kymo_info={}   #to store kymo in case multiple kymograph is analysed
total_track={} # integrate multiple tracks for final plotting 
f_roots=[]  #to store all the directories that contains track csv file
log_file={}

log_file["Directory (tracks):"]=str(dir)
log_file["Dirctory saving:"]=str(save_exe)
log_file["Pixel size:"]=pixel_size
log_file["Frame rate:"]=frame
log_file["Run threshold:"]=run_thr
log_file["Upper pausing threshold:"]=pause_thr
log_file["Lower pausing threshold:"]=pause_thr_l
log_file["Displacement threshold:"]=displacement_thr
log_file["Ovjective:"]=Objective


st_tra=0  #cunting stationary tracks
mo_tra=0  #cunting moving tracks
a_track=0
r_track=0
nu_retr=[]
nu_antr=[]

#######First, go through the entire directory and list all the subdirectories 
#######contain track.csv files
print("Checking files...")
for roots, directory, files in os.walk(dir):
    for f in files:        
        if f.endswith(".csv") and f[0:5]==k_w:
            if not roots in f_roots:
                f_roots.append(roots)

print("No errors detected.")
#######Second, now we have all the subdirectories and we loop through each of
#######to process each track#####################################################################################
print("Start data analysis")
for rs in f_roots:               #loop each directory in f_roots that contains track csv files
    os.chdir(rs)  # go to the directory
    file_list = os.listdir(rs)     # list all the files in the current directory
    kymo_ana[rs]={} #creat a dictionary to store the track info from current subdir
    #print(kymo_ana)
    for fi in file_list:         #loop all the files in current subdirectory
        if fi.endswith(".csv") and fi[0:5]==k_w:  #if file name meets the requirment, open this track file for processing
            #print(fi)
            nu_files+=1  #count file number being processed                                           
            fig_n=str(fi[0:-4]+".svg")  #extract name of current track
            
            print("Directory of the track: ", rs)
            print("Current track being processed", fig_n)
            print("Checking coordinates...")
#####now we open the coordinates of each track and check if there are any artifacts e.g. reverse of time.... 
#####unphysical errors will be removed     
       
            with open(fi, newline='') as csvfile:
                data_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
                next(data_reader)  # Skip header row if there is one
                fi_d = []
                for row in data_reader:
                    x=row[0]                #row[0] is the first row
                    y=row[1]                #row[1] is the second row
                    if not any(y.endswith(suffix) for suffix in ['.250', '.500', '.750']):
                        fi_d.append(row)
            for i in range(len(fi_d)):
                x, y = float(fi_d[i][0]), float(fi_d[i][1])
                if x<=x+0.75:
                    x=math.floor(x)
                else:
                    x=math.floor(x)+1
                fi_d[i][0] = x
                fi_d[i][1] = math.floor(y)
        
            fi_df = {}
            #for i in fi_d:
            #    x, y = i
            #    if y not in fi_df:
            #        fi_df[y] = [x, y]
            
            for i in range(len(fi_d)-1):
                x,y=fi_d[i]
                nx,ny=fi_d[i+1]
                if y not in fi_df and ny>y:
                    fi_df[y]=[x,y]
            lxy=fi_d[len(fi_d)-1]
            last_key=sorted(fi_df.keys())[-1]
            if lxy[1]>fi_df[last_key][1]:
                fi_df[lxy[1]]=lxy        
            fi_df = list(fi_df.values())
            print("No errors detected in the track.")
#############################################################################################################################
#### define pauses, movements & direction of movements in current track

            print("Detecting pauses, movements & direction of movements...")
            anterograde = [] # list to store anterograde movement data
            pausing = [] # list to store pausing data
            retrograde = [] # list to store retrograde movement data

            antero_track  = [] # list of all anterograde movement coords
            pausing_track = [] # list of all stationary coords
            retro_track  = [] # list of all retrograde movement coords
            
            coordinates={} # to store all coordinates of one track

            label = 0   # identifier of retrograde/stationary/anterograde
            #### variable "label" tells you the dx value in the previous step
            #### if the previous dx > 0, then label is defined as 3
            #### elif the previous dx == 0, the for loop sets label to 2
            #### elif the previous dx < 0, the for loop sets label to 1
            for i in range(len(fi_df)-1):
                x1, y1 = fi_df[i]
                x2, y2 = fi_df[i+1]
                dx=x2-x1
                dy=y2-y1
                
                if dx>0:
                    if label == 3 or label == 0:
                        anterograde.append(fi_df[i])
                    elif label == 2:
                        pausing.append(fi_df[i])
                        pausing_track.append(pausing)
                        anterograde.append(fi_df[i])
                        pausing=[]
                    elif label == 1:
                        retrograde.append(fi_df[i])
                        retro_track.append(retrograde)
                        anterograde.append(fi_df[i])
                        retrograde=[]
                    label = 3
        
    
                elif dx == 0:
                    if label == 2 or label == 0:
                        pausing.append(fi_df[i])
                    elif label ==3:
                        anterograde.append(fi_df[i])
                        antero_track.append(anterograde)
                        pausing.append(fi_df[i])
                        anterograde=[]
                    elif label == 1:            
                        retrograde.append(fi_df[i])            
                        retro_track.append(retrograde)
                        pausing.append(fi_df[i])
                        retrograde=[]        
                    label = 2
        
                elif dx<0:
                    if label == 1 or label == 0:
                        retrograde.append(fi_df[i])
                    elif label ==3:
                        anterograde.append(fi_df[i])
                        antero_track.append(anterograde)
                        anterograde=[]
                    elif label ==2:
                        pausing.append(fi_df[i])
                        pausing_track.append(pausing)
                        retrograde.append(fi_df[i])
                        pausing=[]
                    label = 1

####in case the end of the track is continous run/pause add the end of the track############################## 
            pausing_track.append(pausing)
            antero_track.append(anterograde)
            retro_track.append(retrograde)
####to remove the empty lists###########################################
            antero_track=[sublist for sublist in antero_track if sublist]
            retro_track=[sublist for sublist in retro_track if sublist]
            pausing_track=[sublist for sublist in pausing_track if sublist]
            
###adding last coordinate accordingly############################
            x_last, y_last=fi_df[-1]
            x_bl, y_bl=fi_df[-2]
            dx_last=x_last-x_bl
            if dx_last==0:
                if len(pausing_track)>0:
                    pausing_track[len(pausing_track)-1].append([x_last,y_last])
                else:
                    pausing_track.append([[x_bl, y_bl],[x_last,y_last]])
                    
                #print("pausing track",pausing_track)
            elif dx_last<0:
                if len(retro_track)>0:
                    retro_track[len(retro_track)-1].append([x_last,y_last])
                else:
                    retro_track.append([[x_bl, y_bl],[x_last,y_last]])
                    
                #print("retrograde:",retro_track)
            elif dx_last>0:
                if len(antero_track)>0:
                    antero_track[len(antero_track)-1].append([x_last,y_last])
                else:
                    antero_track.append([[x_bl, y_bl],[x_last,y_last]])
###now adding all tracks into the dictionary called coordinates###############################
            print("Storing track information...")
            coordinates["original_track"]=fi_df
            coordinates["anterograde_track"]=antero_track
            coordinates["retrograde_track"]=retro_track
            coordinates["pausing/stationary"]=pausing_track
            #print("anterograde tracks: ", antero_track, "there are ", len(antero_track), "tracks")     #check
            #print("")
            #print("retrograde tracks: ", retro_track, "there are ", len(retro_track), "tracks")       #check
            #print("")
            #print("stationary tracks: ", pausing_track, "there are ", len(pausing_track), "tracks")     #check
            #print("last coordinate: ",x_last, y_last)
            #print("second last coordinate: ", x_bl, y_bl)
######calculating displecement for the current track###############################################################
            print("Calculating displacement...")
            start=int(fi_df[0][0])      #first x coord of the track
            end=int(fi_df[len(fi_df)-1][0])      #last x coord of the track
            dis=end-start               #calculate displacement
            #print("The displacement of this track is (# pixel): ",dis)
            
###################################################################################################################
### Now we calculate individual run length/duration etc
    
            antero_run = 0   #for calculating distance of each anterograde run
            retro_run = 0    #for calculating distance of each retrograde run

            a_run_dur=0      #for calculating anterograde run time
            r_run_dur=0      #for calculating retrograde run time
            pau_dur=0        #for calculating pausing time

            n_ant_r=0
            n_ret_r=0
            n_pause=0
            
            rev_ant=0    #counting reversals to anterograde direction
            rev_ret=0    #counting reversals to retrograde direction
            reversal=0   #counting reversals respect to overall direction of the track
    
            ave_antero_dur=0     #calculate average run time    
            ave_antero_run=0    #calculating average run length    
            ave_antero_v=0
            ave_retro_dur=0     #calculate average run time    
            ave_retro_run=0    #calculating average run length
            ave_retro_v=0
            ave_pausing_time=0
            
            stationary_dur=0
            
            dur_total=0   #for calculating total duration
            Timerunning=0  #for total running time
            Timepausing=0  #for total pausing time
            PercentTimePausing=0  #for percentage of pausing time
            PercentTimeRunning=0
            t_ant_run_dur=0
            t_ret_run_dur=0
            PercentTimeantrun=0
            PercentTimeretrun=0
            T_travel_len=0

            ind_antero_run=[]      #storing distance of each individual anterograde run
            ind_retro_run=[]       #storing distance of each individual retrograde run

            pausing_dur=[]         #store individual pausing time
            antero_dur = []        #store time of each individual anterograde run
            retro_dur = []         #store time of each individual retrograde run

            ind_antero_v=[]        #individual anterograde velocity
            ind_retro_v=[]         #individual retrograde velocity
            
#########for passive movement#####################################################            
            crowling=[]
            n_sub=[]
            c_label=0
            pas_len=0
            passive_dur=0
            passive_dis=0
            passive_vel=0
            PercentTimePassive=0
            
####################################################################################
            
            parameters={}          #creat a dictionary to store all parameters of current track
            
            if not any (len(tracks)>=run_thr for tracks in antero_track+retro_track): #define stationary track
                st_tra+=1 
                stationary=True
                print("No runs detected, this is a stationary track")
                ####calculate stationary time
                stationary_dur=(int(fi_df[len(fi_df)-1][1])-int(fi_df[0][1]))*frame #last coord y - first coord y
                ####store track info
                parameters["Type"]="Stationary"
                parameters["indAntR (um)"]="N/A"
                parameters["indAntD (s)"]="N/A"
                parameters["indAntV (um/s)"]="N/A"
                parameters["aveAntD (s)"]="N/A"
                parameters["aveAntR (um)"]="N/A"
                parameters["aveAntV (um/s)"]="N/A"          
                parameters["# of anterograde run"]="N/A"
                parameters["indRetR (um)"]="N/A"
                parameters["indRetD (s)"]="N/A"
                parameters["indRetV (um/s)"]="N/A"
                parameters["aveRetD (s)"]="N/A"
                parameters["aveRetR (um)"]="N/A"
                parameters["aveRetV (um/s)"]="N/A"
                parameters["# of retrograde run"]="N/A"
                parameters["PauseD (s)"]="N/A"
                parameters["avePauseD (s)"]="N/A" 
                parameters["# of pause"]="N/A"
                parameters["Total_Run_D (s)"]="N/A"
                parameters["Total_Pausing_D (s)"]="N/A"
                parameters["PercTimeAntRun (%)"]="N/A"
                parameters["PercTimeRetRun (%)"]="N/A"
                parameters["total # of run"]="N/A"
                parameters["Stationary_dur (s)"]=stationary_dur
                parameters["Displacement"]=dis
                parameters["Reversals"]="N/A"
                parameters["Direction"]="N/A"
                
                print("Total duration (s)", stationary_dur)
                print("Displacement (pixels):", dis)
                print("")       
            
            ####define moving tracks ############################################################################ 
            else:
                mo_tra+=1
                parameters["Type"]="Mobile"
                for track in range(len(antero_track)):        # loop each track in the list of anterograde tracks
                    #print("current track: ", antero_track[track])    #check
                    #print("# of coords in track: ", len(antero_track[track]))   #check        
                    if len(antero_track[track])>=run_thr:     #movement larger than pixel threshold is defined as run
                        #print("track: ",antero_track[track])    #check
                        #print("# of the coords in track: ", len(antero_track[track]))     #check        
                        for coords in range(len(antero_track[track])-1):   #loop each coords in current track
                            #print(antero_track[track][coords])        #check
                            xa, ya = antero_track[track][coords]      #define x1,y1
                            #print(xa,ya)                             #check
                            xb, yb = antero_track[track][coords+1]    #define x2,y2
                            dx, dy = xb - xa, yb - ya                 #calculate dx and dy
                            unit_length = math.sqrt(dx**2 + dy**2)    #calculate distance between 2 adjcent coords
                            ant_run_time=dy                           #calculate time between 2 adjcent coords
                            antero_run += unit_length*pixel_size      #measure distance of entire current run in microns      
                            a_run_dur +=ant_run_time*frame            #same for time in second
                        #print("length of current run: ", antero_run) #check
                        ind_antero_run.append(antero_run)             #add run length of current track to list
                        antero_dur.append(a_run_dur)                  #add run duration of current track to list
                        antero_v=float(antero_run)/float(a_run_dur)   #calculate velocity of current run
                        ind_antero_v.append(antero_v)                 #add velocity of current run to list           
                        antero_run=0                                  #clear run length of current track for the next one
                        a_run_dur=0                                   #same for time
                        rev_ant+=1                                    #count how many anterograde run calculated 
                        n_ant_r+=1
                        ave_antero_dur=sum(antero_dur)/len(antero_dur)     #calculate average run time    
                        ave_antero_run=sum(ind_antero_run)/len(ind_antero_run)    #calculating average run length    
                        ave_antero_v=sum(ind_antero_v)/len(ind_antero_v)   #calculating average run speed 
                                                                           #number for average will keep updating until the
                                                                           #listing is finished                      
                parameters["indAntR (um)"]=ind_antero_run
                parameters["indAntD (s)"]=antero_dur
                parameters["indAntV (um/s)"]=ind_antero_v
                parameters["aveAntD (s)"]=ave_antero_dur
                parameters["aveAntR (um)"]=ave_antero_run
                parameters["aveAntV (um/s)"]=ave_antero_v          
                parameters["# of anterograde run"]=n_ant_r
################################################################################################################################
### calculating average and individual retrograde run length and duration, same concept and structure for anterograde

                for rtrack in range(len(retro_track)):        # loop each track in the list of anterograde tracks 
                    #print("current track: ", antero_track[track])    #check
                    #print("# of coords in track: ", len(antero_track[track]))   #check    
                    if len(retro_track[rtrack])>=run_thr:
                        #print("track: ",antero_track[track])    #check
                        #print("# of the coords in track: ", len(antero_track[track]))     #check   
                        #re_run_dur=len(retro_track[track])*0.5     #run time. make 0.2 as user input
                        #retro_dur.append(re_run_dur)               #recording duration of each run
                        for coords in range(len(retro_track[rtrack])-1):
                            #print(antero_track[track][coords])        #check
                            xar, yar = retro_track[rtrack][coords]
                            #print(xa,ya)                             #check
                            xbr, ybr = retro_track[rtrack][coords+1]
                            dxr, dyr = xbr - xar, ybr - yar
                            runit_length = math.sqrt(dxr**2 + dyr**2)
                            ret_run_time=dyr
                            retro_run += runit_length*pixel_size  
                            r_run_dur += ret_run_time*frame
                        #print("length of current run: ", antero_run) #check
                        ind_retro_run.append(retro_run)
                        retro_dur.append(r_run_dur)
                        retro_v=float(retro_run)/float(r_run_dur)
                        ind_retro_v.append(retro_v)
                        retro_run=0
                        r_run_dur=0
                        rev_ret+=1
                        ave_retro_dur=sum(retro_dur)/len(retro_dur)     #calculate average run time    
                        ave_retro_run=sum(ind_retro_run)/len(ind_retro_run)    #calculating average run length
                        ave_retro_v=sum(ind_retro_v)/len(ind_retro_v)
                        n_ret_r+=1
            
                parameters["indRetR (um)"]=ind_retro_run
                parameters["indRetD (s)"]=retro_dur
                parameters["indRetV (um/s)"]=ind_retro_v
                parameters["aveRetD (s)"]=ave_retro_dur
                parameters["aveRetR (um)"]=ave_retro_run
                parameters["aveRetV (um/s)"]=ave_retro_v
                parameters["# of retrograde run"]=n_ret_r
################################################################################################################################
# calculating average/individual pausing time

                for p_track in range(len(pausing_track)):           # loop each track in the list of anterograde tracks
                    #print("current track: ", pausing_track[p_track])    #check
                    #print("# of coords in track: ", len(pausing_track[p_track]))   #check
                    if pause_thr_l<len(pausing_track[p_track]) <=pause_thr:
                        for coords in range(len(pausing_track[p_track])-1):
                            #print(antero_track[track][coords])        #check
                            xap, yap = pausing_track[p_track][coords]
                            #print(xa,ya)                             #check
                            xbp, ybp = pausing_track[p_track][coords+1]
                            dxp, dyp = xbp - xap, ybp - yap            
                            pau_time=dyp
                            pau_dur += pau_time*frame
                        #print("length of current run: ", antero_run) #check
                        pausing_dur.append(pau_dur)        
                        pau_dur=0
                        ave_pausing_time=sum(pausing_dur)/len(pausing_dur)
                        n_pause+=1

                parameters["PauseD (s)"]=pausing_dur
                parameters["avePauseD (s)"]=ave_pausing_time 
                parameters["# of pause"]=n_pause
##############################################################################################################################                
                Timerunning=sum(antero_dur)+sum(retro_dur)
                Timepausing=sum(pausing_dur)
                t_ant_run_dur=sum(antero_dur)
                t_ret_run_dur=sum(retro_dur)
                PercentTimeantrun=(t_ant_run_dur/Timerunning)*100
                PercentTimeretrun=(t_ret_run_dur/Timerunning)*100    
###calculate total traveling distance#########################################################################################                
                for i in range(len(fi_df)-1):
                    xtl1, ytl1 = fi_df[i]
                    xtl2, ytl2 = fi_df[i+1]
                    dxtl, dytl = xtl2 - xtl1, ytl2 - ytl1
                    seg_len = math.sqrt(dxtl**2 + dytl**2)
                    T_travel_len += seg_len*pixel_size
                               
                parameters["Total_Run_D (s)"]=Timerunning
                parameters["Total_Pausing_D (s)"]=Timepausing
                parameters["PercTimeAntRun (%)"]=PercentTimeantrun
                parameters["PercTimeRetRun (%)"]=PercentTimeretrun
                parameters["total # of run"]=int(n_ant_r+n_ret_r)
                parameters["Total travel len (um)"]=T_travel_len
                #parameters["% of ant run (%)"]=(n_ant_r/int(n_ant_r+n_ret_r))*100
                #parameters["% of ret run (%)"]=(n_ret_r/int(n_ant_r+n_ret_r))*100
                #nu_retr.append(n_ret_r)
                #nu_antr.append(n_ant_r)                
####passive movement and calculation############################################################################################
                print("Processing passive movement of current track...")
                for r, a, p in zip_longest (retro_track, antero_track, pausing_track, fillvalue='N/A'):
                    #print("anterograde",a, len(a))
                    #print("pausing",p, len(p))
                    #print("retrograde", r, len(r))
                    if 0<len(p)<=pause_thr_l and p!="N/A":
                        #print("below running thre:",p, len(p), "append to crowling for anterograde")
                        crowling.append(p)
                    if 0<len(a)<run_thr and a!="N/A":
                        #print("below pausing thre",a, len(a), "append to crowling for pause")
                        crowling.append(a)
                    if 0<len(r)<run_thr and r != "N/A":
                        #print("below running thre (retro)", r, len(r), "append to crowling for retrograde")
                        crowling.append(r)
                    else:
                        print("Processive movement")
                        
                def get_y_coordinate(coord):
                    return coord[1]
                crowling.sort(key=lambda coord: get_y_coordinate(coord[0]))
                #print("passive movement",crowling)
                #print("length", len(crowling))
                
                if len(crowling)==1:
                    n_sub=crowling
                
                if len(crowling)>1:
                    for cm in range(len(crowling)-1):
                        #print("cm:",crowling[cm],",","cm+1:", crowling[cm+1])
                        #print("cy1",crowling[cm][-1], "cy2",crowling[cm+1][0])
                        cx1,cy1=crowling[cm][-1]
                        cx2,cy2=crowling[cm+1][0]
                        dcx2=cx2-cx1
                        dcy2=cy2-cy1
                        if dcy2==0 and dcx2==0:
                            if c_label==0:
                                #print("same coord detected, append:",crowling[cm],"and",crowling[cm+1][1:len(crowling[cm+1])], "to new sublist")
                                #print("new list will be:", crowling[cm], crowling[cm+1][1:len(crowling[cm+1])])
                                n_sub.append(crowling[cm][0:len(crowling[cm])]+crowling[cm+1][1:len(crowling[cm+1])])
                                #print("new list will be:",n_sub)
                            if c_label!=0:
                                #print("same coord detected, append:", crowling[cm+1][1:len(crowling[cm+1])], "to new sublist")
                                #print("new list will be:", crowling[cm], crowling[cm+1][1:len(crowling[cm+1])])
                                n_sub[len(n_sub)-1]=n_sub[len(n_sub)-1]+crowling[cm+1][1:len(crowling[cm+1])]
                                #print("new list will be:",n_sub)
                            c_label+=1
                    
                        elif dcy2==0 and dcx2!=0:
                            if c_label==0:
                                #print("adjcent coord detected, append:",crowling[cm],"and",crowling[cm+1][1:len(crowling[cm+1])], "to new sublist")
                                #print("new list will be:", crowling[cm], crowling[cm+1][1:len(crowling[cm+1])])
                                n_sub.append(crowling[cm][0:len(crowling[cm])]+crowling[cm+1][1:len(crowling[cm+1])])
                                #print("new list will be:",n_sub)
                            if c_label!=0:
                                #print("adjcent coord detected, append:", crowling[cm+1][1:len(crowling[cm+1])], "to new sublist")
                                #print("new list will be:", crowling[cm], crowling[cm+1][1:len(crowling[cm+1])])
                                n_sub[len(n_sub)-1]=n_sub[len(n_sub)-1]+crowling[cm+1][1:len(crowling[cm+1])]
                                #print("new list will be:",n_sub)
                            c_label+=1
                    
                        elif dcy2>1:
                            if c_label==0:
                                #print("different coord detected, append: ",crowling[cm],"and", crowling[cm+1])
                                n_sub.append(crowling[cm])
                                n_sub.append(crowling[cm+1])
                                #print("the new list will be: ", n_sub)
                            if c_label!=0:
                                #print("different coord detected, append: ", crowling[cm+1])
                                n_sub.append(crowling[cm+1])
                                #print("the new list will be: ", n_sub)
                            c_label+=1
###now calculate the parameters for passive movement
                #print("sorted passive:",n_sub)
                for w in range(len(n_sub)):
                    for cw in range(len(n_sub[w])-1):
                        #print(n_sub[w][cw], n_sub[w][cw+1])
                        cwx1,cwy1=n_sub[w][cw]
                        cwx2,cwy2=n_sub[w][cw+1]
                        dcwx=cwx2-cwx1
                        dcwy=cwy2-cwy1
                        pas_len=math.sqrt(dcwx**2+dcwy**2)
                        passive_dur+=dcwy*frame
                        passive_dis+=pas_len*pixel_size
                        passive_vel+=passive_dis/passive_dur
                
                #print("passive movement duration (sec):", passive_dur)
                #print("passive movement distance (um):", passive_dis)
                #print("passive movement velocity (um/s):", passive_vel)
                
                parameters["psv_dur (sec)"]=passive_dur
                parameters["psv_dis (um)"]=passive_dis
                parameters["psv_vel (um/s)"]=passive_vel
###############################################################################################################################
                #print("Anterograde run")
                #print("There is/are ", len(ind_antero_run), "anterograde run")
                #print("Length of each anterograde run (microns): ", ind_antero_run)
                #print("run duration of each track (s): ", antero_dur,"# of value: ", len(antero_dur))
                #print("average anterograde run time (s): ", ave_antero_dur)
                #print("average anterograde run length (microns): ", ave_antero_run)
                #print("velocity of each run (um/s): ", ind_antero_v)
                #print("average velocity of anterograde run (um/s): ", ave_antero_v)

                #print("")

                #print("Retrograde run")
                #print("There is/are ", len(ind_retro_run), "retrograde run")
                #print("Length of each retrograde run (microns): ", ind_retro_run)
                #print("run duration of each track (s): ", retro_dur, "# of value: ", len(retro_dur))
                #print("average retrograde run time (s): ", ave_retro_dur)
                #print("average retrograde run length (microns): ", ave_retro_run)
                #print("velocity of each run (um/s): ", ind_retro_v)
                #print("average velocity of retrograde run (um/s): ", ave_retro_v)

                #print("")
 
                #print("Pausing")
                #print("There is/are ", len(pausing_dur), "pauses")
                #print("pausing time of individual pausing track (s): ", pausing_dur)
                #print("average pausing (s): ", ave_pausing_time)
                #print("")
                #print("Current track: ", f[0:len(f)-4])
                #print("")
                #print("total travelling length: ", T_travel_len)
#########adding displacement to parameters of each track##############################################
                parameters["Displacement"]=dis
    
###reversal#####################################################################
#first define which direction this vesicle is moving by calculating displacement
#then count how many retrograde/anterograde runs are in this track in respective of the final direction
#print(fi_df[len(fi_df)-1])


                if dis>displacement_thr:  #when the end of the line is 10 pixels away to the anterograde direction, it is considered as anterograde track
                    reversals=rev_ret    #reversals will therefore be the number of retrograde runs
                    print("this is an anterograde track.", "There are/is", reversals, "direction reversals." )
                    parameters["Reversals"]=reversals
                    parameters["Direction"]="Anterograde"
                    r_track+=1
        
                elif dis<displacement_thr:      #when the end of the line is 10 pixels away 
                                                 #to the retrograde direction, it is considered as retrograde track
                    reversals=rev_ant         #the reversals is therefore be the number of anterograde runs
                    print("this is a retrograde track. ", 
                          "There are/is", reversals, "direction reversals." )
                    parameters["Reversals"]=reversals
                    parameters["Direction"]="Retrograde"
                    a_track+=1
        
########## percentage of pausing time and running time ######################################################################

            dur_total=(fi_df[len(fi_df)-1][1]-fi_df[0][1])*frame #use the last y coord of the 
                                                                 #track to calculate entire recording time duration             
            PercentTimePausing=(Timepausing/dur_total)*100
            PercentTimeRunning=(Timerunning/dur_total)*100
    
            parameters["TotalD (s)"]=dur_total
            parameters["PercTimePausing (%)"]=PercentTimePausing
            parameters["PercTimeRunning (%)"]=PercentTimeRunning
#########################################################################################################################            
            PercentTimePassive=(passive_dur/dur_total)*100
            parameters["PercTimepassivem (%)"]=PercentTimePassive    
            #print("")
            #print("total duration of the track is: ", dur_total,"s")
            #print("total duration of pausing: ", Timepausing,"s")
            #print("total duration of running: ", Timerunning,"s")
            #print("percentage of pausing: ", PercentTimePausing,"%")
            #print("percentage of running: ", PercentTimeRunning,"%")
 ############################################################################################################################           
### Now we plot the anterograde/retrograde/pausing tracks against the entire track to check if there are errors in the
##coordinates
            track=np.array(fi_df)
            a,b =track.T
            ret_plot=retro_track
            ant_plot=antero_track
            pas_plot=crowling
            paus_plot=[]
            sta_plot=[]
            
            print("Displaying tracks...")

            for at in ant_plot:        
                #print(t)
                ant_t=np.array(at)
                #print(ant_t)
                ax,ay=ant_t.T
                #print(x,y)
                plt.plot(ax,-ay, color="blue", linewidth=1, label="anterograde")
                if len(at) >= run_thr:
                    ant_r_t=np.array(at)
                    arx,ary=ant_r_t.T
                    plt.plot(arx,-ary, color="pink", linewidth=3, label="ant_run")
    
            for rt in ret_plot:
                #print(rt)
                ret_t=np.array(rt) 
                #print(ret_t)
                rx,ry=ret_t.T
                #print(rx,ry)
                plt.plot(rx,-ry, color="green", linewidth=1, label="retrograde")
                if len(rt) >= run_thr:
                    ret_r_t=np.array(rt)
                    rrx,rry=ret_r_t.T
                    plt.plot(rrx,-rry, color="orange", linewidth=3, label="ret_run")    
            for passive in pas_plot:
                pas_t=np.array(passive) 
                pasx, pasy=pas_t.T
                plt.plot(pasx,-pasy, color="black", linewidth=4, label="passive move")
            
                    
            if not any(len(trs)>=run_thr for trs in ant_plot+ret_plot):
                sta_plot=pausing_track  
                for stt in sta_plot:
                    sta_t=np.array(stt)
                    sx,sy=sta_t.T
                    plt.plot(sx,-sy, color="red", linewidth=1, label="stationary")
                
            else:
                paus_plot=pausing_track
                for pt in paus_plot:
                    pau_t=np.array(pt) 
                    px,py=pau_t.T
                    plt.plot(px,-py, color="red", linewidth=1, label="pausing")
                    if pause_thr_l< len(pt) <= pause_thr:
                        pau_p_t=np.array(pt)
                        ppx,ppy=pau_p_t.T
                        plt.plot(ppx,-ppy, color="grey", linewidth=3, label="pause")
    
            plt.plot(a, -b, color="black", linewidth=0.5)
            Ante = mlines.Line2D([], [], color='blue', label='Anterograde')
            Retr = mlines.Line2D([], [], color='green', label='Retrograde')
            Imb = mlines.Line2D([], [], color='red', label='Stationary')
            Ant_run = mlines.Line2D([], [], color='pink', label='Ant run')
            Ret_run= mlines.Line2D([], [], color='orange', label='Ret run')
            Pause= mlines.Line2D([], [], color='grey', label='Pause')
            passive_m=mlines.Line2D([], [], color='black', label='passive')
            plt.legend(handles=[Ante, Retr, Imb, Ant_run, Ret_run, Pause, passive_m])
            plt.title(fi[0:(len(fi)-4)])
            plt.xlabel("Distance (# of pixels)")
            plt.ylabel("Time (# of frames)")
            plt.savefig(save_fig+"\\"+fig_n)    
            plt.show()
            print("###current track finished")
            print("")
################################################################################################################################
###saving all data into a dictionary contains all tracks information and how many moving or stationary tracks    
    
            kymo_ana[rs][str(fi[0:-4])]=parameters
            total_track[str(fi[0:-4])]=coordinates
            
            
            nu_retr.append(n_ret_r)
            nu_antr.append(n_ant_r)


#print("mo_tra", mo_tra)
#print("st_tra", st_tra)
print("In total ", nu_files,"were processed")
#print(directory)
#kymo["# of moving tracks"]=mo_tra
#kymo["# of stationary tracks"]=st_tra
#kymo["Total # of tracks"]=nu_files
#kymo["# of anterograde tracks"]=a_track
#kymo["# of retrograde tracks"]=r_track
#kymo["# anterograde run"]= sum(nu_antr)
#kymo["# retrograde run"]= sum(nu_retr)
#kymo["Total # of runs"]=sum(nu_antr)+sum(nu_retr)
#kymo["% moving tracks"]=(mo_tra/nu_files)*100
#kymo["% stationary tracks"]=(st_tra/nu_files)*100
#kymo["% anterograde track"]=(a_track/mo_tra)*100
#kymo["% retrograde track"]=(r_track/mo_tra)*100
#kymo["% anterograde run"]= (sum(nu_antr)/(sum(nu_antr)+sum(nu_retr)))*100
#kymo["% retrograde run"]=(sum(nu_retr)/(sum(nu_antr)+sum(nu_retr)))*100
#kymo_info[exp_name]=kymo
log_file["Number of processed tracks: "]=nu_files
####
st_tra=0  #cunting stationary tracks
mo_tra=0  #cunting moving tracks
nu_files=0

#t_df=pd.DataFrame(kymo_ana)
#t_df.to_excel(save_exe+"\\kymo_ana.xlsx", index=[0])
#t_df.to_csv(save_exe+"tra_info.txt", sep='\t',index=[0])
print("Saving results as excel")

writer=pd.ExcelWriter(save_exe+"\\"+exp_name+"_kymo_ana.xlsx", engine='openpyxl')

for sh_n, sh_data in kymo_ana.items():    #d_n is name of directory, di is the sub_dictionary
    #print(sh_n)
    sh_n=sh_n.split("\\")[-2]
    #print(sh_n)
    
    #sheet=writer.book.add_worksheet(sh_n[-4:0])
    #print(sh_n[0:-4])
    df=pd.DataFrame.from_dict(sh_data, orient='index')
    df.to_excel(writer, sheet_name=sh_n)
    #writer.book.create_sheet(sh_n[-4:-1])
#writer.save()   
#writer.close()



#k_df=pd.DataFrame(kymo_info)
#k_df.to_excel(save_exe+"\\kymo_info.xlsx", index=[0])





# In[ ]:


###to plot everything together
print("Generating final plots...")

ks=[]
or_t=[]
aan_t=[]
rre_t=[]
pps_t=[]
sst_t=[]

for t in total_track:
    #print(t)
    #print(total_track[t])
    for k in total_track[t].keys():
        #print(k)
        ks.append(k)
    #print(ks)
    or_t.append(total_track[t][ks[0]])
    aan_t.append(total_track[t][ks[1]])
    rre_t.append(total_track[t][ks[2]])
    if not any (len(l) >= run_thr for l in total_track[t][ks[1]]+total_track[t][ks[2]]):
        sst_t.append(total_track[t][ks[3]])
    else:
        pps_t.append(total_track[t][ks[3]])

for trs in or_t:
    tracks=np.array(trs)
    ta,tb=tracks.T
    plt.plot(ta,-tb, color="black", linewidth=0.5)

for atrs in aan_t:
    for acoords in atrs:
        att=np.array(acoords)
        ata, atb=att.T
        plt.plot(ata,-atb, color="blue", linewidth=0.5, label="anterograde")
        if len(acoords) >= run_thr:
                    att_r=np.array(acoords)
                    atra,atrb=att_r.T
                    plt.plot(atra,-atrb, color="pink", linewidth=1, label="ant_run")

for rtrs in rre_t:
    for rcoords in rtrs:
        rtt=np.array(rcoords)
        rta, rtb=rtt.T
        plt.plot(rta,-rtb, color="green", linewidth=0.5, label="retrograde")
        if len(rcoords) >= run_thr:
                    rtt_r=np.array(rcoords)
                    rtra,rtrb=rtt_r.T
                    plt.plot(rtra,-rtrb, color="orange", linewidth=1, label="ret_run")

for ssts in sst_t:
    for scoords in ssts:
        sstt=np.array(scoords)
        sta, stb=sstt.T
        plt.plot(sta,-stb, color="red", linewidth=0.5, label="retrograde")                   
                    
for ptrs in pps_t:
    for pcoords in ptrs:
        pstt=np.array(pcoords)
        psta, pstb=pstt.T
        plt.plot(psta,-pstb, color="red", linewidth=0.5, label="stationary")
        if pause_thr_l<len(pcoords) <= pause_thr:
                    pstt_r=np.array(pcoords)
                    pstra,pstrb=pstt_r.T
                    plt.plot(pstra,-pstrb, color="black", linewidth=1, label="pause")                    

                    
Ante = mlines.Line2D([], [], color='pink', label='Anterograde run')
Retr = mlines.Line2D([], [], color='orange', label='Retrograde run')
Pause = mlines.Line2D([], [], color='black', label='Pause')
Sta = mlines.Line2D([], [], color='red', label='Stationary')
#plt.legend(handles=[Ante, Retr, Pause, Sta])              
plt.xlabel("Distance (# of pixels)")
plt.ylabel("Time (# of frames)")
plt.savefig(save_fig+"\\"+exp_name+"_summary.svg")
plt.show()            

print("Generating final log file...")
log_f=pd.DataFrame(log_file, index=range(len(log_file)))
log_f.to_csv(save_exe+"\\"+exp_name+"_log.txt", sep='\t',index=False)
    
print("Yuhao is always the best :)")


# In[ ]:


#what do u wanna know after processing all the tracks?
# 4. date of the processing?
# 7. were there any error messages?


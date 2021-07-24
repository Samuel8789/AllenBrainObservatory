# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 14:07:35 2020

@author: sp3660
"""

#Grating analysis and plotting


from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
import os
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from PathStructure import PathStructure


ProjectCodePath, ProjectName, ProjectDataPath, ProjectRAWPath, ProjectTempRAWPath =PathStructure()
boc = BrainObservatoryCache(manifest_file=ProjectDataPath+'/boc/manifest.json')
#os.makedirs(ProjectRAWPath+'/RAWMovies')

dirlist = []

for (root,dirs,files) in os.walk(ProjectDataPath, topdown=True):
    dirlist += dirs

#%%    
SelectedCreLine='Emx1-IRES-Cre'
SelectedStructure='VISp'
SelectedDepth=275
SelectedSession='three_session_A'
#%%  
SelectedDataDir=ProjectDataPath+'\\'+SelectedCreLine+'\\'+SelectedStructure+'\\'+str(SelectedDepth)+'\\'+SelectedSession
filelist=[]
for (root,dirs,files) in os.walk(SelectedDataDir, topdown=True):
    filelist += files 
i=0
j=0
while i<len(filelist):
    if 'events' in filelist[i]:
        filelist.pop(i)
    else:
        i+=1
while j<len(filelist):
    if 'mat' in filelist[j]:
        filelist.pop(j)
        filelist.pop(j-1)
    else:
        j+=1        
    
    
SelectedExpIDs2 = list(map(int, filelist))


SelectedExperiments=boc.get_ophys_experiments(targeted_structures=[SelectedStructure], imaging_depths=[SelectedDepth], cre_lines=[SelectedCreLine], session_types=[SelectedSession])
SelectedContainers = [sub['experiment_container_id']  for sub in SelectedExperiments]
SelectedExpIDs = [sub['id']  for sub in SelectedExperiments] 

NWBexpFiles={}
EventsFiles={}

for i in range(0, len(SelectedExpIDs)):
    FullSelectedDataDir=ProjectDataPath+'\\'+SelectedCreLine+'\\'+SelectedStructure+'\\'+str(SelectedDepth)+'\\'+SelectedSession+'\\ContID_'+str(SelectedContainers[i])+'_ExpID_'+str(SelectedExpIDs[i])
    NWBexpFiles["NWB{0}".format(SelectedExpIDs[i])] = boc.get_ophys_experiment_data(SelectedExpIDs[i],  FullSelectedDataDir+'\\'+str(SelectedExpIDs[i]))
    EventsFiles["Events{0}".format(SelectedExpIDs[i])] = boc.get_ophys_experiment_events(SelectedExpIDs[i],  FullSelectedDataDir+'\\'+str(SelectedExpIDs[i])+'events')

#%%

SelectNWB=0
SelectedExp=SelectedExpIDs[2] 
SelectedNWB= NWBexpFiles["NWB{0}".format(SelectedExp)]
SelectedEvents= EventsFiles["Events{0}".format(SelectedExp)]
    


max_proj=SelectedNWB.get_max_projection()
fig = plt.figure(figsize=(6,6))
plt.imshow(max_proj, cmap='gray')

rois = SelectedNWB.get_roi_mask_array()
print(rois.shape)
print("Number of cells:", rois.shape[0])
plt.figure(figsize=(6,6))
plt.imshow(rois.sum(axis=0))


fig = plt.figure(figsize=(10,8))
for i in range(50):
    plt.plot(2.5*SelectedEvents[i]+(i*2), color='gray')

stim_epoch = SelectedNWB.get_stimulus_epoch_table()


fig = plt.figure(figsize=(10,8))
for i in range(50):
    plt.plot(dff[i,:]+(i*2), color='gray')
    
    
    
    
#%%

selected_frames=stim_epoch[stim_epoch.stimulus=='drifting_gratings']
fig = plt.figure(figsize=(10,8))
for i in range(50):
    plt.plot(2.5*SelectedEvents[i,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]]
+(i*2), color='gray')

grat_only=SelectedEvents[i,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]]



#%%    
#for each stimulus, shade the plot when the stimulus is presented
stim_epoch.stimulus.unique()



fig = plt.figure(figsize=(10,8))
for i in range(50):
    plt.plot(SelectedEvents[i,:]+(i*2), color='gray')
    
#for each stimulus, shade the plot when the stimulus is presented
colors = {'drifting_gratings':   'blue',
          'natural_movie_three':    'orange',
          'spontaneous':       'green',
          'natural_movie_one': 'red'}

for c,stim_name in enumerate(stim_epoch.stimulus.unique()):
    stim = stim_epoch[stim_epoch.stimulus==stim_name]
    for j in range(len(stim)):
        plt.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[stim_name], alpha=0.1)


drifting_gratings_table = SelectedNWB.get_stimulus_table('drifting_gratings')


#%%
dxcm, tsd = SelectedNWB.get_running_speed()
fig = plt.figure(figsize=(10,10))
for i in range(50):
    plt.plot(SelectedEvents[i,:]+(i*2), color='gray')
plt.plot((0.2*dxcm)-20)
    
#for each stimulus, shade the plot when the stimulus is presented
colors = ['blue','orange','green','red']
for c,stim_name in enumerate(stim_epoch.stimulus.unique()):
    stim = stim_epoch[stim_epoch.stimulus==stim_name]
    for j in range(len(stim)):
        plt.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[c], alpha=0.1)


#%%
fig = plt.figure(figsize=(10,8))
for i in range(50):
    plt.plot(    plt.plot(2.5*SelectedEvents[i,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]]
+(i*2), color='gray')+(i*2), color='gray')
    
orient = drifting_gratings_table.orientation.iloc[0]               
#shade traces with the time of each presentation of the above scene
stim_subset = drifting_gratings_table[drifting_gratings_table.orientation==orient]
for j in range(len(stim_subset)):
    plt.axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)


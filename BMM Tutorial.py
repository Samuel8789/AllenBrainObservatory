# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 10:39:18 2020

@author: sp3660
"""
from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PathStructure import PathStructure





ProjectCodePath, ProjectName, ProjectDataPath, ProjectRAWPath, ProjectTempRAWPath =PathStructure()
#%%
# boc = BrainObservatoryCache(manifest_file=ProjectDataPath+'/boc/manifest.json')
boc = BrainObservatoryCache(manifest_file='boc/manifest.json')

#%%
boc.get_all_targeted_structures()
boc.get_all_imaging_depths()
boc.get_all_cre_lines()
boc.get_all_reporter_lines()
boc.get_all_stimuli()
#%%
visual_area = 'VISp'
cre_line ='Cux2-CreERT2'
depth=275
exps = boc.get_experiment_containers(
    targeted_structures=[visual_area], 
    cre_lines=[cre_line])
exps = pd.DataFrame(exps)

#%%
experiment_container_id = 511510736
exp_cont = boc.get_ophys_experiments(
    experiment_container_ids=[experiment_container_id],
)
pd.DataFrame(exp_cont)
#%%
exp_cont_ns = boc.get_ophys_experiments(
    experiment_container_ids=[experiment_container_id], 
    stimuli=['locally_sparse_noise'],
)
pd.DataFrame(exp_cont_ns)

session_id = exp_cont_ns[0]['id']
print(session_id)
#%%
data_set = boc.get_ophys_experiment_data(ophys_experiment_id=session_id)
#%%
max_projection = data_set.get_max_projection()
fig = plt.figure(figsize=(6,6))
plt.imshow(max_projection, cmap='gray')
rois = data_set.get_roi_mask_array()
print(rois.shape)
print("Number of cells:", rois.shape[0])
plt.figure(figsize=(6,6))
plt.imshow(rois.sum(axis=0))
#%%
ts, dff = data_set.get_dff_traces()
dff.shape
fig = plt.figure(figsize=(10,8))
for i in range(50):
    plt.plot(dff[i]+(i*2), color='gray')
#%%
events = boc.get_ophys_experiment_events(session_id)
events.shape
fig = plt.figure(figsize=(10,8))
for i in range(50):
    plt.plot(2.5*events[i]+(i*2), color='gray')
#%%
stim_epoch = data_set.get_stimulus_epoch_table()
stim_epoch
stim_epoch.stimulus.unique()

fig = plt.figure(figsize=(10,8))
for i in range(50):
    plt.plot(dff[i,:]+(i*2), color='gray')
    
#for each stimulus, shade the plot when the stimulus is presented
colors = {'static_gratings':   'blue',
          'natural_scenes':    'orange',
          'spontaneous':       'green',
          'natural_movie_one': 'red'}

for c,stim_name in enumerate(stim_epoch.stimulus.unique()):
    stim = stim_epoch[stim_epoch.stimulus==stim_name]
    for j in range(len(stim)):
        plt.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[stim_name], alpha=0.1)
#%%
natural_scene_table = data_set.get_stimulus_table('natural_scenes')
natural_scene_table.head()
natural_scene_template = data_set.get_stimulus_template('natural_scenes')
natural_scene_template.shape
scene_number = natural_scene_table.frame.iloc[3]
plt.imshow(natural_scene_template[scene_number,:,:], cmap='gray')

fig = plt.figure(figsize=(10,8))
for i in range(50):
    plt.plot(dff[i,:]+(i*2), color='gray')
    
#shade traces with the time of each presentation of the above scene
stim_subset = natural_scene_table[natural_scene_table.frame==scene_number]
for j in range(len(stim_subset)):
    plt.axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)
#%%

dxcm, tsd = data_set.get_running_speed()
fig = plt.figure(figsize=(10,3))
plt.plot(dxcm)
plt.ylabel("Running speed (cm/s)")
plt.xlabel("acquisition frame")
#%%
fig = plt.figure(figsize=(10,10))
for i in range(50):
    plt.plot(dff[i,:]+(i*2), color='gray')
plt.plot((0.2*dxcm)-20)
    
#for each stimulus, shade the plot when the stimulus is presented
colors = ['blue','orange','green','red']
for c,stim_name in enumerate(stim_epoch.stimulus.unique()):
    stim = stim_epoch[stim_epoch.stimulus==stim_name]
    for j in range(len(stim)):
        plt.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[c], alpha=0.1)

#%%
cell_specimens = pd.DataFrame(boc.get_cell_specimens())
cell_specimens.shape
cell_specimens.head()
cell_specimens.keys()
subset = cell_specimens[cell_specimens.experiment_container_id==experiment_container_id]
len(subset)
#%%

plt.imshow(natural_scene_template[22,:,:], cmap='gray')
subset[(subset.p_ns<0.05)&(subset.pref_image_ns==22)].cell_specimen_id


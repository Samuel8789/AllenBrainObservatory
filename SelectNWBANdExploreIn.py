 # -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 09:43:22 2020

@author: sp3660
"""
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
SelectedCreLine='Vip-IRES-Cre'
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



for i in range(0,len(SelectedExpIDs2)):
    Select=i
    SelectExpID=SelectedExpIDs2[Select]
    if SelectExpID==659495103:
        continue
    
    SelectNWB=NWBexpFiles.get('NWB'+str(SelectExpID))
    SelectEvents=EventsFiles["Events{0}".format(SelectExpID)]
    SelectConatiner=SelectNWB.get_metadata()['experiment_container_id']
   

    SelectDataDir=ProjectDataPath+'\\'+SelectedCreLine+'\\'+SelectedStructure+'\\'+str(SelectedDepth)+'\\'+SelectedSession+'\\ContID_'+str(SelectConatiner)+'_ExpID_'+str(SelectedExpIDs2[Select])
    
    
    
    SelectCells=SelectNWB.get_cell_specimen_ids()
    _, raw_traces = SelectNWB.get_fluorescence_traces(SelectCells)
    _, demixed_traces = SelectNWB.get_demixed_traces(SelectCells)
    _, neuropil_traces = SelectNWB.get_neuropil_traces(SelectCells)
    _, corrected_traces = SelectNWB.get_corrected_fluorescence_traces(SelectCells)
    timest, dff_traces = SelectNWB.get_dff_traces(SelectCells)
    
    
    CellIDs=SelectNWB.get_cell_specimen_ids()
    CellIDxs=SelectNWB.get_cell_specimen_indices(CellIDs)
    TimeStmps=SelectNWB.get_fluorescence_timestamps()
    all_roi_masks = SelectNWB.get_roi_mask_array()
    max_projection = SelectNWB.get_max_projection()
    MtData=SelectNWB.get_metadata()
    MtData.pop('session_start_time')
    MotCorr=SelectNWB.get_motion_correction()
    MotCorrNp=MotCorr.to_numpy()
    NeuRP=SelectNWB.get_neuropil_r()
    # while True:
    #     try:
    #         PupTime,Pupil=SelectNWB.get_pupil_location()
    #         break
    #     except NoEyeTrackingException:
    #         pass
    #     try:
    #         PupSIz=SelectNWB.get_pupil_size()
    #         break

        
    ROIIDs=SelectNWB.get_roi_ids()
    ROoiMASks=SelectNWB.get_roi_mask()
    Speed=SelectNWB.get_running_speed()
    #Stim=SelectNWB.get_stimulus(50000)
    StimTab=SelectNWB.get_stimulus_epoch_table()
    StimTabNp=StimTab.to_numpy()
    StimList=SelectNWB.list_stimuli()
    Grats=SelectNWB.get_stimulus_table('drifting_gratings')
    GratsNp=Grats.to_numpy()
    NatMov1Fr=SelectNWB.get_stimulus_table('natural_movie_one')
    NatMov1FrNp=NatMov1Fr.to_numpy()
    NatMov2Fr=SelectNWB.get_stimulus_table('natural_movie_three')
    NatMov2FrNp=NatMov2Fr.to_numpy()
    NatMov1=SelectNWB.get_stimulus_template('natural_movie_one')
    NatMov2=SelectNWB.get_stimulus_template('natural_movie_three')
    
    
    
    
    scipy.io.savemat(SelectDataDir+'\\'+str(SelectExpID)+'.mat', dict(
        TimeStmps=TimeStmps, 
        dff_traces=dff_traces, 
        Datapath=SelectedDataDir, 
        Codepath=ProjectCodePath, 
        Rawpath=ProjectRAWPath,
        TempRawpath=ProjectTempRAWPath,
        max_projection=max_projection,
        MtData=MtData,
        all_roi_masks=all_roi_masks,
        Speed=Speed,
        StimTab=StimTabNp,
        Grats=GratsNp,
        NatMov1Fr=NatMov1FrNp,
        NatMov2Fr=NatMov2FrNp,
        NatMov1=NatMov1,
        NatMov2=NatMov2,
        Events=SelectEvents
        ))

    
      

#%%
# from matplotlib import pyplot as plt


# # plot raw and corrected ROI trace
# plt.figure(figsize=(14,4))
# plt.title("Raw Fluorescence Trace")
# plt.plot(time, raw_traces[0])
# plt.show()

# plt.figure(figsize=(14,4))
# plt.title("Demixed Fluorescence Trace")
# plt.plot(time, demixed_traces[0])
# plt.show()

# plt.figure(figsize=(14,4))
# plt.title("Neuropil-corrected Fluorescence Trace")
# plt.plot(time, corrected_traces[0])
# plt.show()

# plt.figure(figsize=(14,4))
# plt.title("dF/F Trace")
# # warning: dF/F can occasionally be one element longer or shorter 
# # than the time stamps for the original traces.
# plt.plot(time[:len(dff_traces[0])], dff_traces[0])
# plt.show()

# #%%
# import matplotlib.pyplot as plt
# import numpy as np
# %matplotlib inline


# # get masks for specific cells
# roi_mask_list = SelectNWB.get_roi_mask(cell_specimen_ids=SelectCells)

# # plot each mask
# f, axes = plt.subplots(1,2)

# # make a mask of all ROIs in the experiment    
# combined_mask = all_roi_masks.max(axis=0)

# axes[-1].imshow(combined_mask, cmap='gray')
# axes[-1].set_title('all ROIs')

# # show the movie max projection
# axes[-0].imshow(max_projection, cmap='gray')
# axes[-0].set_title('max projection')

# plt.show()

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 10:46:24 2020

@author: sp3660
"""
# This is going to be the code to acces the NWB files for a given line and also acces the id to downlad the raw data
# First import
from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
import os
from collections import Counter
from PathStructure import PathStructure
ProjectCodePath, ProjectName, ProjectDataPath, ProjectRAWPath, ProjectTempRAWPath =PathStructure()
boc = BrainObservatoryCache(manifest_file=ProjectDataPath+'/boc/manifest.json')
#os.makedirs(ProjectRAWPath+'/RAWMovies')

#%% Now select visual area and check what kind of experiments were done, if theyr are from same mice, etc

cre_lines = boc.get_all_cre_lines()
print("all cre lines:\n")
pprint.pprint(cre_lines)
#%%
SelectedCreLine='Emx1-IRES-Cre'
#%%
targeted_structures = boc.get_all_targeted_structures()
print("all targeted structures: " + str(targeted_structures))
#%%
SelectedStructure='VISp'
#%% First Chek Depths available
AllDepths=boc.get_experiment_containers(   
                                  targeted_structures=[SelectedStructure],
                                  cre_lines=[SelectedCreLine],)

AvailableDepths = sorted(dict(Counter([sub['imaging_depth']  for sub in AllDepths ])).items(), key=lambda x: x[1], reverse=True)
print("AvailableDepths=:\n")
for i in AvailableDepths:
	print(i[0], i[1])
    
#%%
SelectedDepth=275
#%%
ImagingPlane = boc.get_experiment_containers(   
                                  targeted_structures=[SelectedStructure],
                                  imaging_depths=[SelectedDepth],
                                  cre_lines=[SelectedCreLine],)

print(SelectedStructure+"," +str(SelectedDepth)+","+ SelectedCreLine+" Mice: %d\n" % len(ImagingPlane))
print("Example experiment container record:")
pprint.pprint(ImagingPlane)
#%%
# So I have multiple  containers from several different mice, now let's check one conatiner. Now I have to select one of the mice

AllSessions=boc.get_all_session_types()
print("all sessions: " + str(AllSessions))
#%%
SelectedSession='three_session_A'
#%%

SelectedExperiments=boc.get_ophys_experiments(targeted_structures=[SelectedStructure], imaging_depths=[SelectedDepth], cre_lines=[SelectedCreLine], session_types=[SelectedSession])
print("Experiments Selected %s: %d\n" %(SelectedSession, len(SelectedExperiments)))
pprint.pprint(SelectedExperiments)




#%%
SelectedContainers = [sub['experiment_container_id']  for sub in SelectedExperiments] 

SelectedExpIDs = [sub['id']  for sub in SelectedExperiments] 
#SelectedExpIDs.sort()
EventsFiles={}
NWBexpFiles={}
for i in range(0, len(SelectedExpIDs)):
    SelectedDataDir=ProjectDataPath+'\\'+SelectedCreLine+'\\'+SelectedStructure+'\\'+str(SelectedDepth)+'\\'+SelectedSession+'\\ContID_'+str(SelectedContainers[i])+'_ExpID_'+str(SelectedExpIDs[i])
    #os.makedirs(SelectedDataDir)
    NWBexpFiles["NWB{0}".format(SelectedExpIDs[i])] = boc.get_ophys_experiment_data(SelectedExpIDs[i],  SelectedDataDir+'\\'+str(SelectedExpIDs[i]))
    EventsFiles["Events{0}".format(SelectedExpIDs[i])] = boc.get_ophys_experiment_events(SelectedExpIDs[i],  SelectedDataDir+'\\'+str(SelectedExpIDs[i])+'events')


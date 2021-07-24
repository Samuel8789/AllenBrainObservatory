# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 11:11:51 2020

@author: Samuel
"""

from allensdk.core.brain_observatory_cache import BrainObservatoryCache
from pathlib import Path


#%% Create Folder Structure
data_folder = Path(r'C:\Users\Samuel\Dropbox\Science\Lab\CodeData\BrainObservatoryData')
project='ExcitatoryUpperV1'
boc = BrainObservatoryCache(manifest_file=data_folder/'boc/manifest.json')

#%% Get Observatory Info

CRE=boc.get_all_cre_lines() #global information about the databse to have it in lists
AREAS=boc.get_all_targeted_structures() #global information about the databse to have it in lists
DEPTHS=boc.get_all_imaging_depths()#global information about the databse to have it in lists
REPORTER=boc.get_all_reporter_lines()#global information about the databse to have it in lists GCAMP used
SESSIONS=boc.get_all_session_types()#global information about the databse to have it in lists
STIMULI=boc.get_all_stimuli()#global information about the databse to have it in lists

#%% Get Container and experiments

visPExc275=boc.get_experiment_containers(imaging_depths=[275],targeted_structures=['VISp'],cre_lines=['Slc17a7-IRES2-Cre'] )
contIDs=[cont['id']for cont in visPExc275]
cont1=boc.get_ophys_experiments(experiment_container_ids=[contIDs[0]])
expIDs=[ids['id'] for ids in cont1]
#%%
exp_sesA = boc.get_ophys_experiment_data(expIDs[1])
exp_sesB = boc.get_ophys_experiment_data(expIDs[0])
exp_sesC2 = boc.get_ophys_experiment_data(expIDs[2])
#%% Get all Traces in matrixes


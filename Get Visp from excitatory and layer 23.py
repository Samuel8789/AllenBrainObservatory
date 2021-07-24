# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 17:07:49 2019

@author: sp3660
"""

from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
import os
#%%
boc = BrainObservatoryCache(manifest_file='boc/manifest.json')
#%%
emx_ecs = boc.get_experiment_containers(targeted_structures=['VISp'], imaging_depths=[275], cre_lines=['Emx1-IRES-Cre'] )
slc_ecs = boc.get_experiment_containers(targeted_structures=['VISp'], imaging_depths=[275], cre_lines=['Slc17a7-IRES2-Cre'])
#%%
# =============================================================================
# emx_ecs_sdonor=sorted(emx_ecs, key=lambda i: i['donor_name'])
# slc_ecs_sdonor=sorted(slc_ecs, key=lambda i: i['donor_name'])
# 
# emx_ecs_donors=[ec['donor_name'] for ec in emx_ecs_sdonor]
# slc_ecs_donors=[ec['donor_name'] for ec in slc_ecs_sdonor]
# #%%
# from collections import Counter
# 
# 
# Counter(emx_ecs_donors).keys()
# Counter(emx_ecs_donors).values()
# 
# emx_multip_donor=dictionary = dict(zip(Counter(emx_ecs_donors).keys(), Counter(emx_ecs_donors).values()))
# 
# Counter(slc_ecs_donors).keys()
# Counter(slc_ecs_donors).values()
# 
# slc_multip_donor=dictionary = dict(zip(Counter(slc_ecs_donors).keys(), Counter(slc_ecs_donors).values()))
# =============================================================================
#%%

emx_ecs_ids=[ec['id'] for ec in emx_ecs]
slc_ecs_ids=[ec['id'] for ec in slc_ecs]


#emx_exps = [[boc.get_ophys_experiments(experiment_container_ids=[x])] for x in emx_ecs_ids]
import allensdk.brain_observatory.stimulus_info as stim_info
emx_exp=boc.get_ophys_experiments(experiment_container_ids=[emx_ecs_ids[0]], stimuli=[stim_info.DRIFTING_GRATINGS])[0]

emx_exp_data=boc.get_ophys_experiment_data(emx_exp['id'])

print("Metadata from NWB file:")
pprint.pprint(emx_exp_data.get_metadata())
#%% Get all cell in the experiment
import matplotlib.pyplot as plt
import numpy as np
#%%
cell_ids = emx_exp_data.get_cell_specimen_ids()
num_cells=len(cell_ids)
cell_idx = emx_exp_data.get_cell_specimen_indices(cell_ids)
#%%
max_p=emx_exp_data.get_max_projection()

#%%
# plot each mask
f, axes = plt.subplots(1,2, figsize=(15, 5))
all_roi_masks = emx_exp_data.get_roi_mask_array()
combined_mask = all_roi_masks.max(axis=0)

axes[0].imshow(combined_mask, cmap='gray')
axes[0].set_title('all ROIs')

axes[1].imshow(max_p, cmap='gray')
axes[1].set_title('max projection')
plt.show()
#%%

time, raw_traces = emx_exp_data.get_fluorescence_traces()
_, demixed_traces = emx_exp_data.get_demixed_traces()
_, neuropil_traces = emx_exp_data.get_neuropil_traces()
_, corrected_traces = emx_exp_data.get_corrected_fluorescence_traces()
_, dff_traces = emx_exp_data.get_dff_traces()

#%%
dxcm, dxtime = emx_exp_data.get_running_speed()
plt.figure(figsize=(14,4))
plt.plot(dxtime, dxcm)
plt.show()

#%%
from allensdk.brain_observatory.brain_observatory_exceptions import NoEyeTrackingException
#%%
try:
    timestamps, locations = emx_exp_data.get_pupil_location()
except NoEyeTrackingException:
    print("No eye tracking for experiment %s." % emx_exp_data.get_metadata()["ophys_experiment_id"])
#%%
    timestamps, locations = emx_exp_data.get_pupil_location()
plt.figure(figsize=(14,4))
plt.plot(timestamps, locations.T[0])
plt.plot(timestamps, locations.T[1])
plt.title("Eye position over time")
plt.xlabel("time (s)")
plt.ylabel("angle (deg)")
plt.legend(['azimuth', 'altitude'])
plt.show()

#pupil size over time
timestamps, area = emx_exp_data.get_pupil_size()
plt.figure(figsize=(14,4))
plt.plot(timestamps, area)
plt.title("Pupil size over time")
plt.xlabel("time (s)")
plt.ylabel("area (px)")
plt.ylim(0, 20000)
plt.show()

# scatter of gaze positions over approximate screen area
plt.figure()
plt.scatter(locations.T[0], locations.T[1], s=2, c="m", edgecolor="")
plt.title("Eye position scatter plot")
plt.xlim(-70, 70)
plt.ylim(-60, 60)
plt.xlabel("azimuth (deg)")
plt.ylabel("altitude (deg)")
plt.show()

#%%

import numpy, scipy.io
scipy.io.savemat('boc/dff.mat', mdict={'dff': dff_traces})
scipy.io.savemat('boc/time.mat', mdict={'time': time})
#%%
# plot raw and corrected ROI trace
plt.figure(figsize=(14,4))
plt.title("Raw Fluorescence Trace")
plt.plot(time, raw_traces[0])
plt.show()

plt.figure(figsize=(14,4))
plt.title("Demixed Fluorescence Trace")
plt.plot(time, demixed_traces[0])
plt.show()

plt.figure(figsize=(14,4))
plt.title("Neuropil-corrected Fluorescence Trace")
plt.plot(time, corrected_traces[0])
plt.show()

plt.figure(figsize=(14,4))
plt.title("dF/F Trace")
# warning: dF/F can occasionally be one element longer or shorter 
# than the time stamps for the original traces.
plt.plot(time[:len(dff_traces[0])], dff_traces[0])
plt.show()
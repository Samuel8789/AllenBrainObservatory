#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 12:57:29 2022

@author: samuel
"""

from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import allensdk.brain_observatory.stimulus_info as si
import scipy.ndimage.interpolation as spndi
from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
import allensdk.brain_observatory.stimulus_info as si
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.interpolate import griddata




#%%
# boc = BrainObservatoryCache(manifest_file=ProjectDataPath+'/boc/manifest.json')
boc = BrainObservatoryCache(manifest_file='boc/manifest.json')

#%%
boc.get_all_stimuli()
#%%
visual_area = 'VISp'
cre_line ='Cux2-CreERT2'
depth=275
= boc.get_experiment_containers(
    targeted_structures=[visual_area], 
    cre_lines=[cre_line])
exps = pd.DataFrame(exp_containers)

#%%
first_experiment_container_id = exp_containers[0]['id']
exps= boc.get_ophys_experiments(
    experiment_container_ids=[first_experiment_container_id],
)
pd.DataFrame(exps)
#%%
stim='locally_sparse_noise'


for stim in boc.get_all_stimuli():
    # get experiment id with stim

    selected_exp= boc.get_ophys_experiments(
        experiment_container_ids=[first_experiment_container_id], 
        stimuli=[stim])

session_id = selected_exp[0]['id']
print(session_id)



VisID=session_id
data_set = boc.get_ophys_experiment_data(VisID)
print("Metadata from NWB file:")
pprint.pprint(data_set.get_metadata())
full_stim_epoch_table=data_set.get_stimulus_epoch_table()
stim_by_frame_table=data_set.get_stimulus_table(stim)


z=data_set.list_stimuli()
y=data_set.get_stimulus_template(stim)
# yy=data_set.get_stimulus_template(stim)[stim_by_frame_table['frame'].values,:,:]


scene_nums = [10, 20]
fig, axes = plt.subplots(1,len(scene_nums))
for ax,scene in zip(axes, scene_nums):
    ax.imshow(y[scene,:,:], cmap='gray')
    ax.set_axis_off()
    ax.set_title('scene %d' % scene)

selected_stim_epoch_table=full_stim_epoch_table.loc[full_stim_epoch_table['stimulus'] == stim]

params, template = data_set.get_stimulus(selected_stim_epoch_table['start'].iloc[0]+15)
pd.Series(params[2])
img = plt.imshow(template, cmap=plt.cm.gray, interpolation='none')
    
#%%]

img = plt.imshow(y[1,:,:], cmap=plt.cm.gray, interpolation='none')
# sio.savemat(stim+'mat',{'y':y})


test= si.all_stimuli()

#%% gratings
test= si.get_spatial_grating(1024, 1.25,90,2,0)
plt.imshow(test)
test2= si.get_spatio_temporal_grating(5,temporal_frequency=2,height=1024, aspect_ratio=1.25, 
                                      ori=90, pix_per_cycle=2, phase=0, 
                                      p2p_amp=2, baseline=0)
plt.imshow(test2)
#%%
griddata
# test3=si.lsn_coordinate_to_monitor_coordinate((2,2),(1024,1280),)

test4=si.make_display_mask((1024,1280))
test5=si.mask_stimulus_template((1024,1280))



params, template = data_set.get_stimulus(selected_stim_epoch_table['start'].iloc[0]+500)
plt.imshow(template)


mym=si.Monitor(1024,1280, 48.2,'cm')
expgeom=si.ExperimentGeometry(15, mym.height, mym.width, (1024,1280), (0.5,0.5))
m=si.BrainObservatoryMonitor()


img_screen=m.lsn_image_to_screen(template)

m.show_image(img_screen, warp=True, origin='upper')
m.show_image(img_screen, warp=False, mask=False, origin='upper')




import scipy.ndimage



test= scipy.ndimage.zoom(template, 60, order=0)





test=spndi.map_coordinates(template, expgeom.warp_coordinates.T).reshape((mym.n_pixels_r,mym.n_pixels_c))
img=mym.show_image(test, mask=True)




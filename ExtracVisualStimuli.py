# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:13:05 2020

@author: sp3660
"""


from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
import allensdk.brain_observatory.stimulus_info as si
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

#%%
boc = BrainObservatoryCache(manifest_file='boc/manifest.json')
#%%
stims = boc.get_all_stimuli()
print("all stimuli:\n")
pprint.pprint(stims)

#%%

VisID=501474098
data_set = boc.get_ophys_experiment_data(VisID)
print("Metadata from NWB file:")
pprint.pprint(data_set.get_metadata())
tab=data_set.get_stimulus_epoch_table()
test=data_set.get_stimulus_table('natural_movie_two')


z=data_set.list_stimuli()
y=data_set.get_stimulus_template('natural_movie_two')
yy=data_set.get_stimulus_template('natural_movie_two')[test['frame'].values,:,:]


scene_nums = [10, 20]
fig, axes = plt.subplots(1,len(scene_nums))
for ax,scene in zip(axes, scene_nums):
    ax.imshow(y[scene,:,:], cmap='gray')
    ax.set_axis_off()
    ax.set_title('scene %d' % scene)


params, template = data_set.get_stimulus(64000)
pd.Series(params[2])
img = plt.imshow(template, cmap=plt.cm.gray, interpolation='none')



#%%

img = plt.imshow(y[1,:,:], cmap=plt.cm.gray, interpolation='none')
sio.savemat('natural_movie_two.mat',{'y':y})




m=si.BrainObservatoryMonitor()
img_screen=m.natural_movie_image_to_screen(template)
m.show_image(img_screen, warp=True)



    
                                
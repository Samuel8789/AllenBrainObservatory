# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 15:49:19 2020

@author: sp3660
"""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


#%%
from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint

# This class uses a 'manifest' to keep track of downloaded data and metadata.  
# All downloaded files will be stored relative to the directory holding the manifest
# file.  If 'manifest_file' is a relative path (as it is below), it will be 
# saved relative to your working directory.  It can also be an absolute path.
boc = BrainObservatoryCache(manifest_file='boc/manifest.json')
#%%
#%%
# Download a list of all targeted areas
targeted_structures = boc.get_all_targeted_structures()
print("all targeted structures: " + str(targeted_structures))

# Download experiment containers for VISp experiments
visp_ecs = boc.get_experiment_containers(targeted_structures=['VISp'])
print("all VISp experiment containers: %d" % len(visp_ecs))

# Download a list of all imaging depths
depths = boc.get_all_imaging_depths()
print("all imaging depths: " + str(depths))

# Download a list of all stimuli
ses = boc.get_all_session_types()
print("all stimuli:\n")
pprint.pprint(ses)

# Download a list of all cre driver lines 
cre_lines = boc.get_all_cre_lines()
print("all cre lines:\n")
pprint.pprint(cre_lines)

#%%

# Download experiment containers for Cux2 experiments
in_ecs = boc.get_experiment_containers(targeted_structures=['VISp'], imaging_depths=depths,cre_lines=['Pvalb-IRES-Cre', 'Sst-IRES-Cre', 'Vip-IRES-Cre'])
pv_ecs = boc.get_experiment_containers(cre_lines=['Pvalb-IRES-Cre'])

print("PV experiments: %d\n" % len(pv_ecs))

print("Example experiment container record:")
pprint.pprint(pv_ecs[0])

#%%

#all=boc.get_ophys_experiments(targeted_structures=['VISp'], imaging_depths=depths,cre_lines=['Pvalb-IRES-Cre', 'Sst-IRES-Cre', 'Vip-IRES-Cre'], session_types=None)

#%%
all=boc.get_ophys_experiments(targeted_structures=['VISp'], imaging_depths=[275], cre_lines=['Sst-IRES-Cre'], session_types=['three_session_A'])

expids=[]    
for exp in all:
    expids.append(exp.get('id'))
    
 #%%   
exp = boc.get_ophys_experiment_data(expids[3])
# print out the metadata available in the NWB file
pprint.pprint(exp.get_metadata())
    #%%
maxp=exp.get_max_projection()
imgplot = plt.imshow(maxp)


time, raw_traces = exp.get_fluorescence_traces()
_, demixed_traces = exp.get_demixed_traces()
_, neuropil_traces = exp.get_neuropil_traces()
_, corrected_traces = exp.get_corrected_fluorescence_traces()
_, dff_traces = exp.get_dff_traces()

#%%
from matplotlib import pyplot as plt


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




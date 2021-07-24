# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 16:32:44 2020

@author: sp3660
"""

#Basic Exploration of NWB files
from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
import os



ProjectCodePath=os.getcwd()
ProjectName=os.path.basename(ProjectCodePath)
ProjectDataPath='G:\\Dropbox\\CodeData\\'+ ProjectName
ProjectRAWPath='F:\\CodeRawData\\'+ ProjectName
boc = BrainObservatoryCache(manifest_file=ProjectDataPath+'/boc/manifest.json')
#%%

from pynwb import NWBFile, TimeSeries, NWBHDF5IO
from datetime import datetime
from dateutil.tz import tzlocal
import numpy as np
import numpy as np
import h5py



#%%
hf = h5py.File(ProjectRAWPath+'\RAWMovies\ophys_experiment_561994407.h5', 'r')
hf.keys()
n1 = hf.get('data')
n1
n1 = np.array(n1)
n1.shape
#%%
import imageio
imageio.mimwrite('output_filename.mp4', n1 , fps = [60])

n2=n1
n2[n2>255]=255
n2[n2<0]=0
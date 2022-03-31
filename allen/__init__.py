# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 12:56:27 2021

@author: sp3660
"""

import glob
from  .allenBrainObservatory import AllenBrainObservatory
import logging
import logging.config


log_dir='\\\\?\\'+r'K:\Projects\LabNY\Full_Mice_Pre_Processed_Data\Logging\AllenModule'

log_files=glob.glob(log_dir+'\\app_errors_allen**.log')
if log_files:
    new_file_number=len(log_files)+1
else:
    new_file_number=1

filename='\\\\?\\'+r'K:\Projects\LabNY\Full_Mice_Pre_Processed_Data\Logging\AllenModule\app_allen'+str(new_file_number)+'.log'
filename_errors='\\\\?\\'+r'K:\Projects\LabNY\Full_Mice_Pre_Processed_Data\Logging\AllenModule\app_errors_allen'+str(new_file_number)+'.log'







logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG) # or whatever


handler = logging.FileHandler(filename,'w', 'utf-8') # or whatever
handler_errors = logging.FileHandler(filename_errors,'w', 'utf-8') # or whatever

handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')) # or whatever
handler_errors.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')) # or whatever


handler.setLevel(logging.DEBUG)

handler_errors.setLevel(logging.ERROR)


logger.addHandler(handler)
logger.addHandler(handler_errors)


logger.info('Starting App')

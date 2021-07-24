# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 09:16:01 2021

@author: sp3660
"""

# Allen main scrfrom
import sys

sys.path.insert(0, r'C:/Users/sp3660/Documents/Github/ProjectManager')
from AllenBrainObservatory import AllenBrainObservatory


#%%
Allen=AllenBrainObservatory()
#%%
Allen.Get_all_mouse_from_trangenic_line('Pvalb-IRES-Cre')

# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 09:13:02 2021

@author: sp3660
"""

import sys

sys.path.insert(0, r'C:/Users/sp3660/Documents/Github/ProjectManager')


from ProjectsCLass import Project


class AllenBrainObservatory(Project):    
    def __init__(self):
        Project.__init__(self, 'AllenBrainObservatory')
        
        self.main_directory='\\\?\\' +self.project_paths['Documents']  

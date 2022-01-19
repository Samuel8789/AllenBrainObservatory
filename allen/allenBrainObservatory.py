# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 09:13:02 2021

@author: sp3660
"""

import sys
import os 

from project_manager.ProjectsCLass import Project
#
import logging 
module_logger = logging.getLogger(__name__)
# sys.path.insert(0, r'C:\Users\sp3660\Documents\Github\ProjectManager\project_manager')
from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
from collections import Counter



class AllenBrainObservatory(Project):    
    def __init__(self,  githubtoken_path, gui=False ):
       Project.__init__(self, 'AllenBrainObservatory',githubtoken_path, Project.computer, Project.platform )
               
       self.main_directory=self.project_paths['G:']  
       self.allen_manifest =BrainObservatoryCache(manifest_file=os.path.join(self.main_directory, 'boc\manifest.json'))
       # self.all_cre_lines = self.allen_manifest.get_all_cre_lines()


    
    def Get_all_mouse_from_trangenic_line(self, transgenic_line):
        all_line_containers=self.allen_manifest.get_experiment_containers(cre_lines=[transgenic_line],)   
        all_line_mouse = sorted(dict(Counter([cont['donor_name']  for cont in all_line_containers ])).items(), key=lambda x: x[1], reverse=True)
        
        if transgenic_line in [self.all_cre_lines[5],self.all_cre_lines[10],self.all_cre_lines[12]]:
            
            if 'Pvalb' in transgenic_line:
                line='PV'
            elif 'Sst' in transgenic_line:
                line='SST'    
            elif 'Vip' in transgenic_line:
                line='VIP'
                                
            intern_line_path=os.path.join(self.main_directory,'Interneurons',line)
            if not os.path.isdir(intern_line_path):
                 os.mkdir(intern_line_path)
            for mouse in  all_line_mouse:  
                 if not os.path.isdir(os.path.join(intern_line_path,mouse[0])):
                     ms_path='Mouse_'+mouse[0]
                     os.mkdir(os.path.join(intern_line_path ,ms_path ))
                     for cont in all_line_containers: 
                         if cont['donor_name']==mouse[0]:
                             nwb_path=os.path.join(intern_line_path ,ms_path)+'\\Container_'+str(cont['id'])
                             SelectedExperiments=self.allen_manifest.get_ophys_experiments(file_name=nwb_path, ids=cont['id'])
                                                                                           
                                 
                                 

                             

          
        

"""
select cre line first
then select  visal area
then select depth
then get the mice

get all videos of that plane

get all data of that plane

    thre sesions

"""
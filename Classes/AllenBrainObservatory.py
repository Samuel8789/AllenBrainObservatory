# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 09:13:02 2021

@author: sp3660
"""

import sys

sys.path.insert(0, r'C:/Users/sp3660/Documents/Github/ProjectManager')
from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import pprint
import os
from collections import Counter
import allensdk.brain_observatory.stimulus_info as stim_info
from ProjectsCLass import Project


class AllenBrainObservatory(Project):    
    def __init__(self):
        Project.__init__(self, 'AllenBrainObservatory')
        
        self.main_directory=self.project_paths['D:']  
        self.allen_manifest = BrainObservatoryCache(manifest_file=os.path.join(self.main_directory, 'boc\manifest.json'))
        self.all_cre_lines = self.allen_manifest.get_all_cre_lines()


    def testing(self):
        # Download a list of all targeted areas
        targeted_structures = self.allen_manifest.get_all_targeted_structures()
        print("all targeted structures: " + str(targeted_structures))
        # %%
        # Download experiment containers for VISp experiments
        visp_ecs = self.allen_manifest.get_experiment_containers(targeted_structures=['VISp'])
        print("all VISp experiment containers: %d" % len(visp_ecs))
        # %%
        # Download a list of all imaging depths
        depths = self.allen_manifest.get_all_imaging_depths()
        print("all imaging depths: " + str(depths))
        # %%
        # Download a list of all stimuli
        stims = self.allen_manifest.get_all_stimuli()
        print("all stimuli:\n")
        pprint.pprint(stims)
        #%%
        # Download a list of all cre driver lines 
        cre_lines = self.allen_manifest.get_all_cre_lines()
        print("all cre lines:\n")
        pprint.pprint(cre_lines)
        #%%
        # Download experiment containers for Cux2 experiments
        cux2_ecs = self.allen_manifest.get_experiment_containers(cre_lines=['Cux2-CreERT2'])
        print("Cux2 experiments: %d\n" % len(cux2_ecs))
        
        print("Example experiment container record:")
        pprint.pprint(cux2_ecs[0])
        #%%
        # Find all of the experiments for an experiment container
        cux2_ec_id = cux2_ecs[0]['id']
        exps = self.allen_manifest.get_ophys_experiments(experiment_container_ids=[cux2_ec_id])
        print("Experiments for experiment_container_id %d: %d\n" % (cux2_ec_id, len(exps)))
        pprint.pprint(exps)
        
        
        cux2_ec_id = cux2_ecs[-1]['id']

        # Find the experiment with the static static gratings stimulus
        exp = self.allen_manifest.get_ophys_experiments(experiment_container_ids=[cux2_ec_id], 
                                        stimuli=[stim_info.STATIC_GRATINGS])[0]
        print("Experiment with static gratings:")
        pprint.pprint(exp)
        #%%
        exp = self.allen_manifest.get_ophys_experiment_data(cux2_ec_id)
        
        # print out the metadata available in the NWB file
        pprint.pprint(exp.get_metadata())
        #%%
        import pandas as pd
        #%%
        # Download cells for a set of experiments and convert to DataFrame
        cells = self.allen_manifest.get_cell_specimens()
        cells = pd.DataFrame.from_records(cells)
        print("total cells: %d" % len(cells))
        #%%
        # find direction selective cells in VISp
        visp_ec_ids = [ ec['id'] for ec in visp_ecs ]
        visp_cells = cells[cells['experiment_container_id'].isin(visp_ec_ids)]
        print("VISp cells: %d" % len(visp_cells))
        #%%
        # significant response to drifting gratings stimulus
        sig_cells = visp_cells[visp_cells['p_dg'] < 0.05]
        print("cells with sig. response to drifting gratings: %d" % len(sig_cells))
        #%%
        # direction selective cells
        dsi_cells = sig_cells[(sig_cells['g_dsi_dg'] > 0.9)]
        print("direction-selective cells: %d" % len(dsi_cells))
        #%%
        import allensdk.brain_observatory.stimulus_info as stim_info
        #%%
        # find experiment containers for those cells
        dsi_ec_ids = dsi_cells['experiment_container_id'].unique()
        print("total dsi experiment containers: %d" % len(dsi_ec_ids))
        
        # Download the ophys experiments containing the drifting gratings stimulus for VISp experiment containers
        dsi_exps = self.allen_manifest.get_ophys_experiments(experiment_container_ids=dsi_ec_ids, stimuli=[stim_info.DRIFTING_GRATINGS])
        print("VISp drifting gratings ophys experiments: %d" % len(dsi_exps))
        
        print("Example ophys experiment:")
        pprint.pprint(dsi_exps[0])
        #%%
        # pick a direction-selective cell and find its NWB file
        dsi_cell = dsi_cells.iloc[0]
        
        # figure out which ophys experiment has the drifting gratings stimulus for that cell
        cell_exp = self.allen_manifest.get_ophys_experiments(cell_specimen_ids=[dsi_cell['cell_specimen_id']],
                                             stimuli=[stim_info.DRIFTING_GRATINGS])[0]
        
        data_set = self.allen_manifest.get_ophys_experiment_data(cell_exp['id'])
        
        print("Metadata from NWB file:")
        pprint.pprint(data_set.get_metadata())
        
        print("stimuli available in this file:")
        print(data_set.list_stimuli())

        
        
#%%

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
                ms_path='Mouse_'+mouse[0]
                if not os.path.isdir(os.path.join(intern_line_path,ms_path)):
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
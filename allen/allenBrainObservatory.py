# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 09:13:02 2021

@author: sp3660
"""
from project_manager.ProjectsCLass import Project

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pickle
import glob

import tkinter as Tkinter
import random
from tkinter import *
import pandas as pd
import tkinter as tk
import numpy as np

import logging 
module_logger = logging.getLogger(__name__)
from collections import Counter

from matplotlib.ticker import MaxNLocator
import matplotlib.patches as patches
import scipy.io as sio
from scipy.interpolate import griddata
import scipy.ndimage.interpolation as spndi
import scipy.ndimage
import scipy as scp
from pprint import pprint

from PIL import Image
import h5py
import caiman as cm
from scipy.ndimage.filters import median_filter

from allensdk.core.brain_observatory_cache import BrainObservatoryCache
from allensdk.brain_observatory.behavior.behavior_project_cache import VisualBehaviorOphysProjectCache
import allensdk.brain_observatory.stimulus_info as si
from allensdk.api.queries.brain_observatory_api import BrainObservatoryApi

from allensdk.brain_observatory.brain_observatory_exceptions import NoEyeTrackingException
from allensdk.brain_observatory.dff import calculate_dff
from allensdk.brain_observatory.r_neuropil import estimate_contamination_ratios
import allensdk.brain_observatory.dff as dff

from allensdk.brain_observatory.static_gratings import StaticGratings
from allensdk.brain_observatory.natural_scenes import NaturalScenes

from allensdk.brain_observatory.natural_movie import NaturalMovie
from allensdk.brain_observatory.locally_sparse_noise import LocallySparseNoise
from allensdk.brain_observatory.drifting_gratings import DriftingGratings
import allensdk.brain_observatory.receptive_field_analysis.visualization as rfvis
import allensdk.brain_observatory.receptive_field_analysis.receptive_field as rf
from allensdk.brain_observatory.locally_sparse_noise import LocallySparseNoise

from pylab import rcParams
rcParams['axes.xmargin'] = 0
rcParams['axes.ymargin'] = 0
import matplotlib.patches as patches



#%%
class AllenBrainObservatory(Project):    
    def __init__(self,  githubtoken_path, gui=False ):
        Project.__init__(self, 'AllenBrainObservatory',githubtoken_path, Project.computer, Project.platform )
        self.main_directory=self.project_paths['Documents']  
        self.allen_manifest =BrainObservatoryCache(manifest_file=os.path.join(self.main_directory, 'boc\manifest.json'))
        

        self.set_up_monitor()
        self.get_selection_options()
        

    def get_selection_options(self):
        
        pprint([self.allen_manifest.get_all_targeted_structures(), 
              self.allen_manifest.get_all_cre_lines(),
              self.allen_manifest.get_all_imaging_depths(),
              self.allen_manifest.get_all_reporter_lines(),
              self.allen_manifest.get_all_session_types(),
              self.allen_manifest.get_all_stimuli()]
              )
    

    def set_up_monitor(self, screen_size=[900,1440], lenght=48.26):
        
        self.screen_size=screen_size
        self.my_monitor=si.Monitor( screen_size[0], screen_size[1], lenght,'cm')
        self.expgeom=si.ExperimentGeometry(15, self.my_monitor.height, self.my_monitor.width, (screen_size[1], screen_size[0]), (0.5,0.5))
        self.obs_monitor=si.BrainObservatoryMonitor()
        
        
     
#%% dataset selection and download
    def select_lines_and_locations(self, area='VISp', line='Vglut1', depth=175, stim=None, sessions=['three_session_A','three_session_B','three_session_C','three_session_C2']):

        parameters_dict = {}

        parameters_dict["area"] = [area]
        parameters_dict["depth"]=[depth]

        if line=='Vglut1':
            parameters_dict["line"] =[self.allen_manifest.get_all_cre_lines()[9]]
            parameters_dict["reporters"]=['Ai93(TITL-GCaMP6f)']

        elif line=='Vip':
            parameters_dict["line"] =[self.allen_manifest.get_all_cre_lines()[-1]]
            parameters_dict["reporters"]=['Ai162(TIT2L-GC6s-ICL-tTA2)']

        elif line=='SST':
            parameters_dict["line"] =[self.allen_manifest.get_all_cre_lines()[10]]  
            parameters_dict["reporters"]=['Ai162(TIT2L-GC6s-ICL-tTA2)']

        elif line=='PV':
            parameters_dict["line"] =[self.allen_manifest.get_all_cre_lines()[5]] 
            parameters_dict["reporters"]=['Ai162(TIT2L-GC6s-ICL-tTA2)']
            
            
        exp_containers= self.allen_manifest.get_experiment_containers(
            targeted_structures=parameters_dict['area'], 
            cre_lines=parameters_dict['line'],
            imaging_depths=parameters_dict['depth'],
            # reporter_lines=self.reporters,
            )
    
        containers = pd.DataFrame(exp_containers)
        
        if stim:
            stim=[stim]
        
        exp_list=self.allen_manifest.get_ophys_experiments(
            targeted_structures=parameters_dict['area'], 
            cre_lines=parameters_dict['line'], 
            imaging_depths=parameters_dict['depth'],
            stimuli=stim,
            transgenic_lines=None,
            session_types=sessions,
            cell_specimen_ids=None
            )

        exps = pd.DataFrame(exp_list)

        return exp_containers, containers, exp_list, exps

    def select_exp_from_imaged_fov(self, exp_containers=None, container_selection_index=0):
        
        # select first cointiner by default
         # by container ID, this ere experiment form single fov from single mouse, however unknown vis stim
        if exp_containers:
            print('There is data')
            selected_container_id=exp_containers[container_selection_index]['id']
            selected_exp_list= self.allen_manifest.get_ophys_experiments(
                    experiment_container_ids=[selected_container_id])
            
            exps = pd.DataFrame(selected_exp_list)

            return selected_exp_list, exps
      
        else:
            print('Not data with this combinaiton')


    def select_experiment(self, selected_exp_list=None, exp_selection_index=0):
        
        if selected_exp_list:
            print('There are exps')
            selected_session_id= selected_exp_list[exp_selection_index]['id']
            selected_container_id= selected_exp_list[exp_selection_index]['experiment_container_id']

            return selected_session_id, selected_container_id
      
        else:
            print('Not data with this combinaiton')
            

    def download_single_imaging_session_nwb(self, selected_session_id=0):
        
        session_stimuli=self.allen_manifest.get_ophys_experiment_stimuli(selected_session_id)
        data_set = self.allen_manifest.get_ophys_experiment_data(selected_session_id)
        deconvolved_spikes=self.allen_manifest.get_ophys_experiment_events(selected_session_id)
        
        return data_set, deconvolved_spikes, session_stimuli
        
#%%       misc 
        
        
        
    def dealing_with_cells(self, selected_container_id):
        # too slow do not use
        all_container_fov_cells=self.allen_manifest.get_cell_specimens(file_name=None, 
                                               ids=None, 
                                               experiment_container_ids=[selected_container_id],
                                               include_failed=False,
                                               simple=True, 
                                               filters=None
                                               )

        return all_container_fov_cells
    #%% dataset loading
    def dataset_exploration(self, data_set):
        
        cell_info=self.get_cell_info(data_set)
        traces=self.get_all_traces(data_set)
        neuropilinfo=self.get_neuropil_info(data_set)
        locomotion_info=self.get_locomotion(data_set)
        metadata=self.dataset_metadata(data_set)
        pupilinfo=self.get_pupil_data(data_set)
        masks=self.mask_lists(data_set, cell_info[1])
        projection=self.max_projection(data_set)
        
        return cell_info, traces, neuropilinfo, locomotion_info, metadata, pupilinfo, masks, projection
 
        
    def get_cell_info(self, data_set):
        
        cell_number=data_set.number_of_cells
        cell_ids=data_set.get_cell_specimen_ids() # global cell ID
        cell_idxs=data_set.get_cell_specimen_indices(cell_ids) # this experiment cell indexes
        
        return cell_number, cell_ids, cell_idxs
       
    def get_all_traces(self, data_set):
        
        timestamps=data_set.get_fluorescence_timestamps()
        time, raw_traces = data_set.get_fluorescence_traces()
        _, demixed_traces = data_set.get_demixed_traces()
        _, neuropil_traces = data_set.get_neuropil_traces()
        _, corrected_traces = data_set.get_corrected_fluorescence_traces()
        _, dff_traces = data_set.get_dff_traces()
        
        return timestamps, time, raw_traces, demixed_traces, neuropil_traces, corrected_traces, dff_traces
        
    def get_neuropil_info(self, data_set):
        
        neuropil_r=data_set.get_neuropil_r()
        _, neuropil_traces=data_set.get_neuropil_traces()
        
        return neuropil_r, neuropil_traces
        
    def dataset_metadata(self, data_set):
        metadata=data_set.get_metadata()
        data_set.nwb_file
        data_set.pipeline_version
        motion_correct_shifts=data_set.get_motion_correction()
        
        return metadata, motion_correct_shifts
        
    def get_locomotion(self, data_set):
        dxcm, dxtime = data_set.get_running_speed()
        # plt.figure(figsize=(14,4))
        # plt.plot(dxtime, dxcm)
        # plt.show()
        
        return dxcm, dxtime
        
        
    def get_pupil_data(self, data_set):
        try:
            pupilloc=data_set.get_pupil_location()
            pupilsize=data_set.get_pupil_size()
            
            return pupilloc, pupilsize
        except:
            print('No pupil')
      
    def mask_lists(self, data_set, cell_ids):
        #%masks

        rois_ids=data_set.get_roi_ids()
        roi_mask_objects=data_set.get_roi_mask(cell_specimen_ids=cell_ids)
        all_roi_masks=data_set.get_roi_mask_array()
        
        return rois_ids, roi_mask_objects, all_roi_masks
    
    def max_projection(self,data_set):
        
        max_projection = data_set.get_max_projection()

        return max_projection

#%% general stim routines
    def exploring_stimulus(self, data_set):
        print(data_set.get_session_type())
        print(data_set.list_stimuli())
        print(data_set.get_stimulus_epoch_table())
       
        
        if data_set.get_session_type()=='three_session_A':
             stim_templ_mov_one=data_set.get_stimulus_template(data_set.list_stimuli()[1])
             stim_templ_mov_three=data_set.get_stimulus_template(data_set.list_stimuli()[2])
             
             stim_tmplates=(stim_templ_mov_one,stim_templ_mov_three)
             
        elif data_set.get_session_type()=='three_session_C2':
             stim_templ_noise_4=data_set.get_stimulus_template(data_set.list_stimuli()[0])
             stim_templ_noise_8=data_set.get_stimulus_template(data_set.list_stimuli()[1])
             stim_templ_mov_one=data_set.get_stimulus_template(data_set.list_stimuli()[2])
             stim_templ_mov_two=data_set.get_stimulus_template(data_set.list_stimuli()[3])
             stim_tmplates=(stim_templ_noise_4,stim_templ_noise_8,stim_templ_mov_one,stim_templ_mov_two)


        elif data_set.get_session_type()=='three_session_B':
             stim_templ_mov_one=data_set.get_stimulus_template(data_set.list_stimuli()[0])
             stim_templ_scenes=data_set.get_stimulus_template(data_set.list_stimuli()[1])
             stim_tmplates=(stim_templ_mov_one,stim_templ_scenes)


        stimtables=[]
        for stim in data_set.list_stimuli():
            stimtables.append(data_set.get_stimulus_table(stim))

        return stimtables , stim_tmplates
        # frame=2000
        # data_set.get_stimulus(frame)# buged in some datasets
  
        
        
    def plotting_traces_and_stim(self, data_set, dff_traces, dxcm):
    
            stim_epoch = data_set.get_stimulus_epoch_table()
            
            color=[ 'y','c','m','tab:orange','tab:purple','tab:pink','r', 'g', 'b']
            colors={}
            
            labels = ['{}'.format(t) for t in stim_epoch.stimulus.unique()]

            for j, stim in enumerate(stim_epoch.stimulus.unique()):
                colors[stim]=color[j]
            

            fig, ax = plt.subplots(1, figsize=(10,10))
            for i in range(dff_traces.shape[0]):
                ax.plot(dff_traces[i,:]+(i*2), color='gray')
     
            for c, stim_name in enumerate(stim_epoch.stimulus.unique()):
                stim = stim_epoch[stim_epoch.stimulus==stim_name]
                for j in range(len(stim)):
                    ax.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[stim_name], alpha=0.3)
                    
            fig, ax = plt.subplots(1, figsize=(16,16))
            for i in range(dff_traces.shape[0]):
                ax.plot(dff_traces[i,:]+(i*2), color='gray')
            ax.plot((0.2*dxcm)-20)
            ax.set_ylim(-30,dff_traces.shape[0]+i+10)
                
            #for each stimulus, shade the plot when the stimulus is presented
            for c, stim_name in enumerate(stim_epoch.stimulus.unique()):
                stim = stim_epoch[stim_epoch.stimulus==stim_name]
                for j in range(len(stim)):
                    ax.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[stim_name], alpha=0.3,  label = stim_name)
                    
            ax.legend(frameon=False, loc='center left',bbox_to_anchor=(-0.1,0.5))
      
                    
    def transform_global_stim_table_to_epoch_specific(self,data_set,  stim):
        
        stim_epoch = data_set.get_stimulus_epoch_table()
        specific=stim_epoch[stim_epoch['stimulus']==stim]
        stim_table = data_set.get_stimulus_table(stim)
        corrected_stim_table=stim_table.copy()

        test1=stim_table['end']<=specific['end'].iloc[0]
        test2=stim_table['end']<=specific['end'].iloc[1]
        first_index_second=test1[::-1].idxmax()+1
        first_index_third=test2[::-1].idxmax()+1

        firstperiodlag=specific['start'].iloc[0]
        secondperiodlag=specific['start'].iloc[1]-specific['end'].iloc[0]
        thirdperiodlqag=specific['start'].iloc[2]-specific['end'].iloc[1]
        
        corrected_stim_table.loc[:, ['start','end']] = corrected_stim_table[['start','end']].apply(lambda x: x - firstperiodlag)
        corrected_stim_table.loc[first_index_second:, ['start','end']] = corrected_stim_table[['start','end']].apply(lambda x: x - secondperiodlag)
        corrected_stim_table.loc[first_index_third:, ['start','end']] = corrected_stim_table[['start','end']].apply(lambda x: x - thirdperiodlqag)
            

        return corrected_stim_table
        
        
    def explore_spontaneous_activity(self, data_set, dff_traces, spikes, selected_cell_index=0, plot=False):
        stim_epoch = data_set.get_stimulus_epoch_table()
        stim_table=stim_epoch[stim_epoch['stimulus']=='spontaneous']
        
        return stim_table
        

   #%% natural movies     
        
    def explore_natural_movie_analysis(self, data_set, dff_traces, spikes, selected_cell_index=0, plot=False):
        print(data_set.list_stimuli())
        cell_index=15
  
        natmov1obj=NaturalMovie(data_set, 'natural_movie_one')
        natmov2obj=NaturalMovie(data_set, 'natural_movie_two')
        natmov3obj=NaturalMovie(data_set, 'natural_movie_three')
        
        if data_set.get_session_type()=='three_session_A':
            peak1=natmov1obj.get_peak()
            sweep1=natmov1obj.get_sweep_response()
            peak3=natmov3obj.get_peak()
            sweep3=natmov3obj.get_sweep_response()
   
            peak=[peak1,peak3]
            sweep=[sweep1,sweep3]

        elif data_set.get_session_type()=='three_session_C2':
            
            peak1=natmov1obj.get_peak()
            sweep1=natmov1obj.get_sweep_response()
            peak2=natmov2obj.get_peak()
            sweep2=natmov2obj.get_sweep_response()
            
            peak=[peak1,peak2]
            sweep=[sweep1,sweep2]

        elif data_set.get_session_type()=='three_session_B':

            peak=natmov1obj.get_peak()
            sweep=natmov1obj.get_sweep_response()
            
            peak=[peak,]
            sweep=[sweep,]

        if plot:
            natmov1obj.open_track_plot(cell_index=cell_index)
            natmov3obj.open_track_plot(cell_index=cell_index)
        
        return natmov1obj,natmov2obj, peak,sweep       
        
       #%% drifting gratings
    def explore_drifting_analysis(self, data_set, dff_traces, spikes, selected_cell_index=0, plot=False):
        
        activity=dff_traces
        activity=spikes
        
        drift_obj=DriftingGratings(data_set)
        # peak_info=drift_obj.get_peak() #gives error
        reponse=drift_obj.get_response() # not proper dimensions
        noise_cor=drift_obj.get_noise_correlation()
        rep_sim=drift_obj.get_representational_similarity()
        sig_cor=drift_obj.get_signal_correlation()
       
        stim_epoch = data_set.get_stimulus_epoch_table()
        selected_frames=stim_epoch[stim_epoch.stimulus=='drifting_gratings']
        drifting_gratings_table = data_set.get_stimulus_table('drifting_gratings')
        new_stim_table=self.transform_global_stim_table_to_epoch_specific(data_set, 'drifting_gratings')

        orient = drift_obj.orivals[peak_info.iloc[selected_cell_index,:]['ori_dg']]
        ops_orient = drift_obj.orivals[peak_info.iloc[selected_cell_index,:]['ori_dg']+4]
        tf = drift_obj.tfvals[peak_info.iloc[selected_cell_index,:]['tf_dg']]
        print(peak_info.iloc[selected_cell_index,:])


        if plot:
            
            fig, ax = plt.subplots(1, figsize=(16,16))
            ax.imshow(reponse[:,1:,selected_cell_index,0])    
            
            plt.imshow(sig_cor[0])
            plt.imshow(sig_cor[1])
            plt.imshow(rep_sim[0])
            plt.imshow(rep_sim[1])



            plt.imshow(noise_cor[0][:,:,0,0])
            plt.imshow(noise_cor[1][:,:,0,0])

            plt.imshow(noise_cor[2])
            plt.imshow(noise_cor[3])
            drift_obj.plot_direction_selectivity()
            drift_obj.plot_orientation_selectivity()
            drift_obj.plot_preferred_direction()
            drift_obj.plot_preferred_temporal_frequency()
            
            drift_obj.open_star_plot(cell_index=selected_cell_index)

 
            fig, ax = plt.subplots(1, figsize=(16,16))
            for i in range(activity.shape[0]):            
                ax.plot(activity[i,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]]+(i*2), color='gray')
            ax.plot((0.2*dxcm[np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]])-10)
            stim_subset = new_stim_table[(new_stim_table['orientation'].isin([orient, ops_orient])) & (new_stim_table.temporal_frequency==tf)]
            for j in range(len(stim_subset)):
                ax.axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)

            fig, ax = plt.subplots(2, figsize=(16,16), sharex=True)
            ax[0].plot(dff_traces[selected_cell_index,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]], color='gray')
            ax[0].plot((0.2*dxcm[np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]])-10)
            for j in range(len(stim_subset)):
                ax[0].axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)
            ax[1].plot(spikes[selected_cell_index,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]], color='gray')
            ax[1].plot((0.2*dxcm[np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]])-10)
            for j in range(len(stim_subset)):
                ax[1].axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)
                

            fig, ax = plt.subplots(1, figsize=(16,16))
            for i in range(activity.shape[0]):
                ax.plot(activity[i,:]+(i*2), color='gray')
            ax.plot((0.2*dxcm)-10)
            stim_subset = drifting_gratings_table[(drifting_gratings_table.orientation==orient)&(drifting_gratings_table.temporal_frequency==tf)]
            for j in range(len(stim_subset)):
                ax.axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)

            fig, ax = plt.subplots(1, figsize=(16,16))
            ax.plot(activity[selected_cell_index,:], color='gray')
            ax.plot((0.2*dxcm)-10)
            for j in range(len(stim_subset)):
                ax.axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)
            
            
        
        return drift_obj, peak_info, reponse, noise_cor, rep_sim, sig_cor, drifting_gratings_table
        
        
        #%% natural scenes
    def explore_natural_scenes(self, data_set, dff_traces, spikes,selected_cell_index=0, plot=False):
        
        print(data_set.list_stimuli())
        stim_epoch = data_set.get_stimulus_epoch_table()
        natural_scenes_table = data_set.get_stimulus_table('natural_scenes')

        new_stim_table=self.transform_global_stim_table_to_epoch_specific(data_set,selected_cell_index, 'natural_scenes')


        scenes_obj=NaturalScenes(data_set)

        peak_info=scenes_obj.get_peak() #gives error
        reponse=scenes_obj.get_response() # not proper dimensions
        noise_cor=scenes_obj.get_noise_correlation()
        rep_sim=scenes_obj.get_representational_similarity()
        sig_cor=scenes_obj.get_signal_correlation()
        
       
        
        
        natural_scene_table = data_set.get_stimulus_table('natural_scenes')
        natural_scene_table.head()
        natural_scene_template = data_set.get_stimulus_template('natural_scenes')
        natural_scene_template.shape



        prefimage=peak_info.loc[selected_cell_index,'scene_ns']
        selected_frames=stim_epoch[stim_epoch.stimulus=='natural_scenes']
        
        if plot:
            scenes_obj.open_corona_plot(cell_index=selected_cell_index)
            scenes_obj.plot_time_to_peak()
            
            plt.imshow(natural_scene_template[scene_number,:,:], cmap='gray')

            
            fig, ax = plt.subplots(1, figsize=(16,16))
            for i in range(dff_traces.shape[0]):
                ax.plot(dff_traces[i,:]+(i*2), color='gray')
            stim_subset = natural_scene_table[natural_scene_table.frame==prefimage]
            for j in range(len(stim_subset)):
                plt.axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)
    
            
            fig, ax = plt.subplots(1, figsize=(16,16))
            for i in range(dff_traces.shape[0]):
                ax.plot(dff_traces[i,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]]+(i*2), color='gray')
            stim_subset = new_stim_table[new_stim_table.frame==prefimage]
            for j in range(len(stim_subset)):
                plt.axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)
                
            
        return scenes_obj, peak_info, reponse, noise_cor, rep_sim, sig_cor, natural_scenes_table   
            
#%%static gratings

    def explore_static_gratings(self, data_set, dff_traces, spikes,selected_cell_index=0, plot=False)  :
        pass
            
 #%%sparse noise       
        
    def explore_locally_sparse_noise_4deg(self, data_set, dff_traces, spikes,selected_cell_index=0, plot=False):
        
        stimname='locally_sparse_noise_4deg'
        noise_4deg_stim_table=data_set.get_stimulus_table(data_set.list_stimuli()[0])
        noise_templ=data_set.get_locally_sparse_noise_stimulus_template(stimname)
        lsn_movie, offscreen_mask =  data_set.get_locally_sparse_noise_stimulus_template(stimname,
                                                                                    mask_off_screen=True)
        
        frame = 200
        plt.imshow(lsn_movie[frame,:,:], cmap='gray', interpolation='nearest')
        plt.axis('off')
        plt.title('frame %d' % frame)
        plt.show()

        # find frames at a given grid location that are 'on'
        loc = (10,15)
        on_frames = np.where(lsn_movie[:,loc[0],loc[1]] == LocallySparseNoise.LSN_ON)[0]
        
        # pull these trials out of the stimulus table
        noise_4deg_stim_table = noise_4deg_stim_table.loc[on_frames]
        
        self.plot_stimulus_table(noise_4deg_stim_table, "loc (%d,%d) " % loc)


#%% misc
    def saving_stuff(self, data_set):
        
        data_set.save_analysis_arrays()
        data_set.save_analysis_dataframes()
        # sio.savemat(os.path.join(self.main_directory,'corrected_traces.mat'),{'corrected_traces':corrected_traces,'dff_traces':dff_traces})   
   
        # movie_path = '/data/allen-brain-observatory/visual-coding-2p/ophys_movies/ophys_experiment_%d.h5' % exp['id']
        # f = h5py.File(movie_path,'r')
        # frames = f["data"][:10,:,:]
        # f.close()

    

    def mask_module(self, data_set, cell_ids, roi_mask_objects,all_roi_masks, max_projection):   
        cell_ids=cell_info[1][0:2]
        roi_mask_objects=masks[1]
      
        f, axes = plt.subplots(1, len(cell_ids)+2, figsize=(15, 3))
        for ax, roi_mask, cid in zip(axes[:-2], roi_mask_objects, cell_ids):
            ax.imshow(roi_mask.get_mask_plane(), cmap='gray')
            ax.set_title('cell %d' % cid)
        combined_mask = masks[2].max(axis=0)
        axes[-2].imshow(combined_mask, cmap='gray')
        axes[-2].set_title('all ROIs')
        axes[-1].imshow(projection, cmap='gray')
        axes[-1].set_title('max projection')
        plt.show()
        
    def dff_module(self,data_set, corrected_traces)  : 
        
        celltoplot=1
        time, raw_traces = data_set.get_fluorescence_traces()
        selected_corrected_trace=        np.expand_dims(corrected_traces[celltoplot],0)

        # get precomputed df/f
        _, dff_traces = self.data_set.get_dff_traces()
        selected_precomputed_dff_trace =dff_traces[celltoplot]
        recaclulated=dff.calculate_dff(selected_corrected_trace).squeeze()
        dff_windowed_mode=dff.calculate_dff(selected_corrected_trace, dff.compute_dff_windowed_mode).squeeze()

        fig, ax=plt.subplots(3,figsize=(14,4), sharex=True)
        # warning: dF/F can occasionally be one element longer or shorter 
        # than the time stamps for the original traces.
        ax[0].plot(time, selected_corrected_trace[0,:])
        ax[1].plot(time, selected_precomputed_dff_trace)
        ax[2].plot(time, recaclulated)
        cor=scp.signal.correlate(selected_precomputed_dff_trace, recaclulated)
        cor /= np.max(cor)
        plt.plot(cor)
        
        cor=scp.signal.correlate(selected_precomputed_dff_trace, dff_windowed_mode)
        cor /= np.max(cor)
        plt.plot(cor)
 


    def plot_stimulus_table(self, stim_table, title):
         fstart = stim_table.start.min()
         fend = stim_table.end.max()
         
         fig = plt.figure(figsize=(15,1))
         ax = fig.gca()
         for i, trial in stim_table.iterrows():    
             x1 = float(trial.start - fstart) / (fend - fstart)
             x2 = float(trial.end - fstart) / (fend - fstart)            
             ax.add_patch(patches.Rectangle((x1, 0.0), x2 - x1, 1.0, color='r'))
         ax.set_xticks((0,1))
         ax.set_xticklabels((fstart, fend))
         ax.set_yticks(())
         ax.set_title(title)
         ax.set_xlabel("frames")


    def Get_all_mouse_from_trangenic_line(self, transgenic_line):
        all_line_containers=self.allen_manifest.get_experiment_containers(cre_lines=[transgenic_line],)   
        all_line_mouse = sorted(dict(Counter([cont['donor_name']  for cont in all_line_containers ])).items(), key=lambda x: x[1], reverse=True)
        
        if transgenic_line in [self.allen_manifest.get_all_cre_lines()[5],self.allen_manifest.get_all_cre_lines()[10],self.allen_manifest.get_all_cre_lines()[12]]:
            
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
                             
                             
                             
 #%% allen visual stim templates                                                                                          
    def get_visual_templates(self):
        self.all_stim =  self.allen_manifest.get_all_stimuli()
        visual_area = 'VISp'
        cre_line ='Cux2-CreERT2'
        self.exp_containers= self.allen_manifest.get_experiment_containers(
            targeted_structures=[visual_area], 
            cre_lines=[cre_line])
        self.exps = pd.DataFrame(self.exp_containers)

          
        first_experiment_container_id = self.exp_containers[0]['id']
        # exps= self.allen_manifest.get_ophys_experiments(
        #     experiment_container_ids=[first_experiment_container_id],
        # )      
        for stim in self.allen_manifest.get_all_stimuli():
            if 'grating' not in stim and 'spont' not in stim and 'deg' not in stim:
                module_logger.info(stim)
                self.selected_exp_list= self.allen_manifest.get_ophys_experiments(
                        experiment_container_ids=[first_experiment_container_id], 
                        stimuli=[stim])

                session_id =  self.selected_exp_list[0]['id']
        
                self.data_set = self.allen_manifest.get_ophys_experiment_data(session_id)
                module_logger.info("Metadata from NWB file:")
                module_logger.info(self.data_set.get_metadata())
                # self.full_stim_epoch_table=self.data_set.get_stimulus_epoch_table()
                self.stim_by_frame_table=self.data_set.get_stimulus_table(stim)
    
                # z=self.data_set.list_stimuli()
                self.all_templates=self.data_set.get_stimulus_template(stim)
                # self.all_template_shown=self.data_set.get_stimulus_template(stim)[self.stim_by_frame_table['frame'].values,:,:]
                # self.selected_stim_epoch_table=self.full_stim_epoch_table.loc[self.full_stim_epoch_table['stimulus'] == stim]

                
                stimfiles=glob.glob(self.main_directory+'\\**'+stim+'**.mat')
                if not stimfiles:
                    self.all_warped_images=np.zeros((self.all_templates.shape[0],self.screen_size[0],self.screen_size[1]))
                    module_logger.info('processing '+stim)
                    for i in range(self.all_templates.shape[0]):
    
                            template=(self.all_templates[i,:,:])
        
                            if 'noise' in stim:
                                # self.obs_m_img_screen=self.obs_monitor.lsn_image_to_screen(template)
                                self.mym_img_screen=self.my_monitor.lsn_image_to_screen(template, stim)
        
                            elif 'scenes' in stim:
                                # self.obs_m_img_screen=self.obs_monitor.natural_scene_image_to_screen(template, origin='upper')
                                # self.obs_monitor.show_image(self.obs_m_img_screen, warp=True, origin='upper') 
                                if template.shape[0]>self.screen_size[0]:
                                    img = Image.fromarray(np.uint8(template))
                                    resized_img = img.resize((template.shape[1], self.screen_size[0]))
                                    template=np.asarray(resized_img)              
                                    si.NATURAL_SCENES_PIXELS=(self.screen_size[0],template.shape[1])
                                self.mym_img_screen=self.my_monitor.natural_scene_image_to_screen(template, origin='upper')
                            elif 'movie' in stim:
                                # self.obs_m_img_screen=self.obs_monitor.natural_movie_image_to_screen(template, origin='upper')
                                si.NATURAL_MOVIE_PIXELS=(700 ,1030)
                                self.mym_img_screen=self.my_monitor.natural_movie_image_to_screen(template, origin='upper')
                            self.mym_img_screen_warped=spndi.map_coordinates(self.mym_img_screen, self.expgeom.warp_coordinates.T).reshape((self.my_monitor.n_pixels_r,self.my_monitor.n_pixels_c))
                            self.all_warped_images[i,:,:]= self.mym_img_screen_warped   
                       
                    if self.all_warped_images[0,:,:].any():
                        
                        if self.all_warped_images.shape[0]==9000:
                            parts=5     
                        elif self.all_warped_images.shape[0]==3600:
                            parts=4              
                        else:
                            parts=1
     
                        self.all_warped_images_parts=np.zeros((int(self.all_warped_images.shape[0]/parts),self.all_warped_images.shape[1],self.all_warped_images.shape[2]))
                        for i in range(parts):
                            
                            self.all_warped_images_parts[:,:,:]=self.all_warped_images[int(self.all_warped_images.shape[0]*i/parts):int(self.all_warped_images.shape[0]*(i+1)/parts),:,:]
                            self.all_warped_images_parts=self.all_warped_images_parts.astype('uint8')
    
                            sio.savemat(os.path.join(self.main_directory,stim+'_' +str(i+1)+'.mat'),{stim+'_all_warped_frames':self.all_warped_images_parts})              
  
    def get_gratings(self):
        module_logger.info("Getting gratings:")
        oris=np.linspace(0,150,6)
        spatials=[0.02,0.04,0.08,0.16,0.32]
        phases=[0,0.25,0.5,0.75]
        self.gratings=np.zeros((self.screen_size[0],self.screen_size[1],6,5,4))        
        for i,ori in enumerate(oris):
            for j,spat in enumerate(spatials):
                for k,ph in enumerate(phases):
                    self.gratings[:,:,i,j,k]=spndi.map_coordinates(self.my_monitor.grating_to_screen(ph, spat, ori, 15) , self.expgeom.warp_coordinates.T).reshape((self.my_monitor.n_pixels_r,self.my_monitor.n_pixels_c))  
        self.gratings=self.gratings.astype('uint8')     
        sio.savemat(os.path.join(self.main_directory,'static_gratings.mat'),{'all_warped_static_gratings':self.gratings})   
      
    def get_drifting_gratings(self):
        module_logger.info("Getting drifting gratings:")

        oris=np.linspace(0,360-(360/8),8)
        spatial=0.04
        temporals=[1,2,4,8,15]
        frequency=60
        self.drifting_gratings=np.zeros((self.screen_size[0],self.screen_size[1],8,5,60))
        for i,ori in enumerate(oris):
            module_logger.info("Processing drifting gratings: "+str(ori))

            for j,temp in enumerate(temporals):
              for k in range(frequency):
                  self.drifting_gratings[:,:,i,j,k]=spndi.map_coordinates(self.my_monitor.grating_to_screen(k*temp/frequency, spatial, ori, 15), self.expgeom.warp_coordinates.T).reshape((self.my_monitor.n_pixels_r,self.my_monitor.n_pixels_c))  
                    
        parts=5
        self.all_warped_images_parts=np.zeros((self.drifting_gratings.shape[0],self.drifting_gratings.shape[1],self.drifting_gratings.shape[2],int(self.drifting_gratings.shape[3]/parts),self.drifting_gratings.shape[4]))
        for i in range(parts):  
            module_logger.info("Saving drifting gratings:")

            self.all_warped_images_parts[:,:,:,:,:]=self.drifting_gratings[:,:,:,int(self.drifting_gratings.shape[3]*i/parts):int(self.drifting_gratings.shape[3]*(i+1)/parts),:]
            self.all_warped_images_parts=self.all_warped_images_parts.astype('uint8')
            sio.savemat(os.path.join(self.main_directory,'drifting_gratings'+'_' +str(i+1)+'.mat'),{'all_warped_driting_gratings':self.all_warped_images_parts})   
            
    def create_video_file(self):
   

        movie_path = r'C:\Users\sp3660\Documents\Projects\AllenBrainObservatory\RawVideos\ophys_experiment_662358233.h5'
        hf = h5py.File(movie_path, 'r')
        frames = hf["data"][:10,:,:]
        hf.keys()
        n1 = hf.get('data')
        n1 = np.array(n1)
        n1.shape
        mov=cm.movie(n1)
        mov.play(fr=300)

        pass
    
    def vis_stim_plotting(self):
        
        
        # scene_nums = [10, 2882]
        # fig, axes = plt.subplots(1,len(scene_nums))
        # for ax,scene in zip(axes, scene_nums):
        #     ax.imshow(self.all_templates[scene,:,:], cmap='gray')
        #     ax.set_axis_off()
        #     ax.set_title('scene %d' % scene)
        
        
        # for i in range(20):
        #     params, template = self.data_set.get_stimulus(self.selected_stim_epoch_table['start'].iloc[2]+i)
        #     params, template = self.data_set.get_stimulus(self.selected_stim_epoch_table['end'].iloc[2]-1)
        #     plt.figure()
        #     img = plt.imshow(template, cmap=plt.cm.gray, interpolation='none')
        #     pd.Series(params[2])

        # img = plt.imshow(self.all_templates[1,:,:], cmap=plt.cm.gray, interpolation='none')
        
        # params, template = self.data_set.get_stimulus(self.selected_stim_epoch_table['start'].iloc[0]+500)
        # plt.imshow(template)


        # self.obs_monitor.show_image(self.obs_m_img_screen, warp=False, mask=False, origin='upper')
        # self.obs_monitor.show_image(self.obs_m_img_screen, warp=False, mask=True, origin='upper')
        # self.obs_monitor.show_image(self.obs_m_img_screen, warp=True, origin='upper')      
                                  
        # self.my_monitor.show_image(self.mym_img_screen, mask=False, origin='upper')
        # self.my_monitor.show_image(self.mym_img_screen, mask=True, origin='upper')
        # self.my_monitor.show_image(self.mym_img_screen_warped, mask=False, origin='upper')

        pass
   
 




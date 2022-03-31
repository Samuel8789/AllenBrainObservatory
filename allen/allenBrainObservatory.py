# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 09:13:02 2021

@author: sp3660
"""
from project_manager.ProjectsCLass import Project

import sys
import os 
import glob
import logging 
module_logger = logging.getLogger(__name__)
from collections import Counter

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as patches
import numpy as np
import scipy.io as sio
from scipy.interpolate import griddata
import scipy.ndimage.interpolation as spndi
import scipy.ndimage
import scipy as scp

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
from allensdk.brain_observatory.locally_sparse_noise import LocallySparseNoise
from allensdk.brain_observatory.drifting_gratings import DriftingGratings
import allensdk.brain_observatory.receptive_field_analysis.visualization as rfvis
import allensdk.brain_observatory.receptive_field_analysis.receptive_field as rf
import matplotlib.pyplot as plt
import numpy as np
from allensdk.brain_observatory.locally_sparse_noise import LocallySparseNoise
import matplotlib.patches as patches
import matplotlib.pyplot as plt


class AllenBrainObservatory(Project):    
    def __init__(self,  githubtoken_path, gui=False ):
        Project.__init__(self, 'AllenBrainObservatory',githubtoken_path, Project.computer, Project.platform )
        self.screen_size=[900,1440]       
        self.main_directory=self.project_paths['Documents']  
        self.allen_manifest =BrainObservatoryCache(manifest_file=os.path.join(self.main_directory, 'boc\manifest.json'))
        
        self.my_monitor=si.Monitor( self.screen_size[0],self.screen_size[1], 48.26,'cm')
        self.expgeom=si.ExperimentGeometry(15, self.my_monitor.height, self.my_monitor.width, (self.screen_size[1],self.screen_size[0]), (0.5,0.5))
        self.obs_monitor=si.BrainObservatoryMonitor()
        # self.get_visual_templates()
        # self.get_gratings()
        # self.get_drifting_gratings()
        
        
    def select_lines_and_locations(self):
        self.targeted_structures =  self.allen_manifest.get_all_targeted_structures()
        self.all_cre_lines = self.allen_manifest.get_all_cre_lines()
        self.all_imaging_depths= self.allen_manifest.get_all_imaging_depths()
        self.all_reporters=self.allen_manifest.get_all_reporter_lines()
        
        line='Vglut1'
        self.visual_areas = ['VISp']
        if line=='Vglut1':
            self.cre_lines =[self.all_cre_lines[9]]
        self.depths=[175]
        self.reporters=['Ai162(TIT2L-GC6s-ICL-tTA2)']
        self.reporters=['Ai93(TITL-GCaMP6f)']

        
    def select_stimulus(self):
        self.allen_manifest.get_all_session_types()
        self.allen_manifest.get_all_stimuli()
        self.stimuli=['locally_sparse_noise_4deg']
        
        
    def get_multiimaged_fov_container(self):
        container_selection_index=0
        # a container correspond to same fov presented with three different stimulus
        self.exp_containers= self.allen_manifest.get_experiment_containers(
            targeted_structures=self.visual_areas, 
            cre_lines=self.cre_lines,
            imaging_depths=self.depths,
            reporter_lines=self.reporters,
            )
       
        self.containers = pd.DataFrame(self.exp_containers)
        selected_container_id=self.exp_containers[container_selection_index]['id']
        return selected_container_id
    
    def get_imaging_session_from_container(self, selected_container_id):
        exp_selection_index=0
        
        
         # by container ID, this ere experiment form single fov from single mouse, however unknown vis stim
        self.selected_exp_list1= self.allen_manifest.get_ophys_experiments(
                experiment_container_ids=[selected_container_id])
        
        # byexp characteristiscs, these includes all mice with these charcateristics, this is to select by visual stim or by cell specimen
        self.selected_exp_list2=self.allen_manifest.get_ophys_experiments(
            targeted_structures=self.visual_areas, 
            cre_lines=self.cre_lines, 
            imaging_depths=self.depths,
            stimuli=self.stimuli,
            transgenic_lines=self.reporters,
            session_types=None,
            cell_specimen_ids=None
            )
        self.selected_exp_list2
        self.exps = pd.DataFrame(self.selected_exp_list2)
        selected_session_id= self.selected_exp_list2[exp_selection_index]['id']
        
        
    def dealing_with_cells(self, selected_container_id):
        # too slow do not use
        self.all_container_fov_cells=self.allen_manifest.get_cell_specimens(file_name=None, 
                                               ids=None, 
                                               experiment_container_ids=[selected_container_id],
                                               include_failed=False,
                                               simple=True, 
                                               filters=None
                                               )



    def download_single_imaging_session_nwb(self, selected_session_id):
        
        self.session_stimuli=self.allen_manifest.get_ophys_experiment_stimuli(selected_session_id)
        self.data_set = self.allen_manifest.get_ophys_experiment_data(selected_session_id)
        self.deconvolved_spikes=self.allen_manifest.get_ophys_experiment_events(selected_session_id)
        
        self.dataset_exploration(self.data_set)
        
    def dataset_exploration(self, data_set):
        data_set=self.data_set
        self.get_cell_info(data_set)
        self.get_all_traces(data_set)
        self.get_neuropil_info(data_set)
        self.get_locomotion(data_set)
        self.dataset_metadata(data_set)
        self.get_pupil_data(data_set)
        self.dff_module(data_set)
        
    def get_cell_info(self, data_set):
        
        cell_number=data_set.number_of_cells
        cell_ids=data_set.get_cell_specimen_ids() # global cell ID
        cell_idxs=data_set.get_cell_specimen_indices(cell_ids) # this experiment cell indexes
        
   
       
    def get_all_traces(self, data_set):
        
        timestamps=data_set.get_fluorescence_timestamps()
        time, raw_traces = data_set.get_fluorescence_traces()
        _, demixed_traces = data_set.get_demixed_traces()
        _, neuropil_traces = data_set.get_neuropil_traces()
        _, corrected_traces = data_set.get_corrected_fluorescence_traces()
        _, dff_traces = data_set.get_dff_traces()
        
    def get_neuropil_info(self, data_set):
        
        neuropil_r=data_set.get_neuropil_r()
        _, neuropil_traces=data_set.get_neuropil_traces()
        
    def dataset_metadata(self, data_set):
        metadata=data_set.get_metadata()
        data_set.nwb_file
        data_set.pipeline_version
        motion_correct_shifts=data_set.get_motion_correction()
        
    def get_locomotion(self, data_set):
        dxcm, dxtime = self.data_set.get_running_speed()
        # plt.figure(figsize=(14,4))
        # plt.plot(dxtime, dxcm)
        # plt.show()
        
        
    def get_pupil_data(self, data_set):
        try:
            data_set.get_pupil_location()
            data_set.get_pupil_size()
        except:
            print('No pupil')
      
    def mask_lists(self, data_set, cell_ids):
        #%%masks

        rois_ids=data_set.get_roi_ids()
        roi_mask_objects=data_set.get_roi_mask(cell_specimen_ids=cell_ids)
        all_roi_masks=data_set.get_roi_mask_array()
        
    def max_projection(self,data_set):
        
        max_projection = data_set.get_max_projection()




    def exploring_stimulus(self, data_set, dff_traces):
        data_set.get_session_type()
        data_set.list_stimuli()
        
        stim_templ=self.data_set.get_stimulus_template(self.data_set.list_stimuli()[0])
        spont_stim_table=self.data_set.get_stimulus_table(self.data_set.list_stimuli()[4])

        
        if 'locally_sparse_noise_4deg' in data_set.list_stimuli():
            self.explore_locally_sparse_noise_4deg(data_set)
            
        frame=2000
        data_set.get_stimulus(frame)# buged in some datasets
        data_set.get_stimulus_epoch_table()# buged in some datasets
        # data_set.stimulus_search# buged in some datasets
        
        
        
    def plotting_traces_and_stim(self, data_set, dff_traces, dxcm):
        
        stim_epoch = data_set.get_stimulus_epoch_table()
        stim_epoch
        stim_epoch.stimulus.unique()

        fig = plt.figure(figsize=(10,8))
        for i in range(50):
            plt.plot(dff_traces[i,:]+(i*2), color='gray')
            
        #for each stimulus, shade the plot when the stimulus is presented
        colors = {'locally_sparse_noise_4deg':   'blue',
                  'natural_movie_two':    'orange',
                  'locally_sparse_noise_8deg':    'yellow',
                  'spontaneous':       'green',
                  'natural_movie_one': 'red'}

        for c, stim_name in enumerate(stim_epoch.stimulus.unique()):
            stim = stim_epoch[stim_epoch.stimulus==stim_name]
            for j in range(len(stim)):
                plt.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[stim_name], alpha=0.1)
                
                
                
        fig = plt.figure(figsize=(10,10))
        for i in range(50):
            plt.plot(dff_traces[i,:]+(i*2), color='gray')
        plt.plot((0.2*dxcm)-20)
            
        #for each stimulus, shade the plot when the stimulus is presented
        colors = ['blue','orange','green','red']
        for c,stim_name in enumerate(stim_epoch.stimulus.unique()):
            stim = stim_epoch[stim_epoch.stimulus==stim_name]
            for j in range(len(stim)):
                plt.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[c], alpha=0.1)
        #%%
       
    def explore_drifitng_analysis(self):
        # to review
        pass
        # selected_frames=stim_epoch[stim_epoch.stimulus=='drifting_gratings']
        # fig = plt.figure(figsize=(10,8))
        # for i in range(50):
        #     plt.plot(2.5*SelectedEvents[i,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]]
        # +(i*2), color='gray')

        # grat_only=SelectedEvents[i,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]]
  
        # #for each stimulus, shade the plot when the stimulus is presented
        # stim_epoch.stimulus.unique()



        # fig = plt.figure(figsize=(10,8))
        # for i in range(50):
        #     plt.plot(SelectedEvents[i,:]+(i*2), color='gray')
            
        # #for each stimulus, shade the plot when the stimulus is presented
        # colors = {'drifting_gratings':   'blue',
        #           'natural_movie_three':    'orange',
        #           'spontaneous':       'green',
        #           'natural_movie_one': 'red'}

        # for c,stim_name in enumerate(stim_epoch.stimulus.unique()):
        #     stim = stim_epoch[stim_epoch.stimulus==stim_name]
        #     for j in range(len(stim)):
        #         plt.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[stim_name], alpha=0.1)


        # drifting_gratings_table = SelectedNWB.get_stimulus_table('drifting_gratings')


        #%%
        # dxcm, tsd = SelectedNWB.get_running_speed()
        # fig = plt.figure(figsize=(10,10))
        # for i in range(50):
        #     plt.plot(SelectedEvents[i,:]+(i*2), color='gray')
        # plt.plot((0.2*dxcm)-20)
            
        # #for each stimulus, shade the plot when the stimulus is presented
        # colors = ['blue','orange','green','red']
        # for c,stim_name in enumerate(stim_epoch.stimulus.unique()):
        #     stim = stim_epoch[stim_epoch.stimulus==stim_name]
        #     for j in range(len(stim)):
        #         plt.axvspan(xmin=stim.start.iloc[j], xmax=stim.end.iloc[j], color=colors[c], alpha=0.1)


        # #%%
        # fig = plt.figure(figsize=(10,8))
        # for i in range(50):
        #     plt.plot(    plt.plot(2.5*SelectedEvents[i,np.r_[selected_frames.iloc[0][1]:selected_frames.iloc[0][2],selected_frames.iloc[1][1]:selected_frames.iloc[1][2],selected_frames.iloc[2][1]:selected_frames.iloc[2][2]]]
        # +(i*2), color='gray')+(i*2), color='gray')
            
        # orient = drifting_gratings_table.orientation.iloc[0]               
        # #shade traces with the time of each presentation of the above scene
        # stim_subset = drifting_gratings_table[drifting_gratings_table.orientation==orient]
        # for j in range(len(stim_subset)):
        #     plt.axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)
        
        
        
        
        
        
    def explore_natural_scenes(self,data_set):
        
        natural_scene_table = data_set.get_stimulus_table('natural_scenes')
        natural_scene_table.head()
        natural_scene_template = data_set.get_stimulus_template('natural_scenes')
        natural_scene_template.shape
        scene_number = natural_scene_table.frame.iloc[3]
        plt.imshow(natural_scene_template[scene_number,:,:], cmap='gray')

        fig = plt.figure(figsize=(10,8))
        for i in range(50):
            plt.plot(dff[i,:]+(i*2), color='gray')
            
        #shade traces with the time of each presentation of the above scene
        stim_subset = natural_scene_table[natural_scene_table.frame==scene_number]
        for j in range(len(stim_subset)):
            plt.axvspan(xmin=stim_subset.start.iloc[j], xmax=stim_subset.end.iloc[j], color='red', alpha=0.4)
        
        
    def explore_locally_sparse_noise_4deg(self, data_set):
        
        stimname='locally_sparse_noise_4deg'
        noise_4deg_stim_table=self.data_set.get_stimulus_table(self.data_set.list_stimuli()[0])
        noise_templ=self.data_set.get_locally_sparse_noise_stimulus_template(stimname)
        lsn_movie, offscreen_mask =  data_set.get_locally_sparse_noise_stimulus_template('locally_sparse_noise_4deg',
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

    def saving_stuff(self, data_set):
        
        self.data_set.save_analysis_arrays()
        self.data_set.save_analysis_dataframes()
        # sio.savemat(os.path.join(self.main_directory,'corrected_traces.mat'),{'corrected_traces':corrected_traces,'dff_traces':dff_traces})   

   
   
        # movie_path = '/data/allen-brain-observatory/visual-coding-2p/ophys_movies/ophys_experiment_%d.h5' % exp['id']
        # f = h5py.File(movie_path,'r')
        # frames = f["data"][:10,:,:]
        # f.close()

    

    def mask_module(self, data_set, cell_ids, roi_mask_objects,all_roi_masks, max_projection):   
        cell_ids=cell_ids[0:2]
      
        f, axes = plt.subplots(1, len(cell_ids)+2, figsize=(15, 3))
        for ax, roi_mask, cid in zip(axes[:-2], roi_mask_objects, cell_ids):
            ax.imshow(roi_mask.get_mask_plane(), cmap='gray')
            ax.set_title('cell %d' % cid)
        combined_mask = all_roi_masks.max(axis=0)
        axes[-2].imshow(combined_mask, cmap='gray')
        axes[-2].set_title('all ROIs')
        axes[-1].imshow(max_projection, cmap='gray')
        axes[-1].set_title('max projection')
        plt.show()
        
    def dff_module(self,data_set, corrected_traces )  : 
        
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
 
      
    def plotting(self):
      
        # plot raw and corrected ROI trace
        pass
        # celltoplot=0
        # plt.figure(figsize=(14,4))
        # plt.title("Raw Fluorescence Trace")
        # plt.plot(time, raw_traces[celltoplot])
        # plt.show()
        
        # plt.figure(figsize=(14,4))
        # plt.title("Demixed Fluorescence Trace")
        # plt.plot(time, demixed_traces[celltoplot])
        # plt.show()
        
        # plt.figure(figsize=(14,4))
        # plt.title("Neuropil-corrected Fluorescence Trace")
        # plt.plot(time, corrected_traces[celltoplot])
        # plt.show()
        
        # plt.figure(figsize=(14,4))
        # plt.title("dF/F Trace")
        # # warning: dF/F can occasionally be one element longer or shorter 
        # # than the time stamps for the original traces.
        # plt.plot(time[:len(dff_traces[celltoplot])], dff_traces[celltoplot])
        # plt.show()
        
        # plt.figure(figsize=(14,4))
        # plt.title("dF/F Trace")
        # # warning: dF/F can occasionally be one element longer or shorter 
        # # than the time stamps for the original traces.
        # plt.plot(time[:len(self.deconvolved_spikes[celltoplot])], self.deconvolved_spikes[celltoplot])
        # plt.show()


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
                                si.NATURAL_MOVIE_PIXELS=(824 ,1030)
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
   
 
#%% yuriys methods  



# test1=self.mym_img_screen
# xv, yv = np.meshgrid(range(1024), range(1280))
# xv=xv-(1024/2)
# yv=yv-(1280/2)

# my_monitor=si.Monitor(1024,1280, 48.2,'cm')

# distance=15


# anglewidth = np.arctan(my_monitor.width / distance)
# angleheight = np.arctan(my_monitor.height / distance)


# xv=xv*angleheight/1024
# yv=yv*anglewidth/1280


# waqrpedx=(np.sin(xv)*distance)/np.cos(xv)
# waqrpedy=(np.sin(yv)*distance)/np.cos(yv)



# test10=[waqrpedx, waqrpedy]
# test11=np.asarray(test10)
# test13=[yv, xv]
# test14=np.asarray(test13)

# test12=test11.reshape(2,test11.shape[1]*test11.shape[2])

# plt.plot(test12[0,:])
# plt.plot(test12[1,:])


# xvec=waqrpedx.T.flatten()
# yvec=waqrpedy.flatten()

# coord=[xvec, yvec]
# coord2=np.asarray(coord)
# cordre=coord2.reshape(2,1024,1280)
# cordre=coord2[1,:].reshape(1024,1280)


# plt.imshow(cordre[:,:])

# mym_img_screen_warped=spndi.map_coordinates(test1, test13).reshape(1024,1280)
# plt.imshow(test1)

# plt.imshow(mym_img_screen_warped)



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the average head position over the patient cohort
Created on Mon Mar 18 13:55:07 2024

@author: verhei, inspired from code by Mikkel Vinding
"""

import os.path as op
import mne
import numpy as np
from mne.transforms import rotation_angles, rotation3d, write_trans, quat_to_rot
from utils import get_subjects, processed_data_dir, get_metadata, get_subjects_ID


# get all the patiets from metadata
metadata = get_metadata()
patients = list(metadata[metadata['group']=='patient'].index)

# Get data directories
task= 'rest_ec' #'empty_room' 
dirs, names = get_subjects(task=task, maxF=True)

name_id = [f_name.replace("NatMEG_", "") for f_name in names] 

#initialize the trans matrix array
trans_matrices = []
rot = []
for subj_dir, name in zip(dirs, name_id):
    
    # determine the avg trans for only patient cohort
    if name in patients:
        
        try:
            raw = mne.io.read_raw_fif(subj_dir,preload=True)
            subj_avg_trans = raw.info['dev_head_t']['trans']
            subj_rot = subj_avg_trans[0:3, 0:3] #take also the rotation part separately
            
            trans_matrices += [subj_avg_trans]
            rot += [subj_rot]
            
        except:
            print(f'Not finding maxfiltered data file from {name}')
            print('Moving on to the next participant.')

                

mean_trans = np.mean(np.stack(trans_matrices), axis=0)

rot_angles = [rotation_angles(r) for r in rot]
mean_rot_trans = rotation3d(*tuple(np.mean(np.stack(rot_angles), axis=0)))  # stack, then average, then make new coords

assert(np.sum(mean_trans[-1]) == 1.0)                   # sanity check result
mean_cohort_trans = raw.info['dev_head_t']              # use the last info as a template
mean_cohort_trans['trans'] = mean_trans                 # replace the transformation
mean_cohort_trans['trans'][:3, :3] = mean_rot_trans     # replace the rotation part
    

mean_trans_file = f'patient_average_{task}-trans.fif'
if op.isfile(mean_trans_file):
    print(f'NOTE: {mean_trans_file} already exists!')
      
else:    
    # Write trans file        
    write_trans(mean_trans_file, mean_cohort_trans)
    print("Wrote "+mean_trans_file)



    
    
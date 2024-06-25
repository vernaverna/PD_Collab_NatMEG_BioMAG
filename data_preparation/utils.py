#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains useful config variables and functions.

Created on Tue Oct 10 15:38:12 2023

@author: verhei
"""

import os
import mne
import numpy as np
import pandas as pd

from sklearn.decomposition import FactorAnalysis
from sklearn.preprocessing import StandardScaler, RobustScaler


# Define directory paths
raw_data_dir = '/archive/20079_parkinsons_longitudinal/MEG/'
processed_data_dir = '/home/verhei/PARKKI/derivatives/'
fs_subjects_dir = '/archive/20079_parkinsons_longitudinal/fs_subjects_dir/'
fsaverage_dir = '/home/verhei/PARKKI/fs_subjects_dir/'
fname_fsaverage_src = fsaverage_dir + '/fsaverage/bem/fsaverage-ico-5-src.fif'

# for creating an fsaverage, the one in analysis was empty/broken 
#mne.datasets.fetch_fsaverage(subjects_dir=fsaverage_dir)

# Read f-bands
freqs_filename = 'f.txt' #should be stored at your wk_dir
log_bands = np.loadtxt(freqs_filename, delimiter = ",") # Load band limits


def get_subjects_ID(): 
    """
    Helper function to get compatible subject IDs


    Returns
    -------
    subs : list of str
        subject fake IDs without prefix

    """
    
    subs = os.listdir(processed_data_dir)
    subs = [sdir.replace("NatMEG_", "") for sdir in subs] 
    
    return subs


def create_labels(subject, parc='aparc.a2009s'):
    """
    
    Parameters
    ----------
    subject : str
        Subject (free-surfer) name. 
    parc : str, optional
        The parcellation to use. The default is 'aparc_sub'.

    Returns
    -------
    
    lables : a list of labels, or None if freesurfer outputs do not exist

    """
    # create the label directory if needed
    labels_dir = os.path.join(processed_data_dir, 'NatMEG_'+subject, f'labels_{parc}')
    
    if not os.path.exists(labels_dir):
        try:
            # sorts the labels by alphabetical order by default
            labels=mne.read_labels_from_annot(subject=subject, parc=parc, subjects_dir=fs_subjects_dir)
            os.makedirs(labels_dir, exist_ok=True)
            
            for i in range(len(labels)):
                labels[i].save(os.path.join(labels_dir, labels[i].name))
        except:
            print("No MRI labels for {}".format(subject))
            return None

    # else just read the labels
    else:
        label_files = sorted([f.path for f in os.scandir(labels_dir)])
        labels = []
        
        for file in label_files:
            labels.append(mne.read_label(file, ''))
        
    return labels



def make_layout_coords(ch_type='all'):
    """
    Helper function to get layout coordinates for plotting purposes

    Parameters
    ----------
    ch_type  : str
               Can be either all, mag, grad or grad_norm (?)

    Returns
    -------
    coords : np Array
            2D numpy array with x, y coords

    """
    layout_name = f'Vectorview-{ch_type}'   
    layout = mne.channels.read_layout(layout_name)
    
    coords = layout.pos[:,0:2]
    
    return coords
    


def get_subjects(task='rest_ec', maxF=True):
    
    """
    Helper function to read out subjectr directories from data archive.
    
    Parameters
    ----------
    task: str. Which task to use (rest_ec, rest_eo, empty_room_before)
    maxF: boolean. Use Maxfiltered data?

    Returns
    -------
    subjects_data: a list of subject data paths
    subject_names: a list of <str> subject names 

    """
    
    #discard _0005, _0245, _0521 (no emptyroom data)
    to_discard = ['NatMEG_0005','NatMEG_0245', 'NatMEG_0521']
    
    dirs = os.listdir(raw_data_dir)
    subjects_data = []
    subject_names = []
    for s_dir in dirs:
        if (s_dir.startswith('NatMEG') and (s_dir not in to_discard)):
            subject_names.append(s_dir)
            
            # get into the subject data (stupid way)
            sub_path = os.path.join(raw_data_dir, s_dir)
            datedir = os.listdir(sub_path)[0]
            
            sub_data_dir = os.path.join(sub_path, datedir)
            
            if 'meg' in os.listdir(sub_data_dir):
                sub_data_dir = os.path.join(sub_data_dir, 'meg')
                task='rs_ec'
            
            if maxF:
                if task != 'empty_room':
                    fname=f'{task}_mc_avgtrans_tsss_corr95.fif'
                else:
                    fname=f'{task}_before_tsss_corr95.fif'
            else:
                fname = f'{task}.fif'
        
            subjects_data.append(os.path.join(sub_data_dir, fname))
    
    #TODO: think if I want to remove the NatMEG prefix?
    
    return subjects_data, subject_names 
        
        
# def get_subject_quat_files(subject):
#     """
#     Get the quarternion files for a single subject
    
    
#     Parameters
#     ----------
#     subject : str. which subject are we


#     Returns
#     -------
#     sub_trans_dir: subject data path

#     """
      
    
#     trans_dir=f'/archive/20079_parkinsons_longitudinal/analysis/meg_data/{sub}/'
#     sub_trans_f = trans_dir + sub + '-trans.fif'
    
#     return sub_trans_f
        
        
def get_subject_trans(subject):
    """
    Get the trans file for a single subject
    
    
    Parameters
    ----------
    subject : str. which subject are we reading the data from


    Returns
    -------
    sub_trans_dir: subject data path

    """
    sub = subject.replace("NatMEG_", "")   
    
    trans_dir=f'/archive/20079_parkinsons_longitudinal/analysis/meg_data/{sub}/'
    sub_trans_f = trans_dir + sub + '-trans.fif'
    
    return sub_trans_f
    


def get_metadata(drop_cols=True, binning=True):
    """
    Parameters
    ----------
    drop_cols : boolean
        Drop unnecessary? columns from csv file. Default is true.
    binning : boolean
        Creating categorigal bins from some of the cont. variables? 
    
    Returns
    -------
    metadata_df : pd.DataFrame 
        Metadata from subjects in an neat format

    """
    
    meta_path = raw_data_dir.replace('MEG', 'analysis')
    metadata = meta_path + 'alldata_subj2.csv'
    
    if drop_cols:
        meta_df = pd.read_csv(metadata, usecols=[1, 8, 9, 10, 11, 12, 13, 14,
                                                 35, 36, 37, 38, 39, 40, 41,
                                                 42, 43, 45])
    else:
        meta_df = pd.read_csv(metadata)
        
    if binning:
        #age_bins = [40, 50, 55, 60, 65, 70, 75, 90]
        #age_labs = ['<50', '50-54', '55-60', '61-65', '66-70', '71-75', '>75']
        #meta_df['age groups']= pd.cut(meta_df['age'], bins=age_bins, labels=age_labs)
        meta_df['age groups']= pd.qcut(meta_df['age'], q=5, precision=2)
        meta_df['MoCA scores'] = pd.qcut(meta_df['MoCA'], q=4, precision=2)
    
    meta_df['subj'] = ['0'+str(sub) for sub in meta_df['subj'] ]
    meta_df = meta_df.set_index('subj')
    
    return meta_df


def compute_latent_cognitive_score():
    """
    Computes a latent 'cognitive score' variable based on cognitive tests, 
    using unidimensional Factor analysis.
    

    Parameters
    ----------
    metadata_df :pd.DataFrame 
        Data frame containing the cognitive scores 

    Returns
    -------
    
    updated_df : an updated DataFrame with a latent 'cognitive score' 
    """
    
    metadata_df = get_metadata()

    X_val = metadata_df[['FAB', 'MoCA', 'BDI', 'UPDRS']]
    X_val = StandardScaler().fit_transform(X_val)
    
    fa = FactorAnalysis(n_components=1, rotation='varimax')
    nan_rows = np.isnan(X_val).any(axis=1)
    nan_idxs = np.where(nan_rows)[0]
    non_nan_idxs = np.where(nan_rows==False)[0]
    
    if len(nan_idxs) > 0:
        X_v = X_val[~nan_rows,:] #remove possible nan-rows
    else:
        X_v = X_val
    cog_score = fa.fit_transform(X_v) #get summary score
    
    cog_vec = np.zeros(shape=(len(X_val)))
    
    if len(nan_idxs) > 0:
        for i in range(len(non_nan_idxs)):
            cog_vec[non_nan_idxs[i]] = cog_score[i]
        for j in nan_idxs:
            cog_vec[j] = None
            
    
    metadata_df['cog_F'] = cog_vec

    
    updated_df = metadata_df
    
    return updated_df



def get_psd_instance(subject, task='rest_ec', fmin=1, fmax=90):
    """
    Helper function to calculate sensor-level psd per subject.
    

    Parameters
    ----------
    subject : str
        The subject filename to process
    task : TYPE, optional
        Which task to use. The default is 'rest_ec'.
    fmin : TYPE, optional
        Lower freqbound. The default is 1.
    fmax : TYPE, optional
        Higher freqbound. The default is 98.

    Returns
    -------
    
    psd : Array conatining the psd
    freqs : Array containing the frequency vector
    ch_names : a list of channel names

    """
    raw_fname = os.path.join(processed_data_dir, subject, task+'.fif')
    raw = mne.io.read_raw_fif(raw_fname, preload=True)
    
    raw0=raw.copy().filter(l_freq=fmin, h_freq=fmax)
    raw0.pick_types(meg=True)
    psd_obj = raw0.compute_psd(method='welch', window='hann', n_fft=2048, fmax=fmax)
    
    return psd_obj

def get_grad_norm_chs():
    layout = mne.channels.read_layout('Vectorview-grad_norm')
    ch_names = layout.names
    ch_names = [ch.replace(' ', '') for ch in ch_names]
    
    return ch_names

def get_roi_indices(roi, ch_names, ch_type='all'):
    """
    

    Parameters
    ----------
    roi : str
        Name of the region of interest. Currently supported ones are
    ch_names : list of str
        list of channel names
    ch_type : str, optional
        Use gradiometers or magnetometer, or 'grad_norm'. The default is 'all'.

    Returns
    -------
    None.

    """
    if roi=='sensorimotor':    
        # define ROI
        ROI_g = ['MEG0213', 'MEG0212', 'MEG0222', 'MEG0223', 'MEG0413', 'MEG0412', 'MEG0422', 'MEG0423', 'MEG0633', 'MEG0632',
                 'MEG0243', 'MEG0242', 'MEG0232', 'MEG0233', 'MEG0443', 'MEG0442', 'MEG0432', 'MEG0433', 'MEG0713', 'MEG0712', 
                 'MEG1613', 'MEG1612', 'MEG1622', 'MEG1623', 'MEG1813', 'MEG1812', 'MEG1822', 'MEG1823', 'MEG0743', 'MEG0742',
                 'MEG1043', 'MEG1042', 'MEG1112', 'MEG1113', 'MEG1123', 'MEG1122', 'MEG1312', 'MEG1313', 'MEG1323', 'MEG1322',
                 'MEG0723', 'MEG0722', 'MEG1142', 'MEG1143', 'MEG1133', 'MEG1132', 'MEG1342', 'MEG1343', 'MEG1333', 'MEG1332', 
                 'MEG0733', 'MEG0732', 'MEG2212', 'MEG2213', 'MEG2223', 'MEG2222', 'MEG2412', 'MEG2413', 'MEG2423', 'MEG2422']
        
        ROI_m = ['MEG0211', 'MEG0221', 'MEG0411', 'MEG0421', 'MEG0631', 
                 'MEG0241', 'MEG0231', 'MEG0441', 'MEG0431', 'MEG0711',
                 'MEG1611', 'MEG1621', 'MEG1811', 'MEG1821', 'MEG0741',
                 'MEG1041', 'MEG1111', 'MEG1121', 'MEG1311', 'MEG1321', 
                 'MEG0721', 'MEG1141', 'MEG1131', 'MEG1341', 'MEG1331',
                 'MEG0731', 'MEG2211', 'MEG2221', 'MEG2411', 'MEG2421']
        
    
    elif roi=='occipital': 
        ROI_g = ['MEG1942', 'MEG1943', 'MEG1922', 'MEG1923', 'MEG2112', 
                 'MEG2113', 'MEG2342', 'MEG2343', 'MEG2322', 'MEG2323',
                 'MEG1733', 'MEG1732', 'MEG1933', 'MEG1932', 'MEG2123', 
                 'MEG2122', 'MEG2333', 'MEG2332', 'MEG2513', 'MEG2512', 
                 'MEG1742', 'MEG1743', 'MEG2142', 'MEG2143', 'MEG2132', 
                 'MEG2133', 'MEG2542', 'MEG2543', 'MEG2532', 'MEG2533']

        ROI_m = ['MEG1941', 'MEG1921', 'MEG2111', 'MEG2341', 'MEG2321',
                 'MEG1731', 'MEG1931', 'MEG2121', 'MEG2331', 'MEG2511',  
                 'MEG1741', 'MEG2141', 'MEG2131', 'MEG2541', 'MEG2531']
            
    elif roi=='temporal':
        ROI_g = ['MEG0112', 'MEG0113', 'MEG0132', 'MEG0133', 
                 'MEG0142', 'MEG0143', 'MEG0212', 'MEG0213', 'MEG0222',
                 'MEG0223', 'MEG0232', 'MEG0233', 'MEG0242', 'MEG0243',
                 'MEG1312', 'MEG1313', 'MEG1322', 'MEG1323', 
                 'MEG1332', 'MEG1333', 'MEG1342', 'MEG1343', 'MEG1422',
                 'MEG1423', 'MEG1432', 'MEG1433', 'MEG1442', 'MEG1443',
                 'MEG1512', 'MEG1513', 'MEG1522', 'MEG1523', 
                 'MEG1532', 'MEG1533', 'MEG1542', 'MEG1543', 'MEG1612',
                 'MEG1613', 'MEG1622', 'MEG1623', 'MEG2412', 'MEG2413',
                 'MEG2422', 'MEG2423', 'MEG2612', 'MEG2613',
                 'MEG2622', 'MEG2623', 'MEG2632', 'MEG2633', 'MEG2642', 'MEG2643']
        
        ROI_m = ['MEG0111', 'MEG0131', 'MEG0141', 'MEG0211', 'MEG0221', 'MEG0231', 
                 'MEG0241', 'MEG1311', 'MEG1321', 'MEG1331', 'MEG1341', 'MEG1421', 
                 'MEG1431', 'MEG1441', 'MEG1511', 'MEG1521', 'MEG1531', 'MEG1541', 
                 'MEG1611', 'MEG1621', 'MEG2411', 'MEG2421', 'MEG2611', 'MEG2621',
                 'MEG2631', 'MEG2641']
        
    elif roi=='parietal':
        ROI_g = ['MEG0412', 'MEG0413', 'MEG0422', 'MEG0423', 'MEG0432', 'MEG0433', 'MEG0442', 
                 'MEG0443', 'MEG0632', 'MEG0633', 'MEG0712', 'MEG0713', 'MEG0722', 'MEG0723',
                 'MEG0732', 'MEG0733', 'MEG0742', 'MEG0743', 'MEG1042', 'MEG1043', 'MEG1112', 
                 'MEG1113', 'MEG1122', 'MEG1123', 'MEG1132', 'MEG1133', 'MEG1142', 'MEG1143',
                 'MEG1632', 'MEG1633', 'MEG1812', 'MEG1813', 'MEG1822', 'MEG1823', 'MEG1832', 
                 'MEG1833', 'MEG1842', 'MEG1843', 'MEG2012', 'MEG2013', 'MEG2022', 'MEG2023',
                 'MEG2212', 'MEG2213', 'MEG2222', 'MEG2223', 'MEG2232', 'MEG2233', 'MEG2242', 
                 'MEG2243', 'MEG2442', 'MEG2443']
        
        ROI_m = ['MEG0411', 'MEG0421', 'MEG0431', 'MEG0441', 'MEG0631', 'MEG0711', 
                 'MEG0721', 'MEG0731', 'MEG0741', 'MEG1041', 'MEG1111', 'MEG1121', 
                 'MEG1131', 'MEG1141', 'MEG1631', 'MEG1811', 'MEG1821', 'MEG1831', 
                 'MEG1841', 'MEG2011', 'MEG2021', 'MEG2211', 'MEG2221', 'MEG2231',
                 'MEG2241', 'MEG2441']

    elif roi=='frontal':
         ROI_g = ['MEG0122', 'MEG0123', 'MEG0342', 'MEG0343', 'MEG0322', 'MEG0323', 
                  'MEG0642', 'MEG0643', 'MEG0622', 'MEG0623', 'MEG1032', 'MEG1033',
                  'MEG1242', 'MEG1243', 'MEG1232', 'MEG1233', 'MEG1222', 'MEG1223', 
                  'MEG1412', 'MEG1413', 'MEG0312', 'MEG0313', 'MEG0542', 'MEG0543',
                  'MEG0612', 'MEG0613', 'MEG1012', 'MEG1013', 'MEG1022', 'MEG1023', 
                  'MEG0932', 'MEG0933', 'MEG1212', 'MEG1213', 'MEG0512', 'MEG0513',
                  'MEG0532', 'MEG0533', 'MEG0822', 'MEG0823', 'MEG0942', 'MEG0943', 
                  'MEG0922', 'MEG0923', 'MEG0522', 'MEG0523', 'MEG0812', 'MEG0813',
                  'MEG0912', 'MEG0913']
         
         ROI_m = ['MEG0121', 'MEG0341', 'MEG0321', 'MEG0641', 'MEG0621', 'MEG1031', 
                  'MEG1241', 'MEG1231', 'MEG1221', 'MEG1411', 'MEG0311', 'MEG0541',
                  'MEG0611', 'MEG1011', 'MEG1021', 'MEG0931', 'MEG1211', 'MEG0511', 
                  'MEG0531', 'MEG0821', 'MEG0941', 'MEG0921', 'MEG0521', 'MEG0811',
                  'MEG0911']

    elif roi=='all':
        ROI = ch_names
    else:
        msg = 'Enter a valid ROI!'
        return msg
    
    if ch_type=='grad':
        ROI=ROI_g   
    elif ch_type=='mag':
        ROI=ROI_m
    elif ch_type=='grad_norm':
        ROI = [ROI_g[i][:-1]+'X' for i in range(len(ROI_g))]
        ROI=list(dict.fromkeys(ROI))
    elif ch_type=='all':
        ROI=ROI_g+ROI_m
    
    roi_ch_idx = [ch_names.index(ch) for ch in ROI] 
    roi_indexes = roi_ch_idx
    
    return roi_indexes
    
    
    
    
    
    
    
    
    
    

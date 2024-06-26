#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example of maxfiltering with python
Created on Tue Jun 18 10:36:14 2024

@author: heikkiv
"""

import os
import subprocess 

# path to bidsified data
#root="/net/theta/fishpool/projects/dbs_pd/neuro-data/dbs_pd/BIDS/"
root="/data/MaxFilter-test-data/"

sub="35P"
ses="01"
run="02"
task="sp"
trans="default" #, "ses1", "avg" or "init"


#input_file = f"{root}sub-{sub}/ses-{ses}/meg/sub-{sub}_ses-{ses}_task-{task}_run-{run}_proc-raw_meg.fif"
#output_file = f"{root}derivatives/sub-{sub}/ses-{ses}/meg/sub-{sub}_ses-{ses}_task-{task}_run-{run}_proc-tsss+mc_meg.fif"
#trans_output_file = f"{root}derivatives/sub-{sub}/ses-{ses}/meg/sub-{sub}_ses-{ses}_task-{task}_run-{run}_proc-tsss+mc_trans-{trans}_meg.fif"
#log_file = f"{root}derivatives/sub-{sub}/ses-{ses}/meg/sub-{sub}_ses-{ses}_task-{task}_run-{run}_proc-tsss+mc_trans-{trans}_meg_log.log"

input_file = "/data/MaxFilter-test-data/RestEyesOpen.fif"
output_file = "/data/MaxFilter-test-data/RestEyesOpen_tsss+mc.fif"
trans_output_file = "/data/MaxFilter-test-data/RestEyesOpen_trans.fif"
log_file = "/data/MaxFilter-test-data/RestEyesOpen_log.log"

#psecif the transformation file / option here 
trans_file= "default" #"/home/heikkiv/parkinsonsdisease/python/verna/KI_patient_average_rest_ec-trans.fif"


# Fit SSS sphere origin to digitization points using fit_sphere_to_isotrak function
#origin_args = ['/neuro/bin/util/fit_sphere_to_isotrak', '-f', input_file]
#output = subprocess.check_output(args=origin_args, text=True).strip()
#origin = output.split(' ')

if trans!="init": #by default, maxfilter transforms the data into intial head position
    args = ['/neuro/bin/util/maxfilter', '-f', input_file, '-o', output_file, \   #Change the path to NatMEG
            '-st', '-corr', '0.80', '-ds', '2', '-lpfilt', '100', \
            '-hpiwin', '1000','-hpistep', '100', '-move-comp', \
            '-v', ‘-frame’, ‘head, ’-origin', '0 0 55', \
            '-force'] #coordinate frame set to head
                                                        
else:                                                    
    args2= ['/neuro/bin/util/mfilter',  '-f', output_file, '-o',  trans_output_file,  '-trans', trans_file, '-autobad', 'off', '-autoflat','off' '-force']                                                        
                                                   
# save the log and errors
log_output = open(log_file, "w") 
        
# run maxfilter
subprocess.run(args=args, stdout=log_output,stderr=log_output)

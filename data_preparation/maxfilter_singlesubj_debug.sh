#!/bin/bash

root="/net/theta/fishpool/projects/dbs_pd/neuro-data/dbs_pd/BIDS/"

sub="35P"
ses="01"
run="02"

trans_f="${root}sub-${sub}/ses-${ses}/meg/sub-${sub}_ses-${ses}_task-sens_run-${run}_proc-raw_meg.fif"
# trans="default"

for task in sp; do
  path="${root}sub-${sub}/ses-${ses}/meg/sub-${sub}_ses-${ses}_task-${task}_run-${run}_proc-raw_meg.fif"
  mc_tsss_path="${root}derivatives/sub-${sub}/ses-${ses}/meg/sub-${sub}_ses-${ses}_task-${task}_run-${run}_proc-tsss+mc_meg.fif"
  nice /neuro/bin/util/mfilter -gui -f $path -o $mc_tsss_path -st -corr 0.80 -ds 2 -lpfilt 100 -movecomp | tee ${root}derivatives/sub-${sub}/ses-${ses}/meg/sub-${sub}_ses-${ses}_task-${task}_run-${run}_proc-tsss+mc_meg.log

done
       

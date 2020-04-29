#!/bin/bash

qsub -v "LAG1=0","LAG2=5","TYPE_OBS=mod","DOMAIN=GULFSTREAM2" submit_fig_xp_all.pbs
qsub -v "LAG1=0","LAG2=5","TYPE_OBS=obs","DOMAIN=GULFSTREAM2" submit_fig_xp_all.pbs


#!/bin/sh

sbatch mono_gpu_OSE.slurm 0 OSMOSIS False
sbatch mono_gpu_OSE.slurm 0 GULFSTREAM False


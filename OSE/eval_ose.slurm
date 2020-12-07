#!/bin/bash
#SBATCH -A yrf@cpu                  # nom du compte
#SBATCH --job-name=cpu_eval_ose     # nom du job
#SBATCH --ntasks=1                  # nombre de tâche (un unique processus ici)
#SBATCH --hint=nomultithread        # on réserve des coeurs physiques et non logiques
#SBATCH --time=00:05:00             # temps exécution maximum demande (HH:MM:SS)
#SBATCH --output=eval_ose%j.out     # nom du fichier de sortie
#SBATCH --error=eval_ose%j.err      # nom du fichier erreur (ici commun avec la sortie)
 
# nettoyage des modules charges en interactif et hérités par défaut
module purge
 
# chargement des modules
module load geos/3.7.3
module load tensorflow-gpu/py3/1.14-openmpi
module load imagemagick/7.0.8-7
module load ffmpeg

export PYTHONPATH=${HOME}/PB_ANDA:${PYTHONPATH}
export QT_QPA_PLATFORM='offscreen'

# echo des commandes lancées
set -x
 
# exécution du code
python -u export_NetCDF.py ${1} ${2} ${3}
python -u analysis_OSE.py ${1} ${2} ${3}

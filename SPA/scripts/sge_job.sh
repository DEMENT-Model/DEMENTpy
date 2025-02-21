#!/bin/bash

#$ -N dmg_data


#$ -q mic

##$ -l mem_free=37G,h_vmem=37G


#$ -j y
#$ -o std_oe/

#$ -m ea


module load anaconda

python3 drought_trait.py base_site site

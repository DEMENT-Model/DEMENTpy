#!/bin/bash
#$ -N dementjupyter
#$ -q mic
#$ -m beas

module load anaconda/3.7-5.3.0

jupyter nbconvert --to notebook --execute 2019110600.ipynb --output 2019110600.ipynb
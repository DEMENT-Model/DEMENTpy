#!/bin/bash
#$ -N dm_
#$ -q pub8i
#$ -m beas


module load anaconda

cd src

python dementpy.py input_ output 20191126


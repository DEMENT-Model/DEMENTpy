#!/bin/bas3
#$ -N bas3
#$ -q pub8i
#$ -m beas


module load anaconda

cd src

python dementpy.py mode1/input mode1/output 201911263

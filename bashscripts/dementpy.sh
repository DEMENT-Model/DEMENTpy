#!/bin/bash


#$ -N dm_
#$ -q mic
#$ -m beas


module load anaconda

cd src

python dementpy.py mode/input_ output 20200120


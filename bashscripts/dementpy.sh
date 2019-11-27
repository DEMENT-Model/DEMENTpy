#!/bin/bash
#$ -N dementpy
#$ -q mic
#$ -m beas


module load anaconda/3.7-5.3.0

cd src

python dementpy.py input1 output1 20191126


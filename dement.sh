#!/bin/bash
#$ -N dementpy
#$ -q pub8i
#$ -m beas


module load anaconda/3.7-5.3.0
python dement.py runtime.txt

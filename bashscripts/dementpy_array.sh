#!/bin/bash

#==========================
#======Set SGE options:
#==========================

#$ -cwd
#$ -S /bin/bash

#$ -N dementpy


#==============================================
#===Queue Resource Request
# -pe smp #:
#          shared memory parallel environment
# -R y:
#      The reserve option allows your job to get a foot in the door on a node,
#      and prevents the node from being constantly loaded by single-core jobs.
# mem_free: 
#          This should be set to be the what you think your job will need (or a little more).
#          This will reserve memory for your job on the node that it is run on.
# h_vmem:
#        This is the “high water mark” for memory for your job.  This should be set to be equal to,
#        your mem_free request.
#==============================================
#$ -q mic
#$ -pe smp 8
#$ -R y
#$ -l mem_free=3G,h_vmem=3G
#$ -l h_rt=30:00:00
#$ -l s_rt=35:00:00

########################
# Merge the standard error into the standard output stream
########################
#$ -j y
##$ -e output/
#$ -o output/

##$ -M wbwenwu@gmail.com
#$ -m ea

## submit an array job
#$ -t 1-8

#==========================
#======Job it self
#==========================

#module load anaconda
#cd src
#python dementpy.py mode/input_ output 20200120


date > output/results.$SGE_TASK_ID

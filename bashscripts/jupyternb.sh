#--------------------------------------------------------#
#----Script submitting jupyter notebook jobs to HPC------#
#--------------------------------------------------------#
#!/bin/bash
#$ -N dementjupyter
#$ -q mic
#$ -m beas

module load anaconda/3.7-5.3.0

jupyter nbconvert --to notebook --execute test_5.ipynb --output test_5.ipynb
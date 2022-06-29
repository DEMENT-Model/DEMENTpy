#--------------------------------------------------------#
#----Script submitting jupyter notebook jobs to HPC------#
#--------------------------------------------------------#

#!/bin/bash
#$ -N dementjupyter
#$ -q mic
#$ -m beas

module load anaconda/3.7-5.3.0

# --ExecutePreprocessor.timeout=180 important to the running of notebooks.
jupyter nbconvert --ExecutePreprocessor.timeout=180 --to notebook --execute test_5.ipynb --output test_5.ipynb
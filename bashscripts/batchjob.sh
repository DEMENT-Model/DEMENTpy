#---------------------------------------------#
#--- Submit jobs from various folers to HPC---#
#---------------------------------------------#

#!/bin/bash

cd dDEMENTpy
qsub dementpy.sh
cd ..

cd dDEMENTpy1
qsub dementpy.sh
cd ..

cd dDEMENTpy2
qsub dementpy.sh
cd ..

cd dDEMENTpy3
qsub dementpy.sh
cd ..
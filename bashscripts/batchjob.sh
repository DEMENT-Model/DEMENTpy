#---------------------------------------------#
#--- Submit jobs from various folers to HPC---#
#---------------------------------------------#

#!/bin/bash

cd dementpy0
qsub dementpy.sh
cd ..

cd dementpy1
qsub dementpy.sh
cd ..

cd dementpy2
qsub dementpy.sh
cd ..
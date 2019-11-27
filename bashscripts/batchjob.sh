#---------------------------------------------#
#--- Submit jobs from various folers to HPC---#
#---------------------------------------------#
# HOW to modify dementpy.sh before being qsubed??
# i.e. supply output file names from this script without need to change one by one.

# version_1

#!/bin/bash
current="20191117"
substitute="20191118"

cd dementpy0
sed -i -e "s/$current/$substitute/" dementpy.sh 
qsub dementpy.sh
cd ..

cd dementpy1
sed -i -e "s/$current/$substitute/" dementpy1.sh
qsub dementpy.sh
cd ..

cd dementpy2
sed -i -e "s/$current/$substitute/" dementpy2.sh
qsub dementpy.sh
cd ..



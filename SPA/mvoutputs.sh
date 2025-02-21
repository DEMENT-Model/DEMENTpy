#----------------------------------------------------------------------------#
#--- Script mving output files from various folers over to the same Folder---#
#----------------------------------------------------------------------------#

#!/bin/bash

cd dementpy/output
mv 2019111800.pickle ../../test_5
cd ../..

cd dementpy1/output
mv 2019111801.pickle ../../test_5
cd ../..

cd dementpy2/output
mv 2019111802.pickle ../../test_5
cd ../..

cd dementpy3/output
mv 2019111803.pickle ../../test_5
cd ../..

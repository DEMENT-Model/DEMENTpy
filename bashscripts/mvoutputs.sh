#----------------------------------------------------------------------------#
#--- Script mving output files from various folers over to the same Folder---#
#----------------------------------------------------------------------------#

#!/bin/bash

cd dDEMENTpy/output
mv 2019110600.pickle ../../test_5
cd ../..

cd dDEMENTpy1/output
mv 2019110601.pickle ../../test_5
cd ../..

cd dDEMENTpy2/output
mv 2019110602.pickle ../../test_5
cd ../..

cd dDEMENTpy3/output
mv 2019110603.pickle ../../test_5
cd ../..
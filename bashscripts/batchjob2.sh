#!/bin/bash

#--------------------------------------------------
# input1 loop
#--------------------------------------------------
for count in 1 2 3
do

  jobname="sev$count"
  sed -i -e "s/sev.*/$jobname/" dementpy.sh

  outname="20191126$count"
  #echo $outname
  sed -i -e "s/20191126.*/$outname/" dementpy.sh

  qsub dementpy.sh

done


#--------------------------------------------------
# input2 loop
#--------------------------------------------------

sed -i -e "s/input/input2/" dementpy.sh
sed -i -e "s/output/output2/" dementpy.sh

for ((count = 2; count = 0; count--))
do
outname = '20191126'+'$count'
sed -i -e "s/20191126/$outname/" dementpy.sh
qsub dementpy.sh
done

#--------------------------------------------------
# input3 loop
#--------------------------------------------------

sed -i -e "s/input/input3/" dementpy.sh
sed -i -e "s/output/output3/" dementpy.sh

for ((count = 2; count = 0; count--))
do
outname = '20191126'+'$count'
sed -i -e "s/20191126/$outname/" dementpy.sh
qsub dementpy.sh
done
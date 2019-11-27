#!/bin/bash

for count in 1 2 3
do
  # job naming
  jobname="sev$count"
  sed -i -e "s/sev.*/$jobname/" dementpy.sh
  
  # output naming
  outname="20191126$count"
  sed -i -e "s/20191126.*/$outname/" dementpy.sh

  qsub dementpy.sh

done
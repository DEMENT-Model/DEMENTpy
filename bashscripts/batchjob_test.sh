#!/bin/bash
for mode in 1 2 3
do
  for count in 1 2 3
  do
    # job naming
    jobname="m_$mode_c_$count"
    sed -i -e "s/m_.*/$jobname/" dementpy.sh

    # output naming
    outname="20191126$count"
    sed -i -e "s/20191126.*/$outname/" dementpy.sh
    
    #folder naming
    input="mode$mode/input$count"

    qsub dementpy.sh

  done
done
#!/bin/bash

# Make changes to dementpy.sh using for loops
for mode in 1
do
  for scenario in bas mid sev
  do
    #input output folder
    folder="mode$mode/input_$scenario output"
    sed -i -e "s%mode.*.output%$folder%" dementpy.sh

    for count in 1 2 3 4 5
    do

      jobname="dm_$mode$scenario$count"
      sed -i -e "s/dm_.*/$jobname/" dementpy.sh

      outname="20200120$count"_"$scenario"
      echo "job name:" $outname
      sed -i -e "s/20200120.*/$outname/" dementpy.sh

      qsub dementpy.sh

    done

  done

done

# restore the dementpy.sh to where it has begin
jobname="dm_"
sed -i -e "s/dm_.*/$jobname/" dementpy.sh

folder="mode/input_ output"
sed -i -e "s%mode.*.output%$folder%" dementpy.sh

outname="20200120"
sed -i -e "s/20200120.*/$outname/" dementpy.sh

#!/bin/bash


# iterate over origin sites
for origin in desert scrubland grassland pineoak subalpine
do

  base_site="$origin"

  sed -i -e "s/base_site/$base_site/" slurm_job.sub

  default="$base_site site"

  # iterate over sites
  for site in desert #scrubland pineoak subalpine grassland
  do
    # change site name
    target="$base_site $site"

    sed -i -e "s/$default/$target/" slurm_job.sub
  
    # submit to HPC
    sbatch slurm_job.sub
  
    # change it back
    sed -i -e "s/$target/$default/" slurm_job.sub

  done

  sed -i -e "s/$base_site/base_site/" slurm_job.sub

done

#!/bin/bash

FOLDER_NAME=off_beam

mkdir -p ${FOLDER_NAME}

while read run; do

    echo ""
    echo "Running script on run: $run"
  
 
    # enter the singularity environment
    singularity shell -B/pnfs:/pnfs,/exp/annie/app/users/doran/temp_directory:/tmp,/exp/annie/data:/exp/annie/data,/exp/annie/app:/exp/annie/app /cvmfs/singularity.opensciencegrid.org/anniesoft/toolanalysis\:latest/ << EOF
   
	root -l -q -b 'extract_data_off_beam.C(${run})'

    exit

EOF

ls -lrth ${FOLDER_NAME}/

done < "FY22_23.list"

echo ""
echo "done"
echo ""

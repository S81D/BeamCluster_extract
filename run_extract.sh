#!/bin/bash

while read run; do

    echo ""
    echo "Running script on run: $run"
  
 
    # enter the singularity environment
    singularity shell -B/pnfs:/pnfs,/exp/annie/app/users/doran/temp_directory:/tmp,/exp/annie/data:/exp/annie/data,/exp/annie/app:/exp/annie/app /cvmfs/singularity.opensciencegrid.org/anniesoft/toolanalysis\:latest/ << EOF
   
	root -l -q -b 'extract_data.C(${run})'

    exit

EOF

ls -lrth 2024_usable_data/

done < "files_2024_usable.list"

echo ""
echo "done"
echo ""

#!/bin/bash

while read run; do

    echo ""
    echo "Running script on run: $run"
  
 
    # enter the singularity environment
    singularity shell -B/pnfs:/pnfs,/exp/annie/app/users/doran/temp_directory:/tmp,/exp/annie/data:/exp/annie/data,/exp/annie/app:/exp/annie/app /cvmfs/singularity.opensciencegrid.org/anniesoft/toolanalysis\:latest/ << EOF
   
	root -l -q -b 'count_triggers.C(${run})'

    exit

EOF

done < "FY22_23.list"

echo ""
echo "done"
echo ""

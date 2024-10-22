#!/bin/bash

#rm -rf output_data.txt

while read run; do

    echo ""
	echo "####################################"
    echo "Running script on run: $run"
  
    # enter the singularity environment
    singularity shell -B/pnfs:/pnfs,/exp/annie/app/users/doran/temp_directory:/tmp,/exp/annie/data:/exp/annie/data,/exp/annie/app:/exp/annie/app /cvmfs/singularity.opensciencegrid.org/anniesoft/toolanalysis\:latest/ << EOF
   
	root -l -q -b 'beamrun_DQ.C(${run})'

    exit

EOF

done < "files_runs.list"

echo ""
echo "done"
echo ""

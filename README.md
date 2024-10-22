# BeamCluster_extract
Extract/Filter information from the PhaseIITreeMaker root trees.

### usage
Populate a ```file.list``` with the runs you wish to loop over. Running ```sh run_extract.sh``` will execute ```extract_data.C``` for each line in the list file, and will create + fill new root files with events that only have the desired cluster information.

#### `extract_data_from_BC/`
Two scripts very similar to the scripts described above. These are extraction scripts that populate an output .txt file with DQ metrics.

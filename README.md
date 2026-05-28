# BeamCluster_extract
Extract/Filter information from the PhaseIITreeMaker root trees. Given the BeamCluster TTrees contain all triggers and device information, they are quite hefty (oftentimes > 1 GB per file), there was a desire to compress/filter the TTrees to only keep good qualty beam information. These scripts exist to filter the beam information to keep only the physics analysis relevant information.

### usage
Populate a ```file.list``` with the runs you wish to loop over. Running ```sh run_extract.sh``` will execute ```extract_data.C``` for each line in the list file, and will create + fill new root files with events that only have the desired cluster information.


## NCQE Cross Section

As part of the NCQE cross section, data was filtered from the BeamCluster ntuples into a condensed format for ease of analysis. The scripts located in `NCQE_Cross_Section/` were used to do this:

- `extract_data.C` -- filter/pull out relevant beam information.
- `extract_data_off_beam.C` -- same thing but only pulls out OFF BEAM triggers, for an off-beam/CIT background assessment.
- `run_extract.sh` -- execute the extract script.
- `count_triggers.C` -- for off-beam subtraction, the total number of triggers must be known for both the beam and off-beam data. This script will count the total number of triggers to make this possible. The filter scripts only save the cluster information, hence it's separate.
- `count_trigs.sh` -- execute the count trigger script.

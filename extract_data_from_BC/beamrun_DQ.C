#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>


// function to calculate statistical errors
float CalculateStatError(float numerator, float denominator);

float CalculateStatError(float numerator, float denominator) {
    // Check if denominator or numerator is zero to avoid division by zero
    if (denominator == 0 || numerator == 0) {
        return 0;
    }

    // Calculate the error: sqrt( (dN/N)^2 + (dD/D)^2 ) * (N/D)
    float relative_error_numerator = sqrt(numerator) / numerator;
    float relative_error_denominator = sqrt(denominator) / denominator;
    float total_relative_error = sqrt((relative_error_numerator * relative_error_numerator) +
                                      (relative_error_denominator * relative_error_denominator));

    // The error in the ratio is: ratio * total_relative_error
    float ratio = 100 * numerator / denominator;
    return ratio * total_relative_error;
}


// function to calculate run counts and rates
void beamrun_DQ(int run_number) {

    // Construct the filename using the run number
    //TString directory = "/pnfs/annie/persistent/processed/BeamClusterTrees/";
    TString directory = "../../root_files/";
    TString fileName = TString::Format("BeamCluster_%d.root", run_number);

    // Open the ROOT file
    TFile *inputFile = TFile::Open(directory + fileName);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Failed to open file: " << fileName << std::endl;
        return;
    }

    // ---------------------------------------------------------

    // Input File initialization
    Int_t TankMRDCoinc, NoVeto, eventNumber, Extended, HasLAPPD;
    std::vector<int>* GroupedTriggerWord = 0;
    std::vector<int>* tracks = 0;
    std::vector<ULong64_t>* GroupedTriggerTime = nullptr;
    ULong64_t eventTimeTank;
    std::vector<double>* clusterTime = nullptr;
    Double_t BRF_fit, beam_downstream, beam_upstream, beam_horn;

            
    TTree *tree = (TTree*)inputFile->Get("Event");

    // Set up variables to read from input tree
    tree->SetBranchAddress("eventNumber", &eventNumber);
	tree->SetBranchAddress("GroupedTriggerWord", &GroupedTriggerWord);
    tree->SetBranchAddress("GroupedTriggerTime", &GroupedTriggerTime);
    tree->SetBranchAddress("beam_E_TOR860", &beam_downstream);
    tree->SetBranchAddress("beam_E_TOR875", &beam_upstream);
    tree->SetBranchAddress("beam_THCURR", &beam_horn);
    tree->SetBranchAddress("NumClusterTracks", &tracks);
    tree->SetBranchAddress("TankMRDCoinc", &TankMRDCoinc);
    tree->SetBranchAddress("NoVeto", &NoVeto);
    tree->SetBranchAddress("eventTimeTank", &eventTimeTank);
    tree->SetBranchAddress("Extended", &Extended);
    tree->SetBranchAddress("HasLAPPD", &HasLAPPD);
    tree->SetBranchAddress("clusterTime", &clusterTime);
    tree->SetBranchAddress("BRFFirstPeakFit", &BRF_fit);

    // ---------------------------------------------------------

    // initialize output print variables
    float_t totalclusters = 0;
    float_t totalevents = 0;
    float_t totalclusters_in_prompt = 0;
    float_t totalclusters_in_ext = 0;

    float_t totalext_rate_1 = 0;
    float_t totalext_rate_2 = 0;
    float_t totalokay_beam = 0;
    float_t totaltmrd_coinc = 0;
    float_t totalveto_hit = 0;
    float_t totalveto_tmrd_coinc = 0;
    float_t totalhas_track = 0;
    float_t totalhas_lappd = 0;
    float_t totalhas_BRF = 0;
    float_t totaltimezero = 0;
    float_t totalclusters_in_spill = 0;

    // ---------------------------------------------------------

    // Loop over events in the tree and get total counts
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        // total beam triggers (UBT only)
        int trigger_index = -1; // Default to -1 if not found
        if (GroupedTriggerWord && !GroupedTriggerWord->empty()) {
            for (size_t k = 0; k < GroupedTriggerWord->size(); k++) {
                if ((*GroupedTriggerWord)[k] == 14) {
                    trigger_index = static_cast<int>(k);
                    break;
                }
            }
        }
        if (trigger_index == -1) {
            continue;
        }

        totalevents++;

        // events with "good" conditions
        // (0.5e12 < pot < 8e12), (172 < horn val < 176), (downstream/upstream within 5%)
        double beam_loss_frac = std::abs(beam_downstream - beam_upstream) / (beam_downstream + beam_upstream);
        if ((beam_downstream > 0.5 && beam_downstream < 8) || (beam_horn > 172 && beam_horn < 176) || (beam_loss_frac < 0.05)) {

            totalokay_beam++;     // all subsequent rates are based on the number of "good" beam triggers

            // events with eventTimeTank == 0 (weird issue where ADV waveforms don't get timestamps)
            if (eventTimeTank == 0){
                totaltimezero++;
            }

            // Extended readout rates (1 = charge-based, 2 = forced, random)
            if (Extended == 1){
                totalext_rate_1++;
            }
            if (Extended == 2){
                totalext_rate_2++;
            }

            // find total number of clusters, and how many of the prompt clusters are in the beam spill
            for (size_t j = 0; j < clusterTime->size(); j++) {                      // loop over all clusters
                totalclusters++;
                if (clusterTime->at(j) < 2000) {                                    // cluster is in prompt window
                    totalclusters_in_prompt++;
                    if (clusterTime->at(j) > 190 && clusterTime->at(j) < 1750) {    // based on expected beam spill structure from data 
                        totalclusters_in_spill++;
                    }  
                }
                else {
                    totalclusters_in_ext++;                                         // clusters in the ext window
                }
            }

			// Does the event have LAPPDs?
            if (HasLAPPD == 1) {
                totalhas_lappd++;
            }

            // coincidences
            if (NoVeto == 0) {
                totalveto_hit++;              // opposite logic to NoVeto
            }
            if (TankMRDCoinc == 1) {
                totaltmrd_coinc++;            // Tank + MRD
                if (NoVeto == 0) {
                    totalveto_tmrd_coinc++;   // Veto + Tank + MRD
                }
                // find MRD tracks
                // Search through the MRD cluster array to see if any of them had a single track. If thats the case (and there was TankMRDCoinc), we can say there was a track
                if (std::find(tracks->begin(), tracks->end(), 1) != tracks->end()) {
                    totalhas_track++;
                } 
            }

            // Does the event have a usable BRF signal (+ fit)
            if (BRF_fit != 0) {
                totalhas_BRF++;
            } 

        }
    }

    // ---------------------------------------------------------

    // Calculate rates to output

    // all % are going to be taken out of total events with an UBT

    float_t events_per_cluster = totalclusters / totalevents;
    float_t ext_rate_1 = totalext_rate_1 * 100 / totalevents;
    float_t ext_rate_2 = totalext_rate_2 * 100 / totalevents;
    float_t okay_beam = totalokay_beam * 100 / totalevents;
    float_t tmrd_coinc = totaltmrd_coinc * 100 / totalevents;
    float_t veto_hit = totalveto_hit * 100 / totalevents;
    float_t veto_tmrd_coinc = totalveto_tmrd_coinc * 100 / totalevents;
    float_t has_track = totalhas_track * 100 / totalevents;
    float_t has_lappd = totalhas_lappd * 100 / totalevents;
    float_t has_BRF = totalhas_BRF * 100 / totalevents;
    float_t timezero = totaltimezero * 100 / totalokay_beam;
    float_t clusters_in_spill = totalclusters_in_spill * 100 / totalclusters;   // except these one of course
    float_t clusters_in_prompt = totalclusters_in_prompt * 100 / totalclusters;
    float_t clusters_in_ext = totalclusters_in_ext * 100 / totalclusters;


    // calculate statistical errors: dA = sqrt( (dB/B)^2  + (dC/C)^2 ) * A
    float_t er_events_per_cluster = CalculateStatError(totalclusters, totalevents);
    float_t er_ext_rate_1 = CalculateStatError(totalext_rate_1, totalevents);
    float_t er_ext_rate_2 = CalculateStatError(totalext_rate_2, totalevents);
    float_t er_okay_beam = CalculateStatError(totalokay_beam, totalevents);
    float_t er_tmrd_coinc = CalculateStatError(totaltmrd_coinc, totalevents);
    float_t er_veto_hit = CalculateStatError(totalveto_hit, totalevents);
    float_t er_veto_tmrd_coinc = CalculateStatError(totalveto_tmrd_coinc, totalevents);
    float_t er_has_track = CalculateStatError(totalhas_track, totalevents);
    float_t er_has_lappd = CalculateStatError(totalhas_lappd, totalevents);
    float_t er_has_BRF = CalculateStatError(totalhas_BRF, totalevents);
    float_t er_timezero = CalculateStatError(totaltimezero, totalevents);
    float_t er_clusters_in_spill = CalculateStatError(totalclusters_in_spill, totalclusters);
    float_t er_clusters_in_prompt = CalculateStatError(totalclusters_in_prompt, totalclusters);
    float_t er_clusters_in_ext = CalculateStatError(totalclusters_in_ext, totalclusters);

    // ---------------------------------------------------------

    // Append results to output files

    int t_s = 5;    // text spacing of output table

    // counts + rates (+ errors)
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    std::ofstream outfile_rates("run_rates.txt", std::ios::app);  // Open in append mode

    // Check if the file opened successfully
    if (!outfile_rates) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // Optionally write headers if needed:
    if (outfile_rates.tellp() == 0) { 

        outfile_rates << std::setw(t_s) << "run_number "
                << std::setw(t_s) << "total_events "
                << std::setw(t_s) << "hasLAPPD "
                << std::setw(t_s) << "hasLAPPD_rate "
                << std::setw(t_s) << "hasLAPPD_error "
                << std::setw(t_s) << "hasBRF "
                << std::setw(t_s) << "hasBRF_rate "
                << std::setw(t_s) << "hasBRF_error "
                << std::setw(t_s) << "eventTimeZero "
                << std::setw(t_s) << "eventTimeZero_rate "
                << std::setw(t_s) << "eventTimeZero_error "
                << std::setw(t_s) << "beam_ok "
                << std::setw(t_s) << "beam_ok_rate "
                << std::setw(t_s) << "beam_ok_error "
                << std::setw(t_s) << "total_clusters "
                << std::setw(t_s) << "clusters_per_event "
                << std::setw(t_s) << "clusters_per_event_error "
                << std::setw(t_s) << "prompt_clusters "
                << std::setw(t_s) << "prompt_clusters_rate "
                << std::setw(t_s) << "prompt_clusters_error "
                << std::setw(t_s) << "spill_clusters "
                << std::setw(t_s) << "spill_clusters_rate "
                << std::setw(t_s) << "spill_clusters_error "
                << std::setw(t_s) << "ext_clusters "
                << std::setw(t_s) << "ext_clusters_rate "
                << std::setw(t_s) << "ext_clusters_error "
                << std::setw(t_s) << "ext_1 "
                << std::setw(t_s) << "ext_1_rate "
                << std::setw(t_s) << "ext_1_error "
                << std::setw(t_s) << "ext_2 "
                << std::setw(t_s) << "ext_2_rate "
                << std::setw(t_s) << "ext_2_error "
                << std::setw(t_s) << "TankMRDCoinc "
                << std::setw(t_s) << "TankMRDCoinc_rate "
                << std::setw(t_s) << "TankMRDCoinc_error "
                << std::setw(t_s) << "veto "
                << std::setw(t_s) << "veto_rate "
                << std::setw(t_s) << "veto_error "
                << std::setw(t_s) << "vetoTankMRDCoinc "
                << std::setw(t_s) << "vetoTankMRDCoinc_rate "
                << std::setw(t_s) << "vetoTankMRDCoinc_error "
                << std::setw(t_s) << "MRDTrack "
                << std::setw(t_s) << "MRDTrack_rate "
                << std::setw(t_s) << "MRDTrack_error "
                << std::endl;

    }

    outfile_rates << std::fixed << std::setprecision(5);    // printed output precision

    // Write data for the current run
    outfile_rates << std::setw(t_s) << run_number << " "
              << std::setw(t_s) << static_cast<int>(totalevents) << " "
              << std::setw(t_s) << static_cast<int>(totalhas_lappd) << " "
              << std::setw(t_s) << has_lappd << " "
              << std::setw(t_s) << er_has_lappd << " "
              << std::setw(t_s) << static_cast<int>(totalhas_BRF) << " "
              << std::setw(t_s) << has_BRF << " "
              << std::setw(t_s) << er_has_BRF << " "
              << std::setw(t_s) << static_cast<int>(totaltimezero) << " "
              << std::setw(t_s) << timezero << " "
              << std::setw(t_s) << er_timezero << " "
              << std::setw(t_s) << static_cast<int>(totalokay_beam) << " "
              << std::setw(t_s) << okay_beam << " "
              << std::setw(t_s) << er_okay_beam << " "
              << std::setw(t_s) << static_cast<int>(totalclusters) << " "
              << std::setw(t_s) << events_per_cluster << " "
              << std::setw(t_s) << er_events_per_cluster << " "
              << std::setw(t_s) << static_cast<int>(totalclusters_in_prompt) << " "
              << std::setw(t_s) << clusters_in_prompt << " "
              << std::setw(t_s) << er_clusters_in_prompt << " "
              << std::setw(t_s) << static_cast<int>(totalclusters_in_spill) << " "
              << std::setw(t_s) << clusters_in_spill << " "
              << std::setw(t_s) << er_clusters_in_spill << " "
              << std::setw(t_s) << static_cast<int>(totalclusters_in_ext) << " "
              << std::setw(t_s) << clusters_in_ext << " "
              << std::setw(t_s) << er_clusters_in_ext << " "
              << std::setw(t_s) << static_cast<int>(totalext_rate_1) << " "
              << std::setw(t_s) << ext_rate_1 << " "
              << std::setw(t_s) << er_ext_rate_1 << " "
              << std::setw(t_s) << static_cast<int>(totalext_rate_2) << " "
              << std::setw(t_s) << ext_rate_2 << " "
              << std::setw(t_s) << er_ext_rate_2 << " "
              << std::setw(t_s) << static_cast<int>(totaltmrd_coinc) << " "
              << std::setw(t_s) << tmrd_coinc << " "
              << std::setw(t_s) << er_tmrd_coinc << " "
              << std::setw(t_s) << static_cast<int>(totalveto_hit) << " "
              << std::setw(t_s) << veto_hit << " "
              << std::setw(t_s) << er_veto_hit << " "
              << std::setw(t_s) << static_cast<int>(totalveto_tmrd_coinc) << " "
              << std::setw(t_s) << veto_tmrd_coinc << " "
              << std::setw(t_s) << er_veto_tmrd_coinc << " "
              << std::setw(t_s) << static_cast<int>(totalhas_track) << " "
              << std::setw(t_s) << has_track << " "
              << std::setw(t_s) << er_has_track << " "
              
              << std::endl;

    outfile_rates.close();  // Close the file after writing

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    // Finalize
    inputFile->Close();
    
}

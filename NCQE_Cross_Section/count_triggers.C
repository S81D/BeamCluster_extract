#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <vector>

void count_triggers(int run_number) {

    // Construct the filename using the run number
    TString directory = "/pnfs/annie/persistent/users/doran/datasets/NCQE_BEAMCLUSTER_DATA/";
    TString fileName = TString::Format("BeamCluster_%d.root", run_number);

    // Open the ROOT file
    TFile *inputFile = TFile::Open(directory + fileName);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Failed to open file: " << fileName << std::endl;
        return;
    }

    // ---------------------------------------------------------

    // Input File initialization
    Int_t beam_ok, BunchRotationOn;
    Double_t beam_pot_875, beam_pot_860, beam_THCURR, BRF_fit;
    std::vector<int>* GroupedTriggerWord = nullptr;
    ULong64_t eventTimeTank;
            
    TTree *tree = (TTree*)inputFile->Get("Event");

    // Set up variables to read from input tree
    tree->SetBranchAddress("GroupedTriggerWord", &GroupedTriggerWord);
    tree->SetBranchAddress("beam_ok", &beam_ok);
    tree->SetBranchAddress("beam_pot_875", &beam_pot_875); 
    tree->SetBranchAddress("beam_E_TOR860", &beam_pot_860);
    tree->SetBranchAddress("beam_THCURR", &beam_THCURR);
    tree->SetBranchAddress("BunchRotationOn", &BunchRotationOn);
    tree->SetBranchAddress("BRFFirstPeakFit", &BRF_fit);
	tree->SetBranchAddress("eventTimeTank", &eventTimeTank);

    Long64_t totalevents = 0;

    // TOR860
    Double_t total_POT = 0.0;


    // Loop over events in the tree
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        // Only concerned with triggers from the UBT (beam triggers)
        int trigger_index = -1; // Default to -1 if not found
        if (GroupedTriggerWord && !GroupedTriggerWord->empty()) {
            for (size_t k = 0; k < GroupedTriggerWord->size(); k++) {
                if ((*GroupedTriggerWord)[k] == 14) {
                    trigger_index = static_cast<int>(k);
                    break;
                }
            }
        }
        // Skip if trigger_index is not found
        if (trigger_index == -1) {
            continue;
        }

        // for runs prior to 4000, horn current > 176 so beam_ok = False --> calculate the beam_ok condition by hand with a higher threshold
        Int_t beam_good = 0;
        if (run_number >= 4000) {
            if (beam_ok == 1) {
                beam_good = 1;
            }
        } else {
            if (beam_THCURR > 172 && beam_THCURR < 178) {   // raise the upper threshold for beam_ok from 176 to 178
                if (beam_pot_860 > 0.5 && beam_pot_860 < 8.0 &&
                    beam_pot_875 > 0.5 && beam_pot_875 < 8.0 &&
                    ((beam_pot_875 - beam_pot_860) / beam_pot_860) < 0.05) {
                    beam_good = 1;
                }
            }
        }

        if (beam_good == 1) {

            // also only want events with 1. a usable BRF signal (rising edge was found) 2. eventTimeTank != 0 and 3. BunchRotation OFF
            if (eventTimeTank == 0 || BRF_fit == 0 || BunchRotationOn == 1){
                continue;
            }

            // only sum POT from recorded triggers with good beam conditions
            total_POT += beam_pot_860;

            // grab the total recorded triggers with good beam conditions
            totalevents++;
            
        }
    }

    inputFile->Close();

    std::cout << std::endl;
    std::cout << "Run " << run_number << std::endl;
    std::cout << "Total accepted triggers: " << totalevents << std::endl;
    std::cout << "Total POT (e12): " << total_POT << std::endl;


    // Finalize output

    TString potCSVFile = "POT_trigger_FY22_23_summary.csv";

    std::ofstream potFile;
    bool fileExists = (gSystem->AccessPathName(potCSVFile) == 0);

    potFile.open(potCSVFile.Data(), std::ios::app);
    if (potFile.is_open()) {

        if (!fileExists) {
            potFile << "RUN,TOR860_POTe12,TOTAL_TRIGGERS" << std::endl;
        }

        potFile << run_number << ","
                << std::scientific << total_POT << ","
                << totalevents << std::endl;

        potFile.close();
    }
    else {
        std::cerr << "Failed to open POT CSV file." << std::endl;
    }
    
}

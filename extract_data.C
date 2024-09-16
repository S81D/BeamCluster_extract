#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <vector>

void extract_data(int run_number) {

    // Construct the filename using the run number
    TString directory = "/pnfs/annie/persistent/processed/BeamClusterTrees/";
    TString fileName = TString::Format("BeamCluster_%d.root", run_number);

    // Open the ROOT file
    TFile *inputFile = TFile::Open(directory + fileName);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Failed to open file: " << fileName << std::endl;
        return;
    }

    // Output file initialization
    // -----------------------------------------------------

    // TFile *outputFile = new TFile("TEST.root", "RECREATE");
    TString outputDir = "2024_usable_data/";
    TString outputFileName = outputDir + TString::Format("R%d_extracted_data.ntuple.root", run_number);
    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    TTree *outputTree = new TTree("data", "Filtered data");

    Float_t CT, CT_minus_BRF, CPE, QB, POT, first_peak_fit;
    std::vector<double> hitT, hitPE, hitID, hitX, hitY, hitZ;
    Int_t isBrightest, MRD_Track, CH, Coinc, NV, part_file, MRD_activity;
    Long64_t event_number, cluster_Number, number_of_clusters;

    // branches of the output tree
    outputTree->Branch("run_number", &run_number, "run_number/I");
    outputTree->Branch("cluster_time", &CT, "cluster_time/F");
    outputTree->Branch("cluster_PE", &CPE, "cluster_PE/F");
    outputTree->Branch("cluster_Qb", &QB, "cluster_Qb/F");
    outputTree->Branch("cluster_Hits", &CH, "cluster_Hits/I");
    outputTree->Branch("cluster_time_BRF", &CT_minus_BRF, "cluster_time_BRF/F");
    outputTree->Branch("isBrightest", &isBrightest, "isBrightest/I");
    outputTree->Branch("MRD_Track", &MRD_Track, "MRD_Track/I");
    outputTree->Branch("TankMRDCoinc", &Coinc, "TankMRDCoinc/I");
    outputTree->Branch("NoVeto", &NV, "NoVeto/I");
    outputTree->Branch("POT", &POT, "POT/F");
    outputTree->Branch("part_file", &part_file, "part_file/I");
    outputTree->Branch("event_number", &event_number, "event_number/L");
    outputTree->Branch("number_of_clusters", &number_of_clusters, "number_of_clusters/L");
    outputTree->Branch("cluster_Number", &cluster_Number, "cluster_Number/L");
    outputTree->Branch("MRD_activity", &MRD_activity, "MRD_activity/I");
    outputTree->Branch("BRFFirstPeakFit", &first_peak_fit, "BRFFirstPeakFit/F");

    outputTree->Branch("hitT", "std::vector<double>", &hitT);
    outputTree->Branch("hitPE", "std::vector<double>", &hitPE);
    outputTree->Branch("hitID", "std::vector<double>", &hitID);
    outputTree->Branch("hitX", "std::vector<double>", &hitX);
    outputTree->Branch("hitY", "std::vector<double>", &hitY);
    outputTree->Branch("hitZ", "std::vector<double>", &hitZ);

    // ---------------------------------------------------------

    // Input File initialization
    std::vector<double>* clusterTime = nullptr;
    std::vector<double>* clusterPE = nullptr;
    std::vector<double>* clusterChargeBalance = nullptr;
    std::vector<int>* clusterHits = nullptr;
    std::vector<std::vector<double>>* Cluster_HitT = nullptr;
    std::vector<std::vector<double>>* Cluster_HitPE = nullptr;
    std::vector<std::vector<double>>* Cluster_HitDetID = nullptr;
    std::vector<std::vector<double>>* Cluster_HitX = nullptr;
    std::vector<std::vector<double>>* Cluster_HitY = nullptr;
    std::vector<std::vector<double>>* Cluster_HitZ = nullptr;
    std::vector<double>* MRDhitDetID = nullptr;
    Int_t partFileNumber, beam_ok, TankMRDCoinc, NoVeto, eventNumber, numberOfClusters;
    Double_t BRF_fit, beam_pot_875;
    std::vector<int>* GroupedTriggerWord = 0;
    std::vector<int>* tracks = 0;
    std::vector<ULong64_t>* GroupedTriggerTime = nullptr;
    ULong64_t eventTimeTank;

            
    TTree *tree = (TTree*)inputFile->Get("Event");

    // Set up variables to read from input tree
    tree->SetBranchAddress("runNumber", &run_number);
    tree->SetBranchAddress("partFileNumber", &partFileNumber);
    tree->SetBranchAddress("eventNumber", &eventNumber);
    tree->SetBranchAddress("clusterTime", &clusterTime);
    tree->SetBranchAddress("clusterPE", &clusterPE);
    tree->SetBranchAddress("clusterChargeBalance", &clusterChargeBalance);
    tree->SetBranchAddress("clusterHits", &clusterHits);
    tree->SetBranchAddress("Cluster_HitT", &Cluster_HitT);
    tree->SetBranchAddress("Cluster_HitPE", &Cluster_HitPE);
    tree->SetBranchAddress("Cluster_HitDetID", &Cluster_HitDetID);
    tree->SetBranchAddress("Cluster_HitX", &Cluster_HitX);
    tree->SetBranchAddress("Cluster_HitY", &Cluster_HitY);
    tree->SetBranchAddress("Cluster_HitZ", &Cluster_HitZ);
    tree->SetBranchAddress("BRFFirstPeakFit", &BRF_fit);
    tree->SetBranchAddress("GroupedTriggerWord", &GroupedTriggerWord);
    tree->SetBranchAddress("GroupedTriggerTime", &GroupedTriggerTime);
    tree->SetBranchAddress("beam_ok", &beam_ok);
    tree->SetBranchAddress("beam_pot_875", &beam_pot_875);
    tree->SetBranchAddress("NumClusterTracks", &tracks);
    tree->SetBranchAddress("numberOfClusters", &numberOfClusters);
    tree->SetBranchAddress("TankMRDCoinc", &TankMRDCoinc);
    tree->SetBranchAddress("NoVeto", &NoVeto);
    tree->SetBranchAddress("MRDhitDetID", &MRDhitDetID);
    tree->SetBranchAddress("eventTimeTank", &eventTimeTank);


    Long64_t totalclusters = 0;
    Long64_t totalevents = 0;

    Long64_t clustercount = 0;
    Long64_t eventcount = 0;

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

        // only want events with "good" beam (0.5e12 < pot < 8e12)
        if (beam_ok == 1) {

            // also only want events with 1. a usable BRF signal (rising edge was found) and 2. eventTimeTank != 0
            if (eventTimeTank == 0 || BRF_fit == 0){
                continue;
            }

            Double_t maxPE = 0;
            size_t brightestIndex = 0;
            for (size_t j = 0; j < clusterPE->size(); j++) {
                Double_t pe = clusterPE->at(j);  // Access PE value directly
                if (pe > maxPE) {
                    maxPE = pe;
                    brightestIndex = j;
                }
            }

            // find MRD tracks
            // Search through the MRD cluster array to see if any of them had a single track. If thats the case (and there was TankMRDCoinc), we can say there was a track
            //size_t theres_a_track = 0;
            if (std::find(tracks->begin(), tracks->end(), 1) != tracks->end()) {
                MRD_Track = 1;
            } else {
                MRD_Track = 0;
            }

            // check if there is any MRD activity (hits vector will be blank if there is none)
            if (MRDhitDetID && MRDhitDetID->empty()) {
                MRD_activity = 0;
            } else if (MRDhitDetID) {
                MRD_activity = 1;
            }

            part_file = partFileNumber;
            event_number = eventNumber;
            number_of_clusters = numberOfClusters;
            Coinc = TankMRDCoinc;
            NV = NoVeto;
            POT = beam_pot_875;

            // BRF first peak fit - shows the beam jitter
            first_peak_fit = BRF_fit;

            totalevents++;
            eventcount++;
            
            // Loop over each cluster
            for (size_t j = 0; j < clusterTime->size(); j++) {

                totalclusters++;
                clustercount++;

                // find the BRF subtracted cluster time
                CT_minus_BRF = clusterTime->at(j) - (BRF_fit / 1000.0);

                CT = clusterTime->at(j);
                CPE = clusterPE->at(j);
                QB = clusterChargeBalance->at(j);
                CH = clusterHits->at(j);
                isBrightest = (j == brightestIndex) ? 1 : 0;
                cluster_Number = j;
                
        
                hitT.clear(); hitPE.clear(); hitID.clear(); hitX.clear(); hitY.clear(); hitZ.clear();
                hitT = Cluster_HitT->at(j);
                hitPE = Cluster_HitPE->at(j);
                hitID = Cluster_HitDetID->at(j);
                hitX = Cluster_HitX->at(j);
                hitY = Cluster_HitY->at(j);
                hitZ = Cluster_HitZ->at(j);

                outputTree->Fill();
                
            }
        }
    }

    std::cout << std::endl;
    std::cout << run_number << " total events (clusters): " << totalevents << " (" << totalclusters << ")" << std::endl;
    std::cout << std::endl;

    // Finalize output
    outputTree->Write();
    outputFile->Close();
    inputFile->Close();
    
}

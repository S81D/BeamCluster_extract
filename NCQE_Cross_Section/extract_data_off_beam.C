#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <vector>

void extract_data_off_beam(int run_number) {

    // Construct the filename using the run number
    TString directory = "/pnfs/annie/persistent/users/doran/datasets/NCQE_BEAMCLUSTER_DATA/";
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
    TString outputDir = "off_beam/";
    TString outputFileName = outputDir + TString::Format("R%d_extracted_off_beam_data.ntuple.root", run_number);
    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    TTree *outputTree = new TTree("data", "Filtered data");

    Float_t CT, CPE, QB, recoVtx, recoVty, recoVtz;
    std::vector<double> hitT, hitPE, hitID, hitX, hitY, hitZ;
    Int_t isBrightest, MRD_Track, CH, Coinc, NV, part_file, MRD_activity, FMV_activity, hadExtended, only_prompt_cluster, throughgoing;
    Long64_t event_number, cluster_Number, number_of_clusters;

    // branches of the output tree
    outputTree->Branch("run_number", &run_number, "run_number/I");
    outputTree->Branch("cluster_time", &CT, "cluster_time/F");
    outputTree->Branch("cluster_PE", &CPE, "cluster_PE/F");
    outputTree->Branch("cluster_Qb", &QB, "cluster_Qb/F");
    outputTree->Branch("cluster_Hits", &CH, "cluster_Hits/I");
    outputTree->Branch("recoVtx", &recoVtx, "recoVtx/F");
    outputTree->Branch("recoVty", &recoVty, "recoVty/F");
    outputTree->Branch("recoVtz", &recoVtz, "recoVtz/F");
    outputTree->Branch("isBrightest", &isBrightest, "isBrightest/I");
    outputTree->Branch("MRD_Track", &MRD_Track, "MRD_Track/I");
    outputTree->Branch("TankMRDCoinc", &Coinc, "TankMRDCoinc/I");
    outputTree->Branch("NoVeto", &NV, "NoVeto/I");
    outputTree->Branch("part_file", &part_file, "part_file/I");
    outputTree->Branch("event_number", &event_number, "event_number/L");
    outputTree->Branch("number_of_clusters", &number_of_clusters, "number_of_clusters/L");
    outputTree->Branch("cluster_Number", &cluster_Number, "cluster_Number/L");
    outputTree->Branch("only_prompt_cluster", &only_prompt_cluster, "only_prompt_cluster/I");
    outputTree->Branch("MRD_activity", &MRD_activity, "MRD_activity/I");
	outputTree->Branch("FMV_activity", &FMV_activity, "FMV_activity/I");
    outputTree->Branch("hadExtended", &hadExtended, "hadExtended/I");
	outputTree->Branch("MRDThrough", &throughgoing, "MRDThrough/I");  

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
    std::vector<double>* MRDhitT = nullptr;
	std::vector<double>* FMVhitT = nullptr;
    std::vector<double>* rVtx = nullptr;
    std::vector<double>* rVty = nullptr;
    std::vector<double>* rVtz = nullptr;
    Int_t partFileNumber, TankMRDCoinc, NoVeto, eventNumber, numberOfClusters, Extended, PrimaryTriggerWord;
    std::vector<int>* tracks = 0;
	std::vector<bool>* Thru = 0;
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
    tree->SetBranchAddress("recoLeastSqVtxX", &rVtx);
    tree->SetBranchAddress("recoLeastSqVtxY", &rVty);
    tree->SetBranchAddress("recoLeastSqVtxZ", &rVtz);
    tree->SetBranchAddress("PrimaryTriggerWord", &PrimaryTriggerWord);
    tree->SetBranchAddress("NumClusterTracks", &tracks);
    tree->SetBranchAddress("numberOfClusters", &numberOfClusters);
    tree->SetBranchAddress("TankMRDCoinc", &TankMRDCoinc);
    tree->SetBranchAddress("NoVeto", &NoVeto);
    tree->SetBranchAddress("MRDhitT", &MRDhitT);
	tree->SetBranchAddress("FMVhitT", &FMVhitT);
    tree->SetBranchAddress("eventTimeTank", &eventTimeTank);
    tree->SetBranchAddress("Extended", &Extended);
	tree->SetBranchAddress("MRDThrough", &Thru);


    Long64_t totalclusters = 0;
    Long64_t totalevents = 0;     // serve as the total number of triggers for the off beam


    // Loop over events in the tree
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        // we want trigword == 31
        if (PrimaryTriggerWord != 31) {
            continue;
        }

        // also only want events with an eventTimeTank != 0
        if (eventTimeTank == 0){
            continue;
        }

        // find brightest cluster and how many are in the prompt window
        Double_t maxPE = 0;
        int clusters_in_prompt = 0;
        size_t brightestIndex = 0;
        for (size_t j = 0; j < clusterPE->size(); j++) {
            Double_t pe = clusterPE->at(j);  // Access PE value directly
            if (clusterTime->at(j) < 2000) { // cluster is in prompt window
                clusters_in_prompt++;
            }
            if (pe > maxPE) {                // brightest cluster
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

        // check if tracks are throughgoing
        if (std::find(Thru->begin(), Thru->end(), true) != Thru->end()) {
            throughgoing = 1;
        } else {
            throughgoing = 0;
        }

        part_file = partFileNumber;
        event_number = eventNumber;
        number_of_clusters = numberOfClusters;
        Coinc = TankMRDCoinc;
        NV = NoVeto;
        hadExtended = Extended;

        totalevents++;
        
        // Loop over each cluster
        for (size_t j = 0; j < clusterTime->size(); j++) {

            totalclusters++;

            CT = clusterTime->at(j);
            CPE = clusterPE->at(j);
            QB = clusterChargeBalance->at(j);
            CH = clusterHits->at(j);
            isBrightest = (j == brightestIndex) ? 1 : 0;
            cluster_Number = j;

            // modified MRD / Veto activity logic
            MRD_activity = 0;
            FMV_activity = 0;
            double cluster_t = clusterTime->at(j);
            // if any MRD or veto hits are within this coincidence window, flag it
            if (MRDhitT) {
                for (const auto& hit_t : *MRDhitT) {
                    double dt = hit_t - cluster_t;
                    if (dt >= 700.0 && dt <= 1000.0) {
                        MRD_activity = 1;
                        break;  // stop once we find one
                    }
                }
            }
            if (FMVhitT) {
                for (const auto& hit_t : *FMVhitT) {
                    double dt = hit_t - cluster_t;
                    if (dt >= 700.0 && dt <= 1000.0) {
                        FMV_activity = 1;
                        break;
                    }
                }
            }

            recoVtx = rVtx->at(j);
            recoVty = rVty->at(j);
            recoVtz = rVtz->at(j);

            only_prompt_cluster = (clusterTime->at(j) < 2000 && clusters_in_prompt == 1) ? 1 : 0;
    
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

	outputTree->Write();
    outputFile->Close();
    inputFile->Close();

    std::cout << std::endl;
    std::cout << "Run " << run_number << std::endl;
    std::cout << "Total accepted triggers: " << totalevents << std::endl;
    std::cout << "Total clusters: " << totalclusters << std::endl;


    // Finalize output

    TString potCSVFile = "offbeam_summary.csv";

    std::ofstream potFile;
    bool fileExists = (gSystem->AccessPathName(potCSVFile) == 0);

    potFile.open(potCSVFile.Data(), std::ios::app);
    if (potFile.is_open()) {

        if (!fileExists) {
            potFile << "RUN,TOTAL_TRIGGERS" << std::endl;
        }

        potFile << run_number << ","
                << totalevents << std::endl;

        potFile.close();
    }
    else {
        std::cerr << "Failed to open POT CSV file." << std::endl;
    }
    
}

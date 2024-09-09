//! HAS to be compiled,
//! root -l macros/PrepUnfolding.cxx+

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TProfile.h>
#include <TLine.h>

#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>

#include "/usr/local/eventStructuredAu/TStarJetVector.h"
#include "/usr/local/eventStructuredAu/TStarJetVectorJet.h"

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <exception>

using namespace std;

struct myJet
{
    TStarJetVectorJet orig;
    double pt;
    double weight;
    myJet(TStarJetVectorJet orig, double pt, double weight) : orig(orig), pt(pt), weight(weight) {};
};

typedef pair<myJet, myJet> matchedJet;

struct InputTreeEntry
{
    InputTreeEntry()
    {
        jets = new TClonesArray("TStarJetVectorJet", 1000);
    }
    int eventid;
    int runid;
    double weight;
    int njets;
    TClonesArray *jets;
};

TTree *initTree(InputTreeEntry &treeEntry, TString name)
{ // //! Set up Geant Chain
    TFile *inputFile = new TFile(name);
    TTree *inputTree = (TTree *)inputFile->Get("ResultTree");
    inputTree->BuildIndex("runid", "eventid");
    inputTree->GetBranch("Jets")->SetAutoDelete(kFALSE);
    inputTree->SetBranchAddress("Jets", &treeEntry.jets);
    inputTree->SetBranchAddress("eventid", &treeEntry.eventid);
    inputTree->SetBranchAddress("runid", &treeEntry.runid);
    inputTree->SetBranchAddress("weight", &treeEntry.weight);
    inputTree->SetBranchAddress("njets", &treeEntry.njets);
    return inputTree;
}

struct OutputTreeEntry
{
    double weight;
    double pt;
    double deltaR;
};

TTree *createTree(OutputTreeEntry &entry, TString name)
{
    TTree *outTree = new TTree(name, name);
    outTree->Branch("weight", &entry.weight);
    outTree->Branch("pt", &entry.pt);
    return outTree;
}

TTree *createTree(OutputTreeEntry &mc, OutputTreeEntry &reco, TString name)
{
    TTree *outTree = new TTree(name, name);
    outTree->Branch("weight", &reco.weight);
    outTree->Branch("ptMc", &mc.pt);
    outTree->Branch("ptReco", &reco.pt);
    outTree->Branch("deltaR", &mc.deltaR);
    return outTree;
}

// The algorithm to match up the jets consists of the following steps that are applied to each event:
// 1. Select particle level jets that satisfy particle level cuts and detector level jets that satisfy detector level cuts.
// 2. Take the Cartesian product of these two sets. For each pair of particle and detector jets from that product, assign a distance deltaR.
// 3. Select the pair with minimum value of deltaR. If deltaR < 0.2, call them “matching”. If not, skip to the next event.
// 4. Remove all pairs that have the same particle or detector jet as a pair from the previous step.
// 5. Go back to step 3 until the list of pairs is empty.
// Main algorithm for matching mcJets and recoJets

void matchJets(const std::vector<myJet> &mcJets, const std::vector<myJet> &recoJets,
               std::vector<std::pair<int, int>> &matches,
               std::vector<int> &misses, std::vector<int> &fakes, double deltaRMax = 0.2)
{
    std::vector<bool> mcMatched(mcJets.size(), false);
    std::vector<bool> recoMatched(recoJets.size(), false);

    // cout << "Matching " << mcJets.size() << " mcJets and " << recoJets.size() << " recoJets" << endl;
    while (true)
    {
        double minDeltaR = 10000;
        int bestMcIndex = -1, bestRecoIndex = -1;

        // Find the pair with the minimum deltaR
        for (size_t i = 0; i < mcJets.size(); ++i)
        {
            if (mcMatched[i])
                continue; // Skip already matched mcJet

            for (size_t j = 0; j < recoJets.size(); ++j)
            {
                if (recoMatched[j])
                    continue; // Skip already matched recoJet

                double deltaR = mcJets[i].orig.DeltaR(recoJets[j].orig);
                // cout << "DeltaR: " << deltaR << endl;
                if (deltaR < minDeltaR)
                {
                    minDeltaR = deltaR;
                    bestMcIndex = i;
                    bestRecoIndex = j;
                }
            }
        }

        // If no valid pair found, or the minimum deltaR is greater than 0.2, exit
        if (bestMcIndex == -1 || bestRecoIndex == -1 || minDeltaR >= deltaRMax)
        {
            break;
        }

        // Mark the jets as matched
        matches.push_back({bestMcIndex, bestRecoIndex});
        mcMatched[bestMcIndex] = true;
        recoMatched[bestRecoIndex] = true;
        // cout << "Matched: " << mcJets[bestMcIndex].orig.Pt() << " " << recoJets[bestRecoIndex].orig.Pt() << endl;
        // cout << "Minimum DeltaR: " << minDeltaR << endl;

        // Collect misses (unmatched mcJets)
        for (size_t i = 0; i < mcJets.size(); ++i)
        {
            if (!mcMatched[i])
            {
                misses.push_back(i); // Collect indices of unmatched mcJets
            }
        }

        // Collect fakes (unmatched recoJets)
        for (size_t j = 0; j < recoJets.size(); ++j)
        {
            if (!recoMatched[j])
            {
                fakes.push_back(j); // Collect indices of unmatched recoJets
            }
        }
    }
    // cout << "=========================   Done Matching   =========================" << endl;
    //
}

int myMatching(float jetRad = 0.4)
{
    TString mcTreeName = "./output/tree_pt-hat2025_41_mc.root";
    TString recoTreeName = "./output/tree_pt-hat2025_41.root";

    float EtaCut = 1.0 - jetRad;
    // Set Minimum Constituent pt
    double ConMinPt = 0.2;
    const double mcJetMinPt = 5;   // GeV
    const double recoJetMinPt = 5; // GeV
    const double deltaRMax = 0.2;

    // Output
    // ------
    TString outFileName = "./output/test_matching.root";

    int MatchNumber = 0;
    int FakeNumber = 0;
    int MissNumber = 0;
    int MissEventNumber = 0;
    int MatchedGeantEventNumber = 0;
    int TotalGeantEventNumber = 0;
    int FakeEventNumber = 0;

    // Set up Trees

    InputTreeEntry inputMc;
    TClonesArray *mcJets = inputMc.jets;

    InputTreeEntry inputReco;
    TClonesArray *recoJets = inputReco.jets;

    TTree *mcTree = initTree(inputMc, mcTreeName);
    TTree *recoTree = initTree(inputReco, recoTreeName);

    // new from Isaac
    const int NUMBEROFPT = 13;
    const static float XSEC[NUMBEROFPT] = {0.00230158, 0.000342755, 0.0000457002, 9.0012, .00000972535, 1.46253, 0.000000469889, 0.354566, 0.0000000269202, 0.00000000143453, 0.151622, 0.0249062, 0.00584527};
    const static float NUMBEROFEVENT[NUMBEROFPT] = {3000000, 3000000, 3000000, 3000000, 2000000, 3000000, 2000000, 3000000, 1000000, 1000000, 3000000, 3000000, 3000000};
    const static float MAXPT[NUMBEROFPT] = {30, 40, 50, 6, 70, 8, 90, 10, 110, 2000, 14, 18, 22};
    const static vector<string> vptbins = {"1115_", "1520_", "2025_", "23_", "2535_", "34_", "3545_", "45_", "4555_", "55999_", "57_", "79_", "911_"};
    // Truth and Measured Jet Bounds are the same for now
    int MeasJetBins = 10;
    double MeasJetBounds[11] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60};
    int TruthJetBins = 10;
    double TruthJetBounds[11] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60};
    // Output and histograms
    TFile *fout = new TFile(outFileName, "RECREATE");

    // SetUp Output Trees
    OutputTreeEntry reco;
    OutputTreeEntry mc;

    TTree *MatchTree = createTree(mc, reco, "MatchTree");
    TTree *FakeTree = createTree(reco, "FakeTree");
    TTree *MissTree = createTree(mc, "MissTree");

    TH1D *hDeltaR = new TH1D("hDeltaR", "#Delta R all; #Delta R", 1000, 0, 5);
    TH1D *hPtMc = new TH1D("hPtMc", "Mc p_{t}; p_{t}, GeV/c", 1000, 0, 100);
    TH1D *hPtReco = new TH1D("hPtReco", "Reco p_{t}; p_{t}, GeV/c", 1000, 0, 100);
    TH2D *hPtMcReco = new TH2D("hPtMcReco", "Mc p_{t} vs Reco p_{t}; Mc p_{t}, GeV/c; Reco p_{t}, GeV/c", 1000, 0, 100, 1000, 0, 100);
    TH1D *hDeltaRMatched = new TH1D("hDeltaRMatched", "#Delta R matched; #Delta R", 100, 0, 1);

    TH2D *hNJets = new TH2D("hNJets", "; Number of Mc jets; Number of Reco jets", 20, 0, 20, 20, 0, 20);

    //================================================================================================
    // Begin Looping over MC events
    //================================================================================================

    int nEventsMc = mcTree->GetEntries();
    vector<int> recoEventNumbers;
    int maxEvents = 1000;

    for (Long64_t iEvent = 0; iEvent < nEventsMc; ++iEvent)
    {
        if (maxEvents > 0 && iEvent > maxEvents)
        {
            break;
        }
        if (!(iEvent % 10000))
            cout << "Working on " << iEvent << " / " << nEventsMc << endl;
        mcTree->GetEntry(iEvent);

        if (mcJets->GetEntries() != inputMc.njets)
        {
            cerr << "mcJets->GetEntries() != inputMc.njets" << endl;
            return -1;
        }

        vector<myJet> myMcJets;
        for (auto j = 0; j < mcJets->GetEntries(); ++j)
        { // get the values for the current jet
            TStarJetVectorJet *jet = (TStarJetVectorJet *)mcJets->At(j);
            //! Fill in jet into pythia result
            if (fabs(jet->Eta()) < EtaCut && jet->Pt() > mcJetMinPt)
            {
                myMcJets.push_back(myJet(*jet, jet->GetSumConstituentsPt(), inputMc.weight));
            }
        }
        // cout << "======================" << endl;
        // cout << "======================" << endl;
        if (myMcJets.size() == 0)
        {
            continue;
        } // skip this event
        // Still in MC level loop, get matching Geant Event
        Long64_t recoEvent = recoTree->GetEntryNumberWithIndex(inputMc.runid, inputMc.eventid);
        if (recoEvent >= 0)
        {
            recoEventNumbers.push_back(recoEvent);
            MatchedGeantEventNumber++;
        }
        //  //bug in new embedding
        //     if(mctotalpT > 23.11003 && mctotalpT < 23.11004){continue;} //11-15
        //     if(mctotalpT > 33.749385 && mctotalpT < 33.749405){continue;} //15-20
        //     if(mctotalpT > 47.09071 && mctotalpT < 47.09072){continue;} //20-25
        //     if(mctotalpT > 9.62226 && mctotalpT < 9.62227){continue;} //2-3
        //     if(mctotalpT > 46.63831 && mctotalpT < 46.63832){continue;} //25-35
        //     if(mctotalpT > 6.90831 && mctotalpT < 6.90832){continue;} //3-4
        //     if(mctotalpT > 82.68752 && mctotalpT < 82.68753){continue;} //35-45
        //     if(mctotalpT > 100.25616 && mctotalpT < 100.25617){continue;} //45-55
        //     if(mctotalpT > 75.10883 && mctotalpT < 75.10884){continue;} //55-999
        //     if(mctotalpT > 3.75004 && mctotalpT < 3.75005){continue;} //5-7
        //     if(mctotalpT > 6.47623 && mctotalpT < 6.47624){continue;} //7-9
        //     if(mctotalpT > 4.22790 && mctotalpT < 4.22791){continue;} //9-11

        // fill in regardless of match/miss to get reference and check which events are candidates for high weight and throw them out

        //  skip events with high weights and high pt jets to avoid bias
        //================================================================================================
        bool isBad = false;
        int weightBin = -1;
        bool old = false;
        for (int i = 0; i < vptbins.size(); ++i)
        {
            if (inputMc.weight == XSEC[i] / NUMBEROFEVENT[i])
            {
                old = true;
                weightBin = i;
            }
        }
        for (auto jet : myMcJets)
        {
            double Jetpt = jet.orig.Pt();
            if (old && Jetpt > MAXPT[weightBin])
            {
                isBad = true;
                break;
            }
        }
        if (isBad)
        {
            continue;
        }
        //================================================================================================

        // matching geant event not found, fill in pythia event as a miss
        if (recoEvent < 0)
        {
            for (auto jet : myMcJets)
            {
                MissEventNumber++;
                mc.pt = jet.orig.Pt();
                mc.weight = jet.weight;
                // if (mc.pt > TruthJetBounds[0] && mc.pt < TruthJetBounds[TruthJetBins])
                MissTree->Fill();
            }
        }
        // Get Geant event if its found
        recoTree->GetEntry(recoEvent);
        vector<myJet> myRecoJets;
        for (int j = 0; j < inputReco.njets; ++j)
        { // get the values for the current jet
            TStarJetVectorJet *jet = (TStarJetVectorJet *)recoJets->At(j);
            //! Fill in jet into pythia result
            if (fabs(jet->Eta()) < EtaCut && jet->Pt() > recoJetMinPt)
            {
                myRecoJets.push_back(myJet(*jet, jet->GetSumConstituentsPt(), inputReco.weight));
            }
        }

        for (auto recoJet : myRecoJets)
        {
            double Jetpt = recoJet.orig.perp();
            if (Jetpt > MAXPT[weightBin])
            {
                isBad = true;
                break;
            }
        }
        if (isBad)
        {
            continue;
        }

        // match and Sort Pythia and Geant jets together
        int nJetsMc = myMcJets.size();
        int nJetsReco = myRecoJets.size();
        int matchedJets = 0;
        int missedJets = 0;
        int fakeJets = 0;

        for (auto mcJetStruct : myMcJets)
        {
            TStarJetVectorJet mcJet = mcJetStruct.orig;
            for (auto recoJetStruct : myRecoJets)
            {
                TStarJetVectorJet recoJet = recoJetStruct.orig;
                hDeltaR->Fill(mcJet.DeltaR(recoJet));
                hPtReco->Fill(recoJet.perp());
            }
            hPtMc->Fill(mcJet.perp());
        }
        hNJets->Fill(nJetsMc, nJetsReco);

        std::vector<std::pair<int, int>> matches;
        std::vector<int> misses; // List of unmatched mcJets (Misses)
        std::vector<int> fakes;  // List of unmatched recoJets (Fakes)

        matchJets(myMcJets, myRecoJets, matches, misses, fakes, deltaRMax);

        for (auto match : matches)
        {
            auto mcJet = myMcJets[match.first];
            auto recoJet = myRecoJets[match.second];
            // delete matched jets from vector
            mc.pt = mcJet.orig.perp();
            mc.weight = mcJet.weight;
            reco.pt = recoJet.orig.perp();
            reco.weight = recoJet.weight;
            mc.deltaR = mcJet.orig.DeltaR(recoJet.orig);

            hDeltaRMatched->Fill(mc.deltaR);
            hPtMcReco->Fill(mc.pt, reco.pt);

            MatchTree->Fill();
            MatchNumber++;
        }

        for (int mcIndex : misses)
        {
            auto mcJet = myMcJets[mcIndex];
            mc.pt = mcJet.orig.perp();
            mc.weight = mcJet.weight;
            MissTree->Fill();
            MissNumber++;
        }

        for (int recoIndex : fakes)
        {
            auto recoJet = myRecoJets[recoIndex];
            reco.pt = recoJet.orig.perp();
            reco.weight = recoJet.weight;
            FakeTree->Fill();
            FakeNumber++;
        }
    }
    fout->cd();
    fout->GetList()->Write();

    fout->Close();

    cout << "Miss Number is " << MissNumber << endl;
    cout << "From Missed Events is " << MissEventNumber << endl;
    cout << "Fake Number is " << FakeNumber << endl;
    cout << "From Fake Events is " << FakeEventNumber << endl;
    cout << "MatchNumber is " << MatchNumber << endl;

    return 0;
}
// end of function
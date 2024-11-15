//! Code to read in geant + pythia output trees and match them

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TProfile.h>
#include <TLine.h>

#include <TROOT.h>
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
#include <TGraph.h>

#include <TStarJetVector.h>
#include <TStarJetVectorJet.h>
#include <TStarJetPicoReader.h>

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <exception>

using namespace std;

struct JetWithInfo
{
    TStarJetVectorJet orig;
    int n_constituents;
    int event_id;
    double weight;
    double multiplicity;
    bool is_rejected;
    double pt;
    double eta;
    double phi;
    double neutral_fraction;
    double vz;

    JetWithInfo(TStarJetVectorJet _orig, int n_constituents,
                int _event_id, double _weight, double _multiplicity,
                bool _is_rejected, double _neutral_fraction, double _vz) : orig(_orig),
                                                                           n_constituents(n_constituents),
                                                                           event_id(_event_id),
                                                                           weight(_weight),
                                                                           multiplicity(_multiplicity),
                                                                           is_rejected(_is_rejected),
                                                                           neutral_fraction(_neutral_fraction),
                                                                           vz(_vz)
    {
        pt = orig.Pt();
        eta = orig.Eta();
        phi = orig.Phi();
    };
};

typedef pair<JetWithInfo, JetWithInfo> MatchedJetWithInfo;

struct InputTreeEntry
{
    InputTreeEntry()
    {
        jets = new TClonesArray("TStarJetVectorJet", 1000);
    }

    ~InputTreeEntry()
    {
        delete jets;
    }

    int runid;
    int runid1;
    int eventid;
    double weight;
    double refmult;
    int njets;
    double vz;
    int mult;
    float event_sum_pt;
    bool is_rejected;
    TClonesArray *jets;
    double neutral_fraction[1000];
    double pt[1000];
    int n_constituents[1000];
    int index[1000];
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
    inputTree->SetBranchAddress("runid1", &treeEntry.runid1);
    inputTree->SetBranchAddress("weight", &treeEntry.weight);
    inputTree->SetBranchAddress("refmult", &treeEntry.refmult);
    inputTree->SetBranchAddress("njets", &treeEntry.njets);
    inputTree->SetBranchAddress("vz", &treeEntry.vz);
    inputTree->SetBranchAddress("mult", &treeEntry.mult);
    inputTree->SetBranchAddress("event_sum_pt", &treeEntry.event_sum_pt);
    inputTree->SetBranchAddress("is_rejected", &treeEntry.is_rejected);
    inputTree->SetBranchAddress("neutral_fraction", treeEntry.neutral_fraction);
    inputTree->SetBranchAddress("pt", treeEntry.pt);
    inputTree->SetBranchAddress("n_constituents", treeEntry.n_constituents);
    inputTree->SetBranchAddress("index", treeEntry.index);
    return inputTree;
}

struct OutputTreeEntry
{
    OutputTreeEntry() {}
    OutputTreeEntry(const JetWithInfo &jet)
    {
        n_constituents = jet.n_constituents;
        event_id = jet.event_id;
        weight = jet.weight;
        multiplicity = jet.multiplicity;
        is_rejected = jet.is_rejected;
        pt = jet.pt;
        eta = jet.eta;
        phi = jet.phi;
        neutral_fraction = jet.neutral_fraction;
        vz = jet.vz;
    }
    int n_constituents;
    int event_id;
    double weight;
    double multiplicity;
    bool is_rejected;
    double pt;
    double eta;
    double phi;
    double neutral_fraction;
    double vz;
};

TTree *createTree(OutputTreeEntry &entry, TString name)
{
    TTree *outTree = new TTree(name, name);
    outTree->Branch("weight", &entry.weight);
    outTree->Branch("pt", &entry.pt);
    outTree->Branch("eta", &entry.eta);
    outTree->Branch("phi", &entry.phi);
    outTree->Branch("n_constituents", &entry.n_constituents);
    outTree->Branch("event_id", &entry.event_id);
    outTree->Branch("multiplicity", &entry.multiplicity);
    outTree->Branch("is_rejected", &entry.is_rejected);
    outTree->Branch("neutral_fraction", &entry.neutral_fraction);
    outTree->Branch("vz", &entry.vz);
    return outTree;
}

TTree *createTree(OutputTreeEntry &mc, OutputTreeEntry &reco, double &deltaR, TString name)
{
    TTree *outTree = new TTree(name, name);
    outTree->Branch("deltaR", &deltaR);
    outTree->Branch("mc_weight", &mc.weight);
    outTree->Branch("reco_weight", &reco.weight);
    outTree->Branch("mc_pt", &mc.pt);
    outTree->Branch("reco_pt", &reco.pt);
    outTree->Branch("mc_n_constituents", &mc.n_constituents);
    outTree->Branch("reco_n_constituents", &reco.n_constituents);
    outTree->Branch("mc_event_id", &mc.event_id);
    outTree->Branch("reco_event_id", &reco.event_id);
    outTree->Branch("mc_multiplicity", &mc.multiplicity);
    outTree->Branch("reco_multiplicity", &reco.multiplicity);
    outTree->Branch("mc_is_rejected", &mc.is_rejected);
    outTree->Branch("reco_is_rejected", &reco.is_rejected);
    outTree->Branch("mc_neutral_fraction", &mc.neutral_fraction);
    outTree->Branch("reco_neutral_fraction", &reco.neutral_fraction);
    outTree->Branch("mc_vz", &mc.vz);
    outTree->Branch("reco_vz", &reco.vz);
    return outTree;
}

// The algorithm to match up the jets consists of the following steps that are applied to each event:
// 1. Select particle level jets that satisfy particle level cuts and detector level jets that satisfy detector level cuts.
// 2. Take the Cartesian product of these two sets. For each pair of particle and detector jets from that product, assign a distance deltaR.
// 3. Select the pair with minimum value of deltaR. If deltaR < 0.2, call them “matching”. If not, skip to the next event.
// 4. Remove all pairs that have the same particle or detector jet as a pair from the previous step.
// 5. Go back to step 3 until the list of pairs is empty.
// Main algorithm for matching mcJets and recoJets

void matchJets(const std::vector<JetWithInfo> &mcJets, const std::vector<JetWithInfo> &recoJets,
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

int MatchGeantToPythia(TString mcTreeName = "output/tree_pt-hat2025_41_mc.root", TString outFileName = "output/matching.root")
{
    const float jetRad = 0.4;
    float EtaCut = 1.0 - jetRad;
    const double mcJetMinPt = 3;   // GeV
    const double recoJetMinPt = 5; // GeV
    const double deltaRMax = 0.4;

    int MatchNumber = 0;
    int FakeNumber = 0;
    int MissNumber = 0;
    int MissEventNumber = 0;
    int MatchedGeantEventNumber = 0;
    int TotalGeantEventNumber = 0;
    int FakeEventNumber = 0;
    // =================================================================================================
    TString mcBaseName = mcTreeName(mcTreeName.Last('/') + 1, mcTreeName.Length());
    TString mcFolder = "/gpfs01/star/pwg/prozorov/jets_pp_2012/output/mc/";
    TString recoFolder = "/gpfs01/star/pwg/prozorov/jets_pp_2012/output/geant/";
    TString recoTreeName = recoFolder + mcBaseName;
    mcTreeName = mcFolder + mcBaseName;

    mcTreeName = "output/test/tree_pt-hat79_19_mc.root";
    recoTreeName = "output/test/tree_pt-hat79_19_geant.root";

    // =================================================================================================
    InputTreeEntry mc;
    InputTreeEntry reco;
    TTree *mcTree = initTree(mc, mcTreeName);
    TTree *recoTree = initTree(reco, recoTreeName);
    // =================================================================================================

    TFile *fout = new TFile(outFileName, "RECREATE");
    OutputTreeEntry reco_result;
    OutputTreeEntry mc_result;
    double deltaR;
    TTree *MatchTree = createTree(mc_result, reco_result, deltaR, "MatchTree");
    TTree *FakeTree = createTree(reco_result, "FakeTree");
    TTree *MissTree = createTree(mc_result, "MissTree");
    // =================================================================================================
    //   histograms
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    TH1D *hDeltaR = new TH1D("hDeltaR", "#Delta R all; #Delta R", 350, 0, 3.5);
    TH1D *hPtMc = new TH1D("hPtMc", "Mc p_{t}; p_{t}, GeV/c", 500, 0, 50);
    TH1D *hPtReco = new TH1D("hPtReco", "Reco p_{t}; p_{t}, GeV/c", 500, 0, 50);
    TH2D *hPtMcReco = new TH2D("hPtMcReco", "Mc p_{t} vs Reco p_{t}; Mc p_{t}, GeV/c; Reco p_{t}, GeV/c", 500, 0, 50, 500, 0, 50);
    TH2D *hNJets = new TH2D("hNJets", "; Number of Mc jets; Number of Reco jets", 20, 0, 20, 20, 0, 20);
    TH1D *hMissedJets = new TH1D("hMissedJets", "Missed Jets; Events", 10, 0, 10);
    TH1D *hFakeJets = new TH1D("hFakeJets", "Fake Jets; Events", 20, 0, 20);
    TH1D *hMatchedJets = new TH1D("hMatchedJets", "Matched Jets; Events", 20, 0, 20);

    TH1D *hMissRate = new TH1D("hMissRate", "Miss Rate; p_{t}, GeV/c", 500, 0, 50);
    TH1D *hFakeRate = new TH1D("hFakeRate", "Fake Rate; p_{t}, GeV/c", 500, 0, 50);

    TH1D *stats = new TH1D("stats", "stats", 4, 0, 4);
    stats->GetXaxis()->SetBinLabel(1, "Match");
    stats->GetXaxis()->SetBinLabel(2, "Miss");
    stats->GetXaxis()->SetBinLabel(3, "Fake");
    stats->GetXaxis()->SetBinLabel(4, "Miss(no geant)");
    //! =================================================================================================
    int N = mcTree->GetEntries();
    cout << "Number of Pythia events: " << N << endl;
    cout << "Number of Geant events:  " << recoTree->GetEntries() << endl;

    set<float> mc_accepted_events_list;
    set<float> reco_list;

    for (Long64_t mcEvent = 0; mcEvent < N; ++mcEvent) // event loop
    {
        if (!(mcEvent % 500000))
            cout << "Working on " << mcEvent << " / " << N << endl;
        mcTree->GetEntry(mcEvent);
        if (mc.is_rejected || abs(mc.vz) > 30) // skip bad events with high weights
            continue;
        if (mc_accepted_events_list.find(mc.event_sum_pt) != mc_accepted_events_list.end())
            continue; // some events have identical total particle pT
        else
            mc_accepted_events_list.insert(mc.event_sum_pt);
        //================================================================================================
        vector<JetWithInfo> mc_jets;
        // Loop over mc jets and select jets that satisfy the cuts
        for (int j = 0; j < mc.njets; ++j)
        {
            TStarJetVectorJet *mc_jet = (TStarJetVectorJet *)mc.jets->At(j);
            if (fabs(mc_jet->Eta()) < EtaCut && mc_jet->Pt() > mcJetMinPt)
                mc_jets.push_back(JetWithInfo(*mc_jet, mc.n_constituents[j], mc.eventid, mc.weight, mc.mult, mc.is_rejected, mc.neutral_fraction[j], mc.vz));
        } // end of mcjet loop

        int recoEvent = recoTree->GetEntryNumberWithIndex(mc.runid, mc.eventid);

        if (recoEvent < 0)
        {
            for (auto mc_jet : mc_jets)
            {
                mc_result = OutputTreeEntry(mc_jet);
                MissTree->Fill();
                MissEventNumber++;
            }
            continue;
        }

        recoTree->GetEntry(recoEvent);

        if (reco.is_rejected || abs(reco.vz) > 30)
            continue; // go to next MC event

        MatchedGeantEventNumber++;

        //================================================================================================
        vector<JetWithInfo> reco_jets;
        for (int j = 0; j < reco.njets; ++j)
        {
            TStarJetVectorJet *reco_jet = (TStarJetVectorJet *)reco.jets->At(j);
            if (fabs(reco_jet->Eta()) < EtaCut && reco_jet->Pt() > recoJetMinPt)
                reco_jets.push_back(JetWithInfo(*reco_jet, reco.n_constituents[j], reco.eventid, reco.weight, reco.mult, reco.is_rejected, reco.neutral_fraction[j], reco.vz));
        }
        //================================================================================================
        hNJets->Fill(mc_jets.size(), reco_jets.size());

        for (auto mc_jet : mc_jets)
        {
            TStarJetVectorJet mcJet = mc_jet.orig;
            hPtMc->Fill(mcJet.Pt());
            for (auto reco_jet : reco_jets)
            {
                TStarJetVectorJet recoJet = reco_jet.orig;
                hDeltaR->Fill(mcJet.DeltaR(recoJet));
            }
        }

        for (auto reco_jet : reco_jets)
        {
            TStarJetVectorJet recoJet = reco_jet.orig;
            hPtReco->Fill(recoJet.Pt());
        }

        //================================================================================================
        std::vector<std::pair<int, int>> matches;
        std::vector<int> misses; // List of unmatched mcJets (Misses)
        std::vector<int> fakes;  // List of unmatched recoJets (Fakes)

        matchJets(mc_jets, reco_jets, matches, misses, fakes, deltaRMax);

        hMatchedJets->Fill(matches.size());

        for (auto match : matches)
        {
            deltaR = mc_jets[match.first].orig.DeltaR(reco_jets[match.second].orig);
            mc_result = OutputTreeEntry(mc_jets[match.first]);
            reco_result = OutputTreeEntry(reco_jets[match.second]);
            hPtMcReco->Fill(mc_result.pt, reco_result.pt);
            MatchTree->Fill();
            MatchNumber++;
        }
        //================================================================================================
        hMissedJets->Fill(misses.size());
        for (auto miss : misses)
        {
            mc_result = OutputTreeEntry(mc_jets[miss]);
            MissTree->Fill();
            hMissRate->Fill(mc_result.pt);
            MissNumber++;
        }
        //================================================================================================
        hFakeJets->Fill(fakes.size());
        for (auto fake : fakes)
        {
            reco_result = OutputTreeEntry(reco_jets[fake]);
            FakeTree->Fill();
            hFakeRate->Fill(reco_result.pt);
            FakeNumber++;
        }
        //================================================================================================
    }
    cout << "Miss Number is " << MissNumber << endl;
    cout << "From Missed Events is " << MissEventNumber << endl;
    cout << "Fake Number is " << FakeNumber << endl;
    cout << "MatchNumber is " << MatchNumber << endl;

    stats->SetBinContent(1, MatchNumber);
    stats->SetBinContent(2, MissNumber);
    stats->SetBinContent(3, FakeNumber);
    stats->SetBinContent(4, MissEventNumber);

    stats->SetTitle(Form("pt_{mc}>%.1f GeV/c, pt_{reco}>%.1f GeV/c, #Delta R < %.1f", mcJetMinPt, recoJetMinPt, deltaRMax));

    fout->cd();
    fout->GetList()->Write();

    fout->Close();

    return 0;
}

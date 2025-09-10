//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar  8 15:05:33 2018 by ROOT version 5.34/34
// from TTree JetTree/ Pico Tree for Jet
// found on file: sum0.root
//////////////////////////////////////////////////////////

#ifndef foo_h
#define foo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxfPrimaryTracks = 73;
   const Int_t kMaxfFtpcPrimaryTracks = 1;
   const Int_t kMaxfTowers = 109;
   const Int_t kMaxfV0s = 1;
   const Int_t kMaxfTrigObjs = 10288;

class foo {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //TStarJetPicoEvent *PicoJetTree;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UInt_t          fEventHeader_fUniqueID;
   UInt_t          fEventHeader_fBits;
   Int_t           fEventHeader_fEventId;
   Int_t           fEventHeader_fRunId;
   Int_t           fEventHeader_fRefMult;
   Int_t           fEventHeader_fGRefMult;
   Int_t           fEventHeader_fRefCent;
   Int_t           fEventHeader_fGRefCent;
   Double_t        fEventHeader_fRefCentWeight;
   Double_t        fEventHeader_fGRefCentWeight;
   Double_t        fEventHeader_fCorRefMult;
   Double_t        fEventHeader_fCorGRefMult;
   Int_t           fEventHeader_fNOfGlobalTracks;
   Float_t         fEventHeader_fReactionPlaneAngle;
   Int_t           fEventHeader_fNOfTriggerIds;
   TArrayI         fEventHeader_fTriggerIdArray;
   Int_t           fEventHeader_fNOfTowerTrackMatched;
   Int_t           fEventHeader_fNOfTowers;
   Int_t           fEventHeader_fNOfPrimaryTracks;
   Int_t           fEventHeader_fNOfMatchedTracks;
   Int_t           fEventHeader_fNOfFtpcPrimaryTracks;
   Int_t           fEventHeader_fNOfV0s;
   Int_t           fEventHeader_fNOfEMCPoints;
   Float_t         fEventHeader_fPVx;
   Float_t         fEventHeader_fPVy;
   Float_t         fEventHeader_fPVz;
   Float_t         fEventHeader_fvpdVz;
   Float_t         fEventHeader_fCTBmult;
   Float_t         fEventHeader_fMeanDip;
   Float_t         fEventHeader_fRank;
   Int_t           fEventHeader_fNOfVertices;
   Int_t           fEventHeader_fNOfTrigObjs;
   Int_t           fEventHeader_fDSMInput;
   Int_t           fEventHeader_fTrigMask;
   Float_t         fEventHeader_fZdcWestRate;
   Float_t         fEventHeader_fZdcEastRate;
   Float_t         fEventHeader_fZdcCoincidenceRate;
   Float_t         fEventHeader_fBbcWestRate;
   Float_t         fEventHeader_fBbcEastRate;
   Float_t         fEventHeader_fBbcCoincidenceRate;
   Float_t         fEventHeader_fBbcBlueBackgroundRate;
   Float_t         fEventHeader_fBbcYellowBackgroundRate;
   Int_t           fEventHeader_fBbcAdcSumEast;
   Int_t           fEventHeader_fBbcOnlineVertex;
   Float_t         fEventHeader_fBbcOfflineVertex;
   Int_t           fEventHeader_fRefMultFTPCE;
   Int_t           fEventHeader_fnumberOfVpdEastHits;
   Int_t           fEventHeader_fnumberOfVpdWestHits;
   Int_t           fPrimaryTracks_;
   UInt_t          fPrimaryTracks_fUniqueID[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   UInt_t          fPrimaryTracks_fBits[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Int_t           fPrimaryTracks_fCharge[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Int_t           fPrimaryTracks_fNFittedHits[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Int_t           fPrimaryTracks_fNHitsPoss[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Int_t           fPrimaryTracks_fKey[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fPx[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fPy[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fPz[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fDCA[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fdEdx[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fdEdxPoints[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fNsigmaPion[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fNsigmaKaon[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fNsigmaProton[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fNsigmaElectron[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fEtaDiffHitProjected[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fPhiDiffHitProjected[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fsDCAxy[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fChi2[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fChi2PV[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Int_t           fPrimaryTracks_fFlag[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Bool_t          fPrimaryTracks_fBemcMatchFlag[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Bool_t          fPrimaryTracks_fTofMatchFlag[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fTofTime[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fTofBeta[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Float_t         fPrimaryTracks_fTofyLocal[kMaxfPrimaryTracks];   //[fPrimaryTracks_]
   Int_t           fFtpcPrimaryTracks_;
   UInt_t          fFtpcPrimaryTracks_fUniqueID[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   UInt_t          fFtpcPrimaryTracks_fBits[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Int_t           fFtpcPrimaryTracks_fCharge[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Int_t           fFtpcPrimaryTracks_fNFittedHits[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Int_t           fFtpcPrimaryTracks_fNHitsPoss[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Int_t           fFtpcPrimaryTracks_fKey[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fPx[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fPy[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fPz[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fDCA[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fdEdx[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fdEdxPoints[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fNsigmaPion[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fNsigmaKaon[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fNsigmaProton[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fNsigmaElectron[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fEtaDiffHitProjected[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fPhiDiffHitProjected[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fsDCAxy[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fChi2[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fChi2PV[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Int_t           fFtpcPrimaryTracks_fFlag[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Bool_t          fFtpcPrimaryTracks_fBemcMatchFlag[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Bool_t          fFtpcPrimaryTracks_fTofMatchFlag[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fTofTime[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fTofBeta[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Float_t         fFtpcPrimaryTracks_fTofyLocal[kMaxfFtpcPrimaryTracks];   //[fFtpcPrimaryTracks_]
   Int_t           fTowers_;
   UInt_t          fTowers_fUniqueID[kMaxfTowers];   //[fTowers_]
   UInt_t          fTowers_fBits[kMaxfTowers];   //[fTowers_]
   Int_t           fTowers_fSMDClusterP[kMaxfTowers];   //[fTowers_]
   Int_t           fTowers_fSMDClusterE[kMaxfTowers];   //[fTowers_]
   Int_t           fTowers_fTowerStatus[kMaxfTowers];   //[fTowers_]
   Int_t           fTowers_fId[kMaxfTowers];   //[fTowers_]
   Float_t         fTowers_fEnergy[kMaxfTowers];   //[fTowers_]
   Float_t         fTowers_fEta[kMaxfTowers];   //[fTowers_]
   Float_t         fTowers_fPhi[kMaxfTowers];   //[fTowers_]
   Int_t           fTowers_fADC[kMaxfTowers];   //[fTowers_]
   Float_t         fTowers_fEtaCorrected[kMaxfTowers];   //[fTowers_]
   Float_t         fTowers_fPhiCorrected[kMaxfTowers];   //[fTowers_]
   Int_t           fTowers_fNAssocTracks[kMaxfTowers];   //[fTowers_]
   TArrayI         fTowers_fMatchedTracks[kMaxfTowers];
   Int_t           fV0s_;
   UInt_t          fV0s_fUniqueID[kMaxfV0s];   //[fV0s_]
   UInt_t          fV0s_fBits[kMaxfV0s];   //[fV0s_]
   Int_t           fV0s_fKeyPos[kMaxfV0s];   //[fV0s_]
   Int_t           fV0s_fKeyNeg[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fPxPos[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fPyPos[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fPzPos[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fPxNeg[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fPyNeg[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fPzNeg[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fDcapn[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fDcaV0[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fDcap[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fDcan[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fDLength[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fDedxPos[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fDedxNeg[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fNSigmaProtonPos[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fNSigmaProtonNeg[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fNSigmaPionPos[kMaxfV0s];   //[fV0s_]
   Float_t         fV0s_fNSigmaPionNeg[kMaxfV0s];   //[fV0s_]
   Int_t           fTrigObjs_;
   UInt_t          fTrigObjs_fUniqueID[kMaxfTrigObjs];   //[fTrigObjs_]
   UInt_t          fTrigObjs_fBits[kMaxfTrigObjs];   //[fTrigObjs_]
   Float_t         fTrigObjs_fEta[kMaxfTrigObjs];   //[fTrigObjs_]
   Float_t         fTrigObjs_fPhi[kMaxfTrigObjs];   //[fTrigObjs_]
   Int_t           fTrigObjs_fTrigFlag[kMaxfTrigObjs];   //[fTrigObjs_]
   Int_t           fTrigObjs_fId[kMaxfTrigObjs];   //[fTrigObjs_]
   Int_t           fTrigObjs_fADC[kMaxfTrigObjs];   //[fTrigObjs_]
   TRef            fLastPrimaryTrack;
   TRef            fLastFtpcPrimaryTrack;
   TRef            fLastTower;
   TRef            fLastV0;
   TRef            fLastTrigObj;
   Double_t        fZdcsmd[2][2][8];

   // List of branches
   TBranch        *b_PicoJetTree_fUniqueID;   //!
   TBranch        *b_PicoJetTree_fBits;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fUniqueID;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fBits;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fEventId;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fRunId;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fRefMult;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fGRefMult;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fRefCent;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fGRefCent;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fRefCentWeight;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fGRefCentWeight;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fCorRefMult;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fCorGRefMult;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfGlobalTracks;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fReactionPlaneAngle;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfTriggerIds;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fTriggerIdArray;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfTowerTrackMatched;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfTowers;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfPrimaryTracks;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfMatchedTracks;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfFtpcPrimaryTracks;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfV0s;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfEMCPoints;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fPVx;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fPVy;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fPVz;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fvpdVz;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fCTBmult;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fMeanDip;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fRank;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfVertices;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fNOfTrigObjs;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fDSMInput;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fTrigMask;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fZdcWestRate;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fZdcEastRate;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fZdcCoincidenceRate;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fBbcWestRate;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fBbcEastRate;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fBbcCoincidenceRate;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fBbcBlueBackgroundRate;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fBbcYellowBackgroundRate;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fBbcAdcSumEast;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fBbcOnlineVertex;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fBbcOfflineVertex;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fRefMultFTPCE;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fnumberOfVpdEastHits;   //!
   TBranch        *b_PicoJetTree_fEventHeader_fnumberOfVpdWestHits;   //!
   TBranch        *b_PicoJetTree_fPrimaryTracks_;   //!
   TBranch        *b_fPrimaryTracks_fUniqueID;   //!
   TBranch        *b_fPrimaryTracks_fBits;   //!
   TBranch        *b_fPrimaryTracks_fCharge;   //!
   TBranch        *b_fPrimaryTracks_fNFittedHits;   //!
   TBranch        *b_fPrimaryTracks_fNHitsPoss;   //!
   TBranch        *b_fPrimaryTracks_fKey;   //!
   TBranch        *b_fPrimaryTracks_fPx;   //!
   TBranch        *b_fPrimaryTracks_fPy;   //!
   TBranch        *b_fPrimaryTracks_fPz;   //!
   TBranch        *b_fPrimaryTracks_fDCA;   //!
   TBranch        *b_fPrimaryTracks_fdEdx;   //!
   TBranch        *b_fPrimaryTracks_fdEdxPoints;   //!
   TBranch        *b_fPrimaryTracks_fNsigmaPion;   //!
   TBranch        *b_fPrimaryTracks_fNsigmaKaon;   //!
   TBranch        *b_fPrimaryTracks_fNsigmaProton;   //!
   TBranch        *b_fPrimaryTracks_fNsigmaElectron;   //!
   TBranch        *b_fPrimaryTracks_fEtaDiffHitProjected;   //!
   TBranch        *b_fPrimaryTracks_fPhiDiffHitProjected;   //!
   TBranch        *b_fPrimaryTracks_fsDCAxy;   //!
   TBranch        *b_fPrimaryTracks_fChi2;   //!
   TBranch        *b_fPrimaryTracks_fChi2PV;   //!
   TBranch        *b_fPrimaryTracks_fFlag;   //!
   TBranch        *b_fPrimaryTracks_fBemcMatchFlag;   //!
   TBranch        *b_fPrimaryTracks_fTofMatchFlag;   //!
   TBranch        *b_fPrimaryTracks_fTofTime;   //!
   TBranch        *b_fPrimaryTracks_fTofBeta;   //!
   TBranch        *b_fPrimaryTracks_fTofyLocal;   //!
   TBranch        *b_PicoJetTree_fFtpcPrimaryTracks_;   //!
   TBranch        *b_fFtpcPrimaryTracks_fUniqueID;   //!
   TBranch        *b_fFtpcPrimaryTracks_fBits;   //!
   TBranch        *b_fFtpcPrimaryTracks_fCharge;   //!
   TBranch        *b_fFtpcPrimaryTracks_fNFittedHits;   //!
   TBranch        *b_fFtpcPrimaryTracks_fNHitsPoss;   //!
   TBranch        *b_fFtpcPrimaryTracks_fKey;   //!
   TBranch        *b_fFtpcPrimaryTracks_fPx;   //!
   TBranch        *b_fFtpcPrimaryTracks_fPy;   //!
   TBranch        *b_fFtpcPrimaryTracks_fPz;   //!
   TBranch        *b_fFtpcPrimaryTracks_fDCA;   //!
   TBranch        *b_fFtpcPrimaryTracks_fdEdx;   //!
   TBranch        *b_fFtpcPrimaryTracks_fdEdxPoints;   //!
   TBranch        *b_fFtpcPrimaryTracks_fNsigmaPion;   //!
   TBranch        *b_fFtpcPrimaryTracks_fNsigmaKaon;   //!
   TBranch        *b_fFtpcPrimaryTracks_fNsigmaProton;   //!
   TBranch        *b_fFtpcPrimaryTracks_fNsigmaElectron;   //!
   TBranch        *b_fFtpcPrimaryTracks_fEtaDiffHitProjected;   //!
   TBranch        *b_fFtpcPrimaryTracks_fPhiDiffHitProjected;   //!
   TBranch        *b_fFtpcPrimaryTracks_fsDCAxy;   //!
   TBranch        *b_fFtpcPrimaryTracks_fChi2;   //!
   TBranch        *b_fFtpcPrimaryTracks_fChi2PV;   //!
   TBranch        *b_fFtpcPrimaryTracks_fFlag;   //!
   TBranch        *b_fFtpcPrimaryTracks_fBemcMatchFlag;   //!
   TBranch        *b_fFtpcPrimaryTracks_fTofMatchFlag;   //!
   TBranch        *b_fFtpcPrimaryTracks_fTofTime;   //!
   TBranch        *b_fFtpcPrimaryTracks_fTofBeta;   //!
   TBranch        *b_fFtpcPrimaryTracks_fTofyLocal;   //!
   TBranch        *b_PicoJetTree_fTowers_;   //!
   TBranch        *b_fTowers_fUniqueID;   //!
   TBranch        *b_fTowers_fBits;   //!
   TBranch        *b_fTowers_fSMDClusterP;   //!
   TBranch        *b_fTowers_fSMDClusterE;   //!
   TBranch        *b_fTowers_fTowerStatus;   //!
   TBranch        *b_fTowers_fId;   //!
   TBranch        *b_fTowers_fEnergy;   //!
   TBranch        *b_fTowers_fEta;   //!
   TBranch        *b_fTowers_fPhi;   //!
   TBranch        *b_fTowers_fADC;   //!
   TBranch        *b_fTowers_fEtaCorrected;   //!
   TBranch        *b_fTowers_fPhiCorrected;   //!
   TBranch        *b_fTowers_fNAssocTracks;   //!
   TBranch        *b_fTowers_fMatchedTracks;   //!
   TBranch        *b_PicoJetTree_fV0s_;   //!
   TBranch        *b_fV0s_fUniqueID;   //!
   TBranch        *b_fV0s_fBits;   //!
   TBranch        *b_fV0s_fKeyPos;   //!
   TBranch        *b_fV0s_fKeyNeg;   //!
   TBranch        *b_fV0s_fPxPos;   //!
   TBranch        *b_fV0s_fPyPos;   //!
   TBranch        *b_fV0s_fPzPos;   //!
   TBranch        *b_fV0s_fPxNeg;   //!
   TBranch        *b_fV0s_fPyNeg;   //!
   TBranch        *b_fV0s_fPzNeg;   //!
   TBranch        *b_fV0s_fDcapn;   //!
   TBranch        *b_fV0s_fDcaV0;   //!
   TBranch        *b_fV0s_fDcap;   //!
   TBranch        *b_fV0s_fDcan;   //!
   TBranch        *b_fV0s_fDLength;   //!
   TBranch        *b_fV0s_fDedxPos;   //!
   TBranch        *b_fV0s_fDedxNeg;   //!
   TBranch        *b_fV0s_fNSigmaProtonPos;   //!
   TBranch        *b_fV0s_fNSigmaProtonNeg;   //!
   TBranch        *b_fV0s_fNSigmaPionPos;   //!
   TBranch        *b_fV0s_fNSigmaPionNeg;   //!
   TBranch        *b_PicoJetTree_fTrigObjs_;   //!
   TBranch        *b_fTrigObjs_fUniqueID;   //!
   TBranch        *b_fTrigObjs_fBits;   //!
   TBranch        *b_fTrigObjs_fEta;   //!
   TBranch        *b_fTrigObjs_fPhi;   //!
   TBranch        *b_fTrigObjs_fTrigFlag;   //!
   TBranch        *b_fTrigObjs_fId;   //!
   TBranch        *b_fTrigObjs_fADC;   //!
   TBranch        *b_PicoJetTree_fLastPrimaryTrack;   //!
   TBranch        *b_PicoJetTree_fLastFtpcPrimaryTrack;   //!
   TBranch        *b_PicoJetTree_fLastTower;   //!
   TBranch        *b_PicoJetTree_fLastV0;   //!
   TBranch        *b_PicoJetTree_fLastTrigObj;   //!
   TBranch        *b_PicoJetTree_fZdcsmd;   //!

   foo(TTree *tree=0);
   virtual ~foo();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef foo_cxx
foo::foo(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sum0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sum0.root");
      }
      f->GetObject("JetTree",tree);

   }
   Init(tree);
}

foo::~foo()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t foo::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t foo::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void foo::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_PicoJetTree_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_PicoJetTree_fBits);
   fChain->SetBranchAddress("fEventHeader.fUniqueID", &fEventHeader_fUniqueID, &b_PicoJetTree_fEventHeader_fUniqueID);
   fChain->SetBranchAddress("fEventHeader.fBits", &fEventHeader_fBits, &b_PicoJetTree_fEventHeader_fBits);
   fChain->SetBranchAddress("fEventHeader.fEventId", &fEventHeader_fEventId, &b_PicoJetTree_fEventHeader_fEventId);
   fChain->SetBranchAddress("fEventHeader.fRunId", &fEventHeader_fRunId, &b_PicoJetTree_fEventHeader_fRunId);
   fChain->SetBranchAddress("fEventHeader.fRefMult", &fEventHeader_fRefMult, &b_PicoJetTree_fEventHeader_fRefMult);
   fChain->SetBranchAddress("fEventHeader.fGRefMult", &fEventHeader_fGRefMult, &b_PicoJetTree_fEventHeader_fGRefMult);
   fChain->SetBranchAddress("fEventHeader.fRefCent", &fEventHeader_fRefCent, &b_PicoJetTree_fEventHeader_fRefCent);
   fChain->SetBranchAddress("fEventHeader.fGRefCent", &fEventHeader_fGRefCent, &b_PicoJetTree_fEventHeader_fGRefCent);
   fChain->SetBranchAddress("fEventHeader.fRefCentWeight", &fEventHeader_fRefCentWeight, &b_PicoJetTree_fEventHeader_fRefCentWeight);
   fChain->SetBranchAddress("fEventHeader.fGRefCentWeight", &fEventHeader_fGRefCentWeight, &b_PicoJetTree_fEventHeader_fGRefCentWeight);
   fChain->SetBranchAddress("fEventHeader.fCorRefMult", &fEventHeader_fCorRefMult, &b_PicoJetTree_fEventHeader_fCorRefMult);
   fChain->SetBranchAddress("fEventHeader.fCorGRefMult", &fEventHeader_fCorGRefMult, &b_PicoJetTree_fEventHeader_fCorGRefMult);
   fChain->SetBranchAddress("fEventHeader.fNOfGlobalTracks", &fEventHeader_fNOfGlobalTracks, &b_PicoJetTree_fEventHeader_fNOfGlobalTracks);
   fChain->SetBranchAddress("fEventHeader.fReactionPlaneAngle", &fEventHeader_fReactionPlaneAngle, &b_PicoJetTree_fEventHeader_fReactionPlaneAngle);
   fChain->SetBranchAddress("fEventHeader.fNOfTriggerIds", &fEventHeader_fNOfTriggerIds, &b_PicoJetTree_fEventHeader_fNOfTriggerIds);
   fChain->SetBranchAddress("fEventHeader.fTriggerIdArray", &fEventHeader_fTriggerIdArray, &b_PicoJetTree_fEventHeader_fTriggerIdArray);
   fChain->SetBranchAddress("fEventHeader.fNOfTowerTrackMatched", &fEventHeader_fNOfTowerTrackMatched, &b_PicoJetTree_fEventHeader_fNOfTowerTrackMatched);
   fChain->SetBranchAddress("fEventHeader.fNOfTowers", &fEventHeader_fNOfTowers, &b_PicoJetTree_fEventHeader_fNOfTowers);
   fChain->SetBranchAddress("fEventHeader.fNOfPrimaryTracks", &fEventHeader_fNOfPrimaryTracks, &b_PicoJetTree_fEventHeader_fNOfPrimaryTracks);
   fChain->SetBranchAddress("fEventHeader.fNOfMatchedTracks", &fEventHeader_fNOfMatchedTracks, &b_PicoJetTree_fEventHeader_fNOfMatchedTracks);
   fChain->SetBranchAddress("fEventHeader.fNOfFtpcPrimaryTracks", &fEventHeader_fNOfFtpcPrimaryTracks, &b_PicoJetTree_fEventHeader_fNOfFtpcPrimaryTracks);
   fChain->SetBranchAddress("fEventHeader.fNOfV0s", &fEventHeader_fNOfV0s, &b_PicoJetTree_fEventHeader_fNOfV0s);
   fChain->SetBranchAddress("fEventHeader.fNOfEMCPoints", &fEventHeader_fNOfEMCPoints, &b_PicoJetTree_fEventHeader_fNOfEMCPoints);
   fChain->SetBranchAddress("fEventHeader.fPVx", &fEventHeader_fPVx, &b_PicoJetTree_fEventHeader_fPVx);
   fChain->SetBranchAddress("fEventHeader.fPVy", &fEventHeader_fPVy, &b_PicoJetTree_fEventHeader_fPVy);
   fChain->SetBranchAddress("fEventHeader.fPVz", &fEventHeader_fPVz, &b_PicoJetTree_fEventHeader_fPVz);
   fChain->SetBranchAddress("fEventHeader.fvpdVz", &fEventHeader_fvpdVz, &b_PicoJetTree_fEventHeader_fvpdVz);
   fChain->SetBranchAddress("fEventHeader.fCTBmult", &fEventHeader_fCTBmult, &b_PicoJetTree_fEventHeader_fCTBmult);
   fChain->SetBranchAddress("fEventHeader.fMeanDip", &fEventHeader_fMeanDip, &b_PicoJetTree_fEventHeader_fMeanDip);
   fChain->SetBranchAddress("fEventHeader.fRank", &fEventHeader_fRank, &b_PicoJetTree_fEventHeader_fRank);
   fChain->SetBranchAddress("fEventHeader.fNOfVertices", &fEventHeader_fNOfVertices, &b_PicoJetTree_fEventHeader_fNOfVertices);
   fChain->SetBranchAddress("fEventHeader.fNOfTrigObjs", &fEventHeader_fNOfTrigObjs, &b_PicoJetTree_fEventHeader_fNOfTrigObjs);
   fChain->SetBranchAddress("fEventHeader.fDSMInput", &fEventHeader_fDSMInput, &b_PicoJetTree_fEventHeader_fDSMInput);
   fChain->SetBranchAddress("fEventHeader.fTrigMask", &fEventHeader_fTrigMask, &b_PicoJetTree_fEventHeader_fTrigMask);
   fChain->SetBranchAddress("fEventHeader.fZdcWestRate", &fEventHeader_fZdcWestRate, &b_PicoJetTree_fEventHeader_fZdcWestRate);
   fChain->SetBranchAddress("fEventHeader.fZdcEastRate", &fEventHeader_fZdcEastRate, &b_PicoJetTree_fEventHeader_fZdcEastRate);
   fChain->SetBranchAddress("fEventHeader.fZdcCoincidenceRate", &fEventHeader_fZdcCoincidenceRate, &b_PicoJetTree_fEventHeader_fZdcCoincidenceRate);
   fChain->SetBranchAddress("fEventHeader.fBbcWestRate", &fEventHeader_fBbcWestRate, &b_PicoJetTree_fEventHeader_fBbcWestRate);
   fChain->SetBranchAddress("fEventHeader.fBbcEastRate", &fEventHeader_fBbcEastRate, &b_PicoJetTree_fEventHeader_fBbcEastRate);
   fChain->SetBranchAddress("fEventHeader.fBbcCoincidenceRate", &fEventHeader_fBbcCoincidenceRate, &b_PicoJetTree_fEventHeader_fBbcCoincidenceRate);
   fChain->SetBranchAddress("fEventHeader.fBbcBlueBackgroundRate", &fEventHeader_fBbcBlueBackgroundRate, &b_PicoJetTree_fEventHeader_fBbcBlueBackgroundRate);
   fChain->SetBranchAddress("fEventHeader.fBbcYellowBackgroundRate", &fEventHeader_fBbcYellowBackgroundRate, &b_PicoJetTree_fEventHeader_fBbcYellowBackgroundRate);
   fChain->SetBranchAddress("fEventHeader.fBbcAdcSumEast", &fEventHeader_fBbcAdcSumEast, &b_PicoJetTree_fEventHeader_fBbcAdcSumEast);
   fChain->SetBranchAddress("fEventHeader.fBbcOnlineVertex", &fEventHeader_fBbcOnlineVertex, &b_PicoJetTree_fEventHeader_fBbcOnlineVertex);
   fChain->SetBranchAddress("fEventHeader.fBbcOfflineVertex", &fEventHeader_fBbcOfflineVertex, &b_PicoJetTree_fEventHeader_fBbcOfflineVertex);
   fChain->SetBranchAddress("fEventHeader.fRefMultFTPCE", &fEventHeader_fRefMultFTPCE, &b_PicoJetTree_fEventHeader_fRefMultFTPCE);
   fChain->SetBranchAddress("fEventHeader.fnumberOfVpdEastHits", &fEventHeader_fnumberOfVpdEastHits, &b_PicoJetTree_fEventHeader_fnumberOfVpdEastHits);
   fChain->SetBranchAddress("fEventHeader.fnumberOfVpdWestHits", &fEventHeader_fnumberOfVpdWestHits, &b_PicoJetTree_fEventHeader_fnumberOfVpdWestHits);
   fChain->SetBranchAddress("fPrimaryTracks", &fPrimaryTracks_, &b_PicoJetTree_fPrimaryTracks_);
   fChain->SetBranchAddress("fPrimaryTracks.fUniqueID", fPrimaryTracks_fUniqueID, &b_fPrimaryTracks_fUniqueID);
   fChain->SetBranchAddress("fPrimaryTracks.fBits", fPrimaryTracks_fBits, &b_fPrimaryTracks_fBits);
   fChain->SetBranchAddress("fPrimaryTracks.fCharge", fPrimaryTracks_fCharge, &b_fPrimaryTracks_fCharge);
   fChain->SetBranchAddress("fPrimaryTracks.fNFittedHits", fPrimaryTracks_fNFittedHits, &b_fPrimaryTracks_fNFittedHits);
   fChain->SetBranchAddress("fPrimaryTracks.fNHitsPoss", fPrimaryTracks_fNHitsPoss, &b_fPrimaryTracks_fNHitsPoss);
   fChain->SetBranchAddress("fPrimaryTracks.fKey", fPrimaryTracks_fKey, &b_fPrimaryTracks_fKey);
   fChain->SetBranchAddress("fPrimaryTracks.fPx", fPrimaryTracks_fPx, &b_fPrimaryTracks_fPx);
   fChain->SetBranchAddress("fPrimaryTracks.fPy", fPrimaryTracks_fPy, &b_fPrimaryTracks_fPy);
   fChain->SetBranchAddress("fPrimaryTracks.fPz", fPrimaryTracks_fPz, &b_fPrimaryTracks_fPz);
   fChain->SetBranchAddress("fPrimaryTracks.fDCA", fPrimaryTracks_fDCA, &b_fPrimaryTracks_fDCA);
   fChain->SetBranchAddress("fPrimaryTracks.fdEdx", fPrimaryTracks_fdEdx, &b_fPrimaryTracks_fdEdx);
   fChain->SetBranchAddress("fPrimaryTracks.fdEdxPoints", fPrimaryTracks_fdEdxPoints, &b_fPrimaryTracks_fdEdxPoints);
   fChain->SetBranchAddress("fPrimaryTracks.fNsigmaPion", fPrimaryTracks_fNsigmaPion, &b_fPrimaryTracks_fNsigmaPion);
   fChain->SetBranchAddress("fPrimaryTracks.fNsigmaKaon", fPrimaryTracks_fNsigmaKaon, &b_fPrimaryTracks_fNsigmaKaon);
   fChain->SetBranchAddress("fPrimaryTracks.fNsigmaProton", fPrimaryTracks_fNsigmaProton, &b_fPrimaryTracks_fNsigmaProton);
   fChain->SetBranchAddress("fPrimaryTracks.fNsigmaElectron", fPrimaryTracks_fNsigmaElectron, &b_fPrimaryTracks_fNsigmaElectron);
   fChain->SetBranchAddress("fPrimaryTracks.fEtaDiffHitProjected", fPrimaryTracks_fEtaDiffHitProjected, &b_fPrimaryTracks_fEtaDiffHitProjected);
   fChain->SetBranchAddress("fPrimaryTracks.fPhiDiffHitProjected", fPrimaryTracks_fPhiDiffHitProjected, &b_fPrimaryTracks_fPhiDiffHitProjected);
   fChain->SetBranchAddress("fPrimaryTracks.fsDCAxy", fPrimaryTracks_fsDCAxy, &b_fPrimaryTracks_fsDCAxy);
   fChain->SetBranchAddress("fPrimaryTracks.fChi2", fPrimaryTracks_fChi2, &b_fPrimaryTracks_fChi2);
   fChain->SetBranchAddress("fPrimaryTracks.fChi2PV", fPrimaryTracks_fChi2PV, &b_fPrimaryTracks_fChi2PV);
   fChain->SetBranchAddress("fPrimaryTracks.fFlag", fPrimaryTracks_fFlag, &b_fPrimaryTracks_fFlag);
   fChain->SetBranchAddress("fPrimaryTracks.fBemcMatchFlag", fPrimaryTracks_fBemcMatchFlag, &b_fPrimaryTracks_fBemcMatchFlag);
   fChain->SetBranchAddress("fPrimaryTracks.fTofMatchFlag", fPrimaryTracks_fTofMatchFlag, &b_fPrimaryTracks_fTofMatchFlag);
   fChain->SetBranchAddress("fPrimaryTracks.fTofTime", fPrimaryTracks_fTofTime, &b_fPrimaryTracks_fTofTime);
   fChain->SetBranchAddress("fPrimaryTracks.fTofBeta", fPrimaryTracks_fTofBeta, &b_fPrimaryTracks_fTofBeta);
   fChain->SetBranchAddress("fPrimaryTracks.fTofyLocal", fPrimaryTracks_fTofyLocal, &b_fPrimaryTracks_fTofyLocal);
   fChain->SetBranchAddress("fFtpcPrimaryTracks", &fFtpcPrimaryTracks_, &b_PicoJetTree_fFtpcPrimaryTracks_);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fUniqueID", &fFtpcPrimaryTracks_fUniqueID, &b_fFtpcPrimaryTracks_fUniqueID);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fBits", &fFtpcPrimaryTracks_fBits, &b_fFtpcPrimaryTracks_fBits);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fCharge", &fFtpcPrimaryTracks_fCharge, &b_fFtpcPrimaryTracks_fCharge);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fNFittedHits", &fFtpcPrimaryTracks_fNFittedHits, &b_fFtpcPrimaryTracks_fNFittedHits);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fNHitsPoss", &fFtpcPrimaryTracks_fNHitsPoss, &b_fFtpcPrimaryTracks_fNHitsPoss);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fKey", &fFtpcPrimaryTracks_fKey, &b_fFtpcPrimaryTracks_fKey);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fPx", &fFtpcPrimaryTracks_fPx, &b_fFtpcPrimaryTracks_fPx);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fPy", &fFtpcPrimaryTracks_fPy, &b_fFtpcPrimaryTracks_fPy);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fPz", &fFtpcPrimaryTracks_fPz, &b_fFtpcPrimaryTracks_fPz);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fDCA", &fFtpcPrimaryTracks_fDCA, &b_fFtpcPrimaryTracks_fDCA);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fdEdx", &fFtpcPrimaryTracks_fdEdx, &b_fFtpcPrimaryTracks_fdEdx);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fdEdxPoints", &fFtpcPrimaryTracks_fdEdxPoints, &b_fFtpcPrimaryTracks_fdEdxPoints);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fNsigmaPion", &fFtpcPrimaryTracks_fNsigmaPion, &b_fFtpcPrimaryTracks_fNsigmaPion);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fNsigmaKaon", &fFtpcPrimaryTracks_fNsigmaKaon, &b_fFtpcPrimaryTracks_fNsigmaKaon);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fNsigmaProton", &fFtpcPrimaryTracks_fNsigmaProton, &b_fFtpcPrimaryTracks_fNsigmaProton);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fNsigmaElectron", &fFtpcPrimaryTracks_fNsigmaElectron, &b_fFtpcPrimaryTracks_fNsigmaElectron);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fEtaDiffHitProjected", &fFtpcPrimaryTracks_fEtaDiffHitProjected, &b_fFtpcPrimaryTracks_fEtaDiffHitProjected);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fPhiDiffHitProjected", &fFtpcPrimaryTracks_fPhiDiffHitProjected, &b_fFtpcPrimaryTracks_fPhiDiffHitProjected);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fsDCAxy", &fFtpcPrimaryTracks_fsDCAxy, &b_fFtpcPrimaryTracks_fsDCAxy);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fChi2", &fFtpcPrimaryTracks_fChi2, &b_fFtpcPrimaryTracks_fChi2);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fChi2PV", &fFtpcPrimaryTracks_fChi2PV, &b_fFtpcPrimaryTracks_fChi2PV);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fFlag", &fFtpcPrimaryTracks_fFlag, &b_fFtpcPrimaryTracks_fFlag);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fBemcMatchFlag", &fFtpcPrimaryTracks_fBemcMatchFlag, &b_fFtpcPrimaryTracks_fBemcMatchFlag);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fTofMatchFlag", &fFtpcPrimaryTracks_fTofMatchFlag, &b_fFtpcPrimaryTracks_fTofMatchFlag);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fTofTime", &fFtpcPrimaryTracks_fTofTime, &b_fFtpcPrimaryTracks_fTofTime);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fTofBeta", &fFtpcPrimaryTracks_fTofBeta, &b_fFtpcPrimaryTracks_fTofBeta);
   fChain->SetBranchAddress("fFtpcPrimaryTracks.fTofyLocal", &fFtpcPrimaryTracks_fTofyLocal, &b_fFtpcPrimaryTracks_fTofyLocal);
   fChain->SetBranchAddress("fTowers", &fTowers_, &b_PicoJetTree_fTowers_);
   fChain->SetBranchAddress("fTowers.fUniqueID", fTowers_fUniqueID, &b_fTowers_fUniqueID);
   fChain->SetBranchAddress("fTowers.fBits", fTowers_fBits, &b_fTowers_fBits);
   fChain->SetBranchAddress("fTowers.fSMDClusterP", fTowers_fSMDClusterP, &b_fTowers_fSMDClusterP);
   fChain->SetBranchAddress("fTowers.fSMDClusterE", fTowers_fSMDClusterE, &b_fTowers_fSMDClusterE);
   fChain->SetBranchAddress("fTowers.fTowerStatus", fTowers_fTowerStatus, &b_fTowers_fTowerStatus);
   fChain->SetBranchAddress("fTowers.fId", fTowers_fId, &b_fTowers_fId);
   fChain->SetBranchAddress("fTowers.fEnergy", fTowers_fEnergy, &b_fTowers_fEnergy);
   fChain->SetBranchAddress("fTowers.fEta", fTowers_fEta, &b_fTowers_fEta);
   fChain->SetBranchAddress("fTowers.fPhi", fTowers_fPhi, &b_fTowers_fPhi);
   fChain->SetBranchAddress("fTowers.fADC", fTowers_fADC, &b_fTowers_fADC);
   fChain->SetBranchAddress("fTowers.fEtaCorrected", fTowers_fEtaCorrected, &b_fTowers_fEtaCorrected);
   fChain->SetBranchAddress("fTowers.fPhiCorrected", fTowers_fPhiCorrected, &b_fTowers_fPhiCorrected);
   fChain->SetBranchAddress("fTowers.fNAssocTracks", fTowers_fNAssocTracks, &b_fTowers_fNAssocTracks);
   fChain->SetBranchAddress("fTowers.fMatchedTracks", fTowers_fMatchedTracks, &b_fTowers_fMatchedTracks);
   fChain->SetBranchAddress("fV0s", &fV0s_, &b_PicoJetTree_fV0s_);
   fChain->SetBranchAddress("fV0s.fUniqueID", &fV0s_fUniqueID, &b_fV0s_fUniqueID);
   fChain->SetBranchAddress("fV0s.fBits", &fV0s_fBits, &b_fV0s_fBits);
   fChain->SetBranchAddress("fV0s.fKeyPos", &fV0s_fKeyPos, &b_fV0s_fKeyPos);
   fChain->SetBranchAddress("fV0s.fKeyNeg", &fV0s_fKeyNeg, &b_fV0s_fKeyNeg);
   fChain->SetBranchAddress("fV0s.fPxPos", &fV0s_fPxPos, &b_fV0s_fPxPos);
   fChain->SetBranchAddress("fV0s.fPyPos", &fV0s_fPyPos, &b_fV0s_fPyPos);
   fChain->SetBranchAddress("fV0s.fPzPos", &fV0s_fPzPos, &b_fV0s_fPzPos);
   fChain->SetBranchAddress("fV0s.fPxNeg", &fV0s_fPxNeg, &b_fV0s_fPxNeg);
   fChain->SetBranchAddress("fV0s.fPyNeg", &fV0s_fPyNeg, &b_fV0s_fPyNeg);
   fChain->SetBranchAddress("fV0s.fPzNeg", &fV0s_fPzNeg, &b_fV0s_fPzNeg);
   fChain->SetBranchAddress("fV0s.fDcapn", &fV0s_fDcapn, &b_fV0s_fDcapn);
   fChain->SetBranchAddress("fV0s.fDcaV0", &fV0s_fDcaV0, &b_fV0s_fDcaV0);
   fChain->SetBranchAddress("fV0s.fDcap", &fV0s_fDcap, &b_fV0s_fDcap);
   fChain->SetBranchAddress("fV0s.fDcan", &fV0s_fDcan, &b_fV0s_fDcan);
   fChain->SetBranchAddress("fV0s.fDLength", &fV0s_fDLength, &b_fV0s_fDLength);
   fChain->SetBranchAddress("fV0s.fDedxPos", &fV0s_fDedxPos, &b_fV0s_fDedxPos);
   fChain->SetBranchAddress("fV0s.fDedxNeg", &fV0s_fDedxNeg, &b_fV0s_fDedxNeg);
   fChain->SetBranchAddress("fV0s.fNSigmaProtonPos", &fV0s_fNSigmaProtonPos, &b_fV0s_fNSigmaProtonPos);
   fChain->SetBranchAddress("fV0s.fNSigmaProtonNeg", &fV0s_fNSigmaProtonNeg, &b_fV0s_fNSigmaProtonNeg);
   fChain->SetBranchAddress("fV0s.fNSigmaPionPos", &fV0s_fNSigmaPionPos, &b_fV0s_fNSigmaPionPos);
   fChain->SetBranchAddress("fV0s.fNSigmaPionNeg", &fV0s_fNSigmaPionNeg, &b_fV0s_fNSigmaPionNeg);
   fChain->SetBranchAddress("fTrigObjs", &fTrigObjs_, &b_PicoJetTree_fTrigObjs_);
   fChain->SetBranchAddress("fTrigObjs.fUniqueID", fTrigObjs_fUniqueID, &b_fTrigObjs_fUniqueID);
   fChain->SetBranchAddress("fTrigObjs.fBits", fTrigObjs_fBits, &b_fTrigObjs_fBits);
   fChain->SetBranchAddress("fTrigObjs.fEta", fTrigObjs_fEta, &b_fTrigObjs_fEta);
   fChain->SetBranchAddress("fTrigObjs.fPhi", fTrigObjs_fPhi, &b_fTrigObjs_fPhi);
   fChain->SetBranchAddress("fTrigObjs.fTrigFlag", fTrigObjs_fTrigFlag, &b_fTrigObjs_fTrigFlag);
   fChain->SetBranchAddress("fTrigObjs.fId", fTrigObjs_fId, &b_fTrigObjs_fId);
   fChain->SetBranchAddress("fTrigObjs.fADC", fTrigObjs_fADC, &b_fTrigObjs_fADC);
   fChain->SetBranchAddress("fLastPrimaryTrack", &fLastPrimaryTrack, &b_PicoJetTree_fLastPrimaryTrack);
   fChain->SetBranchAddress("fLastFtpcPrimaryTrack", &fLastFtpcPrimaryTrack, &b_PicoJetTree_fLastFtpcPrimaryTrack);
   fChain->SetBranchAddress("fLastTower", &fLastTower, &b_PicoJetTree_fLastTower);
   fChain->SetBranchAddress("fLastV0", &fLastV0, &b_PicoJetTree_fLastV0);
   fChain->SetBranchAddress("fLastTrigObj", &fLastTrigObj, &b_PicoJetTree_fLastTrigObj);
   fChain->SetBranchAddress("fZdcsmd[2][2][8]", fZdcsmd, &b_PicoJetTree_fZdcsmd);
   Notify();
}

Bool_t foo::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void foo::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t foo::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef foo_cxx

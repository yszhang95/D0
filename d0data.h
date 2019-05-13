//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 11 15:28:35 2019 by ROOT version 6.06/01
// from TChain d0data/d0data
// found on file: /storage1/users/wl33/D0Trees/Data/Merged_pPbPbpData_MVATree_D0_default_BDTCut03_v1.root
//////////////////////////////////////////////////////////

#ifndef d0data_h
#define d0data_h

#include "d0tree.h"

// Header file for the classes stored in the TChain if any.

class d0data : public d0tree {
public :
   TChain          *fChain;   //!pointer to the analyzed TChain or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TChain if any.

   // Declaration of leaf types
   Float_t         pT;
   Float_t         y;
   Float_t         mass;
   Float_t         mva;
   Int_t           Ntrkoffline;
   Int_t           Npixel;
   Float_t         HFsumET;
   Float_t         bestvtxX;
   Float_t         bestvtxY;
   Float_t         bestvtxZ;
   Float_t         flavor;
   Float_t         eta;
   Float_t         VtxProb;
   Float_t         m3DCosPointingAngle;
   Float_t         m3DPointingAngle;
   Float_t         m2DCosPointingAngle;
   Float_t         m2DPointingAngle;
   Float_t         m3DDecayLengthSignificance;
   Float_t         m3DDecayLength;
   Float_t         m2DDecayLengthSignificance;
   Float_t         m2DDecayLength;
   Float_t         zDCASignificanceDaugther1;
   Float_t         xyDCASignificanceDaugther1;
   Float_t         NHitD1;
   Bool_t          HighPurityDaugther1;
   Float_t         pTD1;
   Float_t         pTerrD1;
   Float_t         EtaD1;
   Float_t         dedxHarmonic2D1;
   Float_t         zDCASignificanceDaugther2;
   Float_t         xyDCASignificanceDaugther2;
   Float_t         NHitD2;
   Bool_t          HighPurityDaugther2;
   Float_t         pTD2;
   Float_t         pTerrD2;
   Float_t         EtaD2;
   Float_t         dedxHarmonic2D2;

   // List of branches
   TBranch        *b_pT;   //!
   TBranch        *b_y;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_mva;   //!
   TBranch        *b_Ntrkoffline;   //!
   TBranch        *b_Npixel;   //!
   TBranch        *b_HFsumET;   //!
   TBranch        *b_bestvtxX;   //!
   TBranch        *b_bestvtxY;   //!
   TBranch        *b_bestvtxZ;   //!
   TBranch        *b_flavor;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_VtxProb;   //!
   TBranch        *b_3DCosPointingAngleF;   //!
   TBranch        *b_3DPointingAngle;   //!
   TBranch        *b_2DCosPointingAngle;   //!
   TBranch        *b_2DPointingAngle;   //!
   TBranch        *b_3DDecayLengthSignificance;   //!
   TBranch        *b_3DDecayLength;   //!
   TBranch        *b_2DDecayLengthSignificance;   //!
   TBranch        *b_2DDecayLength;   //!
   TBranch        *b_zDCASignificanceDaugther1;   //!
   TBranch        *b_xyDCASignificanceDaugther1;   //!
   TBranch        *b_NHitD1;   //!
   TBranch        *b_HighPurityDaugther1;   //!
   TBranch        *b_pTD1;   //!
   TBranch        *b_pTerrD1;   //!
   TBranch        *b_EtaD1;   //!
   TBranch        *b_dedxHarmonic2D1;   //!
   TBranch        *b_zDCASignificanceDaugther2;   //!
   TBranch        *b_xyDCASignificanceDaugther2;   //!
   TBranch        *b_NHitD2;   //!
   TBranch        *b_HighPurityDaugther2;   //!
   TBranch        *b_pTD2;   //!
   TBranch        *b_pTerrD2;   //!
   TBranch        *b_EtaD2;   //!
   TBranch        *b_dedxHarmonic2D2;   //!

   d0data(TChain *tree=0);
   virtual ~d0data();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntries();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual Float_t Pt() const   {return pT;}
   virtual Float_t Y() const  {return y;}
   virtual Float_t Mass() const {return mass;}
   virtual Float_t Mva() const  {return mva;}
   virtual Int_t   nTrkOffline() const  {return Ntrkoffline;}
   virtual Int_t   nPixel() const {return Npixel;}
   virtual Float_t HFSumET() const {return HFsumET;}
   virtual Float_t BestVtxX() const {return bestvtxX;}
   virtual Float_t BestVtxY() const {return bestvtxY;}
   virtual Float_t BestVtxZ() const {return bestvtxZ;}
   virtual Float_t Flavor() const {return flavor;}
   virtual Float_t Eta() const {return eta;}
   virtual Float_t vtxProb() const {return VtxProb;}
   virtual Float_t CosPointingAngle3D() const {return m3DCosPointingAngle;}
   virtual Float_t PointingAngle3D() const {return m3DPointingAngle;}
   virtual Float_t CosPointingAngle2D() const {return m2DCosPointingAngle;}
   virtual Float_t PointingAngle2D() const {return m2DPointingAngle;}
   virtual Float_t DecayLSig3D() const {return m3DDecayLengthSignificance;}
   virtual Float_t DecayL3D() const {return m3DDecayLength;}
   virtual Float_t DecayLSig2D() const {return m2DDecayLengthSignificance;}
   virtual Float_t DecayL2D() const {return m2DDecayLength;}
   virtual Float_t zDCASigD1() const {return zDCASignificanceDaugther1;}
   virtual Float_t xyDCASigD1() const {return xyDCASignificanceDaugther1;}
   virtual Float_t nHitD1() const {return NHitD1;}
   virtual Bool_t  highPurityD1() const {return HighPurityDaugther1;}
   virtual Float_t PtD1() const {return pTD1;}
   virtual Float_t PtErrD1() const {return pTerrD1;}
   virtual Float_t etaD1() const {return EtaD1;}
   virtual Float_t DedxD1() const {return dedxHarmonic2D1;}
   virtual Float_t zDCASigD2() const {return zDCASignificanceDaugther2;}
   virtual Float_t xyDCASigD2() const {return xyDCASignificanceDaugther2;}
   virtual Float_t nHitD2() const {return NHitD2;}
   virtual Bool_t  highPurityD2() const {return HighPurityDaugther2;}
   virtual Float_t PtD2() const {return pTD2;}
   virtual Float_t PtErrD2() const {return pTerrD2;}
   virtual Float_t etaD2() const {return EtaD2;}
   virtual Float_t DedxD2() const {return dedxHarmonic2D2;}
};

#endif

#ifndef d0data_cxx
#define d0data_cxx
d0data::d0data(TChain *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/storage1/users/wl33/D0Trees/Data/Merged_pPbPbpData_MVATree_D0_default_BDTCut03_v1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/storage1/users/wl33/D0Trees/Data/Merged_pPbPbpData_MVATree_D0_default_BDTCut03_v1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/storage1/users/wl33/D0Trees/Data/Merged_pPbPbpData_MVATree_D0_default_BDTCut03_v1.root:/d0ana");
      dir->GetObject("d0data",tree);

   }
   Init(tree);
}

d0data::~d0data()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t d0data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t d0data::GetEntries()
{
   if (!fChain) return 0;
   return fChain->GetEntries();
}

Long64_t d0data::LoadTree(Long64_t entry)
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

void d0data::Init(TChain *tree)
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

   fChain->SetBranchAddress("pT", &pT, &b_pT);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("mva", &mva, &b_mva);
   fChain->SetBranchAddress("Ntrkoffline", &Ntrkoffline, &b_Ntrkoffline);
   fChain->SetBranchAddress("Npixel", &Npixel, &b_Npixel);
   fChain->SetBranchAddress("HFsumET", &HFsumET, &b_HFsumET);
   fChain->SetBranchAddress("bestvtxX", &bestvtxX, &b_bestvtxX);
   fChain->SetBranchAddress("bestvtxY", &bestvtxY, &b_bestvtxY);
   fChain->SetBranchAddress("bestvtxZ", &bestvtxZ, &b_bestvtxZ);
   fChain->SetBranchAddress("flavor", &flavor, &b_flavor);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("VtxProb", &VtxProb, &b_VtxProb);
   fChain->SetBranchAddress("3DCosPointingAngle", &m3DCosPointingAngle, &b_3DCosPointingAngleF);
   fChain->SetBranchAddress("3DPointingAngle", &m3DPointingAngle, &b_3DPointingAngle);
   fChain->SetBranchAddress("2DCosPointingAngle", &m2DCosPointingAngle, &b_2DCosPointingAngle);
   fChain->SetBranchAddress("2DPointingAngle", &m2DPointingAngle, &b_2DPointingAngle);
   fChain->SetBranchAddress("3DDecayLengthSignificance", &m3DDecayLengthSignificance, &b_3DDecayLengthSignificance);
   fChain->SetBranchAddress("3DDecayLength", &m3DDecayLength, &b_3DDecayLength);
   fChain->SetBranchAddress("2DDecayLengthSignificance", &m2DDecayLengthSignificance, &b_2DDecayLengthSignificance);
   fChain->SetBranchAddress("2DDecayLength", &m2DDecayLength, &b_2DDecayLength);
   fChain->SetBranchAddress("zDCASignificanceDaugther1", &zDCASignificanceDaugther1, &b_zDCASignificanceDaugther1);
   fChain->SetBranchAddress("xyDCASignificanceDaugther1", &xyDCASignificanceDaugther1, &b_xyDCASignificanceDaugther1);
   fChain->SetBranchAddress("NHitD1", &NHitD1, &b_NHitD1);
   fChain->SetBranchAddress("HighPurityDaugther1", &HighPurityDaugther1, &b_HighPurityDaugther1);
   fChain->SetBranchAddress("pTD1", &pTD1, &b_pTD1);
   fChain->SetBranchAddress("pTerrD1", &pTerrD1, &b_pTerrD1);
   fChain->SetBranchAddress("EtaD1", &EtaD1, &b_EtaD1);
   fChain->SetBranchAddress("dedxHarmonic2D1", &dedxHarmonic2D1, &b_dedxHarmonic2D1);
   fChain->SetBranchAddress("zDCASignificanceDaugther2", &zDCASignificanceDaugther2, &b_zDCASignificanceDaugther2);
   fChain->SetBranchAddress("xyDCASignificanceDaugther2", &xyDCASignificanceDaugther2, &b_xyDCASignificanceDaugther2);
   fChain->SetBranchAddress("NHitD2", &NHitD2, &b_NHitD2);
   fChain->SetBranchAddress("HighPurityDaugther2", &HighPurityDaugther2, &b_HighPurityDaugther2);
   fChain->SetBranchAddress("pTD2", &pTD2, &b_pTD2);
   fChain->SetBranchAddress("pTerrD2", &pTerrD2, &b_pTerrD2);
   fChain->SetBranchAddress("EtaD2", &EtaD2, &b_EtaD2);
   fChain->SetBranchAddress("dedxHarmonic2D2", &dedxHarmonic2D2, &b_dedxHarmonic2D2);
   Notify();
}

Bool_t d0data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TChain in a TChain or when when a new TChain
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void d0data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t d0data::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifndef d0data_cxx

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 11 15:28:35 2019 by ROOT version 6.06/01
// from TChain d0tree/d0tree
// found on file: /storage1/users/wl33/D0Trees/Data/Merged_pPbPbpData_MVATree_D0_default_BDTCut03_v1.root
//////////////////////////////////////////////////////////

#ifndef d0tree_h
#define d0tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TChain if any.

class d0tree {
public :
   TChain          *fChain;   //!pointer to the analyzed TChain or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TChain if any.

   // Declaration of leaf types

   // List of branches

   d0tree(TChain *tree=0);
   virtual ~d0tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntries();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual Float_t Pt() const   {return -1;}
   virtual Float_t Y() const  {return -1;}
   virtual Float_t Mass() const {return -1;}
   virtual Float_t Mva() const  {return -1;}
   virtual Int_t nTrkOffline() const  {return -1;}
   virtual Int_t nPixel() const {return -1;}
   virtual Float_t HFSumET() const {return -1;}
   virtual Float_t BestVtxX() const {return -1;}
   virtual Float_t BestVtxY() const {return -1;}
   virtual Float_t BestVtxZ() const {return -1;}
   virtual Float_t Flavor() const {return -1;}
   virtual Float_t Eta() const {return -1;}
   virtual Float_t vtxProb() const {return -1;}
   virtual Float_t CosPointingAngle3D() const {return -1;}
   virtual Float_t PointingAngle3D() const {return -1;}
   virtual Float_t CosPointingAngle2D() const {return -1;}
   virtual Float_t PointingAngle2D() const {return -1;}
   virtual Float_t DecayLSig3D() const {return -1;}
   virtual Float_t DecayL3D() const {return -1;}
   virtual Float_t DecayLSig2D() const {return -1;}
   virtual Float_t DecayL2D() const {return -1;}
   virtual Float_t zDCASigD1() const {return -1;}
   virtual Float_t xyDCASigD1() const {return -1;}
   virtual Float_t nHitD1() const {return -1;}
   virtual Bool_t  highPurityD1() const {return -1;}
   virtual Float_t PtD1() const {return -1;}
   virtual Float_t PtErrD1() const {return -1;}
   virtual Float_t etaD1() const {return -1;}
   virtual Float_t DedxD1() const {return -1;}
   virtual Float_t zDCASigD2() const {return -1;}
   virtual Float_t xyDCASigD2() const {return -1;}
   virtual Float_t nHitD2() const {return -1;}
   virtual Bool_t  highPurityD2() const {return -1;}
   virtual Float_t PtD2() const {return -1;}
   virtual Float_t PtErrD2() const {return -1;}
   virtual Float_t etaD2() const {return -1;}
   virtual Float_t DedxD2() const {return -1;}
   virtual Bool_t  IsSwap() const {return false;}
   virtual Int_t   IdMom_Reco() const {return -1;}
   virtual Bool_t  MatchGEN() const {return false;}
};

#endif

#ifndef d0tree_cxx
#define d0tree_cxx
d0tree::d0tree(TChain *tree)  
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
}

d0tree::~d0tree()
{
}

Int_t d0tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
 return 0;
}

Int_t d0tree::GetEntries()
{
   return 0;
}

Long64_t d0tree::LoadTree(Long64_t entry)
{
   return 0;
}

void d0tree::Init(TChain *tree)
{
}

Bool_t d0tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TChain in a TChain or when when a new TChain
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void d0tree::Show(Long64_t entry)
{
}
Int_t d0tree::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifndef d0tree_cxx

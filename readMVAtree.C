#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TMath.h"
#include "TTree.h"

#include "d0data.h"
#include "d0mc.h"
#include "myAnaConsts.h"

void loopTree(d0tree* , std::map<std::string, TH1*>&);
void loopTree(d0tree* , std::map<std::string, TH2*>&, bool isMC);
void loopTree(d0tree* , std::map<std::string, TH3*>&, bool isMC);
bool passKinematicalCut(d0tree*);

void readMVAtree
(int mode = 0, std::string mcpd0file = "", std::string mcnpd0file = "",
 std::string datafile = "")
{
   TH1::SetDefaultSumw2();

   //std::map<std::string, TH1*> prompt;
   //prompt["dca3D"] = new TH1D("hDca3D_pd0", "hDca3D_pd0", 20, 0, 0.04);
   //prompt["pseudo_decayL"] = new TH1D("hDL_pd0", "hpseudoDL_pd0", 40, -0.12, 0.32);
   //std::map<std::string, TH1*> nonprompt;
   //nonprompt["dca3D"] = new TH1D("hDca3D_npd0", "hDca3D_npd0", 20, 0, 0.04);
   //nonprompt["pseudo_decayL"] = new TH1D("hDL_npd0", "hpseudoDL_npd0", 40, -0.12, 0.32);

   /*
   std::map<std::string, TH2*> promptDCA;
   promptDCA["h_match_unswap"] = new TH2D("hDcaVsMassPD0", "hDcaVsMassPD0", 60, 1.7, 2.0,  60, 0, 0.06);
   promptDCA["h_match_all"] = new TH2D("hDcaVsMassPD0_All", "hDcaVsMassPD0_All", 60, 1.7, 2.0,  60, 0, 0.06);

   std::map<std::string, TH2*> nonPromptDCA;
   nonPromptDCA["h_match_unswap"] = new TH2D("hDcaVsMassNPD0", "hDcaVsMassNPD0", 60, 1.7, 2.0,  60, 0, 0.06);
   nonPromptDCA["h_match_all"] = new TH2D("hDcaVsMassNPD0_All", "hDcaVsMassNPD0_All", 60, 1.7, 2.0,  60, 0, 0.06);

   std::map<std::string, TH2*> dataDCA;
   dataDCA["hdata"] = new TH2D("hDcaVsMassDataD0", "hDcaVsMassDataD0", 60, 1.7, 2.0, 60, 0., 0.06);
   */
   std::map<std::string, TH3*> promptDCA;
   promptDCA["h_match_unswap"] = new TH3D("hDcaVsMassAndMvaPD0", "hDcaVsMassAndMvaPD0", 60, 1.7, 2.0, 20, 0.4, 0.8, ana::nDca, ana::dcaMin, ana::dcaMax);
   promptDCA["h_match_all"] = new TH3D("hDcaVsMassAndMvaPD0_All", "hDcaVsMassAndMvaPD0_All", 60, 1.7, 2.0, 20, 0.4, 0.8, ana::nDca, ana::dcaMin, ana::dcaMax);

   std::map<std::string, TH3*> nonPromptDCA;
   nonPromptDCA["h_match_unswap"] = new TH3D("hDcaVsMassAndMvaNPD0", "hDcaVsMassAndMvaNPD0", 60, 1.7, 2.0, 20, 0.4, 0.8, ana::nDca, ana::dcaMin, ana::dcaMax);
   nonPromptDCA["h_match_all"] = new TH3D("hDcaVsMassAndMvaNPD0_All", "hDcaVsMassAndMvaNPD0_All", 60, 1.7, 2.0, 20, 0.4, 0.8, ana::nDca, ana::dcaMin, ana::dcaMax);

   std::map<std::string, TH3*> dataDCA;
   dataDCA["hdata"] = new TH3D("hDcaVsMassAndMvaDataD0", "hDcaVsMassAndMvaDataD0", 60, 1.7, 2.0, 20, 0.4, 0.8, ana::nDca, ana::dcaMin, ana::dcaMax);

   // open mc files
   const std::string dir = "/storage1/users/wl33/D0Trees/";
   mcpd0file = "MC/"
       "Merge_PromptD0_pPbMC_MVATree_signal_combined_BDT_v1.root";
   mcnpd0file = "MC/"
       "Merge_NonPromptD0_pPbMC_MVATree_signal_combined_BDT_v1.root";
   datafile = "Data/"
       "Merged_pPbPbpData_MVATree_D0_default_BDTCut03_v1.root";
   TFile* f1 = new TFile(Form("%s%s", dir.c_str(), 
       mcpd0file.c_str()), "read");
   TFile* f2 = new TFile(Form("%s%s", dir.c_str(), 
       mcnpd0file.c_str()), "read");
   TFile* f3 = new TFile(Form("%s%s", dir.c_str(), 
       datafile.c_str()), "read");
   if(f1->IsOpen()) std::cout << "opened mcpd0 file successfully" << std::endl;
   if(f2->IsOpen()) std::cout << "opened mcnpd0 file successfully" << std::endl;
   if(f3->IsOpen()) std::cout << "opened mcnpd0 file successfully" << std::endl;

   // declare histograms

   // read trees and fill hisotgrams
   f1->cd();
   gROOT->ProcessLine(".pwd");

   std::string mcpd0treeptr = 
                   Form("%s_mc/VertexCompositeNtuple", ana::whichtree[mode].c_str());
   TTree* mcpd0tree = (TTree*) f1->Get(mcpd0treeptr.c_str());
   d0mc *mcpd0 = new d0mc(mcpd0tree);

   std::cout << "started filling histograms of mc promptd0" << std::endl;
   loopTree(mcpd0, promptDCA, true);
   std::cout << "ended filling histograms of mc promptd0" << std::endl;

   delete mcpd0;

   f2->cd();
   gROOT->ProcessLine(".pwd");

   std::string mcnpd0treeptr = 
                   Form("%s_mc/VertexCompositeNtuple", ana::whichtree[mode].c_str());
   TTree* mcnpd0tree = (TTree*) f2->Get(mcnpd0treeptr.c_str());
   d0mc *mcnpd0 = new d0mc(mcnpd0tree);

   std::cout << "started filling histograms of mc non-prompt d0" << std::endl;
   loopTree(mcnpd0, nonPromptDCA, true);
   std::cout << "ended filling histograms of mc non-prompt d0" << std::endl;

   delete mcnpd0;

   f3->cd();
   gROOT->ProcessLine(".pwd");

   std::string datad0treeptr = 
                   Form("%s/VertexCompositeNtuple", ana::whichtree[mode].c_str());
   TTree* datad0tree = (TTree*) f3->Get(datad0treeptr.c_str());
   d0data *datad0 = new d0data(datad0tree);


   std::cout << "started filling histograms of data d0" << std::endl;
   loopTree(datad0, dataDCA, false);
   std::cout << "ended filling histograms of data d0" << std::endl;

   delete datad0;

   
   TFile f4(Form("%s_hists.root", ana::whichtree[mode].c_str()), "recreate");
   std::cout << "ready to write ouput file" << std::endl;
   f4.cd();

   for(const auto& h : promptDCA) {h.second->Write(); delete h.second;}
   for(const auto& h : nonPromptDCA) {h.second->Write(); delete h.second;}
   for(const auto& h : dataDCA) {h.second->Write(); delete h.second;}

   std::cout << "ended writing output file" << std::endl;
  
}

void loopTree(d0tree* d0, std::map<std::string, TH1*>& d0hists)
{
   int nentries = d0->GetEntries();
   std::cout << "total entries: " << nentries << std::endl;
   //nentries = 10000;
   for(int ientry=0; ientry<nentries; ientry++){
      d0->GetEntry(ientry);
      if(!d0->MatchGEN()) continue;
      if(d0->IsSwap()) continue;
      if(!passKinematicalCut(d0)) continue;
      double dca3D = d0->DecayL3D() * std::sin(d0->PointingAngle3D());
      d0hists["dca3D"]->Fill(dca3D);
      double p = d0->Pt() * std::cosh(d0->Eta());
      double scalar_product = d0->DecayL3D() * d0->CosPointingAngle3D() * ana::D0_mass / p;
      d0hists["pseudo_decayL"]->Fill(scalar_product);
   }
}

void loopTree(d0tree* d0, std::map<std::string, TH2*>& d0hists, bool isMC)
{
   int nentries = d0->GetEntries();
   std::cout << "total entries: " << nentries << std::endl;
   //nentries = 10000;
   for(int ientry=0; ientry<nentries; ientry++){
      d0->GetEntry(ientry);
      if(!passKinematicalCut(d0)) continue;
      double dca3D = d0->DecayL3D() * std::sin(d0->PointingAngle3D());
      double mass = d0->Mass();
      if(isMC){
         if(!d0->MatchGEN()) continue;
         if(!d0->IsSwap()) d0hists["h_match_unswap"]->Fill(mass, dca3D);
         d0hists["h_match_all"]->Fill(mass, dca3D);
      } else {

         d0hists["hdata"]->Fill(mass, dca3D);
      }
   }
}

void loopTree(d0tree* d0, std::map<std::string, TH3*>& d0hists, bool isMC)
{
   int nentries = d0->GetEntries();
   std::cout << "total entries: " << nentries << std::endl;
   //nentries = 10000;
   for(int ientry=0; ientry<nentries; ientry++){
      d0->GetEntry(ientry);
      if(!passKinematicalCut(d0)) continue;
      double dca3D = d0->DecayL3D() * std::sin(d0->PointingAngle3D());
      double mass = d0->Mass();
      double mva = d0->Mva();
      if(isMC){
         if(!d0->MatchGEN()) continue;
         if(!d0->IsSwap()) d0hists["h_match_unswap"]->Fill(mass, mva, dca3D);
         d0hists["h_match_all"]->Fill(mass, mva, dca3D);
      } else {
         d0hists["hdata"]->Fill(mass, mva, dca3D);
      }
   }
}

bool passKinematicalCut(d0tree* d0)
{
   //bool passPt = d0->Pt() > 3.5 && d0->Pt() < 4.2;
   bool passPt = d0->Pt() > 4. && d0->Pt() < 5.;
   bool passY = d0->Y() > -1 && d0->Y() < 1;
   bool passPointingAngle = std::fabs(d0->PointingAngle3D()) < 1;
   //passPointingAngle = true;
   bool passTrkEta = std::fabs(d0->etaD1()) < 2.4 && std::fabs(d0->etaD2()) < 2.4;
   bool passTrkPt = d0->PtD1() > 0.7 && d0->PtD2() > 0.7;
   bool passTrkPtErr = d0->PtErrD1()/d0->PtD1() < 0.1 && d0->PtErrD2()/d0->PtD2() < 0.1;
   bool passTrkPurity = d0->highPurityD1() && d0->highPurityD2();
   bool passTrkNhits = d0->nHitD1() >=11 && d0->nHitD2() >=11; 

   //bool passMVA = d0->Mva() > 0.58;
   //passMVA = true;

   if(passPt && passY && passPointingAngle && passTrkEta 
         && passTrkPt && passTrkPtErr && passTrkPurity && passTrkNhits
         //&& passMVA
         ) return true;
   return false;
}

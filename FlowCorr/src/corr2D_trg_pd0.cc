#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <stdlib.h>

#include "TFile.h"
#include "TChain.h"
#include "Event.h"
#include "TTree.h"
#include "TFileCollection.h"
#include "TCollection.h"
#include "THashList.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TVector3.h"
#include "TString.h"

#include "myAnaConsts.h"

using namespace std;

static bool isPromptD0 = true;

void setBranchStatus(Event*);
bool checkBranchStatus(Event*);

bool passGoodTrack(Event*, const unsigned int&);
inline bool passGoodVtx(Event* event);
inline bool passD0Selections(const int&, Event*, const int&, const bool&);
bool passD0PreSelections(Event*, const int&);
bool passD0KinematicCuts(Event*, const int&);
bool passD0MVA(const int&, Event*, const int&, const bool&);
inline bool passNtrkoffline(const double&, const int&);

vector<double> setPtBin(string& dataset);

// par0, main
// par1, datalist
// par2, dataset
// par3, efficiency hists
// par4, output dir

int main(int argc, char** argv)
{
   TH1::SetDefaultSumw2(true);

   if(argc!=5) {
      std::cerr << "The number of arguments is wrong" << std::endl;
      return -1;
   }
   
   string datalist(argv[1]);
   std::cout << datalist << std::endl;

   string dataset(argv[2]);

   const int dataset_trigger = ana::Get_Trigger(dataset);
   if(dataset_trigger<0 || !ana::isHM_PD0_DataSet(dataset)){
      cerr << "wrong dataset name" << endl;
      cout << "name should be:\n" 
         << "PAHM1-6\n"
         << "PPHM_2  //mult > 100\n"
         << "// means comments"
         << endl;
      return -1;
   }

   const vector<double> ptbin = setPtBin(dataset);
   const int nPt = ptbin.size()-1;

   if(!ptbin.size()) {
      cerr << "wrong dataset name" << endl;
      cout << "name should be:\n" 
         << "PAHM1-6\n"
         << "PPHM_2  //mult > 100\n"
         << "// means comments"
         << endl;
      return -1;
   }

   TFile f_eff(argv[3]);
   if(!f_eff.IsOpen()){
      cerr << "cannot find the file of efficiency" << endl;
      return -1;
   }
   TH2D* h_eff = nullptr;
   f_eff.GetObject("hEff", h_eff);
   if(!h_eff) {
      cerr << "cannot find the histogram of efficiency" << endl;
      return -1;
   }

   TChain *chain_d0 = new TChain("d0ana/VertexCompositeNtuple");
   TChain *chain_tracks = new TChain("track_ana/trackTree");

   TFileCollection* fcData = new TFileCollection(datalist.c_str(), "", datalist.c_str());

   chain_d0->AddFileInfoList(fcData->GetList());

   chain_tracks->AddFileInfoList(fcData->GetList());

   Event* evt = new Event(chain_d0, chain_tracks);
   setBranchStatus(evt);
   if(!checkBranchStatus(evt)) return -1;

   // declare hists
   TH1D* hMult;
   TH1D* hMult_ass;

   TH1D* hKET_D0[nPt];
   TH1D* hPt_D0[nPt];
   TH1D* hEta_D0[nPt];
   TH1D* hRapidity_D0[nPt];

   TH1D* hMass_D0[ana::nMass][nPt];

   TH2D* hNtrkofflineVsNtrkgood;

   TH3D* hDcaVsMassAndMva[nPt];

   TH1D* hMult_raw_D0[ana::nMass][nPt];
   TH1D* hMult_eff_D0[ana::nMass][nPt];
   TH2D* hSignal_D0[ana::nMass][nPt];
   TH2D* hBackground_D0[ana::nMass][nPt];

   hMult = new TH1D("hMult", "", 600, 0, 600);
   hMult_ass = new TH1D("hMult_ass", "", 600, 0, 600);

   hNtrkofflineVsNtrkgood = new TH2D("hNtrkofflineVsNtrkgood", "", 300, 0, 300, 300, 0, 300);

   for(int ipt=0; ipt<nPt; ipt++){
      hKET_D0[ipt] = new TH1D(Form("hKET_pt%d", ipt), "", 3000, 0, 30);
      hPt_D0[ipt] = new TH1D(Form("hPt_pt%d", ipt), "", 3000, 0, 30);
      hEta_D0[ipt] = new TH1D(Form("hEta_pt%d", ipt), "", 24, -2.4, 2.4);
      hRapidity_D0[ipt] = new TH1D(Form("hRapidity_pt%d", ipt), "", 24, -2.4, 2.4);
      hDcaVsMassAndMva[ipt] = new TH3D(Form("hDcaVsMassAndMva_pt%d", ipt), "", 60, 1.7, 2.0, 100, -0.3, 0.7, 800, 0, 0.08);
      for(int imass=0; imass<ana::nMass; imass++){
         hMass_D0[imass][ipt] = new TH1D(Form("hMassD0_mass%d_pt%d", imass, ipt),
                  "", 200, 1.5, 2.5);
         hMult_raw_D0[imass][ipt] = new TH1D(Form("hMult_raw_D0_mass%d_pt%d", imass, ipt),
                  "", 50, 0, 50);
         hMult_eff_D0[imass][ipt] = new TH1D(Form("hMult_eff_D0_mass%d_pt%d", imass, ipt),
                  "", 50, 0, 50);
         hSignal_D0[imass][ipt] = new TH2D(Form("hSignal_mass%d_pt%d", imass, ipt),
                  "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
                  ana::nPhiBin, ana::phiBegin, ana::phiEnd);
         hBackground_D0[imass][ipt] = new TH2D(Form("hBackground_mass%d_pt%d", imass, ipt),
                  "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
                  ana::nPhiBin, ana::phiBegin, ana::phiEnd);
      }
   }

   // declare vectors
   vector<TVector3> pVect_trg_d0[ana::nMass][nPt];
   vector<double>   effVect_trg_d0[ana::nMass][nPt];

   vector<TVector3> pVect_dau1_d0[ana::nMass][nPt];

   vector<TVector3> pVect_dau2_d0[ana::nMass][nPt];

   vector<int>      indexVect_d0[ana::nMass][nPt];

   vector<TVector3> pVect_ass;
   list<vector<TVector3>> pVectList_ass[ana::nZ_Vtx_];
   vector<float> effVect_ass;
   list<vector<float>> effVectList_ass[ana::nZ_Vtx_];

   // start timing
   TStopwatch ts;
   ts.Start();

   // loop
   //
   // temporary vectors
   TVector3 p_dau1(0, 0, 0), p_dau2(0, 0, 0), p_d0(0, 0, 0), p_ass(0, 0, 0);

   std::cout << evt->GetEntries() << std::endl;
   long int nentries = evt->GetEntries();
   int percent = 0;
   long int skip = 0;
   for(long int ientry=0; ientry<nentries; ientry++){

      int current = (int)(ientry+1)*100/nentries;
      if( current == 5*percent){
         std::cout << 5*percent++ << "\% percents completed" << std::endl;
      }
      auto bytes = evt->GetEntry(ientry);
      if(bytes == 0 || bytes == -1) {
      // std::cout << "vertex unmatched " << ientry << std::endl;
         skip++;
         continue;
      }

      if(!passGoodVtx(evt)) continue;

      auto iz = ana::findZVtxBin(evt->BestVtxZ());
      if(iz == -1) continue;

      // Ntrkoffline cut
      if(!passNtrkoffline(evt->nTrkOffline(), dataset_trigger)) continue;

      // count number of good tracks per event
      unsigned int nMult_ass_good = 0;
      for(unsigned int itrack=0; itrack<evt->CandSizeTrk(); itrack++){
         // assume all tracks in TTree are good, since error of dz and dxy are not available
         if(passGoodTrack(evt, itrack)) 
            nMult_ass_good++;
      }
      hMult->Fill(evt->nTrkOffline());

      hNtrkofflineVsNtrkgood->Fill(nMult_ass_good, evt->nTrkOffline());

      for(int id0=0; id0<evt->CandSize(); id0++){
         int imass = ana::findMassBin(evt->Mass(id0));
         int ipt = ana::findPtBin(evt->Pt(id0), ptbin);

         if(imass == -1) continue;
         if(ipt == -1) continue;

         if(!passD0Selections(dataset_trigger, evt, id0, isPromptD0)) continue;

         double effks = h_eff->FindBin(evt->Pt(id0), evt->Y(id0));

         hDcaVsMassAndMva[ipt]->Fill(evt->Mass(id0), evt->Mva(id0), 0.);

         p_dau1.SetPtEtaPhi(evt->PtD1(id0), evt->etaD1(id0), evt->phiD1(id0));
         p_dau2.SetPtEtaPhi(evt->PtD2(id0), evt->etaD2(id0), evt->phiD2(id0));
         p_d0 = p_dau1 + p_dau2;

         double KET = sqrt(pow(evt->Mass(id0), 2) + pow(evt->Pt(id0), 2)
                  - evt->Mass(id0));

         hMass_D0[imass][ipt]->Fill(evt->Mass(id0), 1./effks);
         hPt_D0[ipt]->Fill(evt->Pt(id0), 1./effks);
         hEta_D0[ipt]->Fill(p_d0.Eta(), 1./effks);
         hRapidity_D0[ipt]->Fill(evt->Y(id0), 1./effks);
         hKET_D0[ipt]->Fill(KET, 1./effks);

         pVect_trg_d0[imass][ipt].push_back(p_d0);
         effVect_trg_d0[imass][ipt].push_back(effks);
         pVect_dau1_d0[imass][ipt].push_back(p_dau1);
         pVect_dau2_d0[imass][ipt].push_back(p_dau2);
         indexVect_d0[imass][ipt].push_back(id0);
      }

      for(unsigned int itrack=0; itrack<evt->CandSizeTrk(); itrack++){
         // some cuts are not available
         bool passDzErr = true;
         bool passDxyErr = true;
         bool passTrkPurity = true;
         bool passPt = evt->PtTrk(itrack) > ana::ptMin_ass_ && evt->PtTrk(itrack) < ana::ptMax_ass_;
         bool passPtError = true;
         bool passEta = evt->EtaTrk(itrack) > ana::etaMin_ass_ && evt->EtaTrk(itrack) < ana::etaMax_ass_;
         bool pass_ass_ = passDzErr && passDxyErr && passTrkPurity && passPt &&
                           passPtError && passEta;

         if(pass_ass_) {
            p_ass.SetPtEtaPhi(evt->PtTrk(itrack), evt->EtaTrk(itrack), evt->PhiTrk(itrack));
            pVect_ass.push_back(p_ass);
            effVect_ass.push_back(evt->WeightTrk(itrack));
         }
      }

      // calculate signal
      unsigned int nMult_ass = (unsigned int) pVect_ass.size();
      hMult_ass->Fill(nMult_ass);
      
      unsigned int nMult_trg_raw_d0[ana::nMass][nPt]; // Ntrig for mass & pt bins
      double nMult_trg_eff_d0[ana::nMass][nPt]; // eff corrected Ntrig for mass & pt bins

      for(int imass=0; imass<ana::nMass; imass++){
         for(int ipt=0; ipt<nPt; ipt++){
            nMult_trg_raw_d0[imass][ipt] = 0;
            nMult_trg_eff_d0[imass][ipt] = 0;
         }
      }

      for(int imass=0; imass<ana::nMass; imass++){
         for(int ipt=0; ipt<nPt; ipt++){
            unsigned int nMult_trg_d0 = (unsigned int) pVect_trg_d0[imass][ipt].size();
            for(unsigned int id0=0; id0<nMult_trg_d0; id0++){
               if(ipt != ana::findPtBin(pVect_trg_d0[imass][ipt].at(id0).Pt(), ptbin)) 
                  std::cout << "pT bin error" << std::endl;
               double effks = effVect_trg_d0[imass][ipt].at(id0);
               //int index = indexVect_d0[imass][ipt].at(id0);
               nMult_trg_raw_d0[imass][ipt] += 1;
               nMult_trg_eff_d0[imass][ipt] += 1./effks;
            }

            hMult_raw_D0[imass][ipt]->Fill(nMult_trg_raw_d0[imass][ipt]);
            hMult_eff_D0[imass][ipt]->Fill(nMult_trg_eff_d0[imass][ipt]);

            for(unsigned int id0=0; id0<nMult_trg_d0; id0++){
               if(ipt!=ana::findPtBin(pVect_trg_d0[imass][ipt].at(id0).Pt(), ptbin)) 
                  std::cerr << "pT bin error" << std::endl;
               double effks = effVect_trg_d0[imass][ipt].at(id0);
               for(unsigned int iass=0; iass<nMult_ass; iass++){
                  TVector3 pvector_ass = pVect_ass.at(iass);
                  double effweight_ass = effVect_ass.at(iass);
                  if(ana::rejectDaughter_){
                     if(fabs(pvector_ass.Eta() - pVect_dau1_d0[imass][ipt].at(id0).Eta())<0.03
                           && fabs(pvector_ass.DeltaPhi(pVect_dau1_d0[imass][ipt].at(id0)))<0.03)
                        continue;
                     if(fabs(pvector_ass.Eta() - pVect_dau2_d0[imass][ipt].at(id0).Eta())<0.03
                           && fabs(pvector_ass.DeltaPhi(pVect_dau2_d0[imass][ipt].at(id0)))<0.03)
                        continue;
                  }

                  double deltaEta = pvector_ass.Eta() - pVect_trg_d0[imass][ipt].at(id0).Eta();
                  double deltaPhi = pvector_ass.DeltaPhi(pVect_trg_d0[imass][ipt].at(id0));
                  if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

                  hSignal_D0[imass][ipt]->Fill(deltaEta, deltaPhi, 
                           1./nMult_trg_eff_d0[imass][ipt]/effks/effweight_ass);
               }
            }
         }
      }

      // mixed events
      // begin filling background
      // if having enough, erase the first one, fill the current, fill the background
      if(pVectList_ass[iz].size() == ana::nMixedEvts){
         for(int imass=0; imass<ana::nMass; imass++){
            for(int ipt=0; ipt<nPt; ipt++){
               unsigned int nMult_trg_d0 = (unsigned int) pVect_trg_d0[imass][ipt].size();
               for(unsigned int id0=0; id0<nMult_trg_d0; id0++){
                  double effks = effVect_trg_d0[imass][ipt].at(id0);
                  // simultaneously read momentum and efficiency
                  auto ievt_eff_ass = effVectList_ass[iz].begin();
                  for(auto& ievt_p_ass : pVectList_ass[iz]){

                     if(ipt!=ana::findPtBin(pVect_trg_d0[imass][ipt].at(id0).Pt(), ptbin)) 
                        std::cerr << "pT bin error" << std::endl;
                     unsigned int n_ass = ievt_p_ass.size();
                     if(n_ass!=ievt_eff_ass->size()){
                        cout << "wrong" << endl;
                        break;
                     }
                     for(unsigned int iass=0; iass<n_ass; iass++){
                        TVector3 pvector_ass = ievt_p_ass.at(iass);
                        double effweight_ass = ievt_eff_ass->at(iass);

                        double deltaEta = pvector_ass.Eta() - pVect_trg_d0[imass][ipt].at(id0).Eta();
                        double deltaPhi = pvector_ass.DeltaPhi(pVect_trg_d0[imass][ipt].at(id0));
                        if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

                        hBackground_D0[imass][ipt]->Fill(deltaEta, deltaPhi, 
                              1./nMult_trg_eff_d0[imass][ipt]/effks/effweight_ass);
                     }
                     ievt_eff_ass++;
                  }
               }
            }
         }
      }

      // clear all and refresh list/buffer
      if(pVectList_ass[iz].size() == ana::nMixedEvts){
         // update the list/buffer
         pVectList_ass[iz].erase(pVectList_ass[iz].begin());
         pVectList_ass[iz].push_back(pVect_ass);
         effVectList_ass[iz].erase(effVectList_ass[iz].begin());
         effVectList_ass[iz].push_back(effVect_ass);
      }else{
      // if not having enough events, keep filling the list of momentum vector of associated particles
         // fill the list/buffer
         pVectList_ass[iz].push_back(pVect_ass);
         effVectList_ass[iz].push_back(effVect_ass);
      }

      for(int imass=0; imass<ana::nMass; imass++){
         for(int ipt=0; ipt<nPt; ipt++){
            pVect_trg_d0[imass][ipt].clear();
            effVect_trg_d0[imass][ipt].clear();
            pVect_dau1_d0[imass][ipt].clear();
            pVect_dau2_d0[imass][ipt].clear();
            indexVect_d0[imass][ipt].clear();
         }
      }
      pVect_ass.clear();
      effVect_ass.clear();
   }
   std::cout << "completed loop" << std::endl;
   std::cout << skip << " events are skipped" << std::endl;

   ts.Stop();
   ts.Print();

   TString outName; 
   if(isPromptD0){
      string prefix(argv[4]);
      size_t found_slash = datalist.find("/");
      while(found_slash!=std::string::npos) {
         datalist.replace(found_slash, 1, "_");
         found_slash = datalist.find("/");
      }

      if(prefix.size())
         outName = TString::Format("%s/fout_%s_d0ana_HM_%.1f.root", prefix.c_str(), datalist.c_str(), ana::d0_y_max_);
      else
         outName = TString::Format("fout_%s_d0ana_HM_%.1f.root", datalist.c_str(), ana::d0_y_max_);
   } 
   else return -1;

   TFile fout(outName.Data(), "recreate");
   fout.cd();

   // start writing output
   hMult->Write();
   hMult_ass->Write();
   hNtrkofflineVsNtrkgood->Write();
   for(int ipt=0; ipt<nPt; ipt++){
      hKET_D0[ipt]->Write(); 
      hPt_D0[ipt]->Write();
      hEta_D0[ipt]->Write();
      hRapidity_D0[ipt]->Write();
      hDcaVsMassAndMva[ipt]->Write();
      for(int imass=0; imass<ana::nMass; imass++){
         hMass_D0[imass][ipt]->Write();
         hMult_raw_D0[imass][ipt]->Write();
         hMult_eff_D0[imass][ipt]->Write();
         hSignal_D0[imass][ipt]->Write();
         hBackground_D0[imass][ipt]->Write();
      }
   }

   delete hMult;
   delete hMult_ass;
   delete hNtrkofflineVsNtrkgood;
   for(int ipt=0; ipt<nPt; ipt++){
      delete hKET_D0[ipt]; 
      delete hPt_D0[ipt];
      delete hEta_D0[ipt];
      delete hRapidity_D0[ipt];
      delete hDcaVsMassAndMva[ipt];
      for(int imass=0; imass<ana::nMass; imass++){
         delete hMass_D0[imass][ipt];
         delete hMult_raw_D0[imass][ipt];
         delete hMult_eff_D0[imass][ipt];
         delete hSignal_D0[imass][ipt];
         delete hBackground_D0[imass][ipt];
      }
   }

   delete evt;

   return 0;
}

void setBranchStatus(Event* evt)
{
   evt->SetBranchStatus("Ntrkoffline", 1);
   evt->SetBranchStatus("candSize", 1);
   evt->SetBranchStatus("pT", 1);
   evt->SetBranchStatus("mass", 1);
   evt->SetBranchStatus("mva", 1);
   evt->SetBranchStatus("y", 1);
   evt->SetBranchStatus("3DPointingAngle", 1);
   evt->SetBranchStatus("3DDecayLength", 1);
   evt->SetBranchStatus("3DDecayLengthSignificance", 1);
   evt->SetBranchStatus("*D1*", 1);
   evt->SetBranchStatus("*D2*", 1);
   evt->SetBranchStatus("*Daugther1", 1); // mistype daughter... the writer of TTree
   evt->SetBranchStatus("*Daugther2", 1);
   evt->SetBranchStatus("dedx*", 0);

   evt->SetBranchStatus("tracks.candSizeTRK", 1);
   evt->SetBranchStatus("tracks.pTTRK", 1);
   evt->SetBranchStatus("tracks.etaTRK", 1);
   evt->SetBranchStatus("tracks.phiTRK", 1);
   evt->SetBranchStatus("tracks.weightTRK", 1);
}

bool checkBranchStatus(Event* event)
{
   bool check =
      event->GetBranchStatus("candSize")&&

      event->GetBranchStatus("pT");
      event->GetBranchStatus("mass")&&
      event->GetBranchStatus("mva")&&
      event->GetBranchStatus("y")&&
      event->GetBranchStatus("3DPointingAngle")&&
      event->GetBranchStatus("3DDecayLength")&&
      event->GetBranchStatus("3DDecayLengthSignificance")&&

      event->GetBranchStatus("pTD1")&&
      event->GetBranchStatus("pTerrD1")&&
      event->GetBranchStatus("EtaD1")&&
      event->GetBranchStatus("PhiD1")&&
      event->GetBranchStatus("zDCASignificanceDaugther1")&&
      event->GetBranchStatus("xyDCASignificanceDaugther1")&&
      event->GetBranchStatus("NHitD1")&&
      event->GetBranchStatus("HighPurityDaugther1")&&

      event->GetBranchStatus("pTD2")&&
      event->GetBranchStatus("pTerrD2")&&
      event->GetBranchStatus("EtaD2")&&
      event->GetBranchStatus("PhiD2")&&
      event->GetBranchStatus("zDCASignificanceDaugther2")&&
      event->GetBranchStatus("xyDCASignificanceDaugther2")&&
      event->GetBranchStatus("NHitD2")&&
      event->GetBranchStatus("HighPurityDaugther2")&&

      event->GetBranchStatus("tracks.bestvtxX")&&
      event->GetBranchStatus("tracks.bestvtxY")&&
      event->GetBranchStatus("tracks.bestvtxZ")&&
      event->GetBranchStatus("tracks.candSizeTRK")&&
      event->GetBranchStatus("tracks.pTTRK")&&
      event->GetBranchStatus("tracks.etaTRK")&&
      event->GetBranchStatus("tracks.phiTRK")&&
      event->GetBranchStatus("tracks.weightTRK");

   return check;
}

inline bool passD0Selections(const int& trigger, Event* event, const int& icand, const bool& isPromptD0)
{
   if(!passD0PreSelections(event, icand)) return false;
   if(!passD0KinematicCuts(event, icand)) return false;
   if(!passD0MVA(trigger, event, icand, isPromptD0)) return false;
   return true;
}

bool passD0PreSelections(Event* event, const int& icand)
{
   bool passPointingAngle = std::fabs(event->PointingAngle3D(icand)) < 1;
   bool passTrkEta = std::fabs(event->etaD1(icand)) < ana::d0_dau_abs_eta_max_ && std::fabs(event->etaD2(icand)) < ana::d0_dau_abs_eta_max_;
   bool passTrkPt = event->PtD1(icand) > ana::d0_dau_pt_min_ && event->PtD2(icand) > ana::d0_dau_pt_min_;
   bool passTrkPtErr = event->PtErrD1(icand)/event->PtD1(icand) < ana::d0_dau_pterr_max_ && event->PtErrD2(icand)/event->PtD2(icand) < ana::d0_dau_pterr_max_;
   bool passTrkPurity = event->highPurityD1(icand) && event->highPurityD2(icand);
   bool passTrkNhits = event->nHitD1(icand) >=ana::d0_dau_nhit_min_ && event->nHitD2(icand) >= ana::d0_dau_nhit_min_; 
   bool passDeltaEta = std::fabs(event->etaD1(icand) - event->etaD2(icand)) < 1;

   if(passPointingAngle && passTrkEta && passDeltaEta
         && passTrkPt && passTrkPtErr && passTrkPurity && passTrkNhits
         ) return true;
   return false;
}

bool passD0KinematicCuts(Event* event, const int& icand)
{
   // eta is not available in the TTree
   //bool passEta = fabs(event->Eta(icand)<1000.);
//   bool passEta = fabs(event->Eta(icand)<1.5);
   bool passY = event->Y(icand) > ana::d0_y_min_ && event->Y(icand) < ana::d0_y_max_;
   return passY;
}

bool passD0MVA(const int& trigger, Event* event, const int& icand, 
      const bool& isPrompt)
{
   if(isPrompt){
      switch(trigger){
         case 2: return ana::pass_pPb2016_8TeV_PD0_MVA(event->Pt(icand), event->Mva(icand));
                 break;
         case 5: return ana::pass_pp2018_13TeV_PD0_MVA(event->Mva(icand));
                 break;
         default: return false;
      }
   } else{
      // non-prompt would not be measured
      return false;
   }
   return false;
}

inline bool passGoodVtx(Event* event)
{
   if(event->BestVtxZ() < -15. || event->BestVtxZ() > 15.) return false;
   return true;
}

bool passGoodTrack(Event* event, const unsigned int& icand)
{
   bool passHighPurity = true;
   bool passDzErr = true;
   bool passDxyErr = true;
   bool passPt = event->PtTrk(icand) > 0.4;
   bool passPtError = true;
   bool passEta = fabs(event->EtaTrk(icand)) < 2.4;
   return passHighPurity && passDzErr && passDxyErr &&
      passPt && passPtError && passEta;
}

inline bool passNtrkoffline(const double& Ntrkoffline, const int& dataset_trigger)
{
   if(dataset_trigger == ana::dataset_trigger.at("PAHM1-6")) 
      return Ntrkoffline >= ana::multMin_PA_ && Ntrkoffline < ana::multMax_PA_;
   if(dataset_trigger == ana::dataset_trigger.at("PPHM_2")) 
      return Ntrkoffline >= ana::multMin_PP_ && Ntrkoffline < ana::multMax_PP_;
   return false;
}

vector<double> setPtBin(string& dataset)
{
   if(dataset == "PAHM1-6"){
      return vector<double>(ana::ptbin_PD0_pPb, ana::ptbin_PD0_pPb+ana::nPt_PD0_pPb+1);
   }
   if(dataset == "PPHM_2"){
      return vector<double>(ana::ptbin_PD0_pp, ana::ptbin_PD0_pp+ana::nPt_PD0_pp+1);
   }
   return vector<double>();
}

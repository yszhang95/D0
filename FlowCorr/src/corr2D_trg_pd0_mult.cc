#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <stdlib.h>

#include "TFile.h"
#include "TChain.h"
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

#include "Event.h"
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
   const int nTrkBin = ana::Get_N_nTrkBin(dataset);
   if(dataset_trigger<0 || nTrkBin<0){
      cerr << "wrong dataset name" << endl;
      cout << "name should be:\n" 
         << "PAMB\n"
         << "PAHM0\n"
         << "PAHM1-6\n"
         << "PAHM7\n"
         << "PPMB\n"
         << "PPHM_1  //mult 80-100 \n"
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
   std::cout << "tracks ready" << std::endl;

   Event* evt = new Event(chain_d0, chain_tracks);
   setBranchStatus(evt);
   if(!checkBranchStatus(evt)) return -1;

   // declare hists
   TH1D* hMult;
   TH1D* hMult_ass;

   TH1D* hKET_D0[nTrkBin];
   TH1D* hPt_D0[nTrkBin];
   TH1D* hEta_D0[nTrkBin];
   TH1D* hRapidity_D0[nTrkBin];

   TH1D* hMass_D0[ana::nMass][nTrkBin];

   TH2D* hNtrkofflineVsNtrkgood;

   TH3D* hDcaVsMassAndMva[nTrkBin];

   TH1D* hMult_raw_D0[ana::nMass][nTrkBin];
   TH1D* hMult_eff_D0[ana::nMass][nTrkBin];
   TH2D* hSignal_D0[ana::nMass][nTrkBin];
   TH2D* hBackground_D0[ana::nMass][nTrkBin];

   hMult = new TH1D("hMult", "", 600, 0, 600);
   hMult_ass = new TH1D("hMult_ass", "", 600, 0, 600);

   hNtrkofflineVsNtrkgood = new TH2D("hNtrkofflineVsNtrkgood", "", 300, 0, 300, 300, 0, 300);

   for(int iTrkBin=0; iTrkBin<nTrkBin; iTrkBin++){
      hKET_D0[iTrkBin] = new TH1D(Form("hKET_trk%d", iTrkBin), "", 3000, 0, 30);
      hPt_D0[iTrkBin] = new TH1D(Form("hPt_trk%d", iTrkBin), "", 3000, 0, 30);
      hEta_D0[iTrkBin] = new TH1D(Form("hEta_trk%d", iTrkBin), "", 24, -2.4, 2.4);
      hRapidity_D0[iTrkBin] = new TH1D(Form("hRapidity_trk%d", iTrkBin), "", 24, -2.4, 2.4);
      hDcaVsMassAndMva[iTrkBin] = new TH3D(Form("hDcaVsMassAndMva_trk%d", iTrkBin), "", 60, 1.7, 2.0, 100, -0.3, 0.7, 160, 0, 0.08);
      for(int imass=0; imass<ana::nMass; imass++){
         hMass_D0[imass][iTrkBin] = new TH1D(Form("hMassD0_mass%d_trk%d", imass, iTrkBin),
               "", 200, 1.5, 2.5);
         hMult_raw_D0[imass][iTrkBin] = new TH1D(Form("hMult_raw_D0_mass%d_trk%d", imass, iTrkBin),
               "", 50, 0, 50);
         hMult_eff_D0[imass][iTrkBin] = new TH1D(Form("hMult_eff_D0_mass%d_trk%d", imass, iTrkBin),
               "", 50, 0, 50);
         hSignal_D0[imass][iTrkBin] = new TH2D(Form("hSignal_mass%d_trk%d", imass, iTrkBin),
               "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
               ana::nPhiBin, ana::phiBegin, ana::phiEnd);
         hBackground_D0[imass][iTrkBin] = new TH2D(Form("hBackground_mass%d_trk%d", imass, iTrkBin),
               "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
               ana::nPhiBin, ana::phiBegin, ana::phiEnd);
      }
   }

   // declare vectors
   vector<TVector3> pVect_trg_d0[ana::nMass];
   vector<double>   effVect_trg_d0[ana::nMass];
   vector<int>      indexVect_d0[ana::nMass];

   vector<TVector3> pVect_dau1_d0[ana::nMass];
   vector<TVector3> pVect_dau2_d0[ana::nMass];


   vector<TVector3>       pVect_ass;
   list<vector<TVector3>> pVectList_ass[ana::nZ_Vtx_][nTrkBin];
   vector<float>          effVect_ass;
   list<vector<float>>    effVectList_ass[ana::nZ_Vtx_][nTrkBin];

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
      auto iTrkBin = ana::findNtrkBin(evt->nTrkOffline(), dataset_trigger);
      if(iTrkBin<0) continue;

      // count number of good tracks per event
      unsigned int nMult_ass_good = 0;
      for(unsigned int itrack=0; itrack<evt->CandSizeTrk(); itrack++){
         // assume all tracks in TTree are good, since error of dz and dxy are not available
         if(passGoodTrack(evt, itrack)) 
            nMult_ass_good++;
      }

      hMult->Fill(nMult_ass_good);

      hNtrkofflineVsNtrkgood->Fill(nMult_ass_good, evt->nTrkOffline());

      for(int id0=0; id0<evt->CandSize(); id0++){

         int imass = ana::findMassBin(evt->Mass(id0));
         if(imass == -1) continue;

         if(!passD0Selections(dataset_trigger, evt, id0, isPromptD0)) continue;

         double effks = h_eff->FindBin(evt->Pt(id0), evt->Y(id0));
         hDcaVsMassAndMva[iTrkBin]->Fill(evt->Mass(id0), evt->Mva(id0), 0., 1./effks);

         p_dau1.SetPtEtaPhi(evt->PtD1(id0), evt->etaD1(id0), evt->phiD1(id0));
         p_dau2.SetPtEtaPhi(evt->PtD2(id0), evt->etaD2(id0), evt->phiD2(id0));
         p_d0 = p_dau1 + p_dau2;

         hMass_D0[imass][iTrkBin]->Fill(evt->Mass(id0), 1./effks);
         hPt_D0[iTrkBin]->Fill(evt->Pt(id0), 1./effks);
         hEta_D0[iTrkBin]->Fill(p_d0.Eta(), 1./effks);
         hRapidity_D0[iTrkBin]->Fill(evt->Y(id0), 1./effks);
         double KET = sqrt(pow(evt->Mass(id0), 2) + pow(evt->Pt(id0), 2)
               - evt->Mass(id0));
         hKET_D0[iTrkBin]->Fill(KET, 1./effks);

         pVect_trg_d0[imass].push_back(p_d0);
         effVect_trg_d0[imass].push_back(effks);
         indexVect_d0[imass].push_back(id0);
         pVect_dau1_d0[imass].push_back(p_dau1);
         pVect_dau2_d0[imass].push_back(p_dau2);
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
      
      unsigned int nMult_trg_raw_d0[ana::nMass]; // Ntrig for mass & pt bins
      double nMult_trg_eff_d0[ana::nMass]; // eff corrected Ntrig for mass & pt bins

      for(int imass=0; imass<ana::nMass; imass++){
         nMult_trg_raw_d0[imass] = 0;
         nMult_trg_eff_d0[imass] = 0;
      }

      for(int imass=0; imass<ana::nMass; imass++){
         unsigned int nMult_trg_d0 = (unsigned int) pVect_trg_d0[imass].size();
         for(unsigned int id0=0; id0<nMult_trg_d0; id0++){
            double effks = effVect_trg_d0[imass].at(id0);
            //int index = indexVect_d0[imass].at(id0);
            nMult_trg_raw_d0[imass] += 1;
            nMult_trg_eff_d0[imass] += 1./effks;
         }
         hMult_raw_D0[imass][iTrkBin]->Fill(nMult_trg_raw_d0[imass]);
         hMult_eff_D0[imass][iTrkBin]->Fill(nMult_trg_eff_d0[imass]);

         for(unsigned int id0=0; id0<nMult_trg_d0; id0++){
            double effks = effVect_trg_d0[imass].at(id0);
            for(unsigned int iass=0; iass<nMult_ass; iass++){
               TVector3 pvector_ass = pVect_ass.at(iass);
               double effweight_ass = effVect_ass.at(iass);
               if(ana::rejectDaughter_){
                  if(fabs(pvector_ass.Eta() - pVect_dau1_d0[imass].at(id0).Eta())<0.03
                        && fabs(pvector_ass.DeltaPhi(pVect_dau1_d0[imass].at(id0)))<0.03)
                     continue;
                  if(fabs(pvector_ass.Eta() - pVect_dau2_d0[imass].at(id0).Eta())<0.03
                        && fabs(pvector_ass.DeltaPhi(pVect_dau2_d0[imass].at(id0)))<0.03)
                     continue;
               }

               double deltaEta = pvector_ass.Eta() - pVect_trg_d0[imass].at(id0).Eta();
               double deltaPhi = pvector_ass.DeltaPhi(pVect_trg_d0[imass].at(id0));
               if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

               //int index = indexVect_d0[imass].at(id0);
               hSignal_D0[imass][iTrkBin]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0[imass]/effks/effweight_ass);
            }
         }
      }

      // mixed events
      // begin filling background
      // if having enough, erase the first one, fill the current, fill the background
      if(pVectList_ass[iz][iTrkBin].size() == ana::nMixedEvts){
         for(int imass=0; imass<ana::nMass; imass++){
            unsigned int nMult_trg_d0 = (unsigned int) pVect_trg_d0[imass].size();
            for(unsigned int id0=0; id0<nMult_trg_d0; id0++){
               double effks = effVect_trg_d0[imass].at(id0);
               // simultaneously read momentum and efficiency
               auto ievt_eff_ass = effVectList_ass[iz][iTrkBin].begin();
               for(auto& ievt_p_ass : pVectList_ass[iz][iTrkBin]){
                  unsigned int n_ass = ievt_p_ass.size();
                  if(n_ass!=ievt_eff_ass->size()){
                     cout << "wrong" << endl;
                     break;
                  }
                  for(unsigned int iass=0; iass<n_ass; iass++){
                     TVector3 pvector_ass = ievt_p_ass.at(iass);
                     double effweight_ass = ievt_eff_ass->at(iass);

                     double deltaEta = pvector_ass.Eta() - pVect_trg_d0[imass].at(id0).Eta();
                     double deltaPhi = pvector_ass.DeltaPhi(pVect_trg_d0[imass].at(id0));
                     if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

                     //int index = indexVect_d0[imass].at(id0);
                     hBackground_D0[imass][iTrkBin]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0[imass]/effks/effweight_ass);
                  }
                  ievt_eff_ass++;
               }
            }
         }
      }

      // clear all and refresh list/buffer
      if(pVectList_ass[iz][iTrkBin].size() == ana::nMixedEvts){
         // update the list/buffer
         pVectList_ass[iz][iTrkBin].erase(pVectList_ass[iz][iTrkBin].begin());
         pVectList_ass[iz][iTrkBin].push_back(pVect_ass);
         effVectList_ass[iz][iTrkBin].erase(effVectList_ass[iz][iTrkBin].begin());
         effVectList_ass[iz][iTrkBin].push_back(effVect_ass);
      }else{
      // if not having enough events, keep filling the list of momentum vector of associated particles
         // fill the list/buffer
         pVectList_ass[iz][iTrkBin].push_back(pVect_ass);
         effVectList_ass[iz][iTrkBin].push_back(effVect_ass);
      }

      for(int imass=0; imass<ana::nMass; imass++){
         pVect_trg_d0[imass].clear();
         pVect_dau1_d0[imass].clear();
         pVect_dau2_d0[imass].clear();
         indexVect_d0[imass].clear();
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
      size_t found = datalist.find("/");
      if (found!=std::string::npos) datalist.replace(found, 1, "_");
      if(prefix.size())
         outName = TString::Format("%s/fout_%s_d0ana_ntrk_%.1f.root", prefix.c_str(), datalist.c_str(), ana::d0_y_max_);
      else
         outName = TString::Format("fout_%s_d0ana_ntrk_%.1f.root", datalist.c_str(), ana::d0_y_max_);
   }
   else return -1;

   TFile fout(outName.Data(), "recreate");
   fout.cd();

   // start writing output
   hMult->Write();
   hMult_ass->Write();
   hNtrkofflineVsNtrkgood->Write();
   for(int iTrkBin=0; iTrkBin<nTrkBin; iTrkBin++){
      hKET_D0[iTrkBin]->Write(); 
      hPt_D0[iTrkBin]->Write();
      hEta_D0[iTrkBin]->Write();
      hRapidity_D0[iTrkBin]->Write();
      hDcaVsMassAndMva[iTrkBin]->Write();
      for(int imass=0; imass<ana::nMass; imass++){
         hMass_D0[imass][iTrkBin]->Write();
         hMult_raw_D0[imass][iTrkBin]->Write();
         hMult_eff_D0[imass][iTrkBin]->Write();
         hSignal_D0[imass][iTrkBin]->Write();
         hBackground_D0[imass][iTrkBin]->Write();
      }
   }

   delete hMult;
   delete hMult_ass;
   delete hNtrkofflineVsNtrkgood;
   for(int iTrkBin=0; iTrkBin<nTrkBin; iTrkBin++){
      delete hKET_D0[iTrkBin]; 
      delete hPt_D0[iTrkBin];
      delete hEta_D0[iTrkBin];
      delete hRapidity_D0[iTrkBin];
      delete hDcaVsMassAndMva[iTrkBin];
      for(int imass=0; imass<ana::nMass; imass++){
         delete hMass_D0[imass][iTrkBin];
         delete hMult_raw_D0[imass][iTrkBin];
         delete hMult_eff_D0[imass][iTrkBin];
         delete hSignal_D0[imass][iTrkBin];
         delete hBackground_D0[imass][iTrkBin];
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
      //event->GetBranchStatus("zDCASignificanceDaugther1")&&
      ///event->GetBranchStatus("xyDCASignificanceDaugther1")&&
      event->GetBranchStatus("NHitD1")&&
      event->GetBranchStatus("HighPurityDaugther1")&&

      event->GetBranchStatus("pTD2")&&
      event->GetBranchStatus("pTerrD2")&&
      event->GetBranchStatus("EtaD2")&&
      event->GetBranchStatus("PhiD2")&&
      //event->GetBranchStatus("zDCASignificanceDaugther2")&&
      //event->GetBranchStatus("xyDCASignificanceDaugther2")&&
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
   //bool passEta = fabs(event->Eta(icand)<1000.);
   //bool passEta = fabs(event->Eta(icand)<1.5);
   bool passPt = event->Pt(icand) < ana::d0_pt_max_ && event->Pt(icand) >= ana::d0_pt_min_;
   bool passY = event->Y(icand) < ana::d0_y_max_ && event->Y(icand) >= ana::d0_y_min_;
   return //passEta;
      passY &&
      passPt;
}

bool passD0MVA(const int& trigger, Event* event, const int& icand, 
      const bool& isPrompt)
{
   if(isPrompt){
      switch(trigger){
         case 0: return ana::pass_pPb2016_8TeV_PD0_MVA(event->Pt(icand), event->Mva(icand));
                 break;
         case 1: return ana::pass_pPb2016_8TeV_PD0_MVA(event->Pt(icand), event->Mva(icand));
                 break;
         case 2: return ana::pass_pPb2016_8TeV_PD0_MVA(event->Pt(icand), event->Mva(icand));
                 break;
         case 3: return ana::pass_pPb2016_8TeV_PD0_MVA(event->Pt(icand), event->Mva(icand));
                 break;
         case 4: return ana::pass_pp2018_13TeV_PD0_MVA(event->Mva(icand));
                 break;
         case 5: return ana::pass_pp2018_13TeV_PD0_MVA(event->Mva(icand));
                 break;
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

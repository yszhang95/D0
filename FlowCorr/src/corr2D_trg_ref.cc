// ntracks done
// d0 selections done
// d0 daughter rejection done
// correlation same event done
// correlation mixed events

#include <iostream>
#include <vector>
#include <map>
#include <list>

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
#include "TVector3.h"
#include "TString.h"

#include "myAnaConsts.h"

using namespace std;

void setBranchStatus(Event*);
bool checkBranchStatus(Event*);

bool passGoodTrack(Event*, const unsigned int&);
inline bool passGoodVtx(Event* event);

int main(int argc, char** argv)
{

   TH1::SetDefaultSumw2(true);

   if(argc!=3) {
      std::cerr << "The number of arguments is wrong" << std::endl;
      return -1;
   }
   
   string datalist(argv[1]);
   std::cout << datalist << std::endl;


   //bool isPromptD0 = false;
   long int tree = strtol(argv[2], NULL, 10);
   //if(tree == 0)  isPromptD0 = true;

   TChain *chain_d0 = new TChain(TString::Format("%s/VertexCompositeNtuple", ana::treeName[tree].c_str()));
   TChain *chain_tracks = new TChain("track_ana/trackTree");

   TFileCollection* fcData = new TFileCollection(datalist.c_str(), "", datalist.c_str());

   chain_d0->AddFileInfoList(fcData->GetList());
   std::cout << ana::treeName[tree] << " event tree ready" << std::endl;

   chain_tracks->AddFileInfoList(fcData->GetList());
   std::cout << "tracks ready" << std::endl;

   Event* evt = new Event(chain_d0, chain_tracks);
   setBranchStatus(evt);
   if(!checkBranchStatus(evt)) return -1;

   // declare hists
   TH1D* hMult;

   TH2D* hSignal_Ref;
   TH2D* hBackground_Ref;

   hMult = new TH1D("hMult", "", 600, 0, 600);

   hSignal_Ref = new TH2D("hSignal_Ref",
            "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
            ana::nPhiBin, ana::phiBegin, ana::phiEnd);
   hBackground_Ref = new TH2D("hBackground_Ref",
            "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
            ana::nPhiBin, ana::phiBegin, ana::phiEnd);

   // declare vectors
   vector<TVector3> pVect_ref;
   list<vector<TVector3>> pVectList_ref[ana::nZ_Vtx_];
   vector<float> effVect_ref;
   list<vector<float>> effVectList_ref[ana::nZ_Vtx_];

   // start timing
   TStopwatch ts;
   ts.Start();

   // loop
   //
   // temporary vectors
   TVector3 p_ref(0, 0, 0), p_temp(0, 0, 0);

   std::cout << evt->GetEntries() << std::endl;
   long int nentries = evt->GetEntries();
   int percent = 0;
   long int skip = 0;
   for(long int ientry=0; ientry<nentries; ientry++){
      if(pVect_ref.size()!=0){
         cout << "did not clear the vector pVect_ref.size" << endl;
         continue;
      }
      if(effVect_ref.size()!=0){
         cout << "did not clear the vector effVect_ref.size" << endl;
         continue;
      }

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

      // good vertex check
      if(!passGoodVtx(evt)) continue;
      auto iz = ana::findZVtxBin(evt->BestVtxZ());
      if(iz == -1) continue;

      // count number of good tracks per event
      unsigned int nMult_trk_good = 0;
      /*
      for(unsigned int itrack=0; itrack<evt->CandSizeTrk(); itrack++){
         // assume all tracks in TTree are good, since error of dz and dxy are not available
         if(passGoodTrack(evt, itrack)) 
            nMult_trk_good++;
      }
      */
      nMult_trk_good = evt->nTrkOffline();
      hMult->Fill(nMult_trk_good);

      if(nMult_trk_good<ana::multMax_ && nMult_trk_good>=ana::multMin_){
      }else{
         //std::cout << "multiplicity wrong" << std::endl;
         continue;
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
            p_ref.SetPtEtaPhi(evt->PtTrk(itrack), evt->EtaTrk(itrack), evt->PhiTrk(itrack));
            pVect_ref.push_back(p_ref);
            effVect_ref.push_back(evt->WeightTrk(itrack));
         }
      }

      // calculate signal
      unsigned int nMult_ref = (unsigned int) pVect_ref.size();
      double nMult_eff_ref = 0;
      for(unsigned int iref=0; iref<nMult_ref; iref++){
         nMult_eff_ref += 1./effVect_ref.at(iref);
      }

      for(unsigned int itrg=0; itrg<nMult_ref; itrg++){
         for(unsigned int iass=itrg+1; iass<nMult_ref; iass++){

            TVector3 pvector_trg = pVect_ref.at(itrg);
            double effweight_trg = effVect_ref.at(itrg);

            TVector3 pvector_ass = pVect_ref.at(iass);
            double effweight_ass = effVect_ref.at(iass);

            // reflect delta around 0, since one pair is filled twice
            double deltaEta = pvector_ass.Eta() - pvector_trg.Eta();
            double negDeltaEta = -1* deltaEta;
            double deltaPhi = pvector_ass.DeltaPhi(pvector_trg);
            double negDeltaPhi = -1* deltaPhi;
            if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;
            if(negDeltaPhi>-ana::PI && negDeltaPhi<-ana::PI/2.) negDeltaPhi += 2*ana::PI;

            hSignal_Ref->Fill(deltaEta, deltaPhi, 
                  1./nMult_eff_ref/effweight_trg/effweight_ass);
            hSignal_Ref->Fill(negDeltaEta, negDeltaPhi, 
                  1./nMult_eff_ref/effweight_trg/effweight_ass);
         }
      }
      //

      // mixed events
      // begin filling background
      // if having enough, erase the first one, fill the current, fill the background
      if(pVectList_ref[iz].size() == ana::nMixedEvts){
         for(unsigned int itrg=0; itrg<nMult_ref; itrg++){
            // simultaneously read momentum and efficiency
            TVector3 pvector_trg = pVect_ref.at(itrg);
            double effweight_trg = effVect_ref.at(itrg);

            auto ievt_eff_ass = effVectList_ref[iz].begin();
            for(auto ievt_p_ass= pVectList_ref[iz].begin(); ievt_p_ass != pVectList_ref[iz].end(); ievt_p_ass++){

               unsigned int n_ass = ievt_p_ass->size();
               if(n_ass!=ievt_eff_ass->size()){
                  cout << "wrong" << endl;
                  break;
               }
               for(unsigned int iass=0; iass<n_ass; iass++){
                  TVector3 pvector_ass = ievt_p_ass->at(iass);
                  double effweight_ass = ievt_eff_ass->at(iass);
                  double deltaEta = pvector_ass.Eta() - pvector_trg.Eta();
                  double deltaPhi = pvector_ass.DeltaPhi(pvector_trg);
                  if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

                  hBackground_Ref->Fill(deltaEta, deltaPhi, 
                     1./nMult_eff_ref/effweight_trg/effweight_ass);
               }
               ievt_eff_ass++;
            }
         }
      }

      // clear all
      if(pVectList_ref[iz].size() == ana::nMixedEvts){
         // update the list/buffer
         pVectList_ref[iz].erase(pVectList_ref[iz].begin());
         pVectList_ref[iz].push_back(pVect_ref);
         effVectList_ref[iz].erase(effVectList_ref[iz].begin());
         effVectList_ref[iz].push_back(effVect_ref);
      } else{
         // if not having enough events, keep filling the list of momentum vector of associated particles
         // fill the list/buffer
         pVectList_ref[iz].push_back(pVect_ref);
         effVectList_ref[iz].push_back(effVect_ref);
      }
      pVect_ref.clear();
      effVect_ref.clear();
   }
   std::cout << "completed loop" << std::endl;
   std::cout << skip << " events are skipped" << std::endl;

   ts.Stop();
   ts.Print();

   TFile fout(TString::Format("fout_ref_%s_%s.root", datalist.c_str(), ana::treeName[tree].c_str()), "recreate");
   fout.cd();

   // start writing output
   
   hMult->Write();
   hSignal_Ref->Write();
   hBackground_Ref->Write();

   delete hMult;
   delete hSignal_Ref;
   delete hBackground_Ref;

   delete evt;

   return 0;
}

void setBranchStatus(Event* evt)
{
   evt->SetBranchStatus("Ntrkoffline", 1);
   evt->SetBranchStatus("tracks.candSizeTRK", 1);
   evt->SetBranchStatus("tracks.pTTRK", 1);
   evt->SetBranchStatus("tracks.etaTRK", 1);
   evt->SetBranchStatus("tracks.phiTRK", 1);
   evt->SetBranchStatus("tracks.weightTRK", 1);
}

bool checkBranchStatus(Event* event)
{
   bool check =
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

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
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

class corr2D_pd0{
   public:
      corr2D_pd0(std::string name);
      ~corr2D_pd0();

      int write();

      TH1D* histKET() {return hKET_D0;}
      TH1D* histPt() {return hPt_D0;}
      TH1D* histEta() {return hEta_D0;}
      TH1D* histRapidity() {return hRapidity_D0;}
      TH1D* histMass() {return hMass;}

      TH1D** histMult_ass() {return hMult_ass;}
      TH1D** histMass_D0() {return hMass_D0;}
      TH1D** histMult_raw_D0() {return hMult_raw_D0;}
      TH1D** histMult_eff_D0() {return hMult_eff_D0;}
      TH2D** histSignal_D0() {return hSignal_D0;}
      TH2D** histBackground_D0() {return hBackground_D0;};
      TH2D** histSignal_D0_symm() {return hSignal_D0_symm;}
      TH2D** histBackground_D0_symm() {return hBackground_D0_symm;};
      TH2D** histSignal_D0_old() {return hSignal_D0_old;}
      TH2D** histBackground_D0_old() {return hBackground_D0_old;};

      TH3D* histDcaVsMassAndMva() {return hDcaVsMassAndMva;}

   private:
      corr2D_pd0(){};
      corr2D_pd0(const corr2D_pd0& tmp){}
      corr2D_pd0 operator=(const corr2D_pd0& tmp){return *this;}

      TH1D* hMult_ass[ana::nMass];
   
      TH1D* hKET_D0;
      TH1D* hPt_D0;
      TH1D* hEta_D0;
      TH1D* hRapidity_D0;
   
      TH1D* hMass;
      TH1D* hMass_D0[ana::nMass];
      TH1D* hMult_raw_D0[ana::nMass];
      TH1D* hMult_eff_D0[ana::nMass];
   
      TH2D* hSignal_D0[ana::nMass];
      TH2D* hBackground_D0[ana::nMass];
      TH2D* hSignal_D0_symm[ana::nMass];
      TH2D* hBackground_D0_symm[ana::nMass];
      TH2D* hSignal_D0_old[ana::nMass];
      TH2D* hBackground_D0_old[ana::nMass];

      TH3D* hDcaVsMassAndMva;
};

corr2D_pd0::corr2D_pd0(string name)
{
   
   hKET_D0 = new TH1D(Form("hKET_%s", name.c_str()), Form("hKET_%s", name.c_str()), 2000, 0, 20);
   hPt_D0 = new TH1D(Form("hPt_%s", name.c_str()), Form("hPt_%s", name.c_str()), 2000, 0, 20);
   hEta_D0 = new TH1D(Form("hEta_%s", name.c_str()), Form("hEta_%s", name.c_str()), 24, -2.4, 2.4);
   hRapidity_D0 = new TH1D(Form("hRapidity_%s", name.c_str()), Form("hRapidity_%s", name.c_str()), 24, -2.4, 2.4);
   
   hMass = new TH1D(Form("hMass _%s", name.c_str()), Form("hMass _%s", name.c_str()), 60, 1.7, 2.0);
   for(int imass=0; imass<ana::nMass; imass++){
      hMult_ass[imass] = new TH1D(Form("hMult_ass_mass%d_%s", imass, name.c_str()), Form("hMult_ass_mass%d_%s", imass, name.c_str()), 600, 0, 600);
      hMass_D0[imass] = new TH1D(Form("hMass_D0_mass%d_%s", imass, name.c_str()), Form("hMass_D0_mass%d_%s", imass, name.c_str()), 200, 1.5, 2.5);
      hMult_raw_D0[imass] = new TH1D(Form("hMult_raw_D0_mass%d_%s", imass, name.c_str()), Form("hMult_raw_D0_mass%d_%s", imass, name.c_str()), 50, 0, 50);
      hMult_eff_D0[imass] = new TH1D(Form("hMult_eff_D0_mass%d_%s", imass, name.c_str()), Form("hMult_eff_D0_mass%d_%s", imass, name.c_str()), 50, 0, 50);
      hSignal_D0[imass] = new TH2D(Form("hSignal_D0_mass%d_%s", imass, name.c_str()), Form("hSignal_D0_mass%d_%s", imass, name.c_str()), 
            ana::nEtaBin, ana::etaBegin, ana::etaEnd, ana::nPhiBin, ana::phiBegin, ana::phiEnd);
      hBackground_D0[imass] = new TH2D(Form("hBackground_D0_mass%d_%s", imass, name.c_str()), Form("hBackground_D0_mass%d_%s", imass, name.c_str()),
            ana::nEtaBin, ana::etaBegin, ana::etaEnd, ana::nPhiBin, ana::phiBegin, ana::phiEnd);
      hSignal_D0_symm[imass] = new TH2D(Form("hSignal_D0_symm_mass%d_%s", imass, name.c_str()), Form("hSignal_D0_mass%d_%s", imass, name.c_str()), 
            ana::nEtaBin, ana::etaBegin, ana::etaEnd, ana::nPhiBin, ana::phiBegin, ana::phiEnd);
      hBackground_D0_symm[imass] = new TH2D(Form("hBackground_D0_symm_mass%d_%s", imass, name.c_str()), Form("hBackground_D0_mass%d_%s", imass, name.c_str()),
            ana::nEtaBin, ana::etaBegin, ana::etaEnd, ana::nPhiBin, ana::phiBegin, ana::phiEnd);
      hSignal_D0_old[imass] = new TH2D(Form("hSignal_D0_old_mass%d_%s", imass, name.c_str()), Form("hSignal_D0_mass%d_%s", imass, name.c_str()), 
            ana::nEtaBin, ana::etaBegin, ana::etaEnd, ana::nPhiBin, ana::phiBegin_old, ana::phiEnd_old);
      hBackground_D0_old[imass] = new TH2D(Form("hBackground_D0_old_mass%d_%s", imass, name.c_str()), Form("hBackground_D0_mass%d_%s", imass, name.c_str()),
            ana::nEtaBin, ana::etaBegin, ana::etaEnd, ana::nPhiBin, ana::phiBegin_old, ana::phiEnd_old);
   }

   hDcaVsMassAndMva = new TH3D(Form("hDcaVsMassAndMva_%s", name.c_str()), Form("hDcaVsMassAndMva_%s", name.c_str()), 
         60, 1.7, 2.0, 100, -0.3, 0.7, 160, 0, 0.08);
}

corr2D_pd0::~corr2D_pd0()
{
   if(hKET_D0 ) delete hKET_D0 ; 
   if(hPt_D0 ) delete hPt_D0 ; 
   if(hEta_D0 ) delete hEta_D0 ; 
   if(hRapidity_D0 ) delete hRapidity_D0 ; 
   
   if(hMass ) delete hMass ; 
   for(int imass=0; imass<ana::nMass; imass++){
      if(hMult_ass[imass] ) delete hMult_ass[imass]; 
      if(hMass_D0[imass] ) delete hMass_D0[imass] ; 
      if(hMult_raw_D0[imass] ) delete hMult_raw_D0[imass] ; 
      if(hMult_eff_D0[imass] ) delete hMult_eff_D0[imass] ; 
      if(hSignal_D0[imass] ) delete hSignal_D0[imass] ; 
      if(hBackground_D0[imass] ) delete hBackground_D0[imass] ; 
      if(hSignal_D0_symm[imass] ) delete hSignal_D0_symm[imass] ; 
      if(hBackground_D0_symm[imass] ) delete hBackground_D0_symm[imass] ; 
      if(hSignal_D0_old[imass] ) delete hSignal_D0_old[imass] ; 
      if(hBackground_D0_old[imass] ) delete hBackground_D0_old[imass] ; 
   }

   if(hDcaVsMassAndMva ) delete hDcaVsMassAndMva ; 
}

int corr2D_pd0::write()
{
   if(hKET_D0 ) hKET_D0->Write() ; 
   if(hPt_D0 ) hPt_D0->Write() ; 
   if(hEta_D0 ) hEta_D0->Write() ; 
   if(hRapidity_D0 ) hRapidity_D0->Write() ; 
   
   if(hMass ) hMass->Write() ; 
   for(int imass=0; imass<ana::nMass; imass++){
      if(hMult_ass[imass] ) hMult_ass[imass]->Write(); 
      if(hMass_D0[imass] ) hMass_D0[imass]->Write(); 
      if(hMult_raw_D0[imass] ) hMult_raw_D0[imass]->Write(); 
      if(hMult_eff_D0[imass] ) hMult_eff_D0[imass]->Write(); 
      if(hSignal_D0[imass] ) hSignal_D0[imass]->Write(); 
      if(hBackground_D0[imass] ) hBackground_D0[imass]->Write(); 
      if(hSignal_D0_old[imass] ) hSignal_D0_old[imass]->Write(); 
      if(hBackground_D0_old[imass] ) hBackground_D0_old[imass]->Write(); 
      if(hSignal_D0_symm[imass] ) hSignal_D0_symm[imass]->Write(); 
      if(hBackground_D0_symm[imass] ) hBackground_D0_symm[imass]->Write(); 
   }

   if(hDcaVsMassAndMva ) hDcaVsMassAndMva->Write() ; 

   return 0;
}

static bool isPromptD0 = true;

void setBranchStatus(Event*);
bool checkBranchStatus(Event*);

bool passGoodTrack(Event*, const unsigned int&);
inline bool passGoodVtx(Event* event);
inline bool passD0Selections(const int&, Event*, const int&, const bool&);
bool passD0PreSelections(Event*, const int&);
bool passD0KinematicCuts(Event*, const int&);
bool passD0MVA(const int&, Event*, const int&, const bool&, const string&);

// par0, main
// par1, datalist
// par2, dataset
// par3, efficiency hists
// par4, output dir
// par5, pT_Min
// par6, pT_Max
// par7, y_Min
// par8, y_Max

struct kinematicalCuts{
   float pTMin;
   float pTMax;
   float yMin;
   float yMax;
} cuts;

int main(int argc, char** argv)
{
   TH1::SetDefaultSumw2(true);

   vector<double> ptbin = {2., 4., 6., 8.};

   if(argc!=9) {
      std::cerr << "The number of arguments is wrong" << std::endl;
      return -1;
   }
   
   string datalist(argv[1]);
   std::cout << datalist << std::endl;

   string dataset(argv[2]);

   // the kinematical cuts
   cuts.pTMin = stof(argv[5]);
   cuts.pTMax = stof(argv[6]);
   cuts.yMin  = stof(argv[7]);
   cuts.yMax  = stof(argv[8]);

   const int dataset_trigger = ana::Get_Trigger(dataset);
   cout << "dataset number: " << dataset_trigger << endl;
   const int nTrkBin = ana::Get_N_nTrkBin(dataset);
   cout << "number of bins of ntrk: " <<  nTrkBin << endl;
   if(dataset_trigger<0 || nTrkBin<0){
      cerr << "wrong dataset name" << endl;
      cout << "name should be:\n" 
         << "PAMB\n"
         << "PAHM1-6\n"
         << "PAHM7\n"
         << "PPMB\n"
         << "PPHM  //mult 80-100, 100-inf \n"
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
   TChain *chain_evt = new TChain("eventinfoana/EventInfoNtuple");

   TFileCollection* fcData = new TFileCollection(datalist.c_str(), "", datalist.c_str());

   chain_d0->AddFileInfoList(fcData->GetList());

   chain_tracks->AddFileInfoList(fcData->GetList());

   chain_evt->AddFileInfoList(fcData->GetList());

   //std::cout << "tracks ready" << std::endl;

   Event* evt = new Event(chain_d0, chain_tracks, chain_evt);
   setBranchStatus(evt);
   if(!checkBranchStatus(evt)) return -1;

   // declare hists
   const unsigned int nPt = ptbin.size()-1;
   corr2D_pd0* corr2Dptr[nTrkBin][nPt];
   corr2D_pd0* corr2Dptr_loose[nTrkBin][nPt];
   corr2D_pd0* corr2Dptr_tight[nTrkBin][nPt];
   for(int iTrkBin=0; iTrkBin<nTrkBin; iTrkBin++){
      for(int ipt=0; ipt<nPt; ipt++){
         corr2Dptr[iTrkBin][ipt] = new corr2D_pd0(Form("pt%d_trk%d", ipt, iTrkBin));
         corr2Dptr_loose[iTrkBin][ipt] = new corr2D_pd0(Form("pt%d_trk%d_loose", ipt, iTrkBin));
         corr2Dptr_tight[iTrkBin][ipt] = new corr2D_pd0(Form("pt%d_trk%d_tight", ipt, iTrkBin));
      }
   }

   //
   TH1D* hEvt;
   TH1D* hMult;

   TH1D* hNtrk_D0[nTrkBin];

   TH2D* hNtrkofflineVsNtrkgood;

   hEvt  = new TH1D("hEvt", "", 600, 0, 600);
   hMult = new TH1D("hMult", "", 600, 0, 600);

   hNtrkofflineVsNtrkgood = new TH2D("hNtrkofflineVsNtrkgood", "", 300, 0, 300, 300, 0, 300);

   for(int iTrkBin=0; iTrkBin<nTrkBin; iTrkBin++){
      hNtrk_D0[iTrkBin] = new TH1D(Form("hNtrk_trk%d", iTrkBin), "", 3000, 0, 3000);
   }

   // declare vectors
   vector<TVector3> pVect_trg_d0;
   vector<double>   effVect_trg_d0;
   vector<int>      indexVect_d0;
   vector<unordered_map<string, bool>> isd0Vect_d0;
   vector<unordered_map<string, int>> ibinsVect_d0; // [0] ptbin, [1] mass bin

   vector<TVector3> pVect_dau1_d0;
   vector<TVector3> pVect_dau2_d0;

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

   unordered_map<string, int> ibins({{"ipt", -1}, {"imass", -1}});

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

      hEvt->Fill(evt->nTrkOffline());

      if(!passGoodVtx(evt)) continue;

      auto iz = ana::findZVtxBin(evt->BestVtxZ());
      if(iz < 0) continue;
      auto iTrkBin = ana::findNtrkBin(evt->nTrkOffline(), dataset_trigger);
      if(iTrkBin<0) continue;
      hNtrk_D0[iTrkBin]->Fill(evt->nTrkOffline());

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
         if(imass < 0) continue;
         if(ipt < 0) continue;

         if(!passD0PreSelections(evt, id0)) continue;
         if(!passD0KinematicCuts(evt, id0)) continue;

         p_dau1.SetPtEtaPhi(evt->PtD1(id0), evt->etaD1(id0), evt->phiD1(id0));
         p_dau2.SetPtEtaPhi(evt->PtD2(id0), evt->etaD2(id0), evt->phiD2(id0));
         p_d0.SetPtEtaPhi(evt->Pt(id0), evt->Eta(id0), evt->Phi(id0));

         double effks = h_eff->GetBinContent(h_eff->FindBin(evt->Pt(id0), evt->Y(id0)));

         pVect_trg_d0.push_back(p_d0);
         effVect_trg_d0.push_back(effks);
         indexVect_d0.push_back(id0);
         
         ibins.at("ipt") = ipt;
         ibins.at("imass") = imass;

         ibinsVect_d0.push_back(ibins);
         
         pVect_dau1_d0.push_back(p_dau1);
         pVect_dau2_d0.push_back(p_dau2);

         unordered_map<string, bool> passMVA = {
            {
            {"passStd", false},
            {"passLoose", false},
            {"passTight", false}
            }
         };

         if(passD0MVA(dataset_trigger, evt, id0, isPromptD0, "std")){
            passMVA.at("passStd") = true;
         }
         if(passD0MVA(dataset_trigger, evt, id0, isPromptD0, "loose")){
            passMVA.at("passLoose") = true;
         }
         if(passD0MVA(dataset_trigger, evt, id0, isPromptD0, "tight")){
            passMVA.at("passTight") = true;
         }

         isd0Vect_d0.push_back(passMVA);

         double KET = sqrt(pow(evt->Mass(id0), 2) + pow(evt->Pt(id0), 2)
                  - evt->Mass(id0));
         const double DCA = evt->DecayL3D(id0) * sin(evt->PointingAngle3D(id0));
         
         if(passMVA.at("passStd")){
            corr2Dptr[iTrkBin][ipt]->histPt()->Fill(evt->Pt(id0), 1./effks);
            corr2Dptr[iTrkBin][ipt]->histEta()->Fill(evt->Eta(id0), 1./effks);
            corr2Dptr[iTrkBin][ipt]->histRapidity()->Fill(evt->Y(id0), 1./effks);
            corr2Dptr[iTrkBin][ipt]->histKET()->Fill(KET, 1./effks);
            corr2Dptr[iTrkBin][ipt]->histMass_D0()[imass]->Fill(evt->Mass(id0), 1./effks);
            corr2Dptr[iTrkBin][ipt]->histMass()->Fill(evt->Mass(id0), 1./effks);
            corr2Dptr[iTrkBin][ipt]->histDcaVsMassAndMva()->Fill(evt->Mass(id0), evt->Mva(id0), DCA, 1./effks);
         }

         if(passMVA.at("passLoose")){
            corr2Dptr_loose[iTrkBin][ipt]->histPt()->Fill(evt->Pt(id0), 1./effks);
            corr2Dptr_loose[iTrkBin][ipt]->histEta()->Fill(evt->Eta(id0), 1./effks);
            corr2Dptr_loose[iTrkBin][ipt]->histRapidity()->Fill(evt->Y(id0), 1./effks);
            corr2Dptr_loose[iTrkBin][ipt]->histKET()->Fill(KET, 1./effks);
            corr2Dptr_loose[iTrkBin][ipt]->histMass_D0()[imass]->Fill(evt->Mass(id0), 1./effks);
            corr2Dptr_loose[iTrkBin][ipt]->histMass()->Fill(evt->Mass(id0), 1./effks);
            corr2Dptr_loose[iTrkBin][ipt]->histDcaVsMassAndMva()->Fill(evt->Mass(id0), evt->Mva(id0), DCA, 1./effks);
         }

         if(passMVA.at("passTight")){
            corr2Dptr_tight[iTrkBin][ipt]->histPt()->Fill(evt->Pt(id0), 1./effks);
            corr2Dptr_tight[iTrkBin][ipt]->histEta()->Fill(evt->Eta(id0), 1./effks);
            corr2Dptr_tight[iTrkBin][ipt]->histRapidity()->Fill(evt->Y(id0), 1./effks);
            corr2Dptr_tight[iTrkBin][ipt]->histKET()->Fill(KET, 1./effks);
            corr2Dptr_tight[iTrkBin][ipt]->histMass_D0()[imass]->Fill(evt->Mass(id0), 1./effks);
            corr2Dptr_tight[iTrkBin][ipt]->histMass()->Fill(evt->Mass(id0), 1./effks);
            corr2Dptr_tight[iTrkBin][ipt]->histDcaVsMassAndMva()->Fill(evt->Mass(id0), evt->Mva(id0), DCA, 1./effks);
         }
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
      for(int imass=0; imass<ana::nMass; imass++){
         for(int ipt=0; ipt<nPt; ipt++){
            ibins.at("ipt") = ipt;
            ibins.at("imass") = imass;
            bool isFillStd = false;
            bool isFillLoose = false;
            bool isFillTight = false;
            auto it = find(ibinsVect_d0.begin(), ibinsVect_d0.end(), ibins);
            while(it!=ibinsVect_d0.end()){
               auto id0 = distance(ibinsVect_d0.begin(), it);
               const bool& passStd = isd0Vect_d0.at(id0).at("passStd");
               const bool& passLoose = isd0Vect_d0.at(id0).at("passLoose");
               const bool& passTight = isd0Vect_d0.at(id0).at("passTight");

               if(passStd && !isFillStd) {
                  isFillStd = true;
                  corr2Dptr[iTrkBin][ipt]->histMult_ass()[imass]->Fill(nMult_ass);
               }
               if(passLoose && !isFillLoose) {
                  isFillLoose = true;
                  corr2Dptr_loose[iTrkBin][ipt]->histMult_ass()[imass]->Fill(nMult_ass);
               }
               if(passTight && !isFillTight) {
                  isFillTight = true;
                  corr2Dptr_tight[iTrkBin][ipt]->histMult_ass()[imass]->Fill(nMult_ass);
               }
               if(isFillTight && isFillStd && isFillTight) break;
               it++;
            }
         }
      }
      
      unsigned int nMult_trg_raw_d0[nPt][ana::nMass]; // Ntrig for mass & pt bins
      double nMult_trg_eff_d0[nPt][ana::nMass]; // eff corrected Ntrig for mass & pt bins

      unsigned int nMult_trg_raw_d0_loose[nPt][ana::nMass]; // Ntrig for mass & pt bins
      double nMult_trg_eff_d0_loose[nPt][ana::nMass]; // eff corrected Ntrig for mass & pt bins

      unsigned int nMult_trg_raw_d0_tight[nPt][ana::nMass]; // Ntrig for mass & pt bins
      double nMult_trg_eff_d0_tight[nPt][ana::nMass]; // eff corrected Ntrig for mass & pt bins

      for(int ipt=0; ipt<nPt; ipt++){
         for(int imass=0; imass<ana::nMass; imass++){
            nMult_trg_raw_d0[ipt][imass] = 0;
            nMult_trg_eff_d0[ipt][imass] = 0;
            nMult_trg_raw_d0_loose[ipt][imass] = 0;
            nMult_trg_eff_d0_loose[ipt][imass] = 0;
            nMult_trg_raw_d0_tight[ipt][imass] = 0;
            nMult_trg_eff_d0_tight[ipt][imass] = 0;
         }
      }

      for(unsigned int id0=0; id0<indexVect_d0.size(); id0++){
         const double& effks = effVect_trg_d0.at(id0);
         const int& ipt = ibinsVect_d0.at(id0).at("ipt");
         const int& imass = ibinsVect_d0.at(id0).at("imass");
         const bool& passStd = isd0Vect_d0.at(id0).at("passStd");
         const bool& passLoose = isd0Vect_d0.at(id0).at("passLoose");
         const bool& passTight = isd0Vect_d0.at(id0).at("passTight");
         if(passStd){
            nMult_trg_raw_d0[ipt][imass] += 1;
            nMult_trg_eff_d0[ipt][imass] += 1./effks;
         }
         if(passLoose){
            nMult_trg_raw_d0_loose[ipt][imass] += 1;
            nMult_trg_eff_d0_loose[ipt][imass] += 1./effks;
         }
         if(passTight){
            nMult_trg_raw_d0_tight[ipt][imass] += 1;
            nMult_trg_eff_d0_tight[ipt][imass] += 1./effks;
         }
      }

      for(int ipt=0; ipt<nPt; ipt++){
         for(int imass=0; imass<ana::nMass; imass++){
            corr2Dptr[iTrkBin][ipt]->histMult_raw_D0()[imass]->Fill(nMult_trg_raw_d0[ipt][imass]);
            corr2Dptr[iTrkBin][ipt]->histMult_eff_D0()[imass]->Fill(nMult_trg_eff_d0[ipt][imass]);

            corr2Dptr_loose[iTrkBin][ipt]->histMult_raw_D0()[imass]->Fill(nMult_trg_raw_d0_loose[ipt][imass]);
            corr2Dptr_loose[iTrkBin][ipt]->histMult_eff_D0()[imass]->Fill(nMult_trg_eff_d0_loose[ipt][imass]);

            corr2Dptr_tight[iTrkBin][ipt]->histMult_raw_D0()[imass]->Fill(nMult_trg_raw_d0_tight[ipt][imass]);
            corr2Dptr_tight[iTrkBin][ipt]->histMult_eff_D0()[imass]->Fill(nMult_trg_eff_d0_tight[ipt][imass]);
         }
      }

      for(unsigned int id0=0; id0<indexVect_d0.size(); id0++){

         const double& effks = effVect_trg_d0.at(id0);
         const int& imass = ibinsVect_d0.at(id0).at("imass");
         const int& ipt = ibinsVect_d0.at(id0).at("ipt");
         const bool& passStd = isd0Vect_d0.at(id0).at("passStd");
         const bool& passLoose = isd0Vect_d0.at(id0).at("passLoose");
         const bool& passTight = isd0Vect_d0.at(id0).at("passTight");

         for(unsigned int iass=0; iass<nMult_ass; iass++){
            const TVector3& pvector_ass = pVect_ass.at(iass);
            const double& effweight_ass = effVect_ass.at(iass);
            if(ana::rejectDaughter_){
               if(fabs(pvector_ass.Eta() - pVect_dau1_d0.at(id0).Eta())<0.03
                     && fabs(pvector_ass.DeltaPhi(pVect_dau1_d0.at(id0)))<0.03)
                  continue;
               if(fabs(pvector_ass.Eta() - pVect_dau2_d0.at(id0).Eta())<0.03
                     && fabs(pvector_ass.DeltaPhi(pVect_dau2_d0.at(id0)))<0.03)
                  continue;
            }

            double deltaEta = pvector_ass.Eta() - pVect_trg_d0.at(id0).Eta();
            double deltaPhi = pvector_ass.DeltaPhi(pVect_trg_d0.at(id0));

            double absDeltaPhi = fabs(deltaPhi);
            double symAbsDeltaPhi = -100;
            if(absDeltaPhi<ana::PI/2.) symAbsDeltaPhi = -1* absDeltaPhi;
            if(absDeltaPhi>=ana::PI/2.) symAbsDeltaPhi = 2*ana::PI - absDeltaPhi;

            if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

            if(passStd){
               corr2Dptr[iTrkBin][ipt]->histSignal_D0()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr[iTrkBin][ipt]->histSignal_D0_old()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);

               corr2Dptr[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(-1*fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(-1*fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
            }
            if(passLoose){
               corr2Dptr_loose[iTrkBin][ipt]->histSignal_D0()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0_loose[ipt][imass]/effks/effweight_ass);
               corr2Dptr_loose[iTrkBin][ipt]->histSignal_D0_old()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0_loose[ipt][imass]/effks/effweight_ass);

               corr2Dptr_loose[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr_loose[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(-1*fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr_loose[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr_loose[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(-1*fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
            }
            if(passTight){
               corr2Dptr_tight[iTrkBin][ipt]->histSignal_D0()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0_tight[ipt][imass]/effks/effweight_ass);
               corr2Dptr_tight[iTrkBin][ipt]->histSignal_D0_old()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0_tight[ipt][imass]/effks/effweight_ass);

               corr2Dptr_tight[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr_tight[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(-1*fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr_tight[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
               corr2Dptr_tight[iTrkBin][ipt]->histSignal_D0_symm()[imass]->Fill(-1*fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
            }
         }
      }

      // mixed events
      // begin filling background
      // if having enough, erase the first one, fill the current, fill the background
      if(pVectList_ass[iz][iTrkBin].size() == ana::nMixedEvts){
         unsigned int nMult_trg_d0 = (unsigned int) indexVect_d0.size();

         for(unsigned int id0=0; id0<nMult_trg_d0; id0++){

            const double& effks = effVect_trg_d0.at(id0);
            const int& imass = ibinsVect_d0.at(id0).at("imass");
            const int& ipt = ibinsVect_d0.at(id0).at("ipt");
            const bool& passStd = isd0Vect_d0.at(id0).at("passStd");
            const bool& passLoose = isd0Vect_d0.at(id0).at("passLoose");
            const bool& passTight = isd0Vect_d0.at(id0).at("passTight");

            // simultaneously read momentum and efficiency
            auto ievt_eff_ass = effVectList_ass[iz][iTrkBin].begin();
            for(auto& ievt_p_ass : pVectList_ass[iz][iTrkBin]){
               unsigned int n_ass = ievt_p_ass.size();
               if(n_ass!=ievt_eff_ass->size()){
                  cout << "wrong" << endl;
                  break;
               }
               for(unsigned int iass=0; iass<n_ass; iass++){
                  const TVector3& pvector_ass = ievt_p_ass.at(iass);
                  const double& effweight_ass = ievt_eff_ass->at(iass);

                  double deltaEta = pvector_ass.Eta() - pVect_trg_d0.at(id0).Eta();
                  double deltaPhi = pvector_ass.DeltaPhi(pVect_trg_d0.at(id0));
                  double absDeltaPhi = fabs(deltaPhi);
                  double symAbsDeltaPhi = -100;
                  if(absDeltaPhi<ana::PI/2.) symAbsDeltaPhi = -1* absDeltaPhi;
                  if(absDeltaPhi>=ana::PI/2.) symAbsDeltaPhi = 2*ana::PI - absDeltaPhi;

                  if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

                  if(passStd){
                     corr2Dptr[iTrkBin][ipt]->histBackground_D0()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr[iTrkBin][ipt]->histBackground_D0_old()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);

                     corr2Dptr[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(-1*fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(-1*fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                  }

                  if(passLoose){
                     corr2Dptr_loose[iTrkBin][ipt]->histBackground_D0()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0_loose[ipt][imass]/effks/effweight_ass);
                     corr2Dptr_loose[iTrkBin][ipt]->histBackground_D0_old()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0_loose[ipt][imass]/effks/effweight_ass);

                     corr2Dptr_loose[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr_loose[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(-1*fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr_loose[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr_loose[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(-1*fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                  }

                  if(passTight){
                     corr2Dptr_tight[iTrkBin][ipt]->histBackground_D0()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0_tight[ipt][imass]/effks/effweight_ass);
                     corr2Dptr_tight[iTrkBin][ipt]->histBackground_D0_old()[imass]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0_tight[ipt][imass]/effks/effweight_ass);

                     corr2Dptr_tight[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr_tight[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(-1*fabs(deltaEta), absDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr_tight[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                     corr2Dptr_tight[iTrkBin][ipt]->histBackground_D0_symm()[imass]->Fill(-1*fabs(deltaEta), symAbsDeltaPhi, 
                        1./nMult_trg_eff_d0[ipt][imass]/effks/effweight_ass);
                  }
               }
               ievt_eff_ass++;
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

      pVect_trg_d0.clear();
      effVect_trg_d0.clear();
      indexVect_d0.clear();
      isd0Vect_d0.clear();
      ibinsVect_d0.clear();

      pVect_dau1_d0.clear();
      pVect_dau2_d0.clear();

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
         outName = TString::Format("%s/fout_%s_d0ana_ntrk_pT%.1f-%.1f_y%.1f-%.1f.root", prefix.c_str(), 
               datalist.c_str(), cuts.pTMin, cuts.pTMax, cuts.yMin, cuts.yMax);
      else
         outName = TString::Format("fout_%s_d0ana_ntrk_pT%.1f-%.1f_y%.1f-%.1f.root", datalist.c_str(),
               cuts.pTMin, cuts.pTMax, cuts.yMin, cuts.yMax);
   }
   else return -1;

   TFile fout(outName.Data(), "recreate");
   fout.cd();

   // start writing output
   hEvt->Write();
   hMult->Write();
   hNtrkofflineVsNtrkgood->Write();
   for(int iTrkBin=0; iTrkBin<nTrkBin; iTrkBin++){
      hNtrk_D0[iTrkBin]->Write();
   }

   for(int iTrkBin=0; iTrkBin<nTrkBin; iTrkBin++){
      for(int ipt=0; ipt<nPt; ipt++){
         corr2Dptr[iTrkBin][ipt]->write();
         corr2Dptr_loose[iTrkBin][ipt]->write();
         corr2Dptr_tight[iTrkBin][ipt]->write();
      }
   }

   delete hEvt;
   delete hMult;
   delete hNtrkofflineVsNtrkgood;
   for(int iTrkBin=0; iTrkBin<nTrkBin; iTrkBin++){
      delete hNtrk_D0[iTrkBin];
   }

   for(int iTrkBin=0; iTrkBin<nTrkBin; iTrkBin++){
      for(int ipt=0; ipt<nPt; ipt++){
         delete corr2Dptr[iTrkBin][ipt];
         delete corr2Dptr_loose[iTrkBin][ipt];
         delete corr2Dptr_tight[iTrkBin][ipt];
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
   evt->SetBranchStatus("eta", 1);
   evt->SetBranchStatus("phi", 1);
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
      event->GetBranchStatus("eta")&&
      event->GetBranchStatus("phi")&&
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
   bool passPt = event->Pt(icand) < cuts.pTMax && event->Pt(icand) >= cuts.pTMin;
   bool passY = event->Y(icand) < cuts.yMax && event->Y(icand) >= cuts.yMin;
   return //passEta;
      passY &&
      passPt;
}

bool passD0MVA(const int& trigger, Event* event, const int& icand, 
      const bool& isPrompt, const string& type)
{
   if(type == "std"){
      if(isPrompt){
         switch(trigger){
            case 0: return ana::pass_pPb2016_8TeV_PD0_MVA(event->Pt(icand), event->Mva(icand));
                  break;
            case 1: return ana::pass_pPb2016_8TeV_PD0_MVA(event->Pt(icand), event->Mva(icand));
                  break;
            case 2: return ana::pass_pPb2016_8TeV_PD0_MVA(event->Pt(icand), event->Mva(icand));
                  break;
            case 3: return ana::pass_pp2018_13TeV_PD0_MVA(event->Mva(icand));
                  break;
            case 4: return ana::pass_pp2018_13TeV_PD0_MVA(event->Mva(icand));
                  break;
         }
      } else{
         // non-prompt would not be measured
         return false;
      }
   }

   if(type == "loose"){
      if(isPrompt){
         switch(trigger){
            case 0: return ana::pass_pPb2016_8TeV_PD0_looseMVA(event->Pt(icand), event->Mva(icand));
                  break;
            case 1: return ana::pass_pPb2016_8TeV_PD0_looseMVA(event->Pt(icand), event->Mva(icand));
                  break;
            case 2: return ana::pass_pPb2016_8TeV_PD0_looseMVA(event->Pt(icand), event->Mva(icand));
                  break;
            default: return false;
         }
      } else{
         // non-prompt would not be measured
         return false;
      }
   }

   if(type == "tight"){
      if(isPrompt){
         switch(trigger){
            case 0: return ana::pass_pPb2016_8TeV_PD0_tightMVA(event->Pt(icand), event->Mva(icand));
                  break;
            case 1: return ana::pass_pPb2016_8TeV_PD0_tightMVA(event->Pt(icand), event->Mva(icand));
                  break;
            case 2: return ana::pass_pPb2016_8TeV_PD0_tightMVA(event->Pt(icand), event->Mva(icand));
                  break;
            default: return false;
         }
      } else{
         // non-prompt would not be measured
         return false;
      }
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

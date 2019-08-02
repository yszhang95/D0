const double deltaEtaBound = 1;
const int nDca= 2;

const float dcacut[2] = {
   0.012, 0.010
}; // npt

TH1D* proj1D_longrange(TH2*, TH2*, const char*);
std::pair<double, double> proj1D_shortrange_yields(TH2*, TH2*, const char*, TCanvas* c, const int& ipad, TH1D*);
std::pair<double, double> proj1D_longrange_yields(TH2*, TH2*, const char*, TCanvas* c, const int& ipad, TH1D*);

std::pair<double, double> calYields(const pair<double, double>&, const pair<double,double>&);

string ouput_prefix(const char* name);

TF1 draw1D_longrange(TH1*, const char*,
const char*, const char* , const char*, const char*,
TCanvas*, const int& ipad);

void proj1D_HM_NPD0_Process_new(const char* input_d0= "",
      const char* input_ref = "", 
      const string dataset="", 
      const char* input_d0_low = "",
      const char* input_ref_low = "",
      const float pTMin=0.0, const float pTMax=0.0,
      const float yMin =0.0, const float yMax =0.0
      )
{
   string str = ouput_prefix(input_d0);
   //gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   if(dataset!="PAHM1-6" && dataset!="PAMB"){
      std::cerr << "wrong dataset" << std::endl;
      return;
   }

   vector<double> ptbin;
   ptbin.push_back(2.0);
   ptbin.push_back(5.0);
   ptbin.push_back(8.0);

   const int nPt = 2;

   double v2_NPD0[ana::nMass][nPt][nDca];
   double v2_NPD0_err[ana::nMass][nPt][nDca];

   double V2_NPD0[ana::nMass][nPt][nDca];
   double V2_NPD0_err[ana::nMass][nPt][nDca];
   double V2_REF;
   double V2_REF_err;

   double V2_NPD0_low[ana::nMass][nPt][nDca];
   double V2_NPD0_low_err[ana::nMass][nPt][nDca];
   double V2_REF_low;
   double V2_REF_low_err;

   double N_ass[ana::nMass][nPt][nDca];
   double N_ass_err[ana::nMass][nPt][nDca];
   double N_ass_low[ana::nMass][nPt][nDca];
   double N_ass_low_err[ana::nMass][nPt][nDca];

   double N_ass_Ref;
   double N_ass_Ref_err;
   double N_ass_Ref_low;
   double N_ass_Ref_low_err;

   double yields_jet[ana::nMass][nPt][nDca];
   double yields_jet_err[ana::nMass][nPt][nDca];
   double yields_jet_low[ana::nMass][nPt][nDca];
   double yields_jet_low_err[ana::nMass][nPt][nDca];

   double yields_jet_ref;
   double yields_jet_ref_err;
   double yields_jet_ref_low;
   double yields_jet_ref_low_err;

   double V2_Sub_NPD0[ana::nMass][nPt][nDca];
   double V2_Sub_NPD0_err[ana::nMass][nPt][nDca];
   double V2_Sub_REF;
   double V2_Sub_REF_err;

   double v2_Sub_NPD0[ana::nMass][nPt][nDca];
   double v2_Sub_NPD0_err[ana::nMass][nPt][nDca];

   TFile* f1 = new TFile(input_d0);
   TFile* f2 = new TFile(input_ref);
   TFile* f3 = new TFile(input_d0_low);
   TFile* f4 = new TFile(input_ref_low);

   TH3D* hDcaVsMassAndMva[nPt]; 
   
   //high N
   TH2D* h2DSignal_D0[ana::nMass][nPt][nDca];
   TH2D* h2DBackground_D0[ana::nMass][nPt][nDca];
   TH1D* hMult_raw_D0[ana::nMass][nPt][nDca]; // wrong normalized constant
   TH1D* hMass_D0[ana::nMass][nPt][nDca];

   //low N
   TH2D* h2DSignal_D0_low[ana::nMass][nPt][nDca];
   TH2D* h2DBackground_D0_low[ana::nMass][nPt][nDca];
   TH1D* hMult_raw_D0_low[ana::nMass][nPt][nDca]; // wrong normalized constant
   TH1D* hMass_D0_low[ana::nMass][nPt][nDca];

   //high N ref
   TH2D* h2DSignal_Ref = (TH2D*) f2->Get("hSignal_Ref");
   TH2D* h2DBackground_Ref = (TH2D*) f2->Get("hBackground_Ref");
   TH1D* hMult_Ref = (TH1D*) f2->Get("hMult");

   //low N ref
   TH2D* h2DSignal_Ref_low = (TH2D*) f4->Get("hSignal_Ref");
   TH2D* h2DBackground_Ref_low = (TH2D*) f4->Get("hBackground_Ref");
   TH1D* hMult_Ref_low = (TH1D*) f4->Get("hMult");

   // read high N
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            h2DSignal_D0[imass][ipt][idca] = 
               (TH2D*) f1->Get(Form("hSignal_mass%d_pt%d_dca%d", imass, ipt, idca));
            h2DBackground_D0[imass][ipt][idca] = 
               (TH2D*) f1->Get(Form("hBackground_mass%d_pt%d_dca%d", imass, ipt, idca));
            hMult_raw_D0[imass][ipt][idca] = 
               (TH1D*) f1->Get(Form("hMult_raw_D0_mass%d_pt%d_dca%d", imass, ipt, idca));
            hMass_D0[imass][ipt][idca] = 
               (TH1D*) f1->Get(Form("hMassD0_mass%d_pt%d_dca%d", imass, ipt, idca));
         }
      }
   }

   // read low N
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            h2DSignal_D0_low[imass][ipt][idca] = 
               (TH2D*) f3->Get(Form("hSignal_mass%d_pt%d_dca%d", imass, ipt, idca));
            h2DBackground_D0_low[imass][ipt][idca] = 
               (TH2D*) f3->Get(Form("hBackground_mass%d_pt%d_dca%d", imass, ipt, idca));
            hMult_raw_D0_low[imass][ipt][idca] = 
               (TH1D*) f3->Get(Form("hMult_raw_D0_mass%d_pt%d_dca%d", imass, ipt, idca));
            hMass_D0_low[imass][ipt][idca] = 
               (TH1D*) f1->Get(Form("hMassD0_mass%d_pt%d_dca%d", imass, ipt, idca));
         }
      }
   }

   // scaled by number of events
   // high N
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            double nMult_D0= hMult_raw_D0[imass][ipt][idca]->Integral(2, 10000);
            h2DSignal_D0[imass][ipt][idca]->Scale(1./nMult_D0);
            //cout << nMult_D0 << endl;
         }
      }
   }

   // low N
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            double nMult_D0_low= hMult_raw_D0_low[imass][ipt][idca]->Integral(2, 10000);
            h2DSignal_D0_low[imass][ipt][idca]->Scale(1./nMult_D0_low);
         }
      }
   }

   // high N ref and low N ref
   if(true){
      long int nevents_high = hMult_Ref->Integral(3, 100000);
      h2DSignal_Ref->Scale(1./nevents_high);
      long int nevents_low= hMult_Ref_low->Integral(3, 100000);
      h2DSignal_Ref_low->Scale(1./nevents_low);
   }

   // read other histograms
   for(int ipt=0; ipt<nPt; ipt++){
      hDcaVsMassAndMva[ipt] = (TH3D*) f1->Get(Form("hDcaVsMassAndMva_pt%d", ipt));
   }

   TH1D* hPt[nPt][nDca];
   TH1D* hKET[nPt][nDca];
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         f1->GetObject(Form("hPt_pt%d_dca%d", ipt, idca), hPt[ipt][idca]);
         f1->GetObject(Form("hKET_pt%d_dca%d", ipt, idca), hKET[ipt][idca]);
         if(!hPt[ipt][idca] || !hKET[ipt][idca]){
            cout << "cannot find hPt or hKET" << endl;
            return;
         }
      }
   }

   // proj1D long range, high N
   TH1D* hDeltaPhi[ana::nMass][nPt][nDca];
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            hDeltaPhi[imass][ipt][idca] = 
               proj1D_longrange(h2DSignal_D0[imass][ipt][idca], 
                     h2DBackground_D0[imass][ipt][idca], 
                     Form("deltaPhi_mass%d_pt%d_dca%d", imass, ipt, idca));
         }
      }
   }

   // proj1D long range, low N
   TH1D* hDeltaPhi_low[ana::nMass][nPt][nDca];
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            hDeltaPhi_low[imass][ipt][idca] = 
               proj1D_longrange(h2DSignal_D0_low[imass][ipt][idca], 
                  h2DBackground_D0_low[imass][ipt][idca], 
                  Form("deltaPhi_low_mass%d_pt%d_dca%d", imass, ipt, idca));
         }
      }
   }

   // proj1D long range, ref
   TH1D* hDeltaPhi_Ref = proj1D_longrange(h2DSignal_Ref, h2DBackground_Ref, "deltaPhi_Ref");
   TH1D* hDeltaPhi_Ref_low = proj1D_longrange(h2DSignal_Ref_low, h2DBackground_Ref_low, "deltaPhi_Ref_low");
   
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";

   TCanvas* c_deltaPhi[nPt][nDca];
   // high N
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         c_deltaPhi[ipt][idca] = new TCanvas(Form("c_detalPhi_pt%d_dca%d", ipt, idca), "", 550*3, 550*5);
         c_deltaPhi[ipt][idca]->Divide(3, 5);
      }
   }

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            std::string ntrk;
            if(dataset == "PAHM1-6") ntrk = "185< N_{trk}^{offline} < 250";
            string dcacut_str;
            if(idca==0) dcacut_str = string(Form("DCA < %.3f cm", dcacut[ipt]));
            if(idca==1) dcacut_str = string(Form("DCA > %.3f cm", dcacut[ipt]));
            auto func = draw1D_longrange(hDeltaPhi[imass][ipt][idca],
                  Form("../plots/v2vspt/npd0ana1/y%.1f/%s/%s_deltaPhi_mass%d_pt%d_dca%d.png",
                     yMax, dataset.c_str(), str.c_str(), imass, ipt, idca),
                  ntrk.c_str(),
                  Form("%.1f<p_{T}<%.1fGeV, |y|<%.1f", ptbin[ipt], ptbin[ipt+1], yMax),
                  Form("%.3f<mass<%.3f, |#Delta#eta|>1", ana::massbin[imass], ana::massbin[imass+1]),
                  dcacut_str.c_str(),
                  c_deltaPhi[ipt][idca], imass+1
               );
            if(fabs(yMax+yMin)>1e-5) {
                  cout << "|y|<\%.1f is wrong" << endl;
                  return;
                  }
            V2_NPD0[imass][ipt][idca] = func.GetParameter(2);
            V2_NPD0_err[imass][ipt][idca] = func.GetParError(2);
            N_ass[imass][ipt][idca] = func.GetParameter(0);
            N_ass_err[imass][ipt][idca] = func.GetParError(0);
            cout << N_ass[imass][ipt][idca] << endl;
         }
      }
   }
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         c_deltaPhi[ipt][idca]->Print(
            Form("../plots/v2vspt/npd0ana1/y%.1f/%s/%s_deltaPhi_pt%d_dca%d.png",
               yMax, dataset.c_str(), str.c_str(), ipt, idca)
         );
         c_deltaPhi[ipt][idca]->Print(
            Form("../plots/v2vspt/npd0ana1/y%.1f/%s/%s_deltaPhi_pt%d_dca%d.pdf",
               yMax, dataset.c_str(), str.c_str(), ipt, idca)
         );
         delete c_deltaPhi[ipt][idca];
      }
   }

   // low N
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         c_deltaPhi[ipt][idca] = new TCanvas(Form("c_detalPhi_pt%d_dca%d", ipt, idca), "", 550*3, 550*5);
         c_deltaPhi[ipt][idca]->Divide(3, 5);
      }
   }

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            std::string ntrk;
            ntrk = "N_{trk}^{offline} < 35";
            string dcacut_str;
            if(idca==0) dcacut_str = string(Form("DCA < %.3f cm", dcacut[ipt]));
            if(idca==1) dcacut_str = string(Form("DCA > %.3f cm", dcacut[ipt]));
            auto func = draw1D_longrange(hDeltaPhi_low[imass][ipt][idca],
                  Form("../plots/v2vspt/npd0ana1/y%.1f/%s/%s_deltaPhi_mass%d_pt%d_dca%d.png",
                     yMax, dataset.c_str(), str.c_str(), imass, ipt, idca),
                  ntrk.c_str(),
                  Form("%.1f<p_{T}<%.1fGeV, |y|<%.1f", ptbin[ipt], ptbin[ipt+1], yMax),
                  Form("%.3f<mass<%.3f, |#Delta#eta|>1", ana::massbin[imass], ana::massbin[imass+1]),
                  dcacut_str.c_str(),
                  c_deltaPhi[ipt][idca], imass+1
               );
            if(fabs(yMax+yMin)>1e-5) {
                  cout << "|y|<\%.1f is wrong" << endl;
                  return;
                  }
            V2_NPD0_low[imass][ipt][idca] = func.GetParameter(2);
            V2_NPD0_low_err[imass][ipt][idca] = func.GetParError(2);
            N_ass_low[imass][ipt][idca] = func.GetParameter(0);
            N_ass_low_err[imass][ipt][idca] = func.GetParError(0);
            cout << N_ass_low[imass][ipt][idca] << endl;
         }
      }
   }
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         c_deltaPhi[ipt][idca]->Print(
            Form("../plots/v2vspt/npd0ana1/y%.1f/%s/%s_deltaPhi_low_pt%d_dca%d.png",
               yMax, dataset.c_str(), str.c_str(), ipt, idca)
         );
         c_deltaPhi[ipt][idca]->Print(
            Form("../plots/v2vspt/npd0ana1/y%.1f/%s/%s_deltaPhi_low_pt%d_dca%d.pdf",
               yMax, dataset.c_str(), str.c_str(), ipt, idca)
         );
         delete c_deltaPhi[ipt][idca];
      }
   }


   // high N, ref
   if(true){
      TCanvas* c_deltaPhi_Ref = new TCanvas("c_deltaPhi_Ref", "", 550, 550);
      c_deltaPhi_Ref->Divide(1, 1);
      std::string ntrk;
      if(dataset == "PAHM1-6") ntrk = "185< N_{trk}^{offline} < 250";
      auto func_ref = draw1D_longrange(hDeltaPhi_Ref,
            Form("../plots/v2vspt/npd0ana1/%s_ref_deltaPhi.png", 
            str.c_str()),
            ntrk.c_str(),
            "0.3<p_{T}<3.0GeV |#eta|<2.4", 
            "|#Delta#eta|>1", "",
            c_deltaPhi_Ref, 1
         );
      V2_REF = func_ref.GetParameter(2);
      V2_REF_err = func_ref.GetParError(2);
      N_ass_Ref = func_ref.GetParameter(0);
      N_ass_Ref_err = func_ref.GetParError(0);
      c_deltaPhi_Ref->Print(
            Form("../plots/v2vspt/npd0ana1/%s_ref_deltaPhi.png", str.c_str()));
      c_deltaPhi_Ref->Print(
            Form("../plots/v2vspt/npd0ana1/%s_ref_deltaPhi.pdf", str.c_str()));
      cout << N_ass_Ref << endl;
      delete c_deltaPhi_Ref;
   }

   // log N, ref
   if(true){
      TCanvas* c_deltaPhi_Ref = new TCanvas("c_deltaPhi_Ref", "", 550, 550);
      c_deltaPhi_Ref->Divide(1, 1);
      std::string ntrk;
      ntrk = "N_{trk}^{offline} < 35";
      auto func_ref = draw1D_longrange(hDeltaPhi_Ref_low,
            Form("../plots/v2vspt/npd0ana1/%s_ref_deltaPhi.png", 
            str.c_str()),
            ntrk.c_str(),
            "0.3<p_{T}<3.0GeV |#eta|<2.4", 
            "|#Delta#eta|>1", "",
            c_deltaPhi_Ref, 1
         );
      V2_REF_low = func_ref.GetParameter(2);
      V2_REF_low_err = func_ref.GetParError(2);
      N_ass_Ref_low = func_ref.GetParameter(0);
      N_ass_Ref_low_err = func_ref.GetParError(0);
      c_deltaPhi_Ref->Print(
            Form("../plots/v2vspt/npd0ana1/%s_ref_deltaPhi_low.png", str.c_str()));
      c_deltaPhi_Ref->Print(
            Form("../plots/v2vspt/npd0ana1/%s_ref_deltaPhi_low.pdf", str.c_str()));
      cout << N_ass_Ref_low << endl;
      delete c_deltaPhi_Ref;
   }

   TGraphErrors* g_v2[nPt][nDca];
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         g_v2[ipt][idca] = new TGraphErrors(ana::nMass);
         g_v2[ipt][idca]->SetName(Form("g_v2_DCA_pt%d_dca%d", ipt, idca));
         for(int imass=0; imass<ana::nMass; imass++){
            double temp = V2_NPD0[imass][ipt][idca]/ sqrt(V2_REF);
            double temp_err = temp* sqrt(
                  pow(V2_NPD0_err[imass][ipt][idca]/ V2_NPD0[imass][ipt][idca], 2)
               + pow(0.5*V2_REF_err/ V2_REF, 2));

            v2_NPD0[imass][ipt][idca] = temp;
            v2_NPD0_err[imass][ipt][idca] = temp_err;

            g_v2[ipt][idca]->SetPoint(imass, 
                  hMass_D0[imass][ipt][idca]->GetMean(), v2_NPD0[imass][ipt][idca]);
            g_v2[ipt][idca]->SetPointError(imass, 0, v2_NPD0_err[imass][ipt][idca]);
         }
      }
   }

   TH1D* hMass_DCA[nPt][nDca];
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         f1->GetObject(Form("hMass_pt%d_dca%d", ipt, idca), hMass_DCA[ipt][idca]);
      }
   }

   TH1D* hMass_DCA_low[nPt][nDca];
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         f3->GetObject(Form("hMass_pt%d_dca%d", ipt, idca), hMass_DCA_low[ipt][idca]);
         hMass_DCA_low[ipt][idca]->SetName(Form("hMass_low_pt%d_dca%d", ipt, idca));
      }
   }

   TFile fout_v2(Form("%s_v2.root", input_d0), "recreate");

   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         g_v2[ipt][idca]->Write();
         hMass_DCA[ipt][idca]->Write();
         hPt[ipt][idca]->Write();
         hKET[ipt][idca]->Write();
      }
   }

   // fit the yields and then subtract non flow
   // low mult ref
   if(true){
      TCanvas* c_sry_ref_low = new TCanvas("c_sry_ref_low", "", 550, 550);
      c_sry_ref_low->Divide(1,1);
      TCanvas* c_lry_ref_low = new TCanvas("c_lry_ref_low", "", 550, 550);
      c_lry_ref_low->Divide(1,1);
      TH1D* h_sry_ref_low;
      TH1D* h_lry_ref_low;
      auto sr_pair = proj1D_shortrange_yields(h2DSignal_Ref_low, h2DBackground_Ref_low, 
            "", c_sry_ref_low, 1, h_sry_ref_low);
      auto lr_pair = proj1D_longrange_yields(h2DSignal_Ref_low, h2DBackground_Ref_low, 
            "", c_lry_ref_low, 1, h_lry_ref_low);
      auto yields_pair = calYields(sr_pair, lr_pair);
      yields_jet_ref_low = yields_pair.first;
      yields_jet_ref_low_err = yields_pair.second;
      c_sry_ref_low->Print(
            Form("../plots/v2vspt/%s_ref_low_sr.png", str.c_str())
            );
      c_lry_ref_low->Print(
            Form("../plots/v2vspt/%s_ref_low_lr.png", str.c_str())
      );
      delete c_sry_ref_low;
      delete c_lry_ref_low;
   }

   // high mult ref
   if(true){
      TCanvas* c_sry_ref = new TCanvas("c_sry_ref", "", 550, 550);
      c_sry_ref->Divide(1,1);
      TCanvas* c_lry_ref = new TCanvas("c_lry_ref", "", 550, 550);
      c_lry_ref->Divide(1,1);
      TH1D* h_sry_ref;
      TH1D* h_lry_ref;
      auto sr_pair = proj1D_shortrange_yields(h2DSignal_Ref, h2DBackground_Ref, 
            "", c_sry_ref, 1, h_sry_ref);
      auto lr_pair = proj1D_longrange_yields(h2DSignal_Ref, h2DBackground_Ref, 
            "", c_lry_ref, 1, h_lry_ref);
      auto yields_pair = calYields(sr_pair, lr_pair);
      yields_jet_ref = yields_pair.first;
      yields_jet_ref_err = yields_pair.second;
      c_sry_ref->Print(
            Form("../plots/v2vspt/%s_ref_sr.png", str.c_str())
            );
      c_lry_ref->Print(
            Form("../plots/v2vspt/%s_ref_lr.png", str.c_str())
      );
      delete c_sry_ref;
      delete c_lry_ref;
   }
   // higher mult, d0
   TCanvas* c_sry[nPt][nDca];
   TCanvas* c_lry[nPt][nDca];
   TH1D* h_sry[nPt][nDca];
   TH1D* h_lry[nPt][nDca];
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         c_sry[ipt][idca] = new TCanvas(Form("c_sry_%d_%d", ipt, idca), "", 550*3, 550*5);
         c_sry[ipt][idca]->Divide(3, 5);
         c_lry[ipt][idca] = new TCanvas(Form("c_lry_%d_%d", ipt, idca), "", 550*3, 550*5);
         c_lry[ipt][idca]->Divide(3, 5);
      }
   }
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         for(int imass=0; imass<ana::nMass; imass++){
            cout <<ipt << idca << endl;
            string str = ouput_prefix(input_d0);
            string dcacut_str;
            if(idca==0) dcacut_str = string(Form("DCA < %.3f cm", dcacut[ipt]));
            if(idca==1) dcacut_str = string(Form("DCA > %.3f cm", dcacut[ipt]));
            auto sr_pair = proj1D_shortrange_yields(h2DSignal_D0[imass][ipt][idca], h2DBackground_D0[imass][ipt][idca], 
                  dcacut_str.c_str(),
               c_sry[ipt][idca], imass+1, h_sry[ipt][idca]
               );
            auto lr_pair = proj1D_longrange_yields(h2DSignal_D0[imass][ipt][idca], h2DBackground_D0[imass][ipt][idca], 
                  dcacut_str.c_str(),
               c_lry[ipt][idca], imass+1, h_lry[ipt][idca]
                  );
            auto yields_pair = calYields(sr_pair, lr_pair);
            yields_jet[imass][ipt][idca] = yields_pair.first;
            yields_jet_err[imass][ipt][idca] = yields_pair.second;
         }
      }
   }
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         c_sry[ipt][idca]->Print(
               Form("../plots/v2vspt/%s_deltaPhi_pt%d_dca%d_sr.png", str.c_str(), ipt, idca)
            );
         c_lry[ipt][idca]->Print(
               Form("../plots/v2vspt/%s_deltaPhi_pt%d_dca%d_lr.png", str.c_str(), ipt, idca)
            );
         delete c_sry[ipt][idca];
         delete c_lry[ipt][idca];
      }
   }

   // low mult, d0
   TCanvas* c_sry_low[nPt][nDca];
   TCanvas* c_lry_low[nPt][nDca];
   TH1D* h_sry_low[nPt][nDca];
   TH1D* h_lry_low[nPt][nDca];
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         c_sry_low[ipt][idca] = new TCanvas(Form("c_sry_low_%d_%d", ipt, idca), "", 550*3, 550*5);
         c_sry_low[ipt][idca]->Divide(3, 5);
         c_lry_low[ipt][idca] = new TCanvas(Form("c_lry_low_%d_%d", ipt, idca), "", 550*3, 550*5);
         c_lry_low[ipt][idca]->Divide(3, 5);
      }
   }
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         for(int imass=0; imass<ana::nMass; imass++){
            string str = ouput_prefix(input_d0);
            string dcacut_str;
            if(idca==0) dcacut_str = string(Form("DCA < %.3f cm", dcacut[ipt]));
            if(idca==1) dcacut_str = string(Form("DCA > %.3f cm", dcacut[ipt]));
            auto sr_pair = proj1D_shortrange_yields(h2DSignal_D0_low[imass][ipt][idca], h2DBackground_D0_low[imass][ipt][idca], 
                  dcacut_str.c_str(),
               c_sry_low[ipt][idca], imass+1, h_sry_low[ipt][idca]
               );
            auto lr_pair = proj1D_longrange_yields(h2DSignal_D0_low[imass][ipt][idca], h2DBackground_D0_low[imass][ipt][idca], 
                  dcacut_str.c_str(),
               c_lry_low[ipt][idca], imass+1, h_lry_low[ipt][idca]
                  );
            auto yields_pair = calYields(sr_pair, lr_pair);
            yields_jet_low[imass][ipt][idca] = yields_pair.first;
            yields_jet_low_err[imass][ipt][idca] = yields_pair.second;
         }
      }
   }
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         c_sry_low[ipt][idca]->Print(
               Form("../plots/v2vspt/%s_deltaPhi_low_pt%d_dca%d_sr.png", str.c_str(), ipt, idca)
            );
         c_lry_low[ipt][idca]->Print(
               Form("../plots/v2vspt/%s_deltaPhi_low_pt%d_dca%d_lr.png", str.c_str(), ipt, idca)
            );
         delete c_sry_low[ipt][idca];
         delete c_lry_low[ipt][idca];
      }
   }

   TFile fraw(Form("%s_raw.root", input_d0), "recreate");

   TGraphErrors* g_V2[nPt][nDca]; // signal + background, nMass points, need to do fitting, filled
   TGraphErrors* g_V2_Ref; // signal, 1 point, done
   TGraphErrors* g_V2_low[nPt][nDca]; // signal + background, nMass points, need to do fitting, filled, but to define hMass_D0_low
   TGraphErrors* g_V2_Ref_low; // signal, 1 point, done
   TGraphErrors* g_Jets[nPt][nDca]; // signal + background, nMass points, need to do fitting, filled
   TGraphErrors* g_Jets_low[nPt][nDca]; // signal + background, nMass points, need to do fitting, filled
   TGraphErrors* g_Jets_Ref; // signal, 1 point, done
   TGraphErrors* g_Jets_Ref_low; // signal, 1 point, done
   TGraphErrors* g_Nass[nPt][nDca]; // signal, need to do average, nMass points, filled
   TGraphErrors* g_Nass_low[nPt][nDca]; // signal, need to do average, nMass points, filled
   TGraphErrors* g_Nass_Ref; // signal, 1 point, done
   TGraphErrors* g_Nass_Ref_low; // signal, 1 point, done

   // trg d0
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         g_V2[ipt][idca] = new TGraphErrors(ana::nMass);
         g_Jets[ipt][idca] = new TGraphErrors(ana::nMass);
         g_Nass[ipt][idca] = new TGraphErrors(ana::nMass);
         g_V2_low[ipt][idca] = new TGraphErrors(ana::nMass);
         g_Jets_low[ipt][idca] = new TGraphErrors(ana::nMass);
         g_Nass_low[ipt][idca] = new TGraphErrors(ana::nMass);
      }
   }
   // trg ref
   g_V2_Ref = new TGraphErrors(1);
   g_Jets_Ref = new TGraphErrors(1);
   g_Nass_Ref = new TGraphErrors(1);
   g_V2_Ref_low = new TGraphErrors(1);
   g_Jets_Ref_low = new TGraphErrors(1);
   g_Nass_Ref_low = new TGraphErrors(1);

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            g_V2[ipt][idca]->SetPoint(imass, hMass_D0[imass][ipt][idca]->GetMean(), V2_NPD0[imass][ipt][idca]);
            g_V2[ipt][idca]->SetPointError(imass, 0, V2_NPD0_err[imass][ipt][idca]);

            g_Jets[ipt][idca]->SetPoint(imass, hMass_D0[imass][ipt][idca]->GetMean(), yields_jet[imass][ipt][idca]);
            g_Jets[ipt][idca]->SetPointError(imass, 0, yields_jet_err[imass][ipt][idca]);

            g_Nass[ipt][idca]->SetPoint(imass, hMass_D0[imass][ipt][idca]->GetMean(), N_ass[imass][ipt][idca]);
            g_Nass[ipt][idca]->SetPointError(imass, 0, N_ass_err[imass][ipt][idca]);

            g_V2_low[ipt][idca]->SetPoint(imass, hMass_D0_low[imass][ipt][idca]->GetMean(), V2_NPD0_low[imass][ipt][idca]);
            g_V2_low[ipt][idca]->SetPointError(imass, 0, V2_NPD0_low_err[imass][ipt][idca]);

            g_Jets_low[ipt][idca]->SetPoint(imass, hMass_D0_low[imass][ipt][idca]->GetMean(), yields_jet_low[imass][ipt][idca]);
            g_Jets_low[ipt][idca]->SetPointError(imass, 0, yields_jet_low_err[imass][ipt][idca]);

            g_Nass_low[ipt][idca]->SetPoint(imass, hMass_D0_low[imass][ipt][idca]->GetMean(), N_ass_low[imass][ipt][idca]);
            g_Nass_low[ipt][idca]->SetPointError(imass, 0, N_ass_low_err[imass][ipt][idca]);
         }
      }
   }
   // high N, ref
   g_V2_Ref->SetPoint(0, 1, V2_REF);
   g_V2_Ref->SetPointError(0, 1, V2_REF_err);

   g_Jets_Ref->SetPoint(0, 1, yields_jet_ref);
   g_Jets_Ref->SetPointError(0, 1, yields_jet_ref_err);

   g_Nass_Ref->SetPoint(0, 1, N_ass_Ref);
   g_Nass_Ref->SetPointError(0, 1, N_ass_Ref_err);
   // low N, ref
   g_V2_Ref_low->SetPoint(0, 1, V2_REF_low);
   g_V2_Ref_low->SetPointError(0, 1, V2_REF_low_err);

   g_Jets_Ref_low->SetPoint(0, 1, yields_jet_ref_low);
   g_Jets_Ref_low->SetPointError(0, 1, yields_jet_ref_low_err);

   g_Nass_Ref_low->SetPoint(0, 1, N_ass_Ref_low);
   g_Nass_Ref_low->SetPointError(0, 1, N_ass_Ref_low_err);

   cout <<"test" << endl;

   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         g_V2[ipt][idca]->Write(Form("V2vsMass_pt%d_dca%d", ipt, idca)); // signal + background, nMass points, need to do fitting, filled
         g_Jets[ipt][idca]->Write(Form("JetsvsMass_pt%d_dca%d", ipt, idca)); // signal + background, nMass points, need to do fitting, filled
         g_Nass[ipt][idca]->Write(Form("NassvsMass_pt%d_dca%d", ipt, idca)); // signal, need to do average, nMass points, filled

         g_V2_low[ipt][idca]->Write(Form("V2_lowvsMass_pt%d_dca%d", ipt, idca)); // signal + background, nMass points, need to do fitting, filled
         g_Jets_low[ipt][idca]->Write(Form("Jets_lowvsMass_pt%d_dca%d", ipt, idca)); // signal + background, nMass points, need to do fitting, filled
         g_Nass_low[ipt][idca]->Write(Form("Nass_lowvsMass_pt%d_dca%d", ipt, idca)); // signal, need to do average, nMass points, filled

         hMass_DCA[ipt][idca]->Write();
         hPt[ipt][idca]->Write();
      }
   }
   g_V2_Ref_low->Write("V2_Ref_low"); // signal, 1 point, done
   g_Jets_Ref_low->Write("Jets_ref_low"); // signal, 1 point, done
   g_Nass_Ref_low->Write("Nass_ref_low"); // signal, 1 point, done
   g_V2_Ref->Write("V2_Ref"); // signal, 1 point, done
   g_Jets_Ref->Write("Jets_ref"); // signal, 1 point, done
   g_Nass_Ref->Write("Nass_ref"); // signal, 1 point, done

   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         hMass_DCA_low[ipt][idca]->Write();
      }
   }

   fraw.Close();

   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         delete g_V2[ipt][idca]; // signal + background, nMass points, need to do fitting, filled
         delete g_V2_low[ipt][idca]; // signal + background, nMass points, need to do fitting, filled
         delete g_Jets[ipt][idca]; // signal + background, nMass points, need to do fitting, filled
         delete g_Jets_low[ipt][idca]; // signal + background, nMass points, need to do fitting, filled
         delete g_Nass[ipt][idca]; // signal, need to do average, nMass points, filled
         delete g_Nass_low[ipt][idca]; // signal, need to do average, nMass points, filled
      }
   }
   delete g_V2_Ref_low; // signal, 1 point, done
   delete g_Jets_Ref_low; // signal, 1 point, done
   delete g_Nass_Ref_low; // signal, 1 point, done
   delete g_V2_Ref; // signal, 1 point, done
   delete g_Jets_Ref; // signal, 1 point, done
   delete g_Nass_Ref; // signal, 1 point, done
   
   return;
}

TF1 draw1D_longrange(TH1* hDeltaPhi, const char* name, 
      const char* cut1, const char* cut2, const char* cut3, const char* cut4,
      TCanvas* c_detaPhi, const int& ipad)
{
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";
   c_detaPhi->cd(ipad);
   gPad->SetBottomMargin(0.14);
   TF1 func("deltaPhi", function.c_str(), -3.14159*0.5, 3.14159*1.5);
   func.SetParameter(0, hDeltaPhi->GetMaximum());
   func.SetParameter(1, 0.1);
   func.SetParameter(2, 0.1);
   func.SetParameter(3, 0.1);

   hDeltaPhi->SetMarkerStyle(20);

   hDeltaPhi->Fit(&func, "q P");
   hDeltaPhi->Fit(&func, "q P");
   hDeltaPhi->Fit(&func, "m q P");
   hDeltaPhi->Fit(&func, "m q P E");
   auto fitResult = hDeltaPhi->Fit(&func, "m S E q P");

   hDeltaPhi->SetTitle(";#Delta#phi;");
   hDeltaPhi->GetXaxis()->CenterTitle();
   hDeltaPhi->GetXaxis()->SetTitleSize(0.05);

   hDeltaPhi->Draw();
   gPad->Update();
   TPaveStats* pave = (TPaveStats*) hDeltaPhi->FindObject("stats");
   pave->SetX1NDC(0.16);
   pave->SetX2NDC(0.56);
   gPad->Modified();
   gPad->Update();

   TLatex ltx;
   ltx.SetTextSize(0.035);
   ltx.SetTextFont(42);
   ltx.DrawLatexNDC(0.58, 0.42, cut1);
   ltx.DrawLatexNDC(0.56, 0.34, cut2);
   ltx.DrawLatexNDC(0.56, 0.26, cut3);
   ltx.DrawLatexNDC(0.58, 0.18, cut4);

   return func;
}

TH1D* proj1D_longrange(TH2* h2DSignal, TH2* h2DBackground, const char* name)
{
   //int negBinMin = 0;
   int negBinMin = 1;
   //int negBinMax = h2DSignal->GetXaxis()->FindBin(-1.* deltaEtaBound)-1 ;
   int negBinMax = h2DSignal->GetXaxis()->FindBin(-1.* deltaEtaBound) ;
   //int posBinMin = h2DSignal->GetXaxis()->FindBin(1.* deltaEtaBound)+1 ;
   int posBinMin = h2DSignal->GetXaxis()->FindBin(1.* deltaEtaBound) ;
   //int posBinMax = h2DSignal->GetXaxis()->GetNbins()+1;
   int posBinMax = h2DSignal->GetXaxis()->GetNbins();
   TH1D* hNeg = h2DSignal->ProjectionY("hneg", negBinMin, negBinMax);
   TH1D* hPos = h2DSignal->ProjectionY("hpos", posBinMin, posBinMax);
   hNeg->Add(hPos);

   TH1D* temp_neg = h2DBackground->ProjectionY("hneg_bkg", negBinMin, negBinMax);
   TH1D* temp_pos = h2DBackground->ProjectionY("hpos_bkg", posBinMin, posBinMax);
   temp_neg->Add(temp_pos);

   int center = h2DBackground->FindBin(0., 0.);
   hNeg->Divide(temp_neg);
   hNeg->Scale(h2DBackground->GetBinContent(center) / temp_neg->GetBinWidth(1)/ h2DBackground->GetXaxis()->GetBinWidth(1));

   delete hPos;
   delete temp_neg;
   delete temp_pos;

   hNeg->SetName(name);

   return hNeg;
}
std::pair<double, double> proj1D_shortrange_yields(TH2* h2DSignal, TH2* h2DBackground, const char* name, TCanvas* c, const int& ipad, TH1D* hsig)
{
   c->cd(ipad);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111111);
   gPad->SetLeftMargin(0.15);
   gPad->SetBottomMargin(0.15);
   int negBin = h2DSignal->GetXaxis()->FindBin(-1.);
   int posBin = h2DSignal->GetXaxis()->FindBin(1.);
   hsig = h2DSignal->ProjectionY("hsig", negBin, posBin);
   hsig->SetTitle(";#Delta#phi;dN/d(#Delta#phi)");

   TH1D* temp = h2DBackground->ProjectionY("hbkg", negBin, posBin);

   int center = h2DBackground->FindBin(0., 0.);
   hsig->Divide(temp);
   hsig->Scale(h2DBackground->GetBinContent(center) / (TMath::Pi()/16.)/ 0.3);

   TF1 fitter("fitter", "[0]*x^2+[1]*x+[2]", 0.6, 2.2);
   fitter.SetParameters(1, 1, 1);

   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->SetMarkerStyle(20);
   hsig->DrawCopy("E P", Form("pad%d", ipad));

   TLatex ltx;
   ltx.DrawLatexNDC(0.2, 0.8, "#font[42]{|#Delta#eta|<1}");
   ltx.DrawLatexNDC(0.58, 0.2, name);

   double min = fitter.GetMinimum(0.6, 2.2);
   double minX = fitter.GetMinimumX(0.6, 2.2);

   TF1 minfunc("minfunc", "[0]", -0.5*TMath::Pi(), 1.5*TMath::Pi());
   minfunc.SetParameter(0, -min);
   hsig->Add(&minfunc);

   double yields_err;
   double yields = hsig->IntegralAndError(hsig->FindBin(0.), hsig->FindBin(minX), yields_err, "width");
   yields = 2*yields - hsig->GetBinContent(hsig->FindBin(0.)) * TMath::Pi()/16.;

   delete hsig;
   delete temp;

   return std::pair<double, double>(yields, yields_err);
}

std::pair<double, double> proj1D_longrange_yields(TH2* h2DSignal, TH2* h2DBackground, const char* name, TCanvas* c, const int& ipad , TH1D* hsig)
{
   c->cd(ipad);
   //double lw = 0.0;
   double lw = 0.1;
   double up = 2.0;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111111);
   c->SetLeftMargin(0.15);
   c->SetBottomMargin(0.15);

   int negBinMin = 1;
   int negBinMax = 10; // delta eta = -2
   //int negBinMax = 14;
   int posBinMin = 24;// delta eta = 2
   //int posBinMin = 20; 
   int posBinMax = 33;

   TH1D* hNeg = h2DSignal->ProjectionY("hneg", negBinMin, negBinMax);
   TH1D* hPos = h2DSignal->ProjectionY("hpos", posBinMin, posBinMax);
   hNeg->Add(hPos);
   hsig = new TH1D(*hNeg);
   hsig->SetName("hsig");
   delete hNeg;
   delete hPos;

   // symmetrize
   for(int i=1; i<=8; i++){
      hsig->SetBinContent(i+8, hsig->GetBinContent(9-i)+hsig->GetBinContent(i+8));
      hsig->SetBinError(i+8, sqrt(pow(hsig->GetBinError(9-i), 2)+pow(hsig->GetBinError(i+8), 2)));
   }
   for(int i=1; i<=8; i++){
      hsig->SetBinContent(i+16, hsig->GetBinContent(33-i)+hsig->GetBinContent(i+16));
      hsig->SetBinError(i+16, sqrt(pow(hsig->GetBinError(33-i), 2)+pow(hsig->GetBinError(i+16), 2)));
   }
   hsig->GetXaxis()->SetRange(9, 24);

   TH1D* temp_neg = h2DBackground->ProjectionY("hneg_bkg", negBinMin, negBinMax);
   TH1D* temp_pos = h2DBackground->ProjectionY("hpos_bkg", posBinMin, posBinMax);
   temp_neg->Add(temp_pos);

   TH1D* temp = new TH1D(*temp_neg);
   temp->SetName("temp");
   delete temp_neg;
   delete temp_pos;

   int center = h2DBackground->FindBin(0., 0.);
   hsig->Divide(temp);
   hsig->Scale(h2DBackground->GetBinContent(center) / temp->GetBinWidth(1)/h2DBackground->GetXaxis()->GetBinWidth(1));

   TF1 fitter("fitter", "[0]*x^2+[1]*x+[2]", lw, up);
   fitter.SetParameters(1, 1, 1);
   fitter.SetParLimits(0, 0, 10);

   hsig->Fit(&fitter, "R q l");
   hsig->Fit(&fitter, "R q l");
   hsig->Fit(&fitter, "R q l");
   hsig->Fit(&fitter, "R q l");
   hsig->Fit(&fitter, "R q");
   fitter.ReleaseParameter(0);
   hsig->Fit(&fitter, "R q l");
   hsig->Fit(&fitter, "R q l");

   hsig->SetTitle(";#Delta#phi;dN/d(#Delta#phi)");
   hsig->SetMarkerStyle(20);
   hsig->DrawCopy("E P", Form("pad%d", ipad));

   TLatex ltx;
   ltx.DrawLatexNDC(0.2, 0.8, "#font[42]{|#Delta#eta|>1}");
   ltx.DrawLatexNDC(0.58, 0.2, name);

   double min = fitter.GetMinimum(lw, up);
   double minX = fitter.GetMinimumX(lw, up);

   TF1 minfunc("minfunc", "[0]", -0.5*TMath::Pi(), 1.5*TMath::Pi());
   minfunc.SetParameter(0, -min);
   hsig->Add(&minfunc);

   double yields_err;
   double yields = hsig->IntegralAndError(hsig->FindBin(0.), hsig->FindBin(minX), yields_err, "width");
   yields = 2*yields - hsig->GetBinContent(hsig->FindBin(0.)) * TMath::Pi()/16.;

   delete hsig;
   delete temp;

   return std::pair<double, double>(yields, yields_err);
}

std::pair<double, double> calYields(const pair<double, double>& sr, const pair<double,double>& lr)
{
   double yields = sr.first - lr.first;
   double yields_err = sqrt(pow(sr.second, 2) + pow(lr.second, 2));
   return pair<double, double>(yields, yields_err);
}

string ouput_prefix(const char* name)
{
   std::string str = name;
   auto found = str.find("/");
   while(found!=std::string::npos){
      str.replace(found, 1, "_");
      found = str.find("/");
   }
   while(str.at(0) == '.'){
      str.erase(0, 1);
   }
   while(str.at(0) == '_'){
      str.erase(0, 1);
   }
   return str;
}

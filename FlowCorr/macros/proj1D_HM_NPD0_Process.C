const double deltaEtaBound = 1;
const int nDca= 2;

TH1D* proj1D_longrange(TH2*, TH2*, const char*);
TH1D* proj1D_shortrange(TH2*, TH2*);

TF1 draw1D_longrange(TH1*, const char*,
const char*, const char* , const char* );

void proj1D_HM_NPD0_Process(const char* input_d0= "",
      const char* input_ref = "", 
      const string dataset="", const char* input_d0_low_mult = "",
      const float pTMin=0.0, const float pTMax=0.0,
      const float yMin =0.0, const float yMax =0.0
      )
{
   //gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   if(!ana::isHM_PD0_DataSet(dataset)){
      std::cerr << "wrong dataset" << std::endl;
      return;
   }

   vector<double> ptbin;
   int nPtBin;
   if(dataset == "PAHM1-6") {
      ptbin = vector<double>(ana::ptbin_NPD0_pPb, ana::ptbin_NPD0_pPb+ana::nPt_NPD0_pPb+1);
      nPtBin= ptbin.size()-1;
   }

   const int nPt = nPtBin;

   double v2_NPD0[ana::nMass][nPt][nDca];
   double v2_NPD0_err[ana::nMass][nPt][nDca];

   double V2_NPD0[ana::nMass][nPt][nDca];
   double V2_NPD0_err[ana::nMass][nPt][nDca];
   double V2_REF;
   double V2_REF_err;

   double N_ass[ana::nMass][nPt][nDca];
   double N_ass_err[ana::nMass][nPt][nDca];
   double N_ass_low[ana::nMass][nDca];
   double N_ass_low_err[ana::nMass][nDca];

   double yields_jet[ana::nMass][nPt][nDca];
   double yields_jet_err[ana::nMass][nPt][nDca];
   double yields_jet_low[ana::nMass][nDca];

   double V2_NPD0_low[ana::nMass][nDca];
   double V2_NPD0_low_err[ana::nMass][nDca];

   double V2_Sub_NPD0[ana::nMass][nPt][nDca];
   double V2_Sub_NPD0_err[ana::nMass][nPt][nDca];

   double v2_Sub_NPD0[ana::nMass][nPt][nDca];
   double v2_Sub_NPD0_err[ana::nMass][nPt][nDca];

   TFile* f1 = new TFile(input_d0);
   TFile* f2 = new TFile(input_ref);

   TH3D* hDcaVsMassAndMva[nPt]; 
   
   TH2D* h2DSignal_D0[ana::nMass][nPt][nDca];
   TH2D* h2DBackground_D0[ana::nMass][nPt][nDca];
   TH1D* hMult_raw_D0[ana::nMass][nPt][nDca]; // wrong normalized constant
   TH1D* hMass_D0[ana::nMass][nPt][nDca];

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

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            double nMult_D0= hMult_raw_D0[imass][ipt][idca]->GetEntries()
               - hMult_raw_D0[imass][ipt][idca]->GetBinContent(1);
            h2DSignal_D0[imass][ipt][idca]->Scale(1./nMult_D0);
         }
      }
   }

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

   TH2D* h2DSignal_Ref = (TH2D*) f2->Get("hSignal_Ref");
   TH2D* h2DBackground_Ref = (TH2D*) f2->Get("hBackground_Ref");
   TH1D* hMult = (TH1D*) f2->Get("hMult");

   long int nevents = hMult->GetEntries();
   TH1D* hDeltaPhi_Ref = proj1D_longrange(h2DSignal_Ref, h2DBackground_Ref, "deltaPhi_Ref");
   hDeltaPhi_Ref->Scale(1./nevents);
   
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            string str = input_d0;
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
            std::string ntrk;
            if(dataset == "PAHM1-6") ntrk = "185< N_{trk}^{offline} < 250";
            auto func = draw1D_longrange(hDeltaPhi[imass][ipt][idca],
                  Form("../plots/v2vspt/npd0ana1/y%.1f/%s/%s_deltaPhi_mass%d_pt%d_dca%d.png",
                     yMax, dataset.c_str(), str.c_str(), imass, ipt, idca),
                  ntrk.c_str(),
                  Form("%.1f<p_{T}<%.1fGeV, %.1f<y<%.1f", ptbin[ipt], ptbin[ipt+1], yMin, yMax),
                  Form("%.3f<mass<%.3fGeV", ana::massbin[imass], ana::massbin[imass+1])
               );
            V2_NPD0[imass][ipt][idca] = func.GetParameter(2);
            V2_NPD0_err[imass][ipt][idca] = func.GetParError(2);
            N_ass[imass][ipt][idca] = func.GetParameter(0);
            N_ass_err[imass][ipt][idca] = func.GetParError(0);
         }
      }
   }

   if(true){
      std::string str = input_ref;
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
      std::string ntrk;
      if(dataset == "PAHM1-6") ntrk = "185< N_{trk}^{offline} < 250";
      auto func_ref = draw1D_longrange(hDeltaPhi_Ref,
            Form("../plots/v2vspt/npd0ana1/%s_ref_deltaPhi.png", 
            str.c_str()),
            ntrk.c_str(),
            "0.3<p_{T}<3.0GeV |#eta|<2.4", 
            ""
         );
      V2_REF = func_ref.GetParameter(2);
      V2_REF_err = func_ref.GetParError(2);
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

   TFile f3(Form("%s_v2.root", input_d0), "recreate");

   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         g_v2[ipt][idca]->Write();
         hMass_DCA[ipt][idca]->Write();
         hPt[ipt][idca]->Write();
         hKET[ipt][idca]->Write();
      }
   }
   
   return;
}

TF1 draw1D_longrange(TH1* hDeltaPhi, const char* name, 
      const char* cut1, const char* cut2, const char* cut3)
{
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";
   TCanvas* c_deltaPhi = new TCanvas("c_deltaPhi", "", 550, 450);
   c_deltaPhi->SetBottomMargin(0.14);
   TF1 func("deltaPhi", function.c_str(), -3.14159*0.5, 3.14159*1.5);
   func.SetParameter(0, hDeltaPhi->GetMaximum());
   func.SetParameter(1, 0.1);
   func.SetParameter(2, 0.1);
   func.SetParameter(3, 0.1);

   hDeltaPhi->SetMarkerStyle(20);

   hDeltaPhi->Fit(&func, "q");
   hDeltaPhi->Fit(&func, "q");
   hDeltaPhi->Fit(&func, "m q");
   hDeltaPhi->Fit(&func, "m q E");
   auto fitResult = hDeltaPhi->Fit(&func, "m S E q");

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
   ltx.DrawLatexNDC(0.58, 0.34, cut2);
   ltx.DrawLatexNDC(0.58, 0.26, cut3);

   c_deltaPhi->Print(name);

   delete c_deltaPhi;
   return func;
}

TH1D* proj1D_shortrange(TH2* h2DSignal, TH2* h2DBackground)
{
   return nullptr;
}
TH1D* proj1D_longrange(TH2* h2DSignal, TH2* h2DBackground, const char* name)
{
   int negBinMin = 0;
   int negBinMax = h2DSignal->GetXaxis()->FindBin(-1.* deltaEtaBound)-1 ;
   int posBinMin = h2DSignal->GetXaxis()->FindBin(1.* deltaEtaBound)+1 ;
   int posBinMax = h2DSignal->GetXaxis()->GetNbins()+1;
   TH1D* hNeg = h2DSignal->ProjectionY("hneg", negBinMin, negBinMax);
   TH1D* hPos = h2DSignal->ProjectionY("hpos", posBinMin, posBinMax);
   hNeg->Add(hPos);

   TH1D* temp_neg = h2DBackground->ProjectionY("hneg_bkg", negBinMin, negBinMax);
   TH1D* temp_pos = h2DBackground->ProjectionY("hpos_bkg", posBinMin, posBinMax);
   temp_neg->Add(temp_pos);

   //int center = h2DBackground->FindBin(0., 0.);
   int center = temp_neg->FindBin(0.);
   hNeg->Divide(temp_neg);
   //hNeg->Scale(h2DBackground->GetBinContent(center) / temp_neg->GetBinWidth(1));
   hNeg->Scale(temp_neg->GetBinContent(center) / temp_neg->GetBinWidth(1));

   delete hPos;
   delete temp_neg;
   delete temp_pos;

   hNeg->SetName(name);

   return hNeg;
}

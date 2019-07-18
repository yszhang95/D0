const double deltaEtaBound = 1;

TH1D* proj1D_longrange(TH2*, TH2*, const char*);
TH1D* proj1D_shortrange(TH2*, TH2*);

TF1 draw1D_longrange(TH1*, const char*, 
      const char*, const char* , const char* );

void proj1D_HM_PD0_Process(const char* input_d0= "",
      const char* input_ref = "",
      const string dataset="", const char* input_d0_low_mult="",
      const float pTMin=0., const float pTMax=0., 
      const float yMin =0., const float yMax =0.)
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
      ptbin = vector<double>(ana::ptbin_PD0_pPb, ana::ptbin_PD0_pPb+ana::nPt_PD0_pPb+1);
      nPtBin = ptbin.size()-1;
   }
   if(dataset == "PPHM") {
      ptbin = vector<double>(ana::ptbin_PD0_pp, ana::ptbin_PD0_pp+ana::nPt_PD0_pp+1);
      nPtBin = ptbin.size()-1;
   }

   const int nPt = nPtBin;

   double v2_PD0[ana::nMass][nPt];
   double v2_PD0_err[ana::nMass][nPt];

   double V2_PD0[ana::nMass][nPt];
   double V2_PD0_err[ana::nMass][nPt];
   double V2_REF;
   double V2_REF_err;

   double N_ass[ana::nMass][nPt];
   double N_ass_err[ana::nMass][nPt];
   double N_ass_low[ana::nMass];
   double N_ass_low_err[ana::nMass];

   double yields_jet[ana::nMass][nPt];
   double yields_jet_err[ana::nMass][nPt];
   double yields_jet_low[ana::nMass];

   double V2_PD0_low[ana::nMass];
   double V2_PD0_low_err[ana::nMass];

   double V2_Sub_PD0[ana::nMass][nPt];
   double V2_Sub_PD0_err[ana::nMass][nPt];

   double v2_Sub_PD0[ana::nMass][nPt];
   double v2_Sub_PD0_err[ana::nMass][nPt];

   TFile* f1 = new TFile(input_d0);
   TFile* f2 = new TFile(input_ref);

   TH3D* hDcaVsMassAndMva[nPt]; 
   
   TH2D* h2DSignal_D0[ana::nMass][nPt];
   TH2D* h2DBackground_D0[ana::nMass][nPt];
   TH1D* hMult_raw_D0[ana::nMass][nPt]; // wrong normalized constant
   TH1D* hMass_D0[ana::nMass][nPt];

   // read HM 2D
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         h2DSignal_D0[imass][ipt] = (TH2D*) f1->Get(Form("hSignal_mass%d_pt%d", imass, ipt));
         h2DBackground_D0[imass][ipt] = (TH2D*) f1->Get(Form("hBackground_mass%d_pt%d", imass, ipt));
         hMult_raw_D0[imass][ipt] = (TH1D*) f1->Get(Form("hMult_raw_D0_mass%d_pt%d", imass, ipt));
         hMass_D0[imass][ipt] = (TH1D*) f1->Get(Form("hMassD0_mass%d_pt%d", imass, ipt));
      }
   }

   for(int ipt=0; ipt<nPt; ipt++){
      hDcaVsMassAndMva[ipt] = (TH3D*) f1->Get(Form("hDcaVsMassAndMva_pt%d", ipt));
   }

   TH1D* hPt[nPt];
   TH1D* hKET[nPt];
   for(int ipt=0; ipt<nPt; ipt++){
      f1->GetObject(Form("hPt_pt%d", ipt), hPt[ipt]);
      f1->GetObject(Form("hKET_pt%d", ipt), hKET[ipt]);
      if(!hPt[ipt] || !hKET[ipt]){
         cout << "cannot find hPt or hKET" << endl;
         return;
      }
   }

   // scaled by event number
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         double nMult_D0= hMult_raw_D0[imass][ipt]->GetEntries()-hMult_raw_D0[imass][ipt]->GetBinContent(1);
         h2DSignal_D0[imass][ipt]->Scale(1./nMult_D0);
      }
   }

   TH1D* hNeg[ana::nMass][nPt];
   TH1D* hPos[ana::nMass][nPt];
   TH1D* hDeltaPhi[ana::nMass][nPt];

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         hDeltaPhi[imass][ipt] = proj1D_longrange(h2DSignal_D0[imass][ipt], h2DBackground_D0[imass][ipt], Form("deltaPhi_mass%d_pt%d", imass, ipt));
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
         std::string str = input_d0;
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
         if(dataset == "PPHM") ntrk = "100< N_{trk}^{offline} < 250";
         auto func = draw1D_longrange(hDeltaPhi[imass][ipt],
               Form("../plots/v2vspt/d0ana/y%.1f/%s/%s_deltaPhi_mass%d_pt%d.png", 
                  yMax, dataset.c_str(), str.c_str(), imass, ipt),
               ntrk.c_str(),
               Form("%.1f<p_{T}<%.1fGeV, %.1f<y<%.1f", ptbin[ipt], ptbin[ipt+1], yMin, yMax), 
               Form("%.3f<mass<%.3fGeV", ana::massbin[imass], ana::massbin[imass+1])
               );
         V2_PD0[imass][ipt] = func.GetParameter(2);
         V2_PD0_err[imass][ipt] = func.GetParError(2);
         N_ass[imass][ipt] = func.GetParameter(0);
         N_ass_err[imass][ipt] = func.GetParError(0);
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
      if(dataset == "PPHM") ntrk = "100< N_{trk}^{offline} < 250";
      auto func_ref = draw1D_longrange(hDeltaPhi_Ref,
            Form("../plots/v2vspt/d0ana/%s_ref_deltaPhi.png", 
            str.c_str()),
            ntrk.c_str(),
            "0.3<p_{T}<3.0GeV |#eta|<2.4", 
            ""
            );
      V2_REF = func_ref.GetParameter(2);
      V2_REF_err = func_ref.GetParError(2);
   }

   TGraphErrors* g_v2[nPt];
   for(int ipt=0; ipt<nPt; ipt++){
      g_v2[ipt] = new TGraphErrors(ana::nMass);
      g_v2[ipt]->SetName(Form("g_v2_pt%d", ipt));
      for(int imass=0; imass<ana::nMass; imass++){
         double temp = V2_PD0[imass][ipt]/ sqrt(V2_REF);
         double temp_err = temp* sqrt(pow(V2_PD0_err[imass][ipt]/ V2_PD0[imass][ipt], 2) 
               + pow(0.5*V2_REF_err/ V2_REF, 2));

         v2_PD0[imass][ipt] = temp;
         v2_PD0_err[imass][ipt] = temp_err;

         g_v2[ipt]->SetPoint(imass, hMass_D0[imass][ipt]->GetMean(), v2_PD0[imass][ipt]);
         g_v2[ipt]->SetPointError(imass, 0, v2_PD0_err[imass][ipt]);
      }
   }

   TH1D* hMass[nPt];
   for(int ipt=0; ipt<nPt; ipt++){
      hMass[ipt] = hDcaVsMassAndMva[ipt]->ProjectionX(Form("hmass_pt%d", ipt));
   }

   TFile f3(Form("%s_v2.root", input_d0), "recreate");

   for(int ipt=0; ipt<nPt; ipt++){
      g_v2[ipt]->Write();
      hMass[ipt]->Write();
      hPt[ipt]->Write();
      hKET[ipt]->Write();
   }
   return;
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

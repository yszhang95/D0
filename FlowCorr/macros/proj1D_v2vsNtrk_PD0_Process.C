void proj1D_v2vsNtrk_PD0_Process(const char* input_d0= "",
      const char* input_ref = "", 
const float y=0.0, const string dataset="", const char* input_d0_low_mult="")
{
   const double deltaEtaBound = 1;

   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   const int n_trk_bin_ = ana::Get_N_nTrkBin(dataset);
   if(n_trk_bin_<0){
      std::cerr << "wrong dataset naem" << std::endl;
   }

   double v2_PD0[ana::nMass][n_trk_bin_];
   double v2_PD0_err[ana::nMass][n_trk_bin_];

   double V2_PD0[ana::nMass][n_trk_bin_];
   double V2_PD0_err[ana::nMass][n_trk_bin_];
   double V2_REF[n_trk_bin_];
   double V2_REF_err[n_trk_bin_];

   double N_ass[n_trk_bin_];
   double N_ass_low[n_trk_bin_];

   double Upsilon[ana::nMass][n_trk_bin_];
   double Upsilon_low[ana::nMass][n_trk_bin_];

   double V2_Sub_PD0[ana::nMass][n_trk_bin_];
   double V2_Sub_PD0_err[ana::nMass][n_trk_bin_];

   TFile* f1 = new TFile(input_d0);
   TFile* f2 = new TFile(input_ref);

   TH3D* hDcaVsMassAndMva[n_trk_bin_]; 
   
   TH2D* h2DSignal_D0[ana::nMass][n_trk_bin_];
   TH2D* h2DBackground_D0[ana::nMass][n_trk_bin_];
   TH1D* hMult_raw_D0[ana::nMass][n_trk_bin_]; // wrong normalized constant
   TH1D* hMass_D0[ana::nMass][n_trk_bin_];

   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         h2DSignal_D0[imass][i_trk_bin_] = (TH2D*) f1->Get(Form("hSignal_mass%d_trk%d", imass, i_trk_bin_));
         h2DBackground_D0[imass][i_trk_bin_] = (TH2D*) f1->Get(Form("hBackground_mass%d_trk%d", imass, i_trk_bin_));
         hMult_raw_D0[imass][i_trk_bin_] = (TH1D*) f1->Get(Form("hMult_raw_D0_mass%d_trk%d", imass, i_trk_bin_));
         hMass_D0[imass][i_trk_bin_] = (TH1D*) f1->Get(Form("hMassD0_mass%d_trk%d", imass, i_trk_bin_));
      }
   }

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      hDcaVsMassAndMva[i_trk_bin_] = (TH3D*) f1->Get(Form("hDcaVsMassAndMva_trk%d", i_trk_bin_));
   }

   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         double nMult_D0= hMult_raw_D0[imass][i_trk_bin_]->Integral(2, 251);
         h2DSignal_D0[imass][i_trk_bin_]->Scale(1./nMult_D0);
      }
   }

   TH1D* hNeg[ana::nMass][n_trk_bin_];
   TH1D* hPos[ana::nMass][n_trk_bin_];
   TH1D* hDeltaPhi[ana::nMass][n_trk_bin_];

   int negBinMin = 0;
   int negBinMax = h2DSignal_D0[0][0]->GetXaxis()->FindBin(-1.* deltaEtaBound)-1 ;
   int posBinMin = h2DSignal_D0[0][0]->GetXaxis()->FindBin(1.* deltaEtaBound)+1 ;
   int posBinMax = h2DSignal_D0[0][0]->GetXaxis()->GetNbins()+1;

   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         hNeg[imass][i_trk_bin_] = h2DSignal_D0[imass][i_trk_bin_]->ProjectionY(
                  Form("neg_mass%d_trk%d", imass, i_trk_bin_), negBinMin, negBinMax);
         hPos[imass][i_trk_bin_] = h2DSignal_D0[imass][i_trk_bin_]->ProjectionY(
                  Form("pos_mass%d_trk%d", imass, i_trk_bin_), posBinMin, posBinMax);
         hDeltaPhi[imass][i_trk_bin_] = (TH1D*) hNeg[imass][i_trk_bin_]->Clone();
         hDeltaPhi[imass][i_trk_bin_]->SetName(Form("deltaPhi_mass%d_trk%d", imass, i_trk_bin_));
         hDeltaPhi[imass][i_trk_bin_]->Add(hPos[imass][i_trk_bin_]);

         TH1D* temp_neg = h2DBackground_D0[imass][i_trk_bin_]->ProjectionY(
               Form("neg_mass%d_trk%d_temp", imass, i_trk_bin_), negBinMin, negBinMax);
         TH1D* temp_pos = h2DBackground_D0[imass][i_trk_bin_]->ProjectionY(
               Form("pos_mass%d_trk%d__temp", imass, i_trk_bin_), posBinMin, posBinMax);
         TH1D* temp = (TH1D*) temp_neg->Clone();
         temp->SetName(Form("mass%d_trk%d_temp", imass, i_trk_bin_));
         temp->Add(temp_pos);

         int center = h2DBackground_D0[imass][i_trk_bin_]->FindBin(0., 0.);
         hDeltaPhi[imass][i_trk_bin_]->Divide(temp);
         hDeltaPhi[imass][i_trk_bin_]->Scale(h2DBackground_D0[i_trk_bin_][imass]->GetBinContent(center) / temp->GetBinWidth(1));

         delete hNeg[imass][i_trk_bin_];
         delete hPos[imass][i_trk_bin_];
         delete temp_neg;
         delete temp_pos;
         delete temp;
      }
   }

   TH2D* h2DSignal_Ref[n_trk_bin_];
   TH2D* h2DBackground_Ref[n_trk_bin_];
   TH1D* hDeltaPhi_Ref[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      h2DSignal_Ref[i_trk_bin_] = (TH2D*) f2->Get(Form("hSignal_Ref_trk%d", i_trk_bin_));
      h2DBackground_Ref[i_trk_bin_] = (TH2D*) f2->Get(Form("hBackground_Ref_trk%d", i_trk_bin_));
      TH1D* hMult = (TH1D*) f2->Get(Form("hMult_trk%d", i_trk_bin_));

      int negBinMin_Ref = 0;
      int negBinMax_Ref = h2DSignal_Ref[i_trk_bin_]->GetXaxis()->FindBin(-1.*deltaEtaBound)-1 ;
      int posBinMin_Ref = h2DSignal_Ref[i_trk_bin_]->GetXaxis()->FindBin(1.*deltaEtaBound)+1 ;
      int posBinMax_Ref = h2DSignal_Ref[i_trk_bin_]->GetXaxis()->GetNbins()+1;

      TH1D* hNeg_Ref_Signal;
      TH1D* hPos_Ref_Signal;
      hNeg_Ref_Signal = h2DSignal_Ref[i_trk_bin_]->ProjectionY("neg_ref_signal", negBinMin, negBinMax);
      hPos_Ref_Signal = h2DSignal_Ref[i_trk_bin_]->ProjectionY("pos_ref_signal", posBinMin, posBinMax);

      hDeltaPhi_Ref[i_trk_bin_] = (TH1D*) hNeg_Ref_Signal->Clone();
      hDeltaPhi_Ref[i_trk_bin_]->SetName(Form("deltaPhi_Ref_trk%d", i_trk_bin_));
      hDeltaPhi_Ref[i_trk_bin_]->Add(hPos_Ref_Signal);

      TH1D* hNeg_Ref_Background;
      TH1D* hPos_Ref_Background;
      hNeg_Ref_Background = h2DBackground_Ref[i_trk_bin_]->ProjectionY("neg_ref_background", negBinMin, negBinMax);
      hPos_Ref_Background = h2DBackground_Ref[i_trk_bin_]->ProjectionY("pos_ref_background", posBinMin, posBinMax);

      TH1D* hDeltaPhi_Ref_bkg = (TH1D*) hNeg_Ref_Background->Clone();
      hDeltaPhi_Ref_bkg->SetName("deltaPhi_Ref_bkg");
      hDeltaPhi_Ref_bkg->Add(hPos_Ref_Background);

      int center = h2DBackground_Ref[i_trk_bin_]->FindBin(0., 0.);
      hDeltaPhi_Ref[i_trk_bin_]->Divide(hDeltaPhi_Ref_bkg);
      hDeltaPhi_Ref[i_trk_bin_]->Scale(h2DBackground_Ref[i_trk_bin_]->GetBinContent(center)/hDeltaPhi_Ref_bkg->GetBinWidth(1));

      delete hNeg_Ref_Signal;
      delete hPos_Ref_Signal;
      delete hNeg_Ref_Background;
      delete hPos_Ref_Background;
      delete hDeltaPhi_Ref_bkg;
   }
   
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";

   TCanvas* c_deltaPhi[ana::nMass][n_trk_bin_];
   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         c_deltaPhi[imass][i_trk_bin_] = new TCanvas("c_deltaPhi", "", 550, 450);
         TF1 func("deltaPhi", function.c_str(), -3.14159*0.5, 3.14159*1.5);
         func.SetParameter(0, hDeltaPhi[imass][i_trk_bin_]->GetMaximum());
         func.SetParameter(1, 0.1);
         func.SetParameter(2, 0.1);
         func.SetParameter(3, 0.1);

         hDeltaPhi[imass][i_trk_bin_]->SetMarkerStyle(20);

         hDeltaPhi[imass][i_trk_bin_]->Fit(&func, "q");
         hDeltaPhi[imass][i_trk_bin_]->Fit(&func, "q");
         hDeltaPhi[imass][i_trk_bin_]->Fit(&func, "m q");
         hDeltaPhi[imass][i_trk_bin_]->Fit(&func, "m q E");
         auto fitResult = hDeltaPhi[imass][i_trk_bin_]->Fit(&func, "m S E q");

         V2_PD0[imass][i_trk_bin_] = func.GetParameter(2);
         V2_PD0_err[imass][i_trk_bin_] = func.GetParError(2);

         hDeltaPhi[imass][i_trk_bin_]->SetTitle(";#Delta#phi;");

         TLatex ltx;
         ltx.DrawLatexNDC(0.16, 0.85, "1.5<pT<8GeV, |y|<2");
         ltx.DrawLatexNDC(0.16, 0.75, Form("%.3f<mass<%.3fGeV", ana::massbin[imass], ana::massbin[imass+1]));
         auto trkRange = ana::get_Mult_Edges(dataset);
         if(trkRange.size() && trkRange[trkRange.size()]!=std::numeric_limits<unsigned int>::max()) 
            ltx.DrawLatexNDC(0.16, 0.68, Form("%u<Ntrkoffline<%u", trkRange[i_trk_bin_], trkRange[i_trk_bin_+1]));
         if(trkRange.size() && trkRange[trkRange.size()]==std::numeric_limits<unsigned int>::max()) 
            ltx.DrawLatexNDC(0.16, 0.68, Form("Ntrkoffline>%u", trkRange[i_trk_bin_]));

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

         hDeltaPhi[imass][i_trk_bin_]->Draw();
         gPad->Update();
         TPaveStats* pave = (TPaveStats*) hDeltaPhi[imass][i_trk_bin_]->FindObject("stats");
         pave->SetX1NDC(0.16);
         pave->SetX2NDC(0.56);
         gPad->Modified();
         gPad->Update();
         c_deltaPhi[imass][i_trk_bin_]->Print(Form("../plots/v2vsNtrk/y%.1f/%s/%s_deltaPhi_mass%d_trk%d.png", 
                  y, dataset.c_str(), str.c_str(), imass, i_trk_bin_));

         delete c_deltaPhi[imass][i_trk_bin_];
      }
   }

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      TCanvas* cRef = new TCanvas("cRef", "", 550, 450);
      TF1 func_ref("deltaPhi_Ref", function.c_str(), -3.14159*0.5, 3.14159*1.5);
      func_ref.SetParameter(0, hDeltaPhi_Ref[i_trk_bin_]->GetMaximum());
      func_ref.SetParameter(1, 0.1); 
      func_ref.SetParameter(2, 0.1);
      func_ref.SetParameter(3, 0.1);

      V2_REF[i_trk_bin_] = func_ref.GetParameter(2);
      V2_REF_err[i_trk_bin_] = func_ref.GetParError(2);

      hDeltaPhi_Ref[i_trk_bin_]->SetMarkerStyle(20);
      hDeltaPhi_Ref[i_trk_bin_]->Fit(&func_ref, "q");
      hDeltaPhi_Ref[i_trk_bin_]->Fit(&func_ref, "q");
      hDeltaPhi_Ref[i_trk_bin_]->Fit(&func_ref, "m q");
      auto fitResult = hDeltaPhi_Ref[i_trk_bin_]->Fit(&func_ref, "m q S E");
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
      hDeltaPhi_Ref[i_trk_bin_]->Draw();
      gPad->Update();
      TPaveStats* pave = (TPaveStats*) hDeltaPhi_Ref[i_trk_bin_]->FindObject("stats");
      pave->SetX1NDC(0.16);
      pave->SetX2NDC(0.56);
      gPad->Modified();
      gPad->Update();
      cRef->Print(Form("../plots/v2vsNtrk/%s_ref_trk%d_V2.png", str.c_str(), i_trk_bin_));
      delete cRef;
   }

   TGraphErrors* g_v2_[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_v2_[i_trk_bin_] = new TGraphErrors(ana::nMass);
      g_v2_[i_trk_bin_]->SetName(Form("g_v2_trk%d", i_trk_bin_));
      for(int imass=0; imass<ana::nMass; imass++){
         double temp = V2_PD0[imass][i_trk_bin_]/ sqrt(V2_REF[i_trk_bin_]);
         double temp_err = temp* sqrt(pow(V2_PD0_err[imass][i_trk_bin_]/ V2_PD0[imass][i_trk_bin_], 2) 
               + pow(0.5*V2_REF[i_trk_bin_]/ V2_REF[i_trk_bin_], 2));

         v2_PD0[imass][i_trk_bin_] = temp;
         v2_PD0_err[imass][i_trk_bin_] = temp_err;

         g_v2_[i_trk_bin_]->SetPoint(imass, hMass_D0[imass][i_trk_bin_]->GetMean(), v2_PD0[imass][i_trk_bin_]);
         g_v2_[i_trk_bin_]->SetPointError(imass, 0, v2_PD0_err[imass][i_trk_bin_]);
      }
   }

   TH1D* hMass[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      hMass[i_trk_bin_] = hDcaVsMassAndMva[i_trk_bin_]->ProjectionX(Form("hmass_trk%d", i_trk_bin_));
   }

   TH1D* hNtrk[n_trk_bin_];
   TH1D* hKET[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      hNtrk[i_trk_bin_] = (TH1D*) f1->Get(Form("hNtrk_trk%d", i_trk_bin_));
      hKET[i_trk_bin_] = (TH1D*) f1->Get(Form("hKET_trk%d", i_trk_bin_));
   }

   TFile f3(Form("%s_v2.root", input_d0), "recreate");

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_v2_[i_trk_bin_]->Write();
      hMass[i_trk_bin_]->Write();
      hNtrk[i_trk_bin_]->Write();
      hKET[i_trk_bin_]->Write();
   }
   f3.Close();

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      delete g_v2_[i_trk_bin_];
   }
   delete f1;
   delete f2;
   
   return;
}


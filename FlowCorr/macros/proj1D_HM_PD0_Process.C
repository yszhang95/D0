void proj1D_HM_PD0_Process(const char* input_d0= "../data/corr2D_PAHM185-250_d0ana_HM_2.0.root",
      const char* input_ref = "data/corr2D_ref_d0ana.root", 
const float y=2.0, const string dataset="PAHM1-6")
{

   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   if(!ana::isHM_PD0_DataSet(dataset)){
      std::cerr << "wrong dataset" << std::endl;
      return;
   }

   vector<double> ptbin;
   int nPt;
   if(dataset == "PAHM1-6") {
      ptbin = vector<double>(ana::ptbin_PD0_pPb, ana::ptbin_PD0_pPb+ana::nPt_PD0_pPb+1);
      nPt = ptbin.size()-1;
   }
   if(dataset == "PPHM_2") {
      ptbin = vector<double>(ana::ptbin_PD0_pp, ana::ptbin_PD0_pp+ana::nPt_PD0_pp+1);
      nPt = ptbin.size()-1;
   }


   double v2_DCA[ana::nMass][nPt];
   double v2_DCA_err[ana::nMass][nPt];

   TFile* f1 = new TFile(input_d0);
//   TFile* f2 = new TFile(input_ref);

   TH3D* hDcaVsMassAndMva[nPt]; 
   
   TH2D* h2DSignal_D0[ana::nMass][nPt];
   TH2D* h2DBackground_D0[ana::nMass][nPt];
   TH1D* hMult_raw_D0[ana::nMass][nPt]; // wrong normalized constant
   TH1D* hMass_D0[ana::nMass][nPt];

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

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         double nMult_D0= hMult_raw_D0[imass][ipt]->Integral(2, 251);
         h2DSignal_D0[imass][ipt]->Scale(1./nMult_D0);
      }
   }

   TH1D* hNeg[ana::nMass][nPt];
   TH1D* hPos[ana::nMass][nPt];
   TH1D* hDeltaPhi[ana::nMass][nPt];

   int negBinMin = 0;
   int negBinMax = h2DSignal_D0[0][0]->GetXaxis()->FindBin(-1.)-1 ;
   int posBinMin = h2DSignal_D0[0][0]->GetXaxis()->FindBin(1.)+1 ;
   int posBinMax = h2DSignal_D0[0][0]->GetXaxis()->GetNbins()+1;

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         hNeg[imass][ipt] = h2DSignal_D0[imass][ipt]->ProjectionY(
                  Form("neg_mass%d_pt%d", imass, ipt), negBinMin, negBinMax);
         hPos[imass][ipt] = h2DSignal_D0[imass][ipt]->ProjectionY(
                  Form("pos_mass%d_pt%d", imass, ipt), posBinMin, posBinMax);
         hDeltaPhi[imass][ipt] = (TH1D*) hNeg[imass][ipt]->Clone();
         hDeltaPhi[imass][ipt]->SetName(Form("deltaPhi_mass%d_pt%d", imass, ipt));
         hDeltaPhi[imass][ipt]->Add(hPos[imass][ipt]);

         TH1D* temp_neg = h2DBackground_D0[imass][ipt]->ProjectionY(
               Form("neg_mass%d_pt%d_temp", imass, ipt), negBinMin, negBinMax);
         TH1D* temp_pos = h2DBackground_D0[imass][ipt]->ProjectionY(
               Form("pos_mass%d_pt%d__temp", imass, ipt), posBinMin, posBinMax);
         TH1D* temp = (TH1D*) temp_neg->Clone();
         temp->SetName(Form("mass%d_pt%d_temp", imass, ipt));
         temp->Add(temp_pos);

         int center = h2DBackground_D0[imass][ipt]->FindBin(0., 0.);
         hDeltaPhi[imass][ipt]->Divide(temp);
         hDeltaPhi[imass][ipt]->Scale(h2DBackground_D0[ipt][imass]->GetBinContent(center) / temp->GetBinWidth(1));

         delete hNeg[imass][ipt];
         delete hPos[imass][ipt];
         delete temp_neg;
         delete temp_pos;
         delete temp;
      }
   }

   /*
   TH2D* h2DSignal_Ref = (TH2D*) f2->Get("hSignal_Ref");
   TH2D* h2DBackground_Ref = (TH2D*) f2->Get("hBackground_Ref");
   TH1D* hMult = (TH1D*) f2->Get("hMult");


   int negBinMin_Ref = 0;
   int negBinMax_Ref = h2DSignal_Ref->GetXaxis()->FindBin(-1.)-1 ;
   int posBinMin_Ref = h2DSignal_Ref->GetXaxis()->FindBin(1.)+1 ;
   int posBinMax_Ref = h2DSignal_Ref->GetXaxis()->GetNbins()+1;

   TH1D* hNeg_Ref_Signal;
   TH1D* hPos_Ref_Signal;
   hNeg_Ref_Signal = h2DSignal_Ref->ProjectionY("neg_ref_signal", negBinMin, negBinMax);
   hPos_Ref_Signal = h2DSignal_Ref->ProjectionY("pos_ref_signal", posBinMin, posBinMax);

   TH1D* hDeltaPhi_Ref = (TH1D*) hNeg_Ref_Signal->Clone();
   hDeltaPhi_Ref->SetName("deltaPhi_Ref");
   hDeltaPhi_Ref->Add(hPos_Ref_Signal);

   TH1D* hNeg_Ref_Background;
   TH1D* hPos_Ref_Background;
   hNeg_Ref_Background = h2DBackground_Ref->ProjectionY("neg_ref_background", negBinMin, negBinMax);
   hPos_Ref_Background = h2DBackground_Ref->ProjectionY("pos_ref_background", posBinMin, posBinMax);

   TH1D* hDeltaPhi_Ref_bkg = (TH1D*) hNeg_Ref_Background->Clone();
   hDeltaPhi_Ref_bkg->SetName("deltaPhi_Ref_bkg");
   hDeltaPhi_Ref_bkg->Add(hPos_Ref_Background);

   int center = hDeltaPhi_Ref_bkg->FindBin(0.);
   hDeltaPhi_Ref->Divide(hDeltaPhi_Ref_bkg);
   hDeltaPhi_Ref->Scale(hDeltaPhi_Ref_bkg->GetBinContent(center)/hDeltaPhi_Ref_bkg->GetBinWidth(1));

   delete hNeg_Ref_Signal;
   delete hPos_Ref_Signal;
   delete hNeg_Ref_Background;
   delete hPos_Ref_Background;
   delete hDeltaPhi_Ref_bkg;
   */
   
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";

   TCanvas* c_deltaPhi[ana::nMass][nPt];
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<nPt; ipt++){
         c_deltaPhi[imass][ipt] = new TCanvas("c_deltaPhi", "", 550, 450);
         TF1 func("deltaPhi", function.c_str(), -3.14159*0.5, 3.14159*1.5);
         func.SetParameter(0, hDeltaPhi[imass][ipt]->GetMaximum());
         func.SetParameter(1, 0.1);
         func.SetParameter(2, 0.1);
         func.SetParameter(3, 0.1);

         hDeltaPhi[imass][ipt]->Fit(&func, "q");
         hDeltaPhi[imass][ipt]->Fit(&func, "q");
         hDeltaPhi[imass][ipt]->Fit(&func, "m q");
         hDeltaPhi[imass][ipt]->Fit(&func, "m q E");
         auto fitResult = hDeltaPhi[imass][ipt]->Fit(&func, "m S E q");

         v2_DCA[imass][ipt] = func.GetParameter(2);
         v2_DCA_err[imass][ipt] = func.GetParError(2);

         hDeltaPhi[imass][ipt]->SetTitle(";#Delta#phi;");

         TLatex ltx;
         ltx.DrawLatexNDC(0.16, 0.85, Form("%.1f<pT<%.1fGeV", ptbin[ipt], ptbin[ipt+1]));
         ltx.DrawLatexNDC(0.16, 0.75, Form("%.3f<mass<%.3fGeV", ana::massbin[imass], ana::massbin[imass+1]));

         std::string str = input_d0;
         auto found = str.find("/");
         while(found!=std::string::npos){
            str.replace(found, 1, "_");
            found = str.find("/");
         }
         //while(str.find(".") == 0){
         while(str.at(0) == '.'){
            str.erase(0, 1);
         }
         while(str.at(0) == '_'){
            str.erase(0, 1);
         }
         c_deltaPhi[imass][ipt]->Print(Form("../plots/d0ana/y%.1f/%s_deltaPhi_mass%d_pt%d.png", 
                     y, str.c_str(), imass, ipt));

         delete c_deltaPhi[imass][ipt];
      }
   }

   /*
   TCanvas* cRef = new TCanvas("cRef", "", 550, 450);
   TF1 func_ref("deltaPhi_Ref", function.c_str(), -3.14159*0.5, 3.14159*1.5);
   func_ref.SetParameter(0, hDeltaPhi_Ref->GetMaximum());
   func_ref.SetParameter(1, 0.1);
   func_ref.SetParameter(2, 0.1);
   func_ref.SetParameter(3, 0.1);
   hDeltaPhi_Ref->Fit(&func_ref, "q");
   hDeltaPhi_Ref->Fit(&func_ref, "q");
   hDeltaPhi_Ref->Fit(&func_ref, "m q");
   auto fitResult = hDeltaPhi_Ref->Fit(&func_ref, "m S E");
   cRef->Print(Form("plots/%s_ref_V2.png", input_ref));
   delete cRef;

   TGraphErrors* g_v2_DCA[nPt][nDca];
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
      g_v2_DCA[ipt][idca] = new TGraphErrors(ana::nMass);
      g_v2_DCA[ipt][idca]->SetName(Form("g_v2_DCA_pt%d_dca%d", ipt, idca));
         for(int imass=0; imass<ana::nMass; imass++){
            double temp = v2_DCA[imass][ipt][idca]/ sqrt(func_ref.GetParameter(2));
            double temp_err = temp* sqrt(pow(v2_DCA_err[imass][ipt][idca]/ v2_DCA[imass][ipt][idca], 2) 
               + pow(0.5*func_ref.GetParError(2)/ func_ref.GetParameter(2), 2));

            v2_DCA[imass][ipt][idca] = temp;
            v2_DCA_err[imass][ipt][idca] = temp_err;

            g_v2_DCA[ipt][idca]->SetPoint(imass, hMass_D0[imass][ipt][idca]->GetMean(), v2_DCA[imass][ipt][idca]);
            g_v2_DCA[ipt][idca]->SetPointError(imass, 0, v2_DCA_err[imass][ipt][idca]);
         }
      }
   }

   TH1D* hMass_DCA[nPt][nDca];
   int mvaBinMin, mvaBinMax, dcaBinMin, dcaBinMax;
   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         mvaBinMin = hDcaVsMassAndMva[ipt]->GetYaxis()->FindBin(ana::mvaCut_NPD0[ipt]+0.1*hDcaVsMassAndMva[ipt]->GetXaxis()->GetBinWidth(1));
         mvaBinMax = hDcaVsMassAndMva[ipt]->GetYaxis()->GetNbins()+1;
         dcaBinMin = hDcaVsMassAndMva[ipt]->GetZaxis()->FindBin(dcaBin[idca]+0.1*hDcaVsMassAndMva[ipt]->GetZaxis()->GetBinWidth(1));
         dcaBinMax = hDcaVsMassAndMva[ipt]->GetZaxis()->FindBin(dcaBin[idca+1]-0.1*hDcaVsMassAndMva[ipt]->GetZaxis()->GetBinWidth(1));
         //dcaBinMax = hDcaVsMassAndMva[ipt]->GetZaxis()->GetNbins()+1;
         hMass_DCA[ipt][idca] = hDcaVsMassAndMva[ipt]->ProjectionX(Form("hmass_pt%d_dca%d", ipt, idca), mvaBinMin, mvaBinMax, dcaBinMin, dcaBinMax);
      }
   }

   TFile f3(Form("%s_v2.root", input_d0), "recreate");

   for(int ipt=0; ipt<nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         g_v2_DCA[ipt][idca]->Write();
         hMass_DCA[ipt][idca]->Write();
      }
   }

   */
   
   return;
}


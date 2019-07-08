#include "include/myAnaConsts.h"

//const int nDca = 11;
//const double dcaBin[nDca+1] = {0., 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013, 0.014};

const int nDca_NonPrompt = 3;
const double dcaBin_NonPrompt[nDca_NonPrompt+1] = {0., 0.009, 0.014, 10000.};

const int nDca_Prompt = 1;
const double dcaBin_Prompt[nDca_Prompt+1] = {0., 100000.};

//void proj1D_full(const char* input_d0= "data/corr2D_npd0ana1_dcaFull_y2.0.root", 
//      const char* input_ref = "data/corr2D_ref_npd0ana1.root", const float y=2.0, const bool isPrompt = false)
void proj1D_full(const char* input_d0= "data/corr2D_d0ana_y2.0.root",
      const char* input_ref = "data/corr2D_ref_d0ana.root", const float y=2.0, const bool isPrompt = true)
{
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   int nDca;
   double* dcaBin;

   if(isPrompt){
      nDca = nDca_Prompt;
      dcaBin = new double[nDca_Prompt+1];
      std::copy(dcaBin_Prompt, dcaBin_Prompt+nDca_Prompt+1, dcaBin);
   }else{
      nDca = nDca_NonPrompt;
      dcaBin = new double[nDca_NonPrompt+1];
      std::copy(dcaBin_NonPrompt, dcaBin_NonPrompt+nDca_NonPrompt+1, dcaBin);
   }

   double v2_DCA[ana::nMass][ana::nPt][nDca];
   double v2_DCA_err[ana::nMass][ana::nPt][nDca];

   TFile* f1 = new TFile(input_d0);
   TFile* f2 = new TFile(input_ref);

   TH3D* hDcaVsMassAndMva[ana::nPt]; 
   
   TH2D* h2DSignal_D0[ana::nMass][ana::nPt][nDca];
   TH2D* h2DBackground_D0[ana::nMass][ana::nPt][nDca];
   TH1D* hMult_raw_D0[ana::nMass][ana::nPt][nDca]; // wrong normalized constant
   TH1D* hMass_D0[ana::nMass][ana::nPt][nDca];

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<ana::nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            h2DSignal_D0[imass][ipt][idca] = (TH2D*) f1->Get(Form("hSignal_mass%d_pt%d_dca%d", imass, ipt, idca));
            h2DBackground_D0[imass][ipt][idca] = (TH2D*) f1->Get(Form("hBackground_mass%d_pt%d_dca%d", imass, ipt, idca));
            hMult_raw_D0[imass][ipt][idca] = (TH1D*) f1->Get(Form("hMult_raw_D0_mass%d_pt%d_dca%d", imass, ipt, idca));
            hMass_D0[imass][ipt][idca] = (TH1D*) f1->Get(Form("hMassD0_mass%d_pt%d_dca%d", imass, ipt, idca));
         }
      }
   }

   for(int ipt=0; ipt<ana::nPt; ipt++){
      hDcaVsMassAndMva[ipt] = (TH3D*) f1->Get(Form("hDcaVsMassAndMva_pt%d", ipt));
   }

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<ana::nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
   cout << "test" <<endl;
            double nMult_D0= hMult_raw_D0[imass][ipt][idca]->Integral(2, 251);
            h2DSignal_D0[imass][ipt][idca]->Scale(1./nMult_D0);
         }
      }
   }

   TH1D* hNeg[ana::nMass][ana::nPt][nDca];
   TH1D* hPos[ana::nMass][ana::nPt][nDca];
   TH1D* hDeltaPhi[ana::nMass][ana::nPt][nDca];

   int negBinMin = 0;
   int negBinMax = h2DSignal_D0[0][0][0]->GetXaxis()->FindBin(-1.)-1 ;
   int posBinMin = h2DSignal_D0[0][0][0]->GetXaxis()->FindBin(1.)+1 ;
   int posBinMax = h2DSignal_D0[0][0][0]->GetXaxis()->GetNbins()+1;

   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<ana::nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            hNeg[imass][ipt][idca] = h2DSignal_D0[imass][ipt][idca]->ProjectionY(
                  Form("neg_mass%d_pt%d_dca%d", imass, ipt, idca), negBinMin, negBinMax);
            hPos[imass][ipt][idca] = h2DSignal_D0[imass][ipt][idca]->ProjectionY(
                  Form("pos_mass%d_pt%d_dca%d", imass, ipt, idca), posBinMin, posBinMax);
            hDeltaPhi[imass][ipt][idca] = (TH1D*) hNeg[imass][ipt][idca]->Clone();
            hDeltaPhi[imass][ipt][idca]->SetName(Form("deltaPhi_mass%d_pt%d_dca%d", imass, ipt, idca));
            hDeltaPhi[imass][ipt][idca]->Add(hPos[imass][ipt][idca]);

            TH1D* temp_neg = h2DBackground_D0[imass][ipt][idca]->ProjectionY(
                  Form("neg_mass%d_pt%d_dca%d_temp", imass, ipt, idca), negBinMin, negBinMax);
            TH1D* temp_pos = h2DBackground_D0[imass][ipt][idca]->ProjectionY(
                  Form("pos_mass%d_pt%d_dca%d_temp", imass, ipt, idca), posBinMin, posBinMax);
            TH1D* temp = (TH1D*) temp_neg->Clone();
            temp->SetName(Form("mass%d_pt%d_dca%d_temp", imass, ipt, idca));
            temp->Add(temp_pos);

            int center = temp->FindBin(0.);
            hDeltaPhi[imass][ipt][idca]->Divide(temp);
            hDeltaPhi[imass][ipt][idca]->Scale(temp->GetBinContent(center) / temp->GetBinWidth(1));

            delete hNeg[imass][ipt][idca];
            delete hPos[imass][ipt][idca];
            delete temp_neg;
            delete temp_pos;
            delete temp;
         }
      }
   }

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
   
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";

   TCanvas* c_deltaPhi[ana::nMass][ana::nPt][nDca];
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<ana::nPt; ipt++){
         for(int idca=0; idca<nDca; idca++){
            //std::cout << imass << " " << ipt << " " << idca << std::endl;
            c_deltaPhi[imass][ipt][idca] = new TCanvas("c_deltaPhi", "", 550, 450);
            TF1 func("deltaPhi", function.c_str(), -3.14159*0.5, 3.14159*1.5);
            func.SetParameter(0, hDeltaPhi[imass][ipt][idca]->GetMaximum());
            func.SetParameter(1, 0.1);
            func.SetParameter(2, 0.1);
            func.SetParameter(3, 0.1);

            hDeltaPhi[imass][ipt][idca]->Fit(&func, "q");
            hDeltaPhi[imass][ipt][idca]->Fit(&func, "q");
            hDeltaPhi[imass][ipt][idca]->Fit(&func, "m q");
            hDeltaPhi[imass][ipt][idca]->Fit(&func, "m q E");
            auto fitResult = hDeltaPhi[imass][ipt][idca]->Fit(&func, "m S E q");
            //std::cout << fitResult->Status() << std::endl;

            v2_DCA[imass][ipt][idca] = func.GetParameter(2);
            v2_DCA_err[imass][ipt][idca] = func.GetParError(2);

            hDeltaPhi[imass][ipt][idca]->SetTitle(";#Delta#phi;");

            TLatex ltx;
            ltx.DrawLatexNDC(0.16, 0.85, Form("%.1f<pT<%.1fGeV", ana::ptbin[ipt], ana::ptbin[ipt+1]));
            ltx.DrawLatexNDC(0.16, 0.75, Form("%.3f<mass<%.3fGeV", ana::massbin[imass], ana::massbin[imass+1]));
            ltx.DrawLatexNDC(0.16, 0.65, Form("%.2f<DCA<%.2fum", dcaBin[idca]*1000, dcaBin[idca+1]*1000));

            if(!isPrompt) c_deltaPhi[imass][ipt][idca]->Print(Form("plots/dcaFull/y%.1f/%s_deltaPhi_mass%d_pt%d_dca%.02f-%0.2f.png", 
                     y, input_d0, imass, ipt, dcaBin[idca]*1000, dcaBin[idca+1]*1000));
            if(isPrompt) c_deltaPhi[imass][ipt][idca]->Print(Form("plots/d0ana/y%.1f/%s_deltaPhi_mass%d_pt%d_dca%.02f-%0.2f.png", 
                     y, input_d0, imass, ipt, dcaBin[idca]*1000, dcaBin[idca+1]*1000));

            delete c_deltaPhi[imass][ipt][idca];
         }
      }
   }

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

   TGraphErrors* g_v2_DCA[ana::nPt][nDca];
   for(int ipt=0; ipt<ana::nPt; ipt++){
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

   TH1D* hMass_DCA[ana::nPt][nDca];
   int mvaBinMin, mvaBinMax, dcaBinMin, dcaBinMax;
   for(int ipt=0; ipt<ana::nPt; ipt++){
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

   for(int ipt=0; ipt<ana::nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         g_v2_DCA[ipt][idca]->Write();
         hMass_DCA[ipt][idca]->Write();
      }
   }

   delete[] dcaBin;
   
   return;
}

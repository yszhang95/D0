// this macro does not work

#include "myAnaConsts.h"
#include "massfitting.C"
#include <iostream>
#include <memory>
TH1D* hDcaMCAll(TH2*, TH1&, const std::string&);
void dcaHists(int mode = 1)
{
   //TGaxis::SetMaxDigits(3);
   gStyle->SetOptStat(0);
   //TFile* f1 = new TFile(Form("%s_hists.root", ana::whichtree[mode].c_str()));
   std::unique_ptr<TFile> f1 = std::unique_ptr<TFile>(new TFile(Form("%s_hists.root", ana::whichtree[mode].c_str())));
   TH2D* hDcaVsMassPD0 = (TH2D*) f1->Get("hDcaVsMassPD0");
   TH2D* hDcaVsMassPD0_All = (TH2D*) f1->Get("hDcaVsMassPD0_All");
   TH2D* hDcaVsMassNPD0 = (TH2D*) f1->Get("hDcaVsMassNPD0");
   TH2D* hDcaVsMassNPD0_All = (TH2D*) f1->Get("hDcaVsMassNPD0_All");
   TH2D* hDcaVsMassDataD0 = (TH2D*) f1->Get("hDcaVsMassDataD0");

   // DCA distribution of Prompt D0 
  
   TCanvas cPrompt("cPrompt", "", 550, 450);
   cPrompt.SetLogy();
   //std::unique_ptr<TH1> hDcaMCPD0 = std::unique_ptr<TH1>(new TH1D("hDcaMCPD0", "hDCaMCPD0", ana::nuofDca, ana::dcaBin));
   TH1D hDcaMCPD0("hDcaMCPD0", "hDCaMCPD0", ana::nuofDca, ana::dcaBin);
   TH1D* hDcaMCPD0_Proj = hDcaMCAll(hDcaVsMassPD0, hDcaMCPD0, "hprompt");
   hDcaMCPD0.Draw();
   hDcaMCPD0_Proj->Scale(hDcaMCPD0.GetMaximum()/hDcaMCPD0_Proj->GetMaximum());
   hDcaMCPD0_Proj->SetLineColor(kRed);
   hDcaMCPD0_Proj->Draw("same");
   cPrompt.Print("cPrompt.png");
   delete hDcaMCPD0_Proj;

   // DCA distribution of Non-Prompt D0

   TCanvas cNonPrompt("cNonPrompt", "", 550, 450);
   cNonPrompt.SetLogy();
   TH1D hDcaMCNPD0("hDcaMCNPD0", "hDCaMCNPD0", ana::nuofDca, ana::dcaBin);
   TH1D* hDcaMCNPD0_Proj = hDcaMCAll(hDcaVsMassNPD0, hDcaMCNPD0, "hnonprompt");
   hDcaMCNPD0.Draw();
   hDcaMCNPD0_Proj->Scale(hDcaMCNPD0.GetMaximum()/hDcaMCNPD0_Proj->GetMaximum());
   hDcaMCNPD0_Proj->SetLineColor(kRed);
   hDcaMCNPD0_Proj->Draw("same");
   cNonPrompt.Print("cNonPrompt.png");
   delete hDcaMCNPD0_Proj;

   // DCA distribution of signal + swap, by fitting mass;

   TH1D* hMCPromptMass = hDcaVsMassPD0->ProjectionX("hPromptMass");
   TH1D* hMCPromptMassAll = hDcaVsMassPD0_All->ProjectionX("hPromptMassAll");
   TH1D* hMCNonPromptMass = hDcaVsMassNPD0->ProjectionX("hNonPromptMass");
   TH1D* hMCNonPromptMassAll = hDcaVsMassNPD0_All->ProjectionX("hNonPromptMassAll");

   TH1D hDcaDataD0("hDcaDataD0", "hDcaDataD0", ana::nuofDca, ana::dcaBin);
   TH1D hDcaDataD0SideBand("hDcaDataD0SideBand", "hDcaDataD0SideBand", ana::nuofDca, ana::dcaBin);

   TH1D* hTemp;
   TCanvas cTemp("cTemp", "", 550, 450);
   for(int iDca=0; iDca<ana::nuofDca; iDca++){
      int binlw = hDcaVsMassDataD0->GetYaxis()->FindBin(ana::dcaBin[iDca]); // x: mass, y: dca
      int binup = hDcaVsMassDataD0->GetYaxis()->FindBin(ana::dcaBin[iDca+1]) -1; // x: mass, y: dca
      hTemp = hDcaVsMassDataD0->ProjectionX("tmp", binlw, binup); 
      hTemp->GetYaxis()->SetRangeUser(0, hTemp->GetMaximum() * 1.3);
      TF1 f = massfitting(hTemp, hMCNonPromptMass, hMCNonPromptMassAll, "iDca");
      TF1 signal("signal", "[0]* [5] * (" "[4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))"
            "+ (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))" ")", 1.7, 2.0);
      signal.SetLineColor(kOrange-3);
      signal.SetLineWidth(1);
      signal.SetLineStyle(2);
      signal.SetFillColorAlpha(kOrange-3,0.3);
      signal.SetFillStyle(1001);
      for(int ipar=0; ipar<6+1; ipar++){
         signal.FixParameter(ipar, f.GetParameter(ipar));
      }

      TF1 swap("swap", "[0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))"
            "+0 *[1]*[2]*[3]*[4]", 1.7, 2.0);
      swap.SetLineColor(kGreen+4);
      swap.SetLineWidth(1);
      swap.SetLineStyle(1);
      swap.SetFillColorAlpha(kGreen+4,0.3);
      swap.SetFillStyle(1001);
      for(int ipar=0; ipar<8+1; ipar++){
         swap.FixParameter(ipar, f.GetParameter(ipar));
      }

      TF1 bkg("bkg", "[9] + [10]*x + [11]*x*x + [12]*x*x*x"
            "+ 0 *[0]*[1]*[2]*[3]*[4]*[5]*[6]*[7]*[8]", 1.7, 2.0);
      bkg.SetLineColor(4);
      bkg.SetLineWidth(1);
      bkg.SetLineStyle(2);
      for(int ipar=0; ipar<12+1; ipar++){
         bkg.FixParameter(ipar, f.GetParameter(ipar));
      }

      signal.Draw("same FC");
      swap.Draw("same FC");
      bkg.Draw("same L");
      cTemp.Print(Form("temp%d.png", iDca));

      double binWidth = hTemp->GetBinWidth(1) * (ana::dcaBin[iDca+1] - ana::dcaBin[iDca]);
      hDcaDataD0.SetBinContent(iDca+1, f.GetParameter(0)/binWidth);
      hDcaDataD0.SetBinError(iDca+1, f.GetParError(0)/binWidth);

      std::cout << bkg.Integral(1.7, 1.76)/f.Integral(1.7, 1.76) << std::endl;

      hTemp->Delete();
   }

   // extract DCA distribution of signal, by side band method
   TH1D* hMassDataD0 = hDcaVsMassDataD0->ProjectionX();
   TF1 f = massfitting(hMassDataD0, hMCNonPromptMass, hMCNonPromptMassAll, "mass");
   hMassDataD0->GetYaxis()->SetRangeUser(0, 1.3* hMassDataD0->GetMaximum());
   hMassDataD0->Draw();
   TF1 signal("signal", "[0]* [5] * (" "[4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))"
            "+ (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))" ")", 1.7, 2.0);
   signal.SetLineColor(kOrange-3);
   signal.SetLineWidth(1);
   signal.SetLineStyle(2);
   signal.SetFillColorAlpha(kOrange-3,0.3);
   signal.SetFillStyle(1001);
   for(int ipar=0; ipar<6+1; ipar++){
      signal.FixParameter(ipar, f.GetParameter(ipar));
   }

   TF1 swap("swap", "[0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))"
            "+0 *[1]*[2]*[3]*[4]", 1.7, 2.0);
   swap.SetLineColor(kGreen+4);
   swap.SetLineWidth(1);
   swap.SetLineStyle(1);
   swap.SetFillColorAlpha(kGreen+4,0.3);
   swap.SetFillStyle(1001);
   for(int ipar=0; ipar<8+1; ipar++){
      swap.FixParameter(ipar, f.GetParameter(ipar));
   }

   TF1 bkg("bkg", "[9] + [10]*x + [11]*x*x + [12]*x*x*x"
            "+ 0 *[0]*[1]*[2]*[3]*[4]*[5]*[6]*[7]*[8]", 1.7, 2.0);
   bkg.SetLineColor(4);
   bkg.SetLineWidth(1);
   bkg.SetLineStyle(2);
   for(int ipar=0; ipar<12+1; ipar++){
      bkg.FixParameter(ipar, f.GetParameter(ipar));
   }

   signal.Draw("same FC");
   swap.Draw("same FC");
   bkg.Draw("same L");
   cTemp.Print("tempAll.png");

   TF1 s(f);
   s.FixParameter(9, 0);
   s.FixParameter(10, 0);
   s.FixParameter(11, 0);
   s.FixParameter(12, 0);

   std::map<std::string, double> fracLeft, fracPeak, fracRight;
   fracLeft["signal"] = signal.Integral(1.7, 1.8)/ f.Integral(1.7, 1.8);
   fracLeft["swap"] = swap.Integral(1.7, 1.8)/ f.Integral(1.7, 1.8);
   fracLeft["bkg"] = bkg.Integral(1.7, 1.8)/ f.Integral(1.7, 1.8);
   fracPeak["signal"] = signal.Integral(1.85, 1.88)/ f.Integral(1.85, 1.88);
   fracPeak["swap"] = swap.Integral(1.85, 1.88)/ f.Integral(1.85, 1.88);
   fracPeak["bkg"] = bkg.Integral(1.85, 1.88)/ f.Integral(1.85, 1.88);
   fracRight["signal"] = signal.Integral(1.94, 2.0)/ f.Integral(1.94, 2.0);
   fracRight["swap"] = swap.Integral(1.94, 2.0)/ f.Integral(1.94, 2.0);
   fracRight["bkg"] = bkg.Integral(1.94, 2.0)/ f.Integral(1.94, 2.0);

   hMassDataD0->Delete(); 
   hMassDataD0 = nullptr;

   TCanvas cData("cData", "", 550, 450);
   cData.SetLogy();
   hDcaDataD0.Draw();
   cData.Print("data.png");

   // clear pointers

   hMCPromptMass->Delete(); 
   hMCPromptMass = nullptr;
   hMCPromptMassAll->Delete();
   hMCPromptMassAll = nullptr;
   hMCNonPromptMass->Delete();
   hMCNonPromptMass = nullptr;
   hMCNonPromptMassAll->Delete();
   hMCNonPromptMassAll = nullptr;

   TFile f2(Form("%s_dca_hists.root", ana::whichtree[mode].c_str()), "recreate");
   hDcaDataD0.Write();
   hDcaMCPD0.Write();
   hDcaMCNPD0.Write();
   f.Write();
}
TH1D* hDcaMCAll(TH2* hDcaVsMass, TH1& hDcaMC, const  std::string& proj_name){
   TH1D* hDcaMCProj = hDcaVsMass->ProjectionY(proj_name.c_str(), 1, 60); // mass range
   for(int iDca=0; iDca<ana::nuofDca; iDca++){
      int binlw = hDcaMCProj->FindBin(ana::dcaBin[iDca]);
      int binup = hDcaMCProj->FindBin(ana::dcaBin[iDca+1]) - 1;
      double binwidth = ana::dcaBin[iDca+1] - ana::dcaBin[iDca];
      double err = 0;
      double yield = hDcaMCProj->IntegralAndError(binlw, binup, err);
      hDcaMC.SetBinContent(iDca+1, yield/binwidth);
      hDcaMC.SetBinError(iDca+1, err/binwidth);
   }
   
   return hDcaMCProj;
}

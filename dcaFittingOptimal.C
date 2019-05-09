#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TFitter.h"
#include "TFitResult.h"

#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooDataHist.h"

#include "myAnaConsts.h"
#include "massfitting.C"

#include <vector>

using namespace RooFit;

TH1D* hD0DcaMCPSignal;
TH1D* hD0DcaMCNPSignal;
TH1D* hD0DcaData;

void dcaFittingOptimal(int mode=1, int method = 0)
{
   std::unique_ptr<TFile> f1 = std::unique_ptr<TFile>(new TFile(Form("%s_dca_hists.root", ana::whichtree[mode].c_str())));
   double sig[9] = {0};
   //for(int iMva=0; iMva<9; iMva++)
   for(int iMva=0; iMva<9; iMva++)
   {
      int label = 40 + 2*iMva;
      hD0DcaData = (TH1D*) f1->Get(Form("hDcaDataD0mva%d", label));
      hD0DcaMCPSignal = (TH1D*) f1->Get(Form("hDcaMCPD0mva%d", label));
      hD0DcaMCNPSignal = (TH1D*) f1->Get(Form("hDcaMCNPD0mva%d", label));

      double yield_PD0 = hD0DcaMCPSignal->Integral("width");
      double yield_NPD0 = hD0DcaMCNPSignal->Integral("width");
      double yield_Data = hD0DcaData->Integral("width");
      hD0DcaMCPSignal->Scale(1./yield_PD0);
      hD0DcaMCNPSignal->Scale(1./yield_NPD0);
      hD0DcaData->Scale(1./yield_Data);

      TH1D* hMassData = (TH1D*) f1->Get(Form("hMassDataD0mva%d", label));
      TH1D* hMassMCNPD0 = (TH1D*) f1->Get(Form("hMassMCNPD0mva%d", label));
      TH1D* hMassMCNPD0All = (TH1D*) f1->Get(Form("hMassMCNPD0Allmva%d", label));
      TH1D* hMassMCPD0 = (TH1D*) f1->Get(Form("hMassMCPD0mva%d", label));
      TH1D* hMassMCPD0All = (TH1D*) f1->Get(Form("hMassMCPD0Allmva%d", label));

      double massBinWidth = hMassData->GetBinWidth(1);
      hMassData->Scale(1./massBinWidth);
      hMassMCNPD0->Scale(1./massBinWidth);
      hMassMCNPD0All->Scale(1./massBinWidth);
      hMassMCPD0->Scale(1./massBinWidth);
      hMassMCPD0All->Scale(1./massBinWidth);

      // create Pic to store pictures
      const std::string dirPic = "if [ ! -d \"Pic\" ]; then\n"
                                 "    mkdir Pic \n"
                                 "fi";
      gSystem->Exec(dirPic.c_str());

      if(method == 0){
         TLatex ltx;
         ltx.SetTextSize(0.04);

         //create Pic/mva to store pictures
         std::string dirMvaPic = Form( "if [ ! -d \"Pic/mva%d\" ]; then\n"
                                       "    mkdir Pic/mva%d \n"
                                       "fi", iMva, iMva);
         gSystem->Exec(dirMvaPic.c_str());

         TObjArray mc(2);
         mc.Add(hD0DcaMCNPSignal);
         mc.Add(hD0DcaMCPSignal);
         TFractionFitter fit(hD0DcaData, &mc);
         fit.Constrain(0, 0.0, 1.0);
         fit.Constrain(1, 0.0, 1.0);
         fit.SetRangeX(1, 9);
         if(iMva==8) fit.SetRangeX(1, 8);
         int status = fit.Fit();
         std::cout << "fit status: " << status << std::endl;

         TCanvas cDca("cDca", "", 550, 450);
         //cDca.SetLogy();
         cDca.SetLeftMargin(0.16);
         cDca.SetBottomMargin(0.16);

         gStyle->SetErrorX(0);

         TH1D* result = (TH1D*) fit.GetPlot();
         hD0DcaData->GetXaxis()->SetRangeUser(0, 0.036);
         hD0DcaData->SetMarkerStyle(20);
         hD0DcaData->SetTitle(";DCA (cm);Normalized dN / d(D0 Dca) (cm^{-1})");
         hD0DcaData->Draw("E P");
         result->SetLineColor(kBlack);
         result->Draw("SAME");

         double fracNPD0, fracErrNPD0;
         double fracPD0, fracErrPD0;
         fit.GetResult(0, fracNPD0, fracErrNPD0);
         fit.GetResult(1, fracPD0, fracErrPD0);

         /*
         hD0DcaMCNPSignal->Scale(fracNPD0);
         hD0DcaMCPSignal->Scale(fracPD0);
         std::cout << hD0DcaMCNPSignal->GetBinContent(1) << std::endl;
         std::cout << hD0DcaMCPSignal->GetBinContent(1) << std::endl;
         hD0DcaMCNPSignal->SetFillColor(kBlue);
         hD0DcaMCNPSignal->SetFillStyle(1001);
         hD0DcaMCPSignal->SetFillColor(kRed);
         hD0DcaMCPSignal->SetFillStyle(1001);
         */

         TH1D* hMCPredctionNPD0 = (TH1D*) fit.GetMCPrediction(0);
         hMCPredctionNPD0->Scale(fracNPD0);
         hMCPredctionNPD0->SetFillColor(kBlue);
         //hMCPredctionNPD0->SetFillStyle(1001);
         hMCPredctionNPD0->SetFillStyle(3001);
         TH1D* hMCPredctionPD0 = (TH1D*) fit.GetMCPrediction(1);
         hMCPredctionPD0->Scale(fracPD0);
         hMCPredctionPD0->SetFillColor(kRed);
         hMCPredctionPD0->SetFillStyle(1001);

         THStack hs("hs", "");
         //hs.Add(hD0DcaMCPSignal);
         //hs.Add(hD0DcaMCNPSignal);
         hs.Add(hMCPredctionPD0);
         hs.Add(hMCPredctionNPD0);
         hs.Draw("SAME HIST nostack");
         result->Draw("SAME");
         hD0DcaData->Draw("E P SAME ][");

         TLegend lgd(0.55, 0.7, 0.9, 0.9);
         lgd.AddEntry(hD0DcaData, "Data", "p");
         lgd.AddEntry(hMCPredctionNPD0, "Predicted MC B #rightarrow D", "f");
         lgd.AddEntry(hMCPredctionPD0, "Predicted MC Prompt D", "f");
         lgd.AddEntry(result, "Fitting", "l");
         lgd.Draw();

         ltx.DrawLatexNDC(0.35, 0.8, Form("MVA > %.2f", 0.4+0.02*iMva));
         ltx.DrawLatexNDC(0.2, 0.7, "4 < pT < 5GeV  |y|<1");

         double realFracNPD0 = fracNPD0/(fracNPD0+fracPD0);
         // err = errNPD(1/total - 1/total * realFrac) cross product errPD/total * realFrac
         double realFracNPD0Err = std::sqrt(std::pow(fracErrNPD0/(fracNPD0+fracPD0)*(1-realFracNPD0Err), 2) 
               + std::pow(fracErrPD0/(fracNPD0+fracPD0)*realFracNPD0Err, 2));

         ltx.SetTextSize(0.035);
         ltx.DrawLatexNDC(0.55, 0.5, Form("B2D Frac. = %.2f+/-%.2f", realFracNPD0, realFracNPD0Err));

         cDca.Print(Form("Pic/mva%d/DcaFitting.png", iMva));
         cDca.SetLogy();
         cDca.Print(Form("Pic/mva%d/DcaFittingLogScale.png", iMva));

         // fit the invariant mass
         TFitResultPtr fitResultPtr;
         TF1 f = massfitting(hMassData, hMassMCNPD0, hMassMCNPD0All, TString::Format("mva%d", iMva), fitResultPtr);

         // draw the invariant mass
         hMassData->GetYaxis()->SetTitle("Entries /(5 MeV)");
         drawMassFitting(hMassData, f, TString::Format("Pic/mva%d/MassFitting.png", iMva), 
               "4 < pT < 5GeV  |y|<1", TString::Format("MVA > %.2f", 0.4+0.02*iMva));

         // calculate the significance
         TF1 signal("signal", "[0]* [5] * (" "[4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))"
               "+ (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))" ")", 1.7, 2.0);
         for(int ipar=0; ipar<6+1; ipar++){
            signal.FixParameter(ipar, f.GetParameter(ipar));
         }

         double max = signal.GetMaximum(1.7, 2.0);
         double xmin = signal.GetX(max/2., 1.84, 1.865);
         double xmax = signal.GetX(max/2., 1.865, 1.9);

         sig[iMva] = signal.Integral(xmin, xmax) * realFracNPD0 / sqrt(f.Integral(xmin, xmax));
      }
      if(method == 1){
         // following does not work
         TCanvas cTemp("cTemp", "", 550, 450);
         RooRealVar x("x", "x", 0, 0.036);
         RooRealVar fracNPD0("fracNPD0", "fracNPD0", 0.5);
         RooDataHist hNP("hNP", "hNP", RooArgList(x), hD0DcaMCNPSignal);
         RooDataHist hP("hP", "hP", RooArgList(x), hD0DcaMCPSignal);
         RooDataHist hData("hData", "hData", RooArgList(x), hD0DcaData);
         RooHistPdf fNP("fNP", "fNP", x, hNP, 0);
         RooHistPdf fP("fP", "fP", x, hP, 0);
         RooAddPdf sum("sum", "NPD0+PD0", RooArgList(fNP, fP), RooArgList(fracNPD0));
         sum.fitTo(hData, SumW2Error(kTRUE), Extended(kTRUE));
         RooArgSet *parList = sum.getParameters(hData);
         parList->Print("v");
         cTemp.cd();
         RooPlot* xframe = x.frame();
         hData.plotOn(xframe);
         sum.plotOn(xframe);
         xframe->Draw();
         cTemp.Print("test.png");
      }
   }
   for(int iMva=0; iMva<9; iMva++){
      std::cout << sig[iMva] << std::endl;
   }
   double mva[9];
   for(int iMva=0; iMva<9; iMva++){
      mva[iMva] = 0.4 + ana::mvaStep * iMva;
   }
   TCanvas cSig("cSig", "", 550, 450);
   auto hist = cSig.DrawFrame(0.36, 0, 0.6, 10);
   hist->SetTitle(";MVA;S_{B2D}/#sqrt{S_{D}+B}");
   TGraph gSig(9, mva, sig);
   gSig.SetMarkerStyle(20);
   gSig.SetLineColor(kBlue);
   gSig.Draw("PC");
   TLatex ltx;
   ltx.DrawLatexNDC(0.55, 0.55, "4 < pT < 5GeV  |y|<1");
   cSig.Print("Pic/Sig.png");
}

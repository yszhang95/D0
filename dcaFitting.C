// this macro does not work

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

#include "myAnaConsts.h"

#include <vector>

TH1D* hD0DcaMCPSignal;
TH1D* hD0DcaMCNPSignal;
TH1D* hD0DcaData;

Double_t funMix(Double_t* x_, Double_t* para)
{
    float x = x_[0];
    float APrompt = para[0];
    float ANonPrompt = para[1];
    float promptYield = 0;
    float nonPromptYield = 0;
    
    promptYield = hD0DcaMCPSignal->GetBinContent(hD0DcaMCPSignal->GetXaxis()->FindBin(x));
    nonPromptYield = hD0DcaMCNPSignal->GetBinContent(hD0DcaMCNPSignal->GetXaxis()->FindBin(x));
    
    return APrompt*promptYield+ANonPrompt*nonPromptYield;
}

Double_t funNonPrompt(Double_t* x_, Double_t* para)
{
    float x = x_[0];
    float APrompt = para[0];
    float ANonPrompt = para[1];
    float nonPromptYield = 0;
    nonPromptYield = hD0DcaMCNPSignal->GetBinContent(hD0DcaMCNPSignal->GetXaxis()->FindBin(x));
    return ANonPrompt*nonPromptYield;
}

Double_t funPrompt(Double_t* x_, Double_t* para)
{
    float x = x_[0];
    float APrompt = para[0];
    float ANonPrompt = para[1];
    float nonPromptYield = 0;
    promptYield = hD0DcaMCPSignal->GetBinContent(hD0DcaMCPSignal->GetXaxis()->FindBin(x));
    //nonPromptYield = hD0DcaMCNPSignal->GetBinContent(hD0DcaMCNPSignal->GetXaxis()->FindBin(x));
    return APrompt*promptYield;
}

void dcaFitting(int mode=1, int method = 0)
{
   std::unique_ptr<TFile> f1 = std::unique_ptr<TFile>(new TFile(Form("%s_dca_hists.root", ana::whichtree[mode].c_str())));
   hD0DcaData = (TH1D*) f1->Get("hDcaDataD0");
   hD0DcaMCPSignal = (TH1D*) f1->Get("hDcaMCPD0");
   hD0DcaMCNPSignal = (TH1D*) f1->Get("hDcaMCNPD0");
   TF1* f = (TF1*) f1->Get("mass");
   double yield_PD0 = hD0DcaMCPSignal->Integral("width");
   double yield_NPD0 = hD0DcaMCNPSignal->Integral("width");
   double yield_Data = hD0DcaData->Integral("width");
   hD0DcaMCPSignal->Scale(1./yield_PD0);
   hD0DcaMCNPSignal->Scale(1./yield_NPD0);
   hD0DcaData->Scale(1./yield_Data);
   if(method == 0){
      TObjArray mc(2);
      mc.Add(hD0DcaMCNPSignal);
      mc.Add(hD0DcaMCPSignal);
      TFractionFitter fit(hD0DcaData, &mc);
      fit.Constrain(0, 0.0, 1.0);
      fit.Constrain(1, 0.0, 1.0);
      fit.SetRangeX(1, 9);
      int status = fit.Fit();
      std::cout << "fit status: " << status << std::endl;
      TCanvas cDca("cDca", "", 550, 450);
      cDca.SetLogy();
      TH1D* result = (TH1D*) fit.GetPlot();
      hD0DcaData->Draw("ep");
      result->Draw("same");
      cDca.Print("DcaFitting.png");

      double frac, error;
      fit.GetResult(0, frac, error);
      TF1 s(*f);
      s.FixParameter(9, 0);
      s.FixParameter(10, 0);
      s.FixParameter(11, 0);
      s.FixParameter(12, 0);
      double signal = s.Integral(1.7, 2.0)/0.005;
      double all = f->Integral(1.7, 2.0)/0.005;
      std::cout << signal*frac/std::sqrt(signal*frac + all -signal) <<  std::endl;
   }
}

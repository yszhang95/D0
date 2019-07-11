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
#include "TFitResultPtr.h"

#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooDataHist.h"

#include "include/myAnaConsts.h"
#include "macros/massfitting.C"

#include <vector>


const double ptbin[] = {2., 4., 6., 8.};

const bool isUseMassFit = true;

//const int nuofDca = 23;
//const double dcaBin[nuofDca+1] = {0, 0.001, 0.002, 0.003, 0.004, 0.005, 
   //0.006, 0.007, 0.008, 0.009, 0.010,
   //0.011, 0.012, 0.013, 0.014, 0.016, 
   //0.022, 0.026, 0.030, 0.036, 0.042, 
   //0.05, 0.062, 0.08}; // cm

const int nuofDca = 37;
const double dcaBin[nuofDca+1] = {0, 0.0005, 0.001,0.0015, 0.002, 
   0.0025, 0.003, 0.0035, 0.004, 0.0045, 
   0.005, 0.0055, 0.006, 0.0065, 0.007, 
   0.0075,  0.008, 0.0085, 0.009, 0.0095, 
   0.010, 0.0105, 0.011, 0.0115, 0.012, 
   0.0125, 0.013, 0.0135, 0.014, 0.016, 
   0.022, 0.026, 0.030, 0.036, 0.042, 
   0.05, 0.062, 0.08}; // cm

//const int nDca = 11;
//const double dcaCut[nDca] = {0., 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013, 0.014};
const int nDca = 3;
const double dcaCut[nDca+1] = {0., 0.009, 0.014, 1000000.};

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
    float promptYield = 0;
    promptYield = hD0DcaMCPSignal->GetBinContent(hD0DcaMCPSignal->GetXaxis()->FindBin(x));
    return APrompt*promptYield;
}

void dcaFractionFitting(
      const float& fitRangeL,
      const float& fitRangeH,
      TFitResultPtr& fitResult,
      const std::map<std::string, std::string>& picName, 
      const std::string& cut, const std::string& label)
{


   TCanvas cFit("cFit", "", 550, 650);
   cFit.Divide(1, 2, 0, 0);

   TH1D hDrawRatio("hDrawRatio", "", 10000, 0, 1);

   auto pad1 = cFit.cd(1);
   pad1->SetBottomMargin(0);
   pad1->SetLeftMargin(0.16);


   double integralTotalYield = hD0DcaData->Integral(1,hD0DcaData->GetXaxis()->GetNbins(), "width");
   cout<<"Data Total yield: "<<integralTotalYield<<endl;
   TF1 fMix(Form("fMix"),&funMix, -0.5, 0.5, 2);
   fMix.SetParameters(0.5*integralTotalYield,0.5*integralTotalYield);
   fMix.SetParLimits(0, 0, 1.2*integralTotalYield);
   fMix.SetParLimits(1, 0.0*integralTotalYield, 1.2*integralTotalYield);
  
   fMix.SetLineColor(2);
   fMix.SetFillColorAlpha(kRed,0.36);
   fMix.SetFillStyle(1001);
 
   //Set fit range

   int fitStatus = 1;
   double fitPrecision = 1.e-10;
   while(fitStatus)
   {
       TFitter::SetPrecision(fitPrecision);
       fMix.SetParameters(0.5*integralTotalYield,0.5*integralTotalYield);
       fMix.SetParError(0, 0.5*integralTotalYield);
       fMix.SetParError(1, 0.5*integralTotalYield);
       fitResult = hD0DcaData->Fit(&fMix, "I E SNQ0", "", fitRangeL, fitRangeH);
       fitStatus = fitResult->Status();
       if(fitStatus)
          fitPrecision *= 5;
   }
   
   cout<<"============== do main fit ============"<<endl;
   fMix.SetParameters(0.5*integralTotalYield, 0.5*integralTotalYield);
   fMix.SetParError(0, 0.5*integralTotalYield);
   fMix.SetParError(1, 0.5*integralTotalYield);
   fMix.SetNpx(10000);
   fitResult = hD0DcaData->Fit(&fMix,"I E S0", "", fitRangeL, fitRangeH);
   fMix.SetRange(fitRangeL,fitRangeH);
   hD0DcaData->SetMarkerStyle(20);
   hD0DcaData->GetXaxis()->SetRangeUser(fitRangeL, fitRangeH);
   //hD0DcaData->GetYaxis()->SetRangeUser(hD0DcaData->GetMinimum()*0.5, hD0DcaData->GetMaximum()*1.5);
   hD0DcaData->SetTitle(";;dN/dDca (cm^{-1})");
   hD0DcaData->GetYaxis()->SetTitleSize(0.05);
   hD0DcaData->Draw();
   fMix.Draw("same");
   fitStatus = fitResult->Status();
   cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<endl;

   TF1 fP("fP",&funPrompt, 0., 0.5, 2);
   fP.SetParameters(fMix.GetParameter(0),fMix.GetParameter(1));
   fP.SetRange(fitRangeL,fitRangeH);
   fP.SetLineColor(4);
   fP.SetFillStyle(1001);
   fP.SetFillColorAlpha(41,0.8);
   fP.SetNpx(10000);
   fP.Draw("same");

   TF1 fNP("fNP",&funNonPrompt, 0., 0.5, 2);
   fNP.SetParameters(fMix.GetParameter(0),fMix.GetParameter(1));
   fNP.SetRange(fitRangeL,fitRangeH);
   fNP.SetLineColor(4);
   fNP.SetFillStyle(1001);
   fNP.SetFillColorAlpha(kBlue-6,0.8);
   fNP.SetNpx(10000);
   fNP.Draw("same");

   hD0DcaData->Draw("same");

   TLegend lgd(0.6, 0.7, 1, 1);
   lgd.AddEntry(hD0DcaData, "data", "lp");
   lgd.AddEntry(&fMix, "Combination of PD0 and NPD0", "lf");
   lgd.AddEntry(&fNP, "MC NPD0", "lf");
   lgd.AddEntry(&fP, "MC PD0", "lf");
   lgd.Draw();

   TLatex ltx;
   ltx.DrawLatexNDC(0.45, 0.65, Form("Prompt D^{0}: %f +/- %f", fMix.GetParameter(0), fMix.GetParError(0)));
   ltx.DrawLatexNDC(0.45, 0.55, Form("B2D: %f +/- %f", fMix.GetParameter(1), fMix.GetParError(1)));
   ltx.DrawLatexNDC(0.55, 0.45, cut.c_str());
   ltx.DrawLatexNDC(0.55, 0.35, label.c_str());
   ltx.DrawLatexNDC(0.3, 0.8, Form("chi2/ndf: %f/ %d", fitResult->Chi2(), fitResult->Ndf()));
   ltx.DrawLatexNDC(0.5, 0.2, Form("%.3f<=DCA<%.3fmm", fitRangeL*10, fitRangeH*10));

   auto pad2 = cFit.cd(2);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.16);
   pad2->SetLeftMargin(0.16);

   TH1D h("h", "", nuofDca, dcaBin);
   h.SetStats(0);
   h.SetMarkerStyle(20);
   h.GetXaxis()->SetRangeUser(fitRangeL, fitRangeH);
   //h.GetYaxis()->SetRangeUser(0, 5);
   h.GetYaxis()->SetRangeUser(0, 2);
   h.GetXaxis()->SetTitleSize(0.05);
   h.GetYaxis()->SetTitleSize(0.05);
   h.SetTitle(";DCA (cm);data/ fit");

   int dcaMinBin = h.GetXaxis()->FindBin(fitRangeL+0.00001*h.GetXaxis()->GetBinWidth(1));
   int dcaMaxBin = h.GetXaxis()->FindBin(fitRangeH-0.00001*h.GetXaxis()->GetBinWidth(1));
   for(int i=dcaMinBin; i<dcaMaxBin+1; i++){
      if(hD0DcaData->GetBinContent(i)==0) continue;
      h.SetBinContent(i, hD0DcaData->GetBinContent(i)/fMix.Eval(hD0DcaData->GetBinCenter(i)));
      h.SetBinError(i, hD0DcaData->GetBinError(i)/fMix.Eval(hD0DcaData->GetBinCenter(i)));
   }
   h.Draw("E P");
   TLine l(fitRangeL, 1, fitRangeH, 1);
   l.SetLineStyle(10);
   l.Draw("same");
   
   cFit.Print(picName.at("noLog").c_str());

   //hD0DcaData->GetYaxis()->SetRangeUser(hD0DcaData->GetMinimum()*0.5, hD0DcaData->GetMaximum()*10.);
   pad1->SetLogy();
   pad1->Update();
   cFit.Print(picName.at("log").c_str());
}

std::pair<double, double> 
dcaFittingLeastChi2_EachPtAndY(int mode, const float& pTMin, const float& pTMax, const float& yMin, const float& yMax, const float& fitRangeL, const float&fitRangeH)
{
   std::map<std::string, std::string> name;
   if(isUseMassFit) name["noLog"] = Form("plots/dcaFull/%s_pT%.1f-%.1f_y%.1f-%.1f_final_DcaFitting_DCA%.04f-%.04fcm_MassFit.png", 
         ana::treeName[mode].c_str(),pTMin,pTMax,yMin,yMax, fitRangeL, fitRangeH);
   if(isUseMassFit) name["log"] = Form("plots/dcaFull/%s_pT%.1f-%.1f_y%.1f-%.1f_fianl_DcaFitting_DCA%.04f-%.04fcm_MassFit_LogScale.png", 
         ana::treeName[mode].c_str(),pTMin,pTMax,yMin,yMax,  fitRangeL, fitRangeH);

   std::string cut = "";
   std::string label(TString::Format("%.1f<p_{T}<%.1fGeV, %.1f#leq|y|<%.1f",pTMin,pTMax,yMin,yMax));

   TFitResultPtr fitResult;

   dcaFractionFitting(fitRangeL, fitRangeH, fitResult, name, cut, label);

   auto pars = fitResult->GetParams();
   auto errs = fitResult->GetErrors();
   auto covMat = fitResult->GetCovarianceMatrix();

   double yieldNPD0 = pars[1];
   double yieldPD0 = pars[0];
   double yieldErrNPD0 = errs[1];
   double yieldErrPD0 = errs[0];

   double realFracNPD0 = yieldNPD0/(yieldNPD0+yieldPD0);

   double realFracNPD0Err = std::sqrt(std::pow(yieldErrNPD0/(yieldNPD0+yieldPD0)*(1-realFracNPD0), 2)
            + std::pow(yieldErrPD0/(yieldNPD0+yieldPD0)*realFracNPD0, 2)
            + 2*covMat[0][1]*yieldPD0/(yieldNPD0+yieldPD0)/(yieldNPD0+yieldPD0)*(-1)*realFracNPD0/(yieldNPD0+yieldPD0));

   std::cout << "error of non-prompt fraction: " << realFracNPD0Err << std::endl;
   std::cout << "non-prompt: " << realFracNPD0 << " +/- " << realFracNPD0Err << std::endl;

   return std::pair<double, double>(realFracNPD0, realFracNPD0Err);
}

void dcaFittingLeastChi2_npd0(int mode = 2, const float y_Max=2.0, const char* input="data/dcaHists_Proj_y0.0-2.0.root")
{
   gStyle->SetOptStat(0);
   TFile* f1= new TFile(input);
   TFile f2("data/fraction.root", "recreate");
//   for(int idca=0; idca<nDca; idca++){
   for(int idca=0; idca<3; idca++){
   //for(int idca=0; idca<1; idca++){
      TH1D hFractionNPD0(Form("hFractionNPD0_dca%d",idca), "", 3, ptbin);
      for(int ipt=0; ipt<3; ipt++){
         TH1D* hData = (TH1D*) f1->Get(Form("hDcaData_pt%d", ipt));
         hD0DcaData = (TH1D*) hData->Clone();
         hD0DcaData->SetName(Form("hDcaData_pt%d_dca%d", ipt, idca));

         TH1D* hD0DcaMCNP = (TH1D*) f1->Get(Form("hDCA_NPD0_pt%d", ipt));
         TH1D* hD0DcaMCP = (TH1D*) f1->Get(Form("hDCA_PD0_pt%d", ipt));

         hD0DcaMCNPSignal = (TH1D*) hD0DcaMCNP->Clone();
         hD0DcaMCNPSignal->SetName(Form("hDca_NPD0_pt%d_dca%d", ipt, idca));

         hD0DcaMCPSignal = (TH1D*) hD0DcaMCP->Clone();
         hD0DcaMCPSignal->SetName(Form("hDca_PD0_pt%d_dca%d", ipt, idca));

         int dcaCutBinSmall = hD0DcaData->GetXaxis()->FindBin(dcaCut[idca]+0.1*hD0DcaData->GetXaxis()->GetBinWidth(1));
         int dcaCutBinLarge = hD0DcaData->GetXaxis()->FindBin(dcaCut[idca+1]+0.1*hD0DcaData->GetXaxis()->GetBinWidth(1));
         //std::cout << dcaCutBin << std::endl;
         for(int idca=0; idca<dcaCutBinSmall; idca++){
            hD0DcaData->SetBinContent(idca, 0.);
            hD0DcaData->SetBinError(idca, 0.);

            hD0DcaMCNPSignal->SetBinContent(idca,0.);
            hD0DcaMCNPSignal->SetBinError(idca,0.);

            hD0DcaMCPSignal->SetBinContent(idca,0.);
            hD0DcaMCPSignal->SetBinError(idca,0.);
         }
         for(int idca=dcaCutBinLarge; idca<nuofDca+1; idca++){
            hD0DcaData->SetBinContent(idca, 0.);
            hD0DcaData->SetBinError(idca, 0.);

            hD0DcaMCNPSignal->SetBinContent(idca,0.);
            hD0DcaMCNPSignal->SetBinError(idca,0.);

            hD0DcaMCPSignal->SetBinContent(idca,0.);
            hD0DcaMCPSignal->SetBinError(idca,0.);
         }
         hD0DcaData->Scale(1./hD0DcaData->Integral("width"));
         hD0DcaMCNPSignal->Scale(1./hD0DcaMCNPSignal->Integral("width"));
         hD0DcaMCPSignal->Scale(1./hD0DcaMCPSignal->Integral("width"));

         const float yMin=0;
         const float yMax = y_Max;
         const float pTMin = ptbin[ipt];
         const float pTMax = ptbin[ipt+1];
         const float dcaMin = dcaCut[idca]+0.000000001;
         const float dcaMax = dcaCut[idca+1]-0.000000001;
         auto npfrac =dcaFittingLeastChi2_EachPtAndY(mode, pTMin, pTMax, yMin, yMax, dcaMin, dcaMax);

         hFractionNPD0.SetBinContent(ipt+1, npfrac.first);
         hFractionNPD0.SetBinError(ipt+1, npfrac.second);

         delete hD0DcaMCPSignal;
         delete hD0DcaMCNPSignal;
         delete hD0DcaData;
      }
      f2.cd();
      hFractionNPD0.Write();
   }
}

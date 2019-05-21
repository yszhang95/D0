#include "myAnaConsts.h"
#include "massfitting.C"
#include <iostream>
#include <memory>

void massHistsOptimal(int mode = 2)
//void dcaHistsOptimal(int mode = 1)
{
   TH1::SetDefaultSumw2();

   //TGaxis::SetMaxDigits(3);
   gStyle->SetOptStat(0);
   //TFile* f1 = new TFile(Form("%s_hists.root", ana::whichtree[mode].c_str()));
   std::unique_ptr<TFile> f1 = std::unique_ptr<TFile>(new TFile(Form("hists/%s_hists_pT%.1f-%.1f_y%.1f-%.1f.root", ana::whichtree[mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax)));
   std::map<std::string, TH3*> hDcaVsMassAndMvaPD0;
   std::map<std::string, TH3*> hDcaVsMassAndMvaNPD0;
   hDcaVsMassAndMvaPD0["h_match_unswap"] = (TH3D*) f1->Get("hDcaVsMassAndMvaPD0");
   hDcaVsMassAndMvaPD0["h_match_all"] = (TH3D*) f1->Get("hDcaVsMassAndMvaPD0_All");
   hDcaVsMassAndMvaNPD0["h_match_unswap"] = (TH3D*) f1->Get("hDcaVsMassAndMvaNPD0");
   hDcaVsMassAndMvaNPD0["h_match_all"] = (TH3D*) f1->Get("hDcaVsMassAndMvaNPD0_All");
   TH3D* hDcaVsMassAndMvaDataD0 = (TH3D*) f1->Get("hDcaVsMassAndMvaDataD0");

   TFile f2(Form("%s_dca_hists_pT%.1f-%.1f_y%.1f-%.1f.root", ana::whichtree[mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax), "recreate");

   // create Pic to store pictures
   const std::string dirPic = Form("if [ ! -d \"Pic\" ]; then\n"
                              "    mkdir %sPic \n"
                              "fi", ana::whichtree[mode].c_str());
   gSystem->Exec(dirPic.c_str());

   const std::string dirPicPtY = Form("if [ ! -d \"Pic/pT%.1f-%.1f_y%.1f-%.1f\" ]; then\n"
                              "    mkdir %sPic/pT%.1f-%.1f_y%.1f-%.1f \n"
                              "fi", ana::whichtree[mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax,
                                                                 ana::pTMin,ana::pTMax,ana::yMin,ana::yMax);
   gSystem->Exec(dirPicPtY.c_str());

   double significance[100];
   double yields_signal[100];
   double yields_bkg[100];
   double mvaCuts[100];
   for(int iMva=0; iMva<13; iMva++){
   //for(int iMva=2; iMva<4; iMva++){

      //create Pic/mva to store pictures
      std::string dirMvaPic(TString::Format( "if [ ! -d \"%sPic/pT%.1f-%.1f_y%.1f-%.1f/mva%d\" ]; then\n"
                                 "    mkdir %sPic/pT%.1f-%.1f_y%.1f-%.1f/mva%d \n"
                                 "fi", ana::whichtree[mode].c_str(), ana::pTMin,ana::pTMax,ana::yMin,ana::yMax, iMva, 
                                 ana::whichtree[mode].c_str(), ana::pTMin,ana::pTMax,ana::yMin,ana::yMax, iMva));
      gSystem->Exec(dirMvaPic.c_str());

      double mvaCut = (double)iMva*ana::mvaStep + ana::mvaMin;

      // extract invariant mass
      int mvaBinMin = hDcaVsMassAndMvaNPD0["h_match_all"]->GetYaxis()->FindBin(mvaCut+0.1*ana::mvaStep); // 0.1 offset to make sure the correct bin is returned
      int mvaBinMax = hDcaVsMassAndMvaNPD0["h_match_all"]->GetYaxis()->GetNbins()+1;
      int label = mvaCut*100;
      TH1D* hMassNPD0 = hDcaVsMassAndMvaNPD0["h_match_unswap"]->ProjectionX(Form("hMassMCNPD0mva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassNPD0All = hDcaVsMassAndMvaNPD0["h_match_all"]->ProjectionX(Form("hMassMCNPD0Allmva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassPD0 = hDcaVsMassAndMvaPD0["h_match_unswap"]->ProjectionX(Form("hMassMCPD0mva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassPD0All = hDcaVsMassAndMvaPD0["h_match_all"]->ProjectionX(Form("hMassMCPD0Allmva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassData = hDcaVsMassAndMvaDataD0->ProjectionX(Form("hMassDataD0mva%d", label), mvaBinMin, mvaBinMax);

      hMassNPD0->Scale(1.0/hMassNPD0->GetBinWidth(1));
      hMassNPD0All->Scale(1.0/hMassNPD0All->GetBinWidth(1));
      hMassPD0->Scale(1.0/hMassPD0->GetBinWidth(1));
      hMassPD0All->Scale(1.0/hMassPD0All->GetBinWidth(1));
      hMassData->Scale(1.0/hMassData->GetBinWidth(1));

      // TFitResultPtr constructed by std::shared_ptr<TFitResult>, be free of it even if you do not delete it
      // fitting the mass
      // rescale the mass
      
      TFitResultPtr fitResultPtr;
      TF1 fMass;
      if(mode == 0) fMass = massfitting(hMassData, hMassPD0, hMassPD0All, TString::Format("fMassMva%d", iMva), fitResultPtr);
      if(mode == 1) fMass = massfitting(hMassData, hMassNPD0, hMassNPD0All, TString::Format("fMassMva%d", iMva), fitResultPtr);
      if(mode == 2) fMass = massfitting(hMassData, hMassNPD0, hMassNPD0All, TString::Format("fMassMva%d", iMva), fitResultPtr);

      // draw the fitting
      drawMassFitting(hMassData, fMass, TString::Format("%sPic/pT%.1f-%.1f_y%.1f-%.1f/mva%d/MassFittingUnNormalized.png", ana::whichtree[mode].c_str(), ana::pTMin,ana::pTMax,ana::yMin,ana::yMax, iMva), 
            TString::Format("%.1f<p_{T}<%.1fGeV, %.1f<y<%.1f",ana::pTMin,ana::pTMax,ana::yMin,ana::yMax), TString::Format("MVA > %.2f", 0.4+0.02*iMva));

      // https://root.cern/doc/v616/TF1Helper_8cxx_source.html and https://root.cern.ch/doc/master/classTF1.html
      // help you understand how the integral error of TF1 is calculated
      // 1-dimensional function have five parameters, x_low, x_min, array of pars, array of covariance matrix, precision
      // array of pars would replace the original set of pars by last fitting unless no set of pars is passed
      // it is very important to pass in the correct set of pars and covariance matrix
      // explicitly pass the parameters to make sure anything won't go wrong
      // be aware that the fixed parameter would bring effect on the evaluation
      
      // get the covariance matrix
      auto covMat = fitResultPtr->GetCovarianceMatrix();

      // define the background function and then fix or set the parameters
      TF1 bkg("bkg", "[9] + [10]*x + [11]*x*x + [12]*x*x*x"
            "+ 0 *[0]*[1]*[2]*[3]*[4]*[5]*[6]*[7]*[8]", 1.7, 2.0);
      for(int ipar=0; ipar<8+1; ipar++){
         bkg.FixParameter(ipar, fMass.GetParameter(ipar));
      }
      for(int ipar=9; ipar<12+1; ipar++){
         bkg.SetParameter(ipar, fMass.GetParameter(ipar));
      }

      // define the signal function and then fix or set the parameters
      TF1 signal("signal", "[0]*[5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+0 *[7]*[8]*[9]*[10]*[11]*[12]");
      for(int ipar=7; ipar<12+1; ipar++){
         signal.FixParameter(ipar, fMass.GetParameter(ipar));
      }
      for(int ipar=0; ipar<6+1; ipar++){
         signal.SetParameter(ipar, fMass.GetParameter(ipar));
      }

      // calculate the yield significance 
      std::map<std::string, double> yield_signal;
      std::map<std::string, double> yieldError_signal;
      std::map<std::string, double> yield_bkg;
      std::map<std::string, double> yieldError_bkg;

/*      yield_signal["peak"] = signal.Integral(ana::peak_min, ana::peak_max);
      yieldError_signal["peak"] = signal.IntegralError(ana::peak_min, ana::peak_max, 0, covMat.GetMatrixArray());
      yield_bkg["peak"] = bkg.Integral(ana::peak_min, ana::peak_max);
      yieldError_bkg["peak"] = bkg.IntegralError(ana::peak_min, ana::peak_max, 0, covMat.GetMatrixArray());

      yields_signal[iMva] = yield_signal["peak"]; 
      yields_bkg[iMva] = yield_bkg["peak"];
      significance[iMva] = yield_signal["peak"] / sqrt(yield_signal["peak"] + yield_bkg["peak"]);
*/

      significance[iMva] = fMass.GetParameter(0)/fMass.GetParError(0);
      yields_signal[iMva] = fMass.GetParameter(0);
      mvaCuts[iMva] = mvaCut;

      f2.cd();
      fMass.Write();

      hMassNPD0->Write();
      hMassNPD0->Delete();
      hMassNPD0All->Write();
      hMassNPD0All->Delete();
      hMassPD0->Write();
      hMassPD0->Delete();
      hMassPD0All->Write();
      hMassPD0All->Delete();
      hMassData->Write();
      hMassData->Delete();
   }

   TCanvas* cc = new TCanvas("cc","cc",600,600);
   TGraph* gr = new TGraph(12, mvaCuts, significance);
   gr->SetMarkerStyle(20);
   gr->Draw("AP");

   TCanvas* cc1 = new TCanvas("cc1","cc1",600,600);
   TGraph* gr1 = new TGraph(12, mvaCuts, yields_signal);
   gr1->SetMarkerStyle(20);
   gr1->Draw("AP");

   gr1->Write("signalvsmva");
   gr->Write("significancevsmva");
}

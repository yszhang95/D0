#include "include/myAnaConsts.h"
#include "macros/massfitting.C"

const int nuofDca = 37;
const double dcaBin[nuofDca+1] = {0, 0.0005, 0.001,0.0015, 0.002, 
   0.0025, 0.003, 0.0035, 0.004, 0.0045, 
   0.005, 0.0055, 0.006, 0.0065, 0.007, 
   0.0075,  0.008, 0.0085, 0.009, 0.0095, 
   0.010, 0.0105, 0.011, 0.0115, 0.012, 
   0.0125, 0.013, 0.0135, 0.014, 0.016, 
   0.022, 0.026, 0.030, 0.036, 0.042, 
   0.05, 0.062, 0.08}; // cm

TH1D* hDcaMCAll(TH3* hDcaVsMassAndMva, TH1& hDcaMC, const double& mvaCut, const  std::string& proj_name);

void hDcaDataSignal(TH3* hData, std::map<std::string, TH3*> hNP, std::map<std::string, TH3*> hP, 
      TH1& hDataD0, const double& mvaCut, const std::string& tmpName, const char* label, bool isPrompt);

void dcaHists_Proj(const char* input_data="data/corr2D_npd0ana1_dcaFull_y2.0.root", const char* prefix = "npd0ana1", const float yMax=2.0, const bool isPrompt=false)
{
   gStyle->SetOptStat(0);
   TFile* f1 = new TFile(input_data);
   TFile* f2[ana::nPt];
   TFile fout(Form("data/dcaHists_Proj_y%.1f-%.1f.root", 0., yMax), "recreate");
   TH3D* hDcaVsMassAndMva_Data[ana::nPt];
   std::map<std::string, TH3*> hDcaVsMassAndMva_MC_NPD0[ana::nPt];
   std::map<std::string, TH3*> hDcaVsMassAndMva_MC_PD0[ana::nPt];
   for(int ipt=0; ipt<ana::nPt; ipt++){
      hDcaVsMassAndMva_Data[ipt] = (TH3D*) f1->Get(Form("hDcaVsMassAndMva_pt%d", ipt));
      hDcaVsMassAndMva_Data[ipt]->GetEntries();
      f2[ipt] = new TFile(Form("MC/%s_hists_Dca3D_pT%.1f-%.1f_y%.1f-%.1f.root", prefix, ana::ptbin[ipt], ana::ptbin[ipt+1], 0., yMax));

      (hDcaVsMassAndMva_MC_NPD0[ipt])["h_match_unswap"] = (TH3D*)f2[ipt]->Get(Form("hDcaVsMassAndMvaNPD0_pt%d_y0",ipt));
      (hDcaVsMassAndMva_MC_NPD0[ipt])["h_match_all"] = (TH3D*)f2[ipt]->Get(Form("hDcaVsMassAndMvaNPD0_All_pt%d_y0",ipt));

      (hDcaVsMassAndMva_MC_PD0[ipt])["h_match_unswap"] = (TH3D*)f2[ipt]->Get(Form("hDcaVsMassAndMvaPD0_pt%d_y0",ipt));
      (hDcaVsMassAndMva_MC_PD0[ipt])["h_match_all"] = (TH3D*)f2[ipt]->Get(Form("hDcaVsMassAndMvaPD0_All_pt%d_y0",ipt));

      TH1D hDcaMCNP(Form("hDCA_NPD0_pt%d", ipt), "", nuofDca, dcaBin);
      TH1D hDcaMCP(Form("hDCA_PD0_pt%d", ipt), "", nuofDca, dcaBin);
      double mvaCut;
      if(isPrompt){ 
         mvaCut = ana::mvaCut_PD0[ipt];
      } else {
         mvaCut = ana::mvaCut_NPD0[ipt];
      }
      TH1D* hDcaMCNP_proj = hDcaMCAll(hDcaVsMassAndMva_MC_NPD0[ipt].at("h_match_all"), hDcaMCNP, mvaCut, "tempNPD0");
      TH1D* hDcaMCP_proj = hDcaMCAll(hDcaVsMassAndMva_MC_PD0[ipt].at("h_match_all"), hDcaMCP, mvaCut, "tempPD0");

      TH1D hDcaData(Form("hDcaData_pt%d", ipt), "", nuofDca, dcaBin);
      hDcaDataSignal(hDcaVsMassAndMva_Data[ipt], hDcaVsMassAndMva_MC_NPD0[ipt], hDcaVsMassAndMva_MC_PD0[ipt], 
            hDcaData, mvaCut, Form("plots/dcaFull/y%.1f/massPerBin_pt%d_y%.1f-%.1f_",yMax, ipt, 0., yMax), TString::Format("%.1f<p_{T}<%.1f |y|<%.1f", ana::ptbin[ipt], ana::ptbin[ipt+1], yMax), false);
      hDcaData.Scale(1./hDcaData.Integral("width"));

      delete hDcaMCNP_proj;
      delete hDcaMCP_proj;
      delete f2[ipt];
      fout.cd();
      hDcaMCNP.Write();
      hDcaMCP.Write();
      hDcaData.Write();
   }
}
TH1D* hDcaMCAll(TH3* hDcaVsMassAndMva, TH1& hDcaMC, const double& mvaCut, const  std::string& proj_name)
{
   int mvaBinMin = hDcaVsMassAndMva->GetYaxis()->FindBin(mvaCut+0.1*hDcaVsMassAndMva->GetYaxis()->GetBinWidth(1)); // 0.1 offset to make sure the correct bin is returned
   int mvaBinMax = hDcaVsMassAndMva->GetYaxis()->GetNbins()+1;
   TH1D* hDcaMCProj = hDcaVsMassAndMva->ProjectionZ(proj_name.c_str(), 1, 60, mvaBinMin, mvaBinMax); // x: mass, y: dca 
   const double norm = hDcaMCProj->Integral(0, hDcaMCProj->GetXaxis()->GetNbins()+1); // normalized constant
   cout << hDcaMCProj->GetNbinsX()<< endl;
   hDcaMCProj->Scale(1./ (norm * hDcaMCProj->GetBinWidth(1))); // convert histogram to the P.D.F.
   for(int iDca=0; iDca<nuofDca; iDca++){
      int binlw = hDcaMCProj->FindBin(dcaBin[iDca] + hDcaMCProj->GetBinWidth(1)/2.); // binwidth/2 offset to make sure the correct bin is returned
      int binup = hDcaMCProj->FindBin(dcaBin[iDca+1] - hDcaMCProj->GetBinWidth(1)/2.);  // - binwidth/2 offset due to the integral is taken within [binlw, binup]
      // calculate the contents and errors per dca bin
      double binwidth = dcaBin[iDca+1] - dcaBin[iDca];
      double err = 0.;
      double yield = hDcaMCProj->IntegralAndError(binlw, binup, err, "width"); // note the option "width"
      hDcaMC.SetBinContent(iDca+1, yield/binwidth);
      hDcaMC.SetBinError(iDca+1, err/binwidth);
   }
   return hDcaMCProj;
}

void hDcaDataSignal(TH3* hData, std::map<std::string, TH3*> hNP, std::map<std::string, TH3*> hP, 
      TH1& hDataD0, const double& mvaCut, const std::string& tmpName, const char* label, bool isPrompt)
{
   int mvaBinMin = hP["h_match_unswap"]->GetYaxis()->FindBin(
         mvaCut+0.1*hP["h_match_unswap"]->GetYaxis()->GetBinWidth(1)); // 0.1 offset to make sure the correct bin is returned
   int mvaBinMax = hP["h_match_unswap"]->GetYaxis()->GetNbins()+1; // overflow of y axis, mva
   int mcDcaBinMin = 0; // underflow of the z axis, the DCA
   int mcDcaBinMax = hP["h_match_unswap"]->GetZaxis()->GetNbins()+1; // overflow of z axis, the DCA

   // x: mass, y: mva, z: dca
   TH1D* hMCPromptMass = hP["h_match_unswap"]->ProjectionX("hPromptMass", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax); // projection, y-mva:[cut, overflow], z-DCA:[underflow, overflow]
   TH1D* hMCPromptMassAll = hP["h_match_all"]->ProjectionX("hPromptMassAll", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax);
   TH1D* hMCNonPromptMass = hNP["h_match_unswap"]->ProjectionX("hNonPromptMass", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax);
   TH1D* hMCNonPromptMassAll = hNP["h_match_all"]->ProjectionX("hNonPromptMassAll", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax);

   TH1D* hTemp;
   TCanvas cTemp("cTemp", "", 550, 450);
   for(int iDca=0; iDca<nuofDca; iDca++){
      int binlw = hData->GetZaxis()->FindBin(dcaBin[iDca] + hData->GetZaxis()->GetBinWidth(1) * 0.2); // 0.2 offset to make sure the correct bin is returned
      int binup = hData->GetZaxis()->FindBin(dcaBin[iDca+1] - hData->GetZaxis()->GetBinWidth(1) * 0.2); // -0.2 offset due to the integral is taken within [binlw, binup]
      hTemp = hData->ProjectionX("tmp", mvaBinMin, mvaBinMax, binlw, binup); // projection, y-mva:[cut, overflow], z-DCA:[dca[iDca], dca[iDca+1]]
      //hTemp->Scale(1./hTemp->GetBinWidth(1));
      hTemp->GetYaxis()->SetRangeUser(0, hTemp->GetMaximum() * 1.3);
      hTemp->SetTitle(";Mass (GeV);Entries /(5 MeV)");
      TFitResultPtr fitResultPtr;
      TF1 f;
      if(isPrompt) {
         f = massfitting(hTemp, hMCPromptMass, hMCPromptMassAll, TString::Format("iDca%d", iDca), fitResultPtr);
      } else {
         f = massfitting(hTemp, hMCNonPromptMass, hMCNonPromptMassAll, TString::Format("iDca%d", iDca), fitResultPtr);
      }
      drawMassFitting(hTemp, f, TString::Format("%s%d.png", tmpName.c_str(), iDca), false, label,
            TString::Format("%2.f < DCA <%2.fmm", dcaBin[iDca]*1e3, dcaBin[iDca+1]*1e3));

      double binWidth = dcaBin[iDca+1] - dcaBin[iDca];
      hDataD0.SetBinContent(iDca+1, f.GetParameter(0)/binWidth/hTemp->GetBinWidth(1));
      hDataD0.SetBinError(iDca+1, f.GetParError(0)/binWidth/hTemp->GetBinWidth(1));;

      // clear pointers
      hTemp->Delete();
   }
   hMCPromptMass->Delete(); 
   hMCPromptMass = nullptr;
   hMCPromptMassAll->Delete();
   hMCPromptMassAll = nullptr;
   hMCNonPromptMass->Delete();
   hMCNonPromptMass = nullptr;
   hMCNonPromptMassAll->Delete();
   hMCNonPromptMassAll = nullptr;
}

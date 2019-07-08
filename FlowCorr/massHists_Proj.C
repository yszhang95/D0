#include "include/myAnaConsts.h"

const int nDca = 11;
const double dcaBin[nDca] = {0., 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013, 0.014};

void massHists_Proj(const char* output="MC/hMC_mass_npd0ana1_y2.0.root", const float yMax=2.0)
{
   TFile* f1[ana::nPt];
   TH3D* hDcaVsMassAndMva[ana::nPt];
   TH3D* hDcaVsMassAndMva_All[ana::nPt];
   for(int ipt=0; ipt<ana::nPt; ipt++){
      f1[ipt] = new TFile(Form("MC/npd0ana1_hists_Dca3D_pT%.1f-%.1f_y%.1f-%.1f.root", ana::ptbin[ipt], ana::ptbin[ipt+1], 0., yMax));
      hDcaVsMassAndMva[ipt] = (TH3D*)f1[ipt]->Get(Form("hDcaVsMassAndMvaNPD0_pt%d_y0", ipt));
      hDcaVsMassAndMva_All[ipt] = (TH3D*)f1[ipt]->Get(Form("hDcaVsMassAndMvaNPD0_All_pt%d_y0", ipt));
   }
   TH1D* hMass_DCA[ana::nPt][nDca];
   TH1D* hMass_DCA_all[ana::nPt][nDca];
   int mvaBinMin, mvaBinMax, dcaBinMin, dcaBinMax;
   for(int ipt=0; ipt<ana::nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         mvaBinMin = hDcaVsMassAndMva[ipt]->GetYaxis()->FindBin(ana::mvaCut_NPD0[ipt]+0.1*hDcaVsMassAndMva[ipt]->GetXaxis()->GetBinWidth(1));
         mvaBinMax = hDcaVsMassAndMva[ipt]->GetYaxis()->GetNbins()+1;
         dcaBinMin = hDcaVsMassAndMva[ipt]->GetZaxis()->FindBin(dcaBin[idca]+0.1*hDcaVsMassAndMva[ipt]->GetZaxis()->GetBinWidth(1));
         dcaBinMax = hDcaVsMassAndMva[ipt]->GetZaxis()->GetNbins()+1;
         hMass_DCA[ipt][idca] = hDcaVsMassAndMva[ipt]->ProjectionX(Form("hMC_mass_pt%d_dca%d", ipt, idca), 
               mvaBinMin, mvaBinMax, dcaBinMin, dcaBinMax);
         hMass_DCA_all[ipt][idca] = hDcaVsMassAndMva_All[ipt]->ProjectionX(Form("hMC_mass_all_pt%d_dca%d", ipt, idca), 
               mvaBinMin, mvaBinMax, dcaBinMin, dcaBinMax);
      }
   }
   TFile f2(output, "recreate");
   for(int ipt=0; ipt<ana::nPt; ipt++){
      for(int idca=0; idca<nDca; idca++){
         hMass_DCA[ipt][idca]->Write();
         hMass_DCA_all[ipt][idca]->Write();
      }
   }
}

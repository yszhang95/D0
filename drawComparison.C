#include "myAnaConsts.h"
void drawComparison(int mode = 0)
{
   TGaxis::SetMaxDigits(3);
   gStyle->SetOptStat(0);
   //TFile* f1 = new TFile("d0ana_hists.root");
   TFile* f1 = new TFile(Form("%s_hists.root", ana::whichtree[mode].c_str()));
   TH1D* hDca3D_pd0 = (TH1D*) f1->Get("hDca3D_pd0");
   TH1D* hDca3D_npd0 = (TH1D*) f1->Get("hDca3D_npd0");
   TH1D* hDL_pd0 = (TH1D*) f1->Get("hDL_pd0");
   TH1D* hDL_npd0 = (TH1D*) f1->Get("hDL_npd0");
   TCanvas* c1 = new TCanvas("c1", "", 550, 450);
   c1->SetBottomMargin(0.16);
   c1->SetLogy();
   hDca3D_pd0->SetTitle(";3D DCA (cm);counts");
   hDca3D_pd0->SetLineColor(kRed);
   hDca3D_npd0->SetLineColor(kBlue);
   hDca3D_pd0->Draw();
   hDca3D_npd0->Draw("same");
   TLegend* lgdDca = new TLegend(0.6, 0.75, 0.9, 0.9);
   lgdDca->AddEntry(hDca3D_pd0, "prompt D0", "pl");
   lgdDca->AddEntry(hDca3D_npd0, "non-prompt D0", "pl");
   lgdDca->Draw();

   TCanvas* c2 = new TCanvas("c2", "", 550, 450);
   c2->SetBottomMargin(0.16);
   c2->SetLogy();
   hDL_pd0->SetTitle(";pseudo decay length (cm);counts");
   hDL_pd0->SetLineColor(kRed);
   hDL_npd0->SetLineColor(kBlue);
   hDL_pd0->Draw();
   hDL_npd0->Draw("same");
   TLegend* lgdDL = new TLegend(0.6, 0.75, 0.9, 0.9);
   lgdDL->AddEntry(hDL_pd0, "prompt D0", "pl");
   lgdDL->AddEntry(hDL_npd0, "non-prompt D0", "pl");
   lgdDL->Draw();
}

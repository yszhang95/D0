#include "TStyle.h"
string files[] = {
"../data/Jets_lowvspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0",
"../data/Jetsvspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0",
"../data/V2_lowvspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0",
"../data/V2vspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0",
"../data/nass_npd0_pT2.0-8.0_y-1.0-1.0_new_binning.root",
"../data/nass_npd0_pT2.0-8.0_y-1.0-1.0_new_binning.root"
   };
   string graphs[] = {
      "Jetsplot",
      "Jetsplot",
      "V2plot",
      "V2plot",
      "LM",
         "HM"
   };
   string outputs[] = {
      "Jets_low",
      "Jets",
      "V2_low",
      "V2",
      "Nass_low",
      "Nass"
   };

double ptbin[] = {2., 5., 8.};

/*
//  old, wrong template
double smallDca[2] ={ // nPt
0.079946,
0.0500907
};
double largeDca[2] ={ // nPt
0.686508,
0.769427
};
double smallDca_e[2] = { // nPt
0.000958518,
0.0106871
};
double largeDca_e[2] = { // nPt
0.0207466,
0.00681814
};
*/

double smallDca[2] ={ // nPt
   0.105, 0.110
};
double smallDca_e[2] = {0.006, 0.005};
double largeDca[2] = {0.624, 0.735};
double largeDca_e[2] = {0.03, 0.03};

/*
// fit range 0.006 - 0.08 // std
double smallDca[2] ={ // nPt
0.098428,
0.103564 
};
double largeDca[2] ={ // nPt
0.577233,
0.691316
};
double smallDca_e[2] = { // nPt
0.00470071,
0.00453739
};
double largeDca_e[2] = { // nPt
0.0256403,
0.0242854
};
*/

/*
// fit range 0.004 - 0.08 // syst errors
double smallDca[2] ={ // nPt
   0.0898474,
   0.0933014
};
double largeDca[2] ={ // nPt
   0.552491,
   0.666089
};
double smallDca_e[2] = { // nPt
   0.000981349,
   0.00353246
};
double largeDca_e[2] = { // nPt
   0.00296709,
   0.0199069
};
*/

/*
// fit range 0.008 - 0.08 // syst errors
double smallDca[2] ={ // nPt
   0.102076,
   0.112225
};
double largeDca[2] ={ // nPt
   0.587071,
   0.71019
};
double smallDca_e[2] = { // nPt
   0.00564512,
   0.00690441
};
double largeDca_e[2] = { // nPt
   0.0315394,
   0.0351376
};
*/

TGraphErrors* fitNPD0_Each(
      const char* large,
      const char* small,
      const char* g_large,
      const char* g_small,
      string outputs
      )
{
   TFile* flarge = new TFile(large);
   TFile* fsmall = new TFile(small);
   TGraphErrors* glarge;
   TGraphErrors* gsmall;
   flarge->GetObject(g_large, glarge);
   fsmall->GetObject(g_small, gsmall);
   TGraphErrors* var[2]; // nPt
   for(int i=0; i<2; i++){// npt
      var[i] = new TGraphErrors(2);
      var[i]->SetPoint(0, smallDca[i], gsmall->GetY()[i]);
      var[i]->SetPointError(0, smallDca_e[i], gsmall->GetEY()[i]);
      var[i]->SetPoint(1, largeDca[i], glarge->GetY()[i]);
      var[i]->SetPointError(1, largeDca_e[i], glarge->GetEY()[i]);
   }
   TGraphErrors* ret = new TGraphErrors(2);
   TCanvas* c[2];
   TF1* fl[2];
   for(int i=0; i<2; i++){
      //TGaxis::SetMaxDigits(2);
      setTDRStyle();
      gStyle->SetOptFit(0);
      ///gStyle->SetPadTickX(1);
      //gStyle->SetPadTickY(1);
      c[i] = new TCanvas(Form("c_%s_%d",outputs.c_str(), i), Form("pt%d", i), 500, 500);
      c[i]->SetLeftMargin(0.2);
      //c[i]->SetLeftMargin(0.16);
      fl[i] = new TF1(Form("f1_%d", i), "[0]*x+(1-x)*[1]", 0, 1);
      fl[i]->SetParameters(var[i]->GetY()[1], var[i]->GetY()[0]);
      var[i]->Fit(fl[i], "QRE");
      var[i]->Fit(fl[i], "QRE");
      var[i]->Fit(fl[i], "QRE");
      var[i]->Fit(fl[i], "RE");

      auto h = var[i]->GetHistogram();
      h->SetTitle(Form(";Nonprompt D^{0} Fraction;%s", outputs.c_str()));
      h->GetXaxis()->SetLimits(0, 1);
      h->SetMinimum(0.0000001);
      h->GetXaxis()->CenterTitle();
      h->GetYaxis()->CenterTitle();
      double ymax = max(fl[i]->GetParameter(0), fl[i]->GetParameter(1));
      cout << gsmall->GetY()[i] << endl;
      h->SetMaximum(ymax*2);

      if(outputs == "V2") h->SetMaximum(0.016);
      if(outputs == "V2_low") h->SetMaximum(0.03);

      h->GetYaxis()->SetTitleOffset(1.5);
      h->Draw();
      var[i]->SetMarkerStyle(20);
      var[i]->SetMarkerSize(1.2);
      var[i]->Draw("E P A");
      ret->SetPoint(i, i+1, fl[i]->GetParameter(0));
      ret->SetPointError(i, 0, fl[i]->GetParError(0));
      TLatex *ltx = new TLatex;
      ltx->SetTextSize(0.04);
      ltx->DrawLatexNDC(0.27, 0.35, Form("%.1f #leq p_{T} < %.1fGeV", ptbin[i], ptbin[i+1]));
      ltx->DrawLatexNDC(0.33, 0.27, "|y|<1");
      //ltx->DrawLatexNDC(0.2, 0.88, "#scale[1.5]{CMS}");
      ltx->DrawLatexNDC(0.2, 0.88, "#scale[1.5]{CMS} #it{#scale[1.4]{Preliminary}}");
      ltx->DrawLatexNDC(0.67, 0.88, "#scale[1.3]{pPb 8.16TeV}");
      if(outputs == "Jets_low") ltx->DrawLatexNDC(0.3, 0.2, "N_{trk}^{offline} < 35");
      if(outputs == "Jets") ltx->DrawLatexNDC(0.3, 0.2, " 185#leq N_{trk}^{offline} < 250");
      if(outputs == "V2_low") ltx->DrawLatexNDC(0.3, 0.2, "N_{trk}^{offline} < 35");
      if(outputs == "V2") ltx->DrawLatexNDC(0.3, 0.2, " 185#leq N_{trk}^{offline} < 250");
      if(outputs == "Nass_low") ltx->DrawLatexNDC(0.3, 0.2, "N_{trk}^{offline} < 35");
      if(outputs == "Nass") ltx->DrawLatexNDC(0.3, 0.2, " 185#leq N_{trk}^{offline} < 250");
      if(outputs == "Jets_low") h->GetYaxis()->SetTitle("Y_{jets}");
      if(outputs == "Jets") h->GetYaxis()->SetTitle("Y_{jets}");
      if(outputs == "V2_low") h->GetYaxis()->SetTitle("V_{2}");
      if(outputs == "V2") h->GetYaxis()->SetTitle("V_{2}");
      if(outputs == "Nass_low") h->GetYaxis()->SetTitle("N_{assoc}");
      if(outputs == "Nass") h->GetYaxis()->SetTitle("N_{assoc}");
      if(outputs == "Jets_low") {
         //ltx->DrawLatexNDC(0.2, 0.73, Form("Y_{jets}^{signal} = Y_{jets}^{nonprompt D^{0}}*Frac. + Y_{jets}^{prompt D^{0}}*(1-Frac.)"));
         //ltx->DrawLatexNDC(0.2, 0.65, Form("Y_{jets}^{nonprompt D^{0}} = %g+/-%g", fl[i]->GetParameter(0), fl[i]->GetParError(0) ));
      }
      if(outputs == "Jets") {
         //ltx->DrawLatexNDC(0.2, 0.73, Form("Y_{jets}^{signal} = Y_{jets}^{nonprompt D^{0}}*Frac. + Y_{jets}^{prompt D^{0}}*(1-Frac.)"));
         //ltx->DrawLatexNDC(0.2, 0.65, Form("Y_{jets}^{nonprompt D^{0}} = %g+/-%g", fl[i]->GetParameter(0), fl[i]->GetParError(0) ));
      }
      if(outputs == "V2_low") {
         //ltx->DrawLatexNDC(0.2, 0.73, Form("V_{2}^{signal} = V_{2}^{nonprompt D^{0}}*Frac. + V_{2}^{prompt D^{0}}*(1-Frac.)"));
         //ltx->DrawLatexNDC(0.2, 0.65, Form("V_{2}^{nonprompt D^{0}} = %g+/-%g", fl[i]->GetParameter(0), fl[i]->GetParError(0) ));
      }
      if(outputs == "V2") {
         //ltx->DrawLatexNDC(0.2, 0.73, Form("V_{2}^{signal} = V_{2}^{nonprompt D^{0}}*Frac. + V_{2}^{prompt D^{0}}*(1-Frac.)"));
         //ltx->DrawLatexNDC(0.2, 0.65, Form("V_{2}^{nonprompt D^{0}} = %g+/-%g", fl[i]->GetParameter(0), fl[i]->GetParError(0) ));
      }
      if(outputs == "Nass_low"){
         //ltx->DrawLatexNDC(0.2, 0.73, Form("N_{assoc}^{signal} = N_{assoc}^{nonprompt D^{0}}*Frac. + N_{assoc}^{prompt D^{0}}*(1-Frac.)"));
         //ltx->DrawLatexNDC(0.2, 0.65, Form("N_{assoc}^{nonprompt D^{0}} = %g+/-%g", fl[i]->GetParameter(0), fl[i]->GetParError(0) ));
      }
      if(outputs == "Nass"){
         //ltx->DrawLatexNDC(0.2, 0.73, Form("N_{assoc}^{signal} = N_{assoc}^{nonprompt D^{0}}*Frac. + N_{assoc}^{prompt D^{0}}*(1-Frac.)"));
         //ltx->DrawLatexNDC(0.2, 0.65, Form("N_{assoc}^{nonprompt D^{0}} = %g+/-%g", fl[i]->GetParameter(0), fl[i]->GetParError(0) ));
      }
      c[i]->Print(Form("c_%s_%d.pdf",outputs.c_str(), i));
   }
   return ret;
}

void fitNPD0(){
   TGraphErrors* ret[6];
   TFile ofile("ofile_npd0.root", "recreate");
   for(int i=0; i<6; i++){
   //for(int i=2; i<3; i++){
   //for(int i=3; i<4; i++){
      if(i<4) {
      ret[i] = fitNPD0_Each(Form("%s_dca1_new_binning_raw.root", files[i].c_str()),
            Form("%s_dca0_new_binning_raw.root", files[i].c_str()),
            Form("%s_dca1", graphs[i].c_str()),
            Form("%s_dca0", graphs[i].c_str()),
            outputs[i].c_str()
         );
      }
      else {
      ret[i] = fitNPD0_Each(files[i].c_str(),
            files[i].c_str(),
            Form("%s_dca1", graphs[i].c_str()),
            Form("%s_dca0", graphs[i].c_str()),
            outputs[i].c_str()
            );

      }
      ofile.cd();
      ret[i]->Write(outputs[i].c_str());
   }
}

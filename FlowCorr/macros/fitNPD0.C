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
      c[i] = new TCanvas(Form("c_%s_%d",outputs.c_str(), i), Form("pt%d", i), 500, 500);
      c[i]->SetLeftMargin(0.16);
      fl[i] = new TF1(Form("f1_%d", i), "[0]*x+(1-x)*[1]", 0, 1);
      fl[i]->SetParameters(var[i]->GetY()[1], var[i]->GetY()[0]);
      var[i]->Fit(fl[i], "QRE");
      var[i]->Fit(fl[i], "QRE");
      var[i]->Fit(fl[i], "QRE");
      var[i]->Fit(fl[i], "RE");

      auto h = var[i]->GetHistogram();
      h->SetTitle(Form(";NPD0 fraction;%s", outputs.c_str()));
      h->GetXaxis()->SetLimits(0, 1);
      h->SetMinimum(0.0000001);
      double ymax = max(fl[i]->GetParameter(0), fl[i]->GetParameter(1));
      cout << gsmall->GetY()[i] << endl;
      h->SetMaximum(ymax*2);
      h->Draw();
      var[i]->Draw("E P A");
      ret->SetPoint(i, i+1, fl[i]->GetParameter(0));
      ret->SetPointError(i, 0, fl[i]->GetParError(0));
      TLatex *ltx = new TLatex;
      ltx->DrawLatexNDC(0.3, 0.7, Form("%.1f<pT<%.1f", ptbin[i], ptbin[i+1]));
      c[i]->Print(Form("c_%s_%d.pdf",outputs.c_str(), i));
   }
   return ret;
}

void fitNPD0(){
   TGraphErrors* ret[4];
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

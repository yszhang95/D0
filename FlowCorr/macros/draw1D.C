#ifndef __MY_INCLUDE__
#define __MY_INCLUDE__
R__ADD_INCLUDE_PATH(../include);
R__ADD_INCLUDE_PATH(../src);
#include "myAnaConsts.h"
#include "functions.cxx"
#endif

using namespace std;

typedef TH1* TH1ptr;

bool isMB;
bool isHM1to6;

Double_t fpoly2(Double_t *x, Double_t *par)
{
//    if (reject && x[0] < 0.6 && x[0] > -0.6) {
   if (x[0] < 0.6 && x[0] > -0.6) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*fabs(x[0]) + par[2]*x[0]*x[0];
}

void proj1D_longrange(TH2* h2DSignal, TH2D* h2DBackground, TH1ptr& h, const char* name)
{
   const double deltaEtaBound = 1.;
   int negBinMin = 0;
   int negBinMax = -1;
   int posBinMin = h2DSignal->GetXaxis()->GetNbins()+1;
   if(isMB) negBinMax = h2DSignal->GetXaxis()->FindBin(-1.* deltaEtaBound) - 1;
   if(isMB) posBinMin = h2DSignal->GetXaxis()->FindBin(1.* deltaEtaBound) +1 ;
   if(!isMB) negBinMax = h2DSignal->GetXaxis()->FindBin(-1.* deltaEtaBound);
   if(!isMB) posBinMin = h2DSignal->GetXaxis()->FindBin(1.* deltaEtaBound);
   int posBinMax = h2DSignal->GetXaxis()->GetNbins()+1;
   TH1D* hNeg = h2DSignal->ProjectionY("hneg", negBinMin, negBinMax);
   TH1D* hPos = h2DSignal->ProjectionY("hpos", posBinMin, posBinMax);
   hNeg->Add(hPos);

   TH1D* temp_neg = h2DBackground->ProjectionY("hneg_bkg", negBinMin, negBinMax);
   TH1D* temp_pos = h2DBackground->ProjectionY("hpos_bkg", posBinMin, posBinMax);
   temp_neg->Add(temp_pos);

   int center = h2DBackground->FindBin(0., 0.);
   hNeg->Divide(temp_neg);
   hNeg->Scale(h2DBackground->GetBinContent(center) / temp_neg->GetBinWidth(1)/h2DBackground->GetXaxis()->GetBinWidth(1));

   delete hPos;
   delete temp_neg;
   delete temp_pos;

   hNeg->SetName(name);
   h = hNeg;
}

void proj1D_sr(TH2* h2DSignal, TH2D* h2DBackground, TH1ptr& h, const char* name)
{
   //h = new TH1D(name, ";#Delta#phi;dN/d(#Delta#phi)", ana::nPhiBin+2, -1*ana::PI-1./16*ana::PI, ana::PI + 1./16*ana::PI);
   h = new TH1D(name, ";#Delta#phi;dN/d(#Delta#phi)", ana::nPhiBin, -1*ana::PI, ana::PI);

   int negBin = h2DSignal->GetXaxis()->FindBin(-1.);
   int posBin = h2DSignal->GetXaxis()->FindBin(1.);
   TH1D* hsig = h2DSignal->ProjectionY("hsig", negBin, posBin);

   TH1D* temp = h2DBackground->ProjectionY("hbkg", negBin, posBin);

   int center = h2DBackground->FindBin(0., 0.);
   hsig->Divide(temp);
   hsig->Scale(h2DBackground->GetBinContent(center) / temp->GetBinWidth(1)/h2DBackground->GetXaxis()->GetBinWidth(1));

   for(int ibin=0; ibin<ana::nPhiBin; ibin++){
      const double offset = 0.01* hsig->GetBinWidth(1);
      const double center = -1*ana::PI +  ibin* hsig->GetBinWidth(1);
      int bin = h->FindBin(center+offset);
      int orginal_bin = -1;
      if(center <(-0.5*ana::PI-offset)) orginal_bin = hsig->FindBin(center+offset + 2*ana::PI);
      else orginal_bin = hsig->FindBin(center+offset);
      h->SetBinContent(bin, hsig->GetBinContent(orginal_bin));
      h->SetBinError(bin, hsig->GetBinError(orginal_bin));
   }

   delete hsig;
}

void proj1D_lr(TH2* h2DSignal, TH2D* h2DBackground, TH1ptr& h, const char* name)
{
   //h = new TH1D(name, ";#Delta#phi;dN/d(#Delta#phi)", ana::nPhiBin+2, -1*ana::PI-1./16*ana::PI, ana::PI + 1./16*ana::PI);
   h = new TH1D(name, ";#Delta#phi;dN/d(#Delta#phi)", ana::nPhiBin, -1*ana::PI, ana::PI);

   int negBinMin = 0;
   int negBinMax = h2DSignal->GetXaxis()->FindBin(-2.);
   int posBinMin = h2DSignal->GetXaxis()->FindBin(2.);
   int posBinMax = 33+1;

   TH1D* hNeg = h2DSignal->ProjectionY("hneg", negBinMin, negBinMax);
   TH1D* hPos = h2DSignal->ProjectionY("hpos", posBinMin, posBinMax);
   hNeg->Add(hPos);
   TH1D* hsig = new TH1D(*hNeg);
   hsig->SetName("hsig");
   delete hNeg;
   delete hPos;

   TH1D* temp_neg = h2DBackground->ProjectionY("hneg_bkg", negBinMin, negBinMax);
   TH1D* temp_pos = h2DBackground->ProjectionY("hpos_bkg", posBinMin, posBinMax);
   temp_neg->Add(temp_pos);

   TH1D* temp = new TH1D(*temp_neg);
   temp->SetName("temp");
   delete temp_neg;
   delete temp_pos;

   int center = h2DBackground->FindBin(0., 0.);
   hsig->Divide(temp);
   hsig->Scale(h2DBackground->GetBinContent(center) / temp->GetBinWidth(1)/h2DBackground->GetXaxis()->GetBinWidth(1));

   for(int ibin=0; ibin<ana::nPhiBin; ibin++){
      const double offset = 0.01* hsig->GetBinWidth(1);
      const double center = -1*ana::PI +  ibin* hsig->GetBinWidth(1);
      int bin = h->FindBin(center+offset);
      int orginal_bin = -1;
      if(center <(-0.5*ana::PI-offset)) orginal_bin = hsig->FindBin(center+offset + 2*ana::PI);
      else orginal_bin = hsig->FindBin(center+offset);
      h->SetBinContent(bin, hsig->GetBinContent(orginal_bin));
      h->SetBinError(bin, hsig->GetBinError(orginal_bin));
   }

   delete hsig;
}

void draw_proj1D(string dataset = "", string input = "", string output = "",
      string type = "", int itrk=-1, int ipt=-1, 
      const vector<float>* ptbin=nullptr, const vector<unsigned int>* trkbin=nullptr, 
      unordered_map<string, TGraphErrors*>* graph_map = nullptr)
{
   //gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   TFile* fin = TFile::Open(input.c_str());
   if(isFailed(fin, "input files")) return;

   TH2D* h2DSignal_D0[ana::nMass];
   TH2D* h2DBackground_D0[ana::nMass];
   TH1D* hMult_raw_D0[ana::nMass];
   TH1D* hMass_D0[ana::nMass];

   for(int imass=0; imass<ana::nMass; imass++){
      fin->GetObject(Form("hSignal_D0_mass%d_pt%d_trk%d%s", imass, ipt, itrk, type.c_str()), h2DSignal_D0[imass]);
      if(isFailed(h2DSignal_D0[imass], "hsignal")) return;
      fin->GetObject(Form("hBackground_D0_mass%d_pt%d_trk%d%s", imass, ipt, itrk, type.c_str()), h2DBackground_D0[imass]);
      if(isFailed(h2DBackground_D0[imass], "hbkg")) return;
      fin->GetObject(Form("hMult_raw_D0_mass%d_pt%d_trk%d%s", imass, ipt, itrk, type.c_str()), hMult_raw_D0[imass]);
      if(isFailed(hMult_raw_D0[imass], "hmult")) return;
      fin->GetObject(Form("hMass_D0_mass%d_pt%d_trk%d%s", imass, ipt, itrk, type.c_str()), hMass_D0[imass]);
      if(isFailed(hMass_D0[imass], "hmass")) return;
   }

   // scaled by event number
   for(int imass=0; imass<ana::nMass; imass++){
      double nMult_D0= hMult_raw_D0[imass]->Integral(2, 100000);
      h2DSignal_D0[imass]->Scale(1./nMult_D0);
   }

   // projection
   TH1ptr hDeltaPhi_D0[ana::nMass];
   for(int imass=0; imass<ana::nMass; imass++){
      proj1D_longrange(h2DSignal_D0[imass], h2DBackground_D0[imass], hDeltaPhi_D0[imass], 
            Form("h1D_mass%d_pt%d_trk%d_%s", imass, ipt, itrk, dataset.c_str()));
      //if(isFailed(hDeltaPhi_D0[imass]), "deltaPhi D0") return;
   }

   TCanvas c("c", "", 550*3, 500*5);
   c.Divide(3, 5);
   TLatex ltx;

   // fit Vn
   double Vn[ana::nMass];
   double Vn_err[ana::nMass];
   double nass[ana::nMass];
   double nass_err[ana::nMass];
   for(int imass=0; imass<ana::nMass; imass++){
      c.cd(imass+1);
      std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";

      gPad->SetBottomMargin(0.14);
      TF1 func("deltaPhi", function.c_str(), -3.14159*0.5, 3.14159*1.5);
      func.SetParameter(0, hDeltaPhi_D0[imass]->GetMaximum());
      func.SetParameter(1, 0.1);
      func.SetParameter(2, 0.1);
      func.SetParameter(3, 0.1);

      hDeltaPhi_D0[imass]->SetMarkerStyle(20);

      hDeltaPhi_D0[imass]->Fit(&func, "q");
      hDeltaPhi_D0[imass]->Fit(&func, "q");
      hDeltaPhi_D0[imass]->Fit(&func, "m q");
      hDeltaPhi_D0[imass]->Fit(&func, "m q E");
      auto fitResult = hDeltaPhi_D0[imass]->Fit(&func, "m S E q");

      hDeltaPhi_D0[imass]->SetTitle(";#Delta#phi;");
      hDeltaPhi_D0[imass]->GetXaxis()->CenterTitle();
      hDeltaPhi_D0[imass]->GetXaxis()->SetTitleSize(0.05);

      auto h = hDeltaPhi_D0[imass]->DrawCopy();
      gPad->Update();
      TPaveStats* pave = (TPaveStats*) h->FindObject("stats");
      pave->SetX1NDC(0.16);
      pave->SetX2NDC(0.56);
      gPad->Modified();
      gPad->Update();

      TString ntrk_str;
      if(trkbin->at(itrk+1) == numeric_limits<unsigned int>::max()) ntrk_str = Form("N_{trk}^{offline} #geq %u", trkbin->at(itrk));
      else ntrk_str = Form("%u #leq N_{trk}^{offline} < %u", trkbin->at(itrk), trkbin->at(itrk+1));
      ltx.DrawLatexNDC(0.65, 0.9, ntrk_str.Data());
      ltx.DrawLatexNDC(0.65, 0.38, "|#Delta#eta|>1");
      ltx.DrawLatexNDC(0.55, 0.31, Form("%.1f #leq p_{T} < %.1f GeV", ptbin->at(ipt), ptbin->at(ipt+1)));
      ltx.DrawLatexNDC(0.65, 0.24, "|y|<1");
      ltx.DrawLatexNDC(0.45, 0.17, Form("%.4f #leq mass < %.4f GeV", ana::massbin[imass], ana::massbin[imass+1]));

      Vn[imass] = func.GetParameter(2);
      Vn_err[imass] = func.GetParError(2);
      nass[imass] = func.GetParameter(0);
      nass_err[imass] = func.GetParError(0);
   }
   c.Print(Form("%s.pdf", output.c_str()));

   for(int imass=0; imass<ana::nMass; imass++){
      delete hDeltaPhi_D0[imass];
   }

   // jet calculation
   // short range
   TH1ptr hsr[ana::nMass];

   double min_sr[ana::nMass];
   double minX_sr[ana::nMass];

   double yields_sr[ana::nMass];
   double yields_sr_err[ana::nMass];

   for(int imass=0; imass<ana::nMass; imass++){
      proj1D_sr(h2DSignal_D0[imass], h2DBackground_D0[imass], hsr[imass],
            Form("h1Dsr_mass%d_pt%d_trk%d_%s", imass, ipt, itrk, dataset.c_str())
            );
   }

   TCanvas csr("csr", "", 550*3, 500*5);
   csr.Divide(3, 5);
   for(int imass=0; imass<ana::nMass; imass++){
      csr.cd(imass+1);
      TF1 func("func", fpoly2, -2.0, 2.0, 3);
      func.GetNpar();
      func.SetParameters(1, 1, 1);
      hsr[imass]->Fit(&func, "Q R 0", "", -2.0, 2.0);
      hsr[imass]->Fit(&func, "Q R 0", "", -2.0, 2.0);
      hsr[imass]->Fit(&func, "Q E R 0", "", -2.0, 2.0);
      auto fitResult = hsr[imass]->Fit(&func, "Q E R S 0", "", -2.0, 2.0);
      //fitResult->Print("v");
      cout << "sr fit result: " << fitResult->Status() << endl;
      TF1 fleft("fleft", fpoly2, -2.0, -0.6, 3);
      TF1 fright("fright", fpoly2, 0.6, 2.0, 3);
      fleft.SetParameters(func.GetParameters());
      fright.SetParameters(func.GetParameters());

      hsr[imass]->SetMarkerStyle(20); 
      hsr[imass]->Draw(); 
      fleft.DrawCopy("SAME");
      fright.DrawCopy("SAME");

      TString ntrk_str;
      if(trkbin->at(itrk+1) == numeric_limits<unsigned int>::max()) ntrk_str = Form("N_{trk}^{offline} #geq %u", trkbin->at(itrk));
      else ntrk_str = Form("%u #leq N_{trk}^{offline} < %u", trkbin->at(itrk), trkbin->at(itrk+1));
      ltx.DrawLatexNDC(0.15, 0.9, ntrk_str.Data());
      ltx.DrawLatexNDC(0.15, 0.8, Form("%.1f #leq p_{T} < %.1f GeV", ptbin->at(ipt), ptbin->at(ipt+1)));
      ltx.DrawLatexNDC(0.15, 0.73, "|y|<1");
      ltx.DrawLatexNDC(0.15, 0.66, "|#Delta#eta|<1");
      ltx.DrawLatexNDC(0.12, 0.59, Form("%.4f #leq mass < %.4f GeV", ana::massbin[imass], ana::massbin[imass+1]));

      min_sr[imass] = func.GetMinimum(0.6, 2.0);
      minX_sr[imass] = func.GetMinimumX(0.6, 2.0);
   }
   //csr.Print(Form("%s_sr.png", output.c_str()));
   csr.Print(Form("%s_sr.pdf", output.c_str()));
   
   for(int imass=0; imass<ana::nMass; imass++){
      TF1 fmin("fmin", "[0]", -1.* ana::PI, 1.*ana::PI);
      fmin.SetParameter(0, min_sr[imass]);
      hsr[imass]->Add(&fmin, -1);

      int binlw =hsr[imass]->GetXaxis()->FindBin(-1.*minX_sr[imass]);
      int binup =hsr[imass]->GetXaxis()->FindBin( 1.*minX_sr[imass]);

      yields_sr[imass] = hsr[imass]->IntegralAndError(binlw, binup, yields_sr_err[imass], "width");
   }

   for(int imass=0; imass<ana::nMass; imass++){
      delete hsr[imass];
   }

   // long range
   TH1ptr hlr[ana::nMass];

   double min_lr[ana::nMass];
   double minX_lr[ana::nMass];

   double yields_lr[ana::nMass];
   double yields_lr_err[ana::nMass];

   for(int imass=0; imass<ana::nMass; imass++){
      proj1D_lr(h2DSignal_D0[imass], h2DBackground_D0[imass], hlr[imass],
            Form("h1Dlr_mass%d_pt%d_trk%d_%s", imass, ipt, itrk, dataset.c_str())
            );
   }

   TCanvas clr("clr", "", 550*3, 500*5);
   clr.Divide(3, 5);
   for(int imass=0; imass<ana::nMass; imass++){
      clr.cd(imass+1);
      TF1 func("func", "[0]+[1]*fabs(x)+[2]*x*x", -2.0, 2.0);
      func.GetNpar();
      func.SetParameters(1, 1, 1);
      if(isHM1to6){
         hlr[imass]->Fit(&func, "Q R", "", 0.1, 2.0);
         hlr[imass]->Fit(&func, "Q R", "", 0.1, 2.0);
         hlr[imass]->Fit(&func, "Q E R", "", 0.1, 2.0);
      }
      else{
         hlr[imass]->Fit(&func, "Q R", "", -2.0, 2.0);
         hlr[imass]->Fit(&func, "Q R", "", -2.0, 2.0);
         hlr[imass]->Fit(&func, "Q E R", "", -2.0, 2.0);
      }
      auto fitResult = hlr[imass]->Fit(&func, "Q E R S", "", -2.0, 2.0);
      //fitResult->Print("v");
      cout << "lr fit result: " << fitResult->Status() << endl;

      hlr[imass]->SetMarkerStyle(20);
      hlr[imass]->Draw();

      TString ntrk_str;
      if(trkbin->at(itrk+1) == numeric_limits<unsigned int>::max()) ntrk_str = Form("N_{trk}^{offline} #geq %u", trkbin->at(itrk));
      else ntrk_str = Form("%u #leq N_{trk}^{offline} < %u", trkbin->at(itrk), trkbin->at(itrk+1));
      ltx.DrawLatexNDC(0.15, 0.9, ntrk_str.Data());
      ltx.DrawLatexNDC(0.2, 0.8, Form("%.1f #leq p_{T} < %.1f GeV", ptbin->at(ipt), ptbin->at(ipt+1)));
      ltx.DrawLatexNDC(0.35, 0.73, "|y|<1");
      ltx.DrawLatexNDC(0.35, 0.66, "|#Delta#eta|>2");
      ltx.DrawLatexNDC(0.2, 0.59, Form("%.4f #leq mass < %.4f GeV", ana::massbin[imass], ana::massbin[imass+1]));

      min_lr[imass] = func.GetMinimum(0.0, 2.0);
      minX_lr[imass] = func.GetMinimumX(0.0, 2.0);
   }
   //clr.Print(Form("%s_lr.png", output.c_str()));
   clr.Print(Form("%s_lr.pdf", output.c_str()));
   //
   for(int imass=0; imass<ana::nMass; imass++){
      TF1 fmin("fmin", "[0]", -1.* ana::PI, 1.*ana::PI);
      fmin.SetParameter(0, min_lr[imass]);
      hlr[imass]->Add(&fmin, -1);

      int binlw =hlr[imass]->GetXaxis()->FindBin(-1.*minX_lr[imass]);
      int binup =hlr[imass]->GetXaxis()->FindBin( 1.*minX_lr[imass]);

      yields_lr[imass] = hlr[imass]->IntegralAndError(binlw, binup, yields_lr_err[imass], "width");
   }

   for(int imass=0; imass<ana::nMass; imass++){
      delete hlr[imass];
   }

   double yields_jet[ana::nMass];
   double yields_jet_err[ana::nMass];

   for(int imass=0; imass<ana::nMass; imass++){
      yields_jet[imass] = yields_sr[imass] - yields_lr[imass];
      yields_jet_err[imass] = sqrt(
               pow(yields_sr_err[imass], 2) + pow(yields_lr_err[imass], 2)
            );
   }

   // mean values in each mass bin
   double mean[ana::nMass];
   double mean_err[ana::nMass];
   for(int imass=0; imass<ana::nMass; imass++){
      mean[imass] = hMass_D0[imass]->GetMean();
      mean_err[imass] = 0.;
   }

   // create TGraphErrors
   (*graph_map)["Vn"] = new TGraphErrors(ana::nMass, mean, Vn, mean_err, Vn_err);
   (*graph_map)["jets"] = new TGraphErrors(ana::nMass, mean, yields_jet, mean_err, yields_jet_err);
   (*graph_map)["nass"] = new TGraphErrors(ana::nMass, mean, nass, mean_err, nass_err);
}

void draw1D()
{
   vector<float> ptbin = {2., 4., 6., 8.};
   const int nPt = ptbin.size()-1;

   string dataset[] = {"PAMB", "PAHM1-6", "PAHM7"};
   string dataMult[] = {"PAMB0-185", "PAHM185-250", "PAHM250-inf"};
   string dataTrigger[] = {"PAMB", "PAHM", "PAHM"};

   string type[] = {"", "_loose", "_tight"};
   int itype = 0;

   string cmdDirPlots = 
   Form("if [ ! -d \"../plots\" ]; then\n"
         "    mkdir ../plots\n"
         "fi");
   string cmdDirPlotsCorr2D = 
   Form("if [ ! -d \"../plots/proj1D\" ]; then\n"
         "    mkdir ../plots/proj1D\n"
         "fi");
   gSystem->Exec(cmdDirPlots.c_str());
   gSystem->Exec(cmdDirPlotsCorr2D.c_str());


   for(int iset=0; iset<3; iset++){

      cout << dataset[iset] << endl;

      if(iset == 0) isMB = true;
      else isMB = false;
      if(iset == 1) isHM1to6 = true;
      else isHM1to6 = false;

      auto trkbin = ana::get_Mult_Edges(dataset[iset]);
      int ntrk = trkbin.size() - 1;

      unordered_map<string, TGraphErrors*> graph_map[ntrk][nPt];

      for(int itrk=0; itrk<ntrk; itrk++){
         for(int ipt=0; ipt<nPt; ipt++){
            TString output;
            if(trkbin.at(itrk+1) == numeric_limits<unsigned int>::max()) 
               output= Form("../plots/proj1D/%s%u-inf_pT%.0f-%.0f%s", dataTrigger[iset].c_str(), 
                     trkbin[itrk], ptbin[ipt], ptbin[ipt+1], type[itype].c_str());
            else 
               output = Form("../plots/proj1D/%s%u-%u_pT%.0f-%.0f%s", dataTrigger[iset].c_str(), 
                     trkbin[itrk], trkbin[itrk+1], ptbin[ipt], ptbin[ipt+1], type[itype].c_str());

            draw_proj1D(dataset[iset], Form("../data/corr2D_trg_pd0_%s.root", dataMult[iset].c_str()), 
                  output.Data(), type[itype], itrk, ipt, &ptbin, &trkbin, &graph_map[itrk][ipt]);
         }
      }

      TFile fout(Form("../data/%s_VarVsMass%s.root", dataMult[iset].c_str(), 
                     type[itype].c_str()), "recreate");
      for(int itrk=0; itrk<ntrk; itrk++){
         for(int ipt=0; ipt<nPt; ipt++){
            for(auto& graph : graph_map[itrk][ipt]) graph.second->Write(Form("%s_pt%d_trk%d", graph.first.c_str(), ipt, itrk));
         }
      }
   }
}

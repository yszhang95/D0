#ifndef __MY_INCLUDE__
#define __MY_INCLUDE__
R__ADD_INCLUDE_PATH(../include);
R__ADD_INCLUDE_PATH(../src);
#include "myAnaConsts.h"
#include "functions.cxx"
#endif

template <class T>
bool isFailed(T* ptr, const char* name = "")
{
   if(ptr) return false;
   cerr << name << " is nullptr!" << endl;
   return true;
}

using namespace std;
void draw_corr2D(string dataset = "", string input = "", string output = "",
      string type = "", int itrk=-1, int ipt=-1, 
      const vector<float>* ptbin=nullptr, const vector<unsigned int>* trkbin=nullptr)
{
   gErrorIgnoreLevel = kWarning;
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

   // scaled by event number, bin width, efficiency correction
   for(int imass=0; imass<ana::nMass; imass++){
      double nMult_D0= hMult_raw_D0[imass]->Integral(2, 100000);
      h2DSignal_D0[imass]->Scale(1./nMult_D0);
      const double bz = 
         h2DBackground_D0[imass]->GetBinContent(h2DBackground_D0[imass]->FindBin(0., 0.));
      const double binwidth_eta = h2DBackground_D0[imass]->GetXaxis()->GetBinWidth(1);
      const double binwidth_phi = h2DBackground_D0[imass]->GetYaxis()->GetBinWidth(1);
      h2DSignal_D0[imass]->Scale(bz/binwidth_eta/binwidth_phi);
      h2DSignal_D0[imass]->Divide(h2DBackground_D0[imass]);
   }
   TCanvas c("c", "", 550*3, 500*5);
   c.Divide(3, 5);
   TLatex ltx;
   for(int imass=0; imass<ana::nMass; imass++){
      c.cd(imass+1);
      h2DSignal_D0[imass]->SetTitle(";#Delta#eta;#Delta#phi");
      h2DSignal_D0[imass]->GetXaxis()->SetRangeUser(-2.4, 2.4);
      h2DSignal_D0[imass]->Draw("SURF1");
      TString ntrk_str;
      if(trkbin->at(itrk+1) == numeric_limits<unsigned int>::max()) ntrk_str = Form("N_{trk}^{offline} #geq %u", trkbin->at(itrk));
      else ntrk_str = Form("%u #leq N_{trk}^{offline} < %u", trkbin->at(itrk), trkbin->at(itrk+1));
      ltx.DrawLatexNDC(0.15, 0.88, ntrk_str.Data());
      ltx.DrawLatexNDC(0.15, 0.81, Form("%.1f #leq p_{T} < %.1f GeV", ptbin->at(ipt), ptbin->at(ipt+1)));
      ltx.DrawLatexNDC(0.15, 0.74, "|y|<1");
      ltx.DrawLatexNDC(0.15, 0.67, Form("%.4f #leq mass < %.4f GeV", ana::massbin[imass], ana::massbin[imass+1]));
   }
   c.Print(Form("%s.pdf", output.c_str()));
}

void draw2D()
{
   vector<float> ptbin = {2., 4., 6., 8.};
   const int nPt = ptbin.size()-1;

   string dataset[] = {"PAMB", "PAHM1-6", "PAHM7"};
   string dataMult[] = {"PAMB0-185", "PAHM185-250", "PAHM250-inf"};
   string dataTrigger[] = {"PAMB", "PAHM", "PAHM"};

   string type[] = {"", "_loose", "_tight"};
   int itype = 1;

   string cmdDirPlots = 
   Form("if [ ! -d \"../plots\" ]; then\n"
         "    mkdir ../plots\n"
         "fi");
   string cmdDirPlotsCorr2D = 
   Form("if [ ! -d \"../plots/corr2D\" ]; then\n"
         "    mkdir ../plots/corr2D\n"
         "fi");
   gSystem->Exec(cmdDirPlots.c_str());
   gSystem->Exec(cmdDirPlotsCorr2D.c_str());


   for(int iset=0; iset<3; iset++){
      auto trkbin = ana::get_Mult_Edges(dataset[iset]);
      int ntrk = trkbin.size() - 1;
      for(int itrk=0; itrk<ntrk; itrk++){
         for(int ipt=0; ipt<nPt; ipt++){
            TString output;
            if(trkbin.at(itrk+1) == numeric_limits<unsigned int>::max()) 
               output= Form("../plots/corr2D/%s%u-inf_pT%.0f-%.0f%s", dataTrigger[iset].c_str(), 
                     trkbin[itrk], ptbin[ipt], ptbin[ipt+1], type[itype].c_str());
            else 
               output = Form("../plots/corr2D/%s%u-%u_pT%.0f-%.0f%s", dataTrigger[iset].c_str(), 
                     trkbin[itrk], trkbin[itrk+1], ptbin[ipt], ptbin[ipt+1], type[itype].c_str());
            draw_corr2D(dataset[iset], Form("../data/corr2D_trg_pd0_%s.root", dataMult[iset].c_str()), 
                  output.Data(), "", itrk, ipt, &ptbin, &trkbin);
         }
      }
   }
}

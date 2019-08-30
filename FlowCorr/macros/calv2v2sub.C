#ifndef __MY_INCLUDE__
#define __MY_INCLUDE__
R__ADD_INCLUDE_PATH(../include);
R__ADD_INCLUDE_PATH(../src);
#include "myAnaConsts.h"
#include "functions.cxx"
#endif

int Offset(string dataset){
   if(dataset == "PAMB") return 0;
   if(dataset == "PAHM1-6") return 4;
   if(dataset == "PAHM7") return 5;
   return 0;
}

void calv2v2sub()
{
   vector<float> ptbin = {2., 4., 6., 8.};
   const int nPt = ptbin.size()-1;

   string dataset[] = {"PAMB", "PAHM1-6", "PAHM7"};
   string dataMult[] = {"PAMB0-185", "PAHM185-250", "PAHM250-inf"};
   string dataTrigger[] = {"PAMB", "PAHM", "PAHM"};

   string type[] = {"", "_loose", "_tight"};
   int itype = 1;

   const int n_v2 = 6;
   const int n_v2sub = 5;

   //double mean[n_v2];
   double *mean;

   double v2[nPt][n_v2];
   double v2sub[nPt][n_v2sub];
   double v2_err[nPt][n_v2];
   double v2sub_err[nPt][n_v2sub];

   double v2_upper[nPt][n_v2];
   double v2sub_upper[nPt][n_v2sub];
   double v2_err_upper[nPt][n_v2];
   double v2sub_err_upper[nPt][n_v2sub];

   double v2_lower[nPt][n_v2];
   double v2sub_lower[nPt][n_v2sub];
   double v2_err_lower[nPt][n_v2];
   double v2sub_err_lower[nPt][n_v2sub];

   double Vn_ref[n_v2];
   double Vn_ref_err[n_v2];
   double jet_ref[n_v2];
   double jet_ref_err[n_v2];
   double nass_ref[n_v2];
   double nass_ref_err[n_v2];

   double* Vn[nPt];
   double* Vn_err[nPt];
   double* jet[nPt];
   double* jet_err[nPt];
   double* nass[nPt];
   double* nass_err[nPt];

   for(int ipt=0; ipt<nPt; ipt++){
      Vn[ipt] = new double[n_v2];
      Vn_err[ipt] = new double[n_v2];
      jet[ipt] = new double[n_v2];
      jet_err[ipt] = new double[n_v2];
      nass[ipt] = new double[n_v2];
      nass_err[ipt] = new double[n_v2];
   }

   /*
   double Vn[nPt][n_v2];
   double Vn_err[nPt][n_v2];
   double jet[nPt][n_v2];
   double jet_err[nPt][n_v2];
   double nass[nPt][n_v2];
   double nass_err[nPt][n_v2];
   */

   double jet_ratio_d0[nPt][n_v2];
   double jet_ratio_d0_err[nPt][n_v2];

   double jet_ratio_ref[n_v2];
   double jet_ratio_ref_err[n_v2];

   vector<double> tmp_ref_Vn;
   vector<double> tmp_ref_Vn_err;
   vector<double> tmp_ref_jet;
   vector<double> tmp_ref_jet_err;
   vector<double> tmp_ref_nass;
   vector<double> tmp_ref_nass_err;

   // all mult
   for(int iset=0; iset<3; iset++){

      TString input_ref = 
         Form("../data/%s_ref.root", dataMult[iset].c_str());

      TFile* f_ref = TFile::Open(input_ref.Data());

      auto trkbin = ana::get_Mult_Edges(dataset[iset]);
      int ntrk = trkbin.size() - 1;

      //cout << ntrk << endl;

      for(int itrk=0; itrk<ntrk; itrk++){
         int offset = Offset(dataset[iset]);
         TGraphErrors* gVn = (TGraphErrors*) f_ref->Get(Form("Vn_ref_trk%d", itrk));
         tmp_ref_Vn.push_back( gVn->GetY()[0] );
         tmp_ref_Vn_err.push_back( gVn->GetEY()[0] );
         TGraphErrors* gjet = (TGraphErrors*) f_ref->Get(Form("jets_ref_trk%d", itrk));
         tmp_ref_jet.push_back( gjet->GetY()[0] );
         tmp_ref_jet_err.push_back( gjet->GetEY()[0] );
         TGraphErrors* gnass = (TGraphErrors*) f_ref->Get(Form("nass_ref_trk%d", itrk));
         tmp_ref_nass.push_back( gnass->GetY()[0] );
         tmp_ref_nass_err.push_back( gnass->GetEY()[0] );
      }
   }

   copy(tmp_ref_Vn.begin(), tmp_ref_Vn.end(), Vn_ref);
   copy(tmp_ref_Vn_err.begin(), tmp_ref_Vn_err.end(), Vn_ref_err);
   copy(tmp_ref_jet.begin(), tmp_ref_jet.end(), jet_ref);
   copy(tmp_ref_jet_err.begin(), tmp_ref_jet_err.end(), jet_ref_err);
   copy(tmp_ref_nass.begin(), tmp_ref_nass.end(), nass_ref);
   copy(tmp_ref_nass_err.begin(), tmp_ref_nass_err.end(), nass_ref_err);

   for(int ipt=0; ipt<nPt; ipt++){
      TFile* fin_d0 = TFile::Open(Form("../data/PA_vars_%s.root", type[itype].c_str()));
      TGraphErrors* gVn = (TGraphErrors*) fin_d0->Get(Form("gVn_pt%d", ipt));
      TGraphErrors* gjet = (TGraphErrors*) fin_d0->Get(Form("gjet_pt%d", ipt));
      TGraphErrors* gnass = (TGraphErrors*) fin_d0->Get(Form("gnass_pt%d", ipt));
      Vn[ipt] = gVn->GetY();
      Vn_err[ipt] = gVn->GetEY();
      jet[ipt] = gjet->GetY();
      jet_err[ipt] = gjet->GetEY();
      nass[ipt] = gnass->GetY();
      nass_err[ipt] = gnass->GetEY();

      mean = gVn->GetX();
   }

   TCanvas* cjetratio = new TCanvas("cjetratio", "", 450*3, 550);
   cjetratio->Divide(3, 1);
   TGraphErrors* gjet_ratio[3];
   TGraphErrors* gconf[3];
   TGraphErrors* gratio_limits[3];
   TF1* fratio[3];
   for(int ipt=0; ipt<3; ipt++){
      cjetratio->cd(ipt+1);
      gjet_ratio[ipt] = new TGraphErrors(n_v2);
      gratio_limits[ipt] = new TGraphErrors(n_v2);
      gconf[ipt] = new TGraphErrors(300);
      fratio[ipt] = new TF1(Form("fratio_pt%d", ipt), "1.+exp([0]+[1]/x)", 0, 300);
      for(int ip=0; ip<n_v2; ip++){
         double ratio = jet[ipt][ip]/jet[ipt][0];
         double ratio_err = ratio*sqrt(
               pow(jet_err[ipt][ip]/jet[ipt][ip], 2) + pow(jet_err[ipt][0]/jet[ipt][0], 2)
               );
         if(ip==0) ratio_err = 0.;
         cout << ratio_err << endl;
         gjet_ratio[ipt]->SetPoint(ip, mean[ip], ratio);
         gjet_ratio[ipt]->SetPointError(ip, 0, ratio_err);
         gratio_limits[ipt]->SetPoint(ip, mean[ip], 0);
      }
      for(int ip=0; ip<300; ip++){
         gconf[ipt]->SetPoint(ip, ip, 0);
      }
      fratio[ipt]->SetParameters(0.01, 0.001);
      gjet_ratio[ipt]->Fit(fratio[ipt], "R", 0, 300);
      gjet_ratio[ipt]->Fit(fratio[ipt], "R", 0, 300);
      gjet_ratio[ipt]->Fit(fratio[ipt], "R", 0, 300);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gconf[ipt], 0.683);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gratio_limits[ipt], 0.683);
      gjet_ratio[ipt]->SetMarkerStyle(20);
      gjet_ratio[ipt]->SetMarkerColor(kRed);
      gjet_ratio[ipt]->SetTitle(";N_{trk}^{offline};Y_{jet} ratio");
      auto h = gjet_ratio[ipt]->GetHistogram();
      h->SetMaximum(2.5);
      h->SetMinimum(0.);
      gjet_ratio[ipt]->Draw("A P");
      gconf[ipt]->SetFillStyle(1001);
      gconf[ipt]->SetFillColorAlpha(kBlue-7, 0.3);
      gconf[ipt]->Draw("E3 SAME");
   }
   cjetratio->Print(Form("jet_ratio%s.pdf", type[itype].c_str()));
   cjetratio->Print(Form("jet_ratio%s.svg", type[itype].c_str()));

   for(int ipt = 0; ipt < nPt; ipt++){
      for(int i=0; i<n_v2; i++){
         v2[ipt][i] = Vn[ipt][i]/sqrt(Vn_ref[i]);
         v2_err[ipt][i] = v2[ipt][i]* sqrt(pow(Vn_err[ipt][i]/Vn[ipt][i], 2) 
                  + pow(Vn_ref_err[i]/Vn_ref[i]/2., 2));
      }
      //note, i=1
      for(int i=1; i<n_v2; i++){
         //double scale_d0 = nass[ipt][0]/nass[ipt][i] * jet[ipt][i]/jet[ipt][0];
         double scale_d0 = nass[ipt][0]/nass[ipt][i] * fratio[ipt]->Eval(mean[i]);
         double scale_d0_upper = nass[ipt][0]/nass[ipt][i] * (gratio_limits[ipt]->GetY()[i] + gratio_limits[ipt]->GetEY()[i]) ;
         double scale_d0_lower = nass[ipt][0]/nass[ipt][i] * (gratio_limits[ipt]->GetY()[i] - gratio_limits[ipt]->GetEY()[i]) ;
         //cout << mean[i] << "  " << fratio[ipt]->Eval(mean[i]) << "  " << jet[ipt][i]/jet[ipt][0] << endl;
         double V2_sub = Vn[ipt][i] - Vn[ipt][0] * scale_d0;
         double V2_sub_err = sqrt(pow(Vn_err[ipt][i], 2) + pow(Vn_err[ipt][0]*scale_d0, 2));

         double V2_sub_upper = Vn[ipt][i] - Vn[ipt][0] * scale_d0_upper;
         double V2_sub_err_upper = sqrt(pow(Vn_err[ipt][i], 2) + pow(Vn_err[ipt][0]*scale_d0_upper, 2));

         double V2_sub_lower = Vn[ipt][i] - Vn[ipt][0] * scale_d0_lower;
         double V2_sub_err_lower = sqrt(pow(Vn_err[ipt][i], 2) + pow(Vn_err[ipt][0]*scale_d0_lower, 2));

         //cout << V2_sub << endl;

         double scale_ref = nass_ref[0]/nass_ref[i] * jet_ref[i] / jet_ref[0];
         //cout << nass_ref[0] / nass_ref[i] << endl;
         //cout << scale_ref << endl;

         double V2_ref_sub = Vn_ref[i] - Vn_ref[0]* scale_ref;
         double V2_ref_sub_err = sqrt(
                     pow(Vn_ref_err[i], 2) + pow(Vn_ref_err[0]*scale_ref, 2)
                  );
         //cout << V2_ref_sub << endl;

         v2sub[ipt][i-1] = V2_sub/ sqrt(V2_ref_sub);
         v2sub_err[ipt][i-1] = fabs(v2sub[ipt][i-1])* sqrt(pow(V2_sub_err/V2_sub, 2) 
                  + pow(V2_ref_sub_err/V2_ref_sub/2., 2));

         v2sub_upper[ipt][i-1] = V2_sub_upper/ sqrt(V2_ref_sub);
         v2sub_err_upper[ipt][i-1] = fabs(v2sub_upper[ipt][i-1])* sqrt(pow(V2_sub_err_upper/V2_sub_upper, 2) 
                  + pow(V2_ref_sub_err/V2_ref_sub/2., 2));

         v2sub_lower[ipt][i-1] = V2_sub_lower/ sqrt(V2_ref_sub);
         v2sub_err_lower[ipt][i-1] = fabs(v2sub_lower[ipt][i-1])* sqrt(pow(V2_sub_err_lower/V2_sub_lower, 2) 
                  + pow(V2_ref_sub_err/V2_ref_sub/2., 2));
      }
   }

   for(int ipt=0; ipt<nPt; ipt++){
      cout << ipt << endl;
      for(int i=0; i<n_v2sub; i++){
         cout << Form("trk%d: ", i) 
            << v2sub[ipt][i] << ", " << v2sub_err[ipt][i] << ", "
            << Vn[ipt][i+1]<<", " << Vn_err[ipt][i+1] << ", "
            << Vn[ipt][0] << ", " << Vn_err[ipt][0] << ", "
            << Vn_ref[i+1] << ", " << Vn_ref_err[i+1]  << ", "
            << Vn_ref[0] << ", " << Vn_ref_err[0]  << ", "
            << endl;
      }
   }

   double e[n_v2] = {0.};
   TFile fout(Form("output%s.root", type[itype].c_str()), "recreate");
   for(int ipt=0; ipt<nPt; ipt++){
      //TGraphErrors* jetratio = new TGraphErrors(n_v2, mean, jet_ratio_d0[ipt], e, jet_ratio_d0_err[ipt]);
      //jetratio->Write(Form("jets_ratio_pt%d", ipt));
      TGraphErrors* gv2 = new TGraphErrors(n_v2, mean, v2[ipt], e, v2_err[ipt]);
      gv2->Write(Form("v2_pt%d", ipt));
      TGraphErrors* gv2sub = new TGraphErrors(n_v2sub, &mean[1], v2sub[ipt], e, v2sub_err[ipt]);
      gv2sub->Write(Form("v2sub_pt%d", ipt));

      TGraphErrors* gv2sub_upper = new TGraphErrors(n_v2sub, &mean[1], v2sub_upper[ipt], e, v2sub_err_upper[ipt]);
      gv2sub_upper->Write(Form("v2sub_upper_pt%d", ipt));

      TGraphErrors* gv2sub_lower = new TGraphErrors(n_v2sub, &mean[1], v2sub_lower[ipt], e, v2sub_err_lower[ipt]);
      gv2sub_lower->Write(Form("v2sub_lower_pt%d", ipt));
   }
}

void fitNass_npd0(const char* input, const char* output)
{
   TFile* file0 = new TFile(input);
   TFile fout(output, "recreate");
   TGraphErrors* gNass[2][2];
   TGraphErrors* gNass_low[2][2];
   for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
         file0->GetObject(Form("NassvsMass_pt%d_dca%d", i, j), gNass[i][j]);
         file0->GetObject(Form("Nass_lowvsMass_pt%d_dca%d", i, j), gNass_low[i][j]);
      }
   }

   double Nass[2][2];
   double Nass_low[2][2];
   
   for(int ipt=0; ipt<2; ipt++){
      for(int j=0; j<2; j++){
         TF1 f("f", "[0]", 1.7, 2.0);
         gNass[ipt][j]->Fit(&f, "Q R");
         gNass[ipt][j]->Fit(&f, "Q R");
         gNass[ipt][j]->Fit(&f, "");
         Nass[ipt][j] = f.GetParameter(0);
      }
   }

   for(int ipt=0; ipt<2; ipt++){
      for(int j=0; j<2; j++){
         TF1 f("f", "[0]", 1.7, 2.0);
         gNass_low[ipt][j]->Fit(&f, "Q R");
         gNass_low[ipt][j]->Fit(&f, "Q R");
         gNass_low[ipt][j]->Fit(&f, "");
         Nass_low[ipt][j] = f.GetParameter(0);
      }
   }
   TGraphErrors* g[2] ;
   TGraphErrors* g_low[2] ;
   for(int j=0; j<2; j++){
      g[j] = new TGraphErrors(2);
      for(int ipt=0; ipt<2; ipt++){
         g[j]->SetPoint(ipt, ipt, Nass[ipt][j]);
      }
   }

   for(int j=0; j<2; j++){
      g_low[j] = new TGraphErrors(2);
      for(int ipt=0; ipt<2; ipt++){
         g_low[j]->SetPoint(ipt, ipt, Nass_low[ipt][j]);
      }
   }
   fout.cd();
   g[0]->Write("HM_dca0");
   g[1]->Write("HM_dca1");
   g_low[0]->Write("LM_dca0");
   g_low[1]->Write("LM_dca1");
}

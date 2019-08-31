void calv2()
{
   TFile* fHMD0 = new TFile("../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root");
   TFile* fMBD0 = new TFile("../data/corr2D_trg_npd0_MB0-185_pT2.0-8.0_y-1.0-1.0_new.root");

   TH1D* hPt_pt0_HM_dca_0 = (TH1D*) fHMD0->Get("hPt_pt0_dca0");
   TH1D* hPt_pt0_HM_dca_1 = (TH1D*) fHMD0->Get("hPt_pt0_dca1");
   TH1D* hPt_pt0_MB_dca_0 = (TH1D*) fMBD0->Get("hPt_pt0_dca0");
   TH1D* hPt_pt0_MB_dca_1 = (TH1D*) fMBD0->Get("hPt_pt0_dca1");

   TH1D* hPt_pt1_HM_dca_0 = (TH1D*) fHMD0->Get("hPt_pt1_dca0");
   TH1D* hPt_pt1_HM_dca_1 = (TH1D*) fHMD0->Get("hPt_pt1_dca1");
   TH1D* hPt_pt1_MB_dca_0 = (TH1D*) fMBD0->Get("hPt_pt1_dca0");
   TH1D* hPt_pt1_MB_dca_1 = (TH1D*) fMBD0->Get("hPt_pt1_dca1");

   TH1D* hPt_pt0_HM = new TH1D(*hPt_pt0_HM_dca_0);
   hPt_pt0_HM->SetName("hPt_pt0_HM");
   TH1D* hPt_pt0_MB = new TH1D(*hPt_pt0_MB_dca_0);
   hPt_pt0_MB->SetName("hPt_pt0_MB");

   TH1D* hPt_pt1_HM = new TH1D(*hPt_pt1_HM_dca_0);
   hPt_pt1_HM->SetName("hPt_pt1_HM");
   TH1D* hPt_pt1_MB = new TH1D(*hPt_pt1_MB_dca_0);
   hPt_pt1_MB->SetName("hPt_pt1_MB");

   hPt_pt0_HM->Add(hPt_pt0_HM_dca_1);
   hPt_pt0_MB->Add(hPt_pt0_MB_dca_1);

   hPt_pt1_HM->Add(hPt_pt1_HM_dca_1);
   hPt_pt1_MB->Add(hPt_pt1_MB_dca_1);

   cout << hPt_pt0_HM->GetMean() << endl;
   cout << hPt_pt1_HM->GetMean() << endl;

   cout << hPt_pt0_MB->GetMean() << endl;
   cout << hPt_pt1_MB->GetMean() << endl;

   double pt_HM[2];
   double pt_MB[2];
   pt_HM[0] = hPt_pt0_HM->GetMean();
   pt_HM[1] = hPt_pt1_HM->GetMean();
   pt_MB[0] = hPt_pt0_MB->GetMean();
   pt_MB[1] = hPt_pt1_MB->GetMean();

   TFile* file0 = new TFile("ofile_npd0.root");
   //TFile* file0 = new TFile("ofile_pd0.root");
   TFile* file1 = new TFile("../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_raw.root");
   TGraphErrors*	V2_Ref_low  = (TGraphErrors*)file1->Get("V2_Ref_low"); 
   TGraphErrors*	Jets_ref_low  = (TGraphErrors* )file1->Get("Jets_ref_low"); 
   TGraphErrors*	Nass_ref_low  = (TGraphErrors*) file1->Get("Nass_ref_low"); 
   TGraphErrors*	V2_Ref= (TGraphErrors*) file1->Get("V2_Ref"); 
   TGraphErrors*	Jets_ref= (TGraphErrors*) file1->Get("Jets_ref"); 
   TGraphErrors*	Nass_ref= (TGraphErrors*) file1->Get("Nass_ref"); 
   cout << "REF jets ratio: " << Jets_ref->GetY()[0]/Jets_ref_low->GetY()[0] << endl;
   cout << Nass_ref_low->GetY()[0]<< ", "
      << Nass_ref->GetY()[0] << ", "
      << Jets_ref_low->GetY()[0] << ", "
      << Jets_ref->GetY()[0] << ", "
      << V2_Ref_low->GetY()[0] << ", "
      << V2_Ref->GetY()[0] << ","
      << V2_Ref_low->GetEY()[0] << ","
      << V2_Ref->GetEY()[0] << endl;

   double scale_ref = 
      Nass_ref_low->GetY()[0]/Nass_ref->GetY()[0] * Jets_ref->GetY()[0]/Jets_ref_low->GetY()[0]; 
   double V2_sub_REF = V2_Ref->GetY()[0] - V2_Ref_low->GetY()[0] *scale_ref;
   double V2_sub_REF_err = sqrt( pow(V2_Ref->GetEY()[0],2) + pow(V2_Ref_low->GetEY()[0] * scale_ref, 2));
   cout << "V2_Sub_REF: " <<  V2_sub_REF << "+/-" <<V2_sub_REF_err <<endl;

   TGraphErrors* Jets_low = (TGraphErrors*) file0->Get("Jets_low");
   TGraphErrors* Jets   = (TGraphErrors*) file0->Get("Jets");
   TGraphErrors* V2_low = (TGraphErrors*) file0->Get("V2_low");
   TGraphErrors* V2   = (TGraphErrors*) file0->Get("V2");
   TGraphErrors* Nass_low = (TGraphErrors*) file0->Get("Nass_low");
   TGraphErrors* Nass   = (TGraphErrors*) file0->Get("Nass");

   TGraphErrors* g_v2sub = new TGraphErrors(2);
   TGraphErrors* g_v2sub_lower = new TGraphErrors(2);
   TGraphErrors* g_v2sub_upper = new TGraphErrors(2);
   TGraphErrors* g_v2 = new TGraphErrors(2);
   TGraphErrors* g_v2_low = new TGraphErrors(2);

   for(int i=0; i<2; i++){
      cout << "pt" << i << endl;
      cout << "jets ratio: " << Jets->GetY()[i]/Jets_low->GetY()[i]  << endl;
      cout << "jets ratio err: " << 
         sqrt(pow(Jets->GetEY()[i]/Jets_low->GetY()[i], 2)+
            pow(Jets->GetY()[i]*Jets_low->GetEY()[i]/Jets->GetY()[i]/Jets->GetY()[i], 2)
            ) 
         << endl;
      double scale = 
         Nass_low->GetY()[i]/Nass->GetY()[i] * Jets->GetY()[i]/Jets_low->GetY()[i]; 
      cout << Nass_low->GetY()[i] << ","
       << Nass->GetY()[i] << ","
       << Jets_low->GetY()[i] << ","
       << Jets->GetY()[i] << ","
       << V2_low->GetY()[i] <<","
       << V2->GetY()[i] << ","
       << V2_low->GetEY()[i] <<","
       << V2->GetEY()[i] << endl;
      cout << scale << endl;
      cout << V2_low->GetY()[i] * scale << endl;
      double V2_sub_D0 = V2->GetY()[i] - V2_low->GetY()[i] *scale;
      double V2_sub_D0_err = sqrt( pow(V2->GetEY()[i],2) + pow(V2_low->GetEY()[i] * scale, 2));
      cout << "V2_sub_D0: " <<  V2_sub_D0 << "+/-" <<V2_sub_D0_err <<endl;
      double v2_sub_D0 = V2_sub_D0/ sqrt(V2_sub_REF);
      double v2_sub_D0_err = fabs(v2_sub_D0)*sqrt(pow(V2_sub_D0_err/V2_sub_D0, 2)+ pow(V2_sub_REF_err/V2_sub_REF/2., 2));
      cout << "v2_sub_D0: " <<  v2_sub_D0 << "+/-" <<v2_sub_D0_err <<endl;

      g_v2sub->SetPoint(i, pt_HM[i], v2_sub_D0);
      g_v2sub->SetPointError(i, 0, v2_sub_D0_err);

      double scale_upper = 
         Nass_low->GetY()[i]/Nass->GetY()[i] * (Jets->GetY()[i]/Jets_low->GetY()[i] + 
            sqrt(pow(Jets->GetEY()[i]/Jets_low->GetY()[i], 2)+
               pow(Jets->GetY()[i]*Jets_low->GetEY()[i]/Jets->GetY()[i]/Jets->GetY()[i], 2)
            ) 
               ); 

      double V2_sub_D0_upper = V2->GetY()[i] - V2_low->GetY()[i] *scale_upper;
      double V2_sub_D0_upper_err = sqrt( pow(V2->GetEY()[i],2) + pow(V2_low->GetEY()[i] * scale_upper, 2));

      double v2_sub_D0_upper = V2_sub_D0_upper/ sqrt(V2_sub_REF);
      double v2_sub_D0_upper_err = fabs(v2_sub_D0_upper)*sqrt(pow(V2_sub_D0_upper_err/V2_sub_D0_upper, 2)+ pow(V2_sub_REF_err/V2_sub_REF/2., 2));

      g_v2sub_upper->SetPoint(i, pt_HM[i], v2_sub_D0_upper);
      g_v2sub_upper->SetPointError(i, 0, v2_sub_D0_upper_err);

      double scale_lower = 
         Nass_low->GetY()[i]/Nass->GetY()[i] * (Jets->GetY()[i]/Jets_low->GetY()[i] -
            sqrt(pow(Jets->GetEY()[i]/Jets_low->GetY()[i], 2)+
               pow(Jets->GetY()[i]*Jets_low->GetEY()[i]/Jets->GetY()[i]/Jets->GetY()[i], 2)
            ) 
               ); 

      double V2_sub_D0_lower = V2->GetY()[i] - V2_low->GetY()[i] *scale_lower;
      double V2_sub_D0_lower_err = sqrt( pow(V2->GetEY()[i],2) + pow(V2_low->GetEY()[i] * scale_lower, 2));

      double v2_sub_D0_lower = V2_sub_D0_lower/ sqrt(V2_sub_REF);
      double v2_sub_D0_lower_err = fabs(v2_sub_D0_lower)*sqrt(pow(V2_sub_D0_lower_err/V2_sub_D0_lower, 2)+ pow(V2_sub_REF_err/V2_sub_REF/2., 2));

      g_v2sub_lower->SetPoint(i, pt_HM[i], v2_sub_D0_lower);
      g_v2sub_lower->SetPointError(i, 0, v2_sub_D0_lower_err);
   }

   for(int ipt=0; ipt<2; ipt++){
      double v2 = V2->GetY()[ipt]/sqrt(V2_Ref->GetY()[0]);
      double v2_err = fabs(v2)*sqrt(pow(V2->GetEY()[ipt]/V2->GetY()[ipt], 2)
            + pow(V2_Ref->GetEY()[0]/V2_Ref->GetY()[0]/2., 2)
            );
      g_v2->SetPoint(ipt, pt_HM[ipt], v2);
      g_v2->SetPointError(ipt, 0, v2_err);
   }

   for(int ipt=0; ipt<2; ipt++){
      double v2 = V2_low->GetY()[ipt]/sqrt(V2_Ref_low->GetY()[0]);
      double v2_err = fabs(v2)*sqrt(pow(V2_low->GetEY()[ipt]/V2_low->GetY()[ipt], 2)
            + pow(V2_Ref_low->GetEY()[0]/V2_Ref_low->GetY()[0]/2., 2)
            );
      g_v2_low->SetPoint(ipt, pt_MB[ipt], v2);
      g_v2_low->SetPointError(ipt, 0, v2_err);
   }

   TFile ofile("nonprompt.root", "recreate");
   g_v2sub->Write("v2sub");
   g_v2sub_upper->Write("v2sub_upper");
   g_v2sub_lower->Write("v2sub_lower");
   g_v2->Write("v2");
   g_v2_low->Write("v2_low");
}

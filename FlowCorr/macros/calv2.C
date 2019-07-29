void calv2()
{
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
   for(int i=0; i<2; i++){
      cout << "pt" << i << endl;
      cout << "jets ratio: " << Jets->GetY()[i]/Jets_low->GetY()[i]  << endl;
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
   }
}

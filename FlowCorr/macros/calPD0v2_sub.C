inline double sub(double, double, double);
inline double sub_err(double, double, double);
void calPD0v2_sub(
      const char* input_Vn,
      const char* input_Vn_low,
      const char* input_Jets,
      const char* input_Jets_low,
      const char* nass,
      const char* ref,
      const char* dataset
      )
{
   TFile* fVn = new TFile(input_Vn);
   TFile* fVn_low = new TFile(input_Vn_low);
   TFile* fJets = new TFile(input_Jets);
   TFile* fJets_low = new TFile(input_Jets_low);
   TFile* fnass = new TFile(nass);
   TFile* fref = new TFile(ref);

   auto vec_trkbin = ana::get_Mult_Edges(dataset);
   if(!vec_trkbin.size()) {
      std::cerr << "wrong dataset" << std::endl;
      return;
   }
   const unsigned int n_trk_bin_ = vec_trkbin.size() - 1;

   TGraphErrors* Vn;
   TGraphErrors* Jets;
   TGraphErrors* Nass;
   TGraphErrors* Vn_low;
   TGraphErrors* Jets_low;
   TGraphErrors* Nass_low;

   TGraphErrors* Vn_ref[n_trk_bin_];
   TGraphErrors* Jets_ref[n_trk_bin_];
   TGraphErrors* Nass_ref[n_trk_bin_];
   TGraphErrors* Vn_ref_low;
   TGraphErrors* Jets_ref_low;
   TGraphErrors* Nass_ref_low;

   TH1D* hNtrk[n_trk_bin_];

   // ref, high N and low N
   double Vn_Ref[n_trk_bin_];
   double Vn_Ref_err[n_trk_bin_];
   double Vn_Ref_low;
   double Vn_Ref_low_err;
   double Jets_Ref[n_trk_bin_];
   double Jets_Ref_low;
   double Nass_Ref[n_trk_bin_];
   double Nass_Ref_low;

   double* Vn_D0;
   double* Vn_D0_err;
   double* Jets_D0;
   double* Nass_D0;
   double Vn_low_D0;
   double Vn_low_D0_err;
   double Jets_low_D0;
   double Nass_low_D0;

   double* Jets_D0_err;
   double Jets_low_D0_err;

   fVn->GetObject("v2vsNtrk", Vn);
   if(!Vn) {
      cerr<< "not found Vn" << endl;
      return;
   }

   fVn_low->GetObject("v2vsNtrk", Vn_low);
   if(!Vn_low) {
      cerr<< "not found Vn_low" << endl;
      return;
   }

   fJets->GetObject("jetvsNtrk", Jets);
   if(!Jets) {
      cerr<< "not found Jets" << endl;
      return;
   }

   fJets_low->GetObject("jetvsNtrk", Jets_low);
   if(!Jets_low) {
      cerr<< "not found Jets_low" << endl;
      return;
   }

   fnass->GetObject("nass", Nass);
   if(!Nass) {
      cerr<< "not found nass" << endl;
      return;
   }

   fnass->GetObject("nass_low", Nass_low);
   if(!Nass) {
      cerr<< "not found nass" << endl;
      return;
   }

   for(int i=0; i<n_trk_bin_; i++){
      fref->GetObject(Form("V2_Ref_trk%d", i), Vn_ref[i]);
      fref->GetObject(Form("Jets_Ref_trk%d", i), Jets_ref[i]);
      fref->GetObject(Form("Nass_Ref_trk%d", i), Nass_ref[i]);
      Vn_Ref[i] = Vn_ref[i]->GetY()[0];
      Jets_Ref[i] = Jets_ref[i]->GetY()[0];
      Nass_Ref[i] = Nass_ref[i]->GetY()[0];
   }
   fref->GetObject("V2_Ref_low", Vn_ref_low);
   fref->GetObject("Jets_ref_low", Jets_ref_low);
   fref->GetObject("Nass_ref_low", Nass_ref_low);
   Vn_Ref_low = Vn_ref_low->GetY()[0]; // 1
   Jets_Ref_low = Jets_ref_low->GetY()[0]; // 1
   Nass_Ref_low = Nass_ref_low->GetY()[0]; // 1

   Vn_D0 = Vn->GetY(); // n_trk_bin_
   Vn_D0_err = Vn->GetEY(); // n_trk_bin_
   Jets_D0 = Jets->GetY(); // n_trk_bin_
   Nass_D0 = Nass->GetY(); // n_trk_bin_
   Vn_low_D0 = Vn_low->GetY()[0]; // 1
   Vn_low_D0_err = Vn_low->GetEY()[0]; // 1
   Jets_low_D0 = Jets_low->GetY()[0]; // 1 
   Nass_low_D0 = Nass_low->GetY()[0]; // 1

   Jets_D0_err = Jets->GetEY();
   Jets_low_D0_err = Jets_low->GetEY()[0];

   for(int i=0; i<n_trk_bin_; i++){
      //cout << Jets_D0[i]/ Jets_low_D0 << endl;
   }

   TGraphErrors gratio(n_trk_bin_);

   ofstream ofile(Form("%s.csv", ref));
   cout << ref << endl;
   ofile << "#"
         << "vn_sub," << "vn_sub_err,"
         << "Vn_D0, " << "Vn_D0_err,"
         << "Vn_low_D0," << "Vn_low_D0_err,"
         << "Vn_Ref,"  << "Vn_Ref_err,"
         << "Vn_Ref_low,"  << "Vn_Ref_low_err,"
         << "Nass_D0," << "Nass_low_D0,"
         << "Nass_Ref," << "Nass_low_Ref,"
         << "Jets," << "Jets_low,"
         << "Jets_ref," << "Jets_ref_low,"
         << "Jets_ratio," << "Jets_ratio_err,"
         << endl;
   for(int i=0; i<n_trk_bin_; i++){
      double scale_ref = Jets_Ref[i]/Jets_Ref_low * Nass_Ref_low/Nass_Ref[i];
      //cout << scale_ref << endl;
      double scale_d0 = Jets_D0[i]/Jets_low_D0 * Nass_low_D0/Nass_D0[i];
      //cout << scale_d0 << endl;
      double Vn_Ref_sub = sub(Vn_Ref[i], Vn_Ref_low, scale_ref);
      double Vn_D0_sub = sub(Vn_D0[i], Vn_low_D0, scale_d0);
      double Vn_Ref_sub_err = sub_err(Vn_Ref_err[i], Vn_Ref_low_err, scale_ref);
      double Vn_D0_sub_err = sub_err(Vn_D0_err[i], Vn_low_D0_err, scale_ref);
      double vn_sub = Vn_D0_sub / sqrt(Vn_Ref_sub);
      double vn_sub_err = fabs(vn_sub)*sqrt(pow(Vn_D0_sub_err/Vn_D0_sub, 2) + pow(Vn_Ref_sub_err/Vn_Ref_sub/2., 2));

      double ratio = Jets_D0[i]/Jets_low_D0;
      double ratio_err = ratio * sqrt(pow(Jets_D0_err[i]/Jets_D0[i], 2) + pow(Jets_low_D0_err/Jets_low_D0, 2));

      cout << Vn_D0_sub / sqrt(Vn_Ref_sub) << "+/-" << vn_sub_err << endl;
      ofile << vn_sub << "," << vn_sub_err <<  ","
         << Vn_D0[i]<< "," << Vn_D0_err[i] << ","
         << Vn_low_D0<< "," << Vn_low_D0_err << ","
         << Vn_Ref[i]<< "," << Vn_Ref_err[i] << ","
         << Vn_Ref_low<< "," << Vn_Ref_low_err << ","
         << Nass_D0[i] << "," << Nass_low_D0 << ","
         << Nass_Ref[i] << "," << Nass_Ref_low << ","
         << Jets_D0[i] << "," << Jets_low_D0<< ","
         << Jets_Ref[i] << "," << Jets_Ref_low << ","
         << ratio << "," << ratio_err << ","
         << endl;

      gratio.SetPoint(i, i+1, ratio);
      gratio.SetPointError(i, 0, ratio_err);
   }
   TFile fout(Form("%s_jet_ratio.root", input_Vn), "recreate");
   gratio.Write();
}

inline double sub(double Vn, double Vn_low, double scale)
{
   return Vn - Vn_low*scale;
}

inline double sub_err(double Vn_err, double Vn_low_err, double scale)
{
   return sqrt(Vn_err*Vn_err + Vn_low_err*scale * Vn_low_err*scale);
}

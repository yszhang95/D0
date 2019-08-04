void fitNass(const char* input, const char* output, const char* dataset)
{
   TFile* file0 = new TFile(input);
   TFile fout(output, "recreate");
   auto vec_trkbin = ana::get_Mult_Edges(dataset);
   if(!vec_trkbin.size()) {
      std::cerr << "wrong dataset" << std::endl;
      return;
   }

   const unsigned int n_trk_bin_ = vec_trkbin.size() - 1;

   TGraphErrors* gNass[n_trk_bin_];
   TGraphErrors* gNass_low;
   for(int i=0; i<n_trk_bin_; i++){
      file0->GetObject(Form("NassvsMass_trk%d", i), gNass[i]);
   }
   file0->GetObject("Nass_lowvsMass", gNass_low);

   double Nass[n_trk_bin_];
   double Nass_low;
   
   for(int i=0; i<n_trk_bin_; i++){
      TF1 f("f", "[0]", 1.7, 2.0);
      gNass[i]->Fit(&f, "Q R");
      gNass[i]->Fit(&f, "Q R");
      gNass[i]->Fit(&f, "");
      Nass[i] = f.GetParameter(0);
   }

   if(true){
      TF1 f("f", "[0]", 1.7, 2.0);
      gNass_low->Fit(&f, "Q R");
      gNass_low->Fit(&f, "Q R");
      gNass_low->Fit(&f, "");
      Nass_low = f.GetParameter(0);
   }

   TGraphErrors* g = new TGraphErrors(n_trk_bin_);
   TGraphErrors* g_low = new TGraphErrors(1);
   for(int j=0; j<n_trk_bin_; j++){
      g->SetPoint(j, j, Nass[j]);
   }

   g_low->SetPoint(0, 1, Nass_low);
   fout.cd();
   g->Write("nass");
   g_low->Write("nass_low");
   delete g;
   delete g_low;
   delete file0;
}

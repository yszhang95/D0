const double deltaEtaBound = 1;

TH1D* proj1D_longrange(TH2*, TH2*, const char*);
TH1D* proj1D_shortrange(TH2*, TH2*);

TF1 draw1D_longrange(TH1*, const char*, 
      const char*, const char* , const char* );

void proj1D_v2vsNtrk_PD0_Process(const char* input_d0= "",
      const char* input_ref = "", 
const string dataset="", const char* input_d0_low_mult="", 
const float pTMin=0., const float pTMax=0., 
const float yMin =0., const float yMax =0.)
{

   //gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   const int n_trk_bin_ = ana::Get_N_nTrkBin(dataset);
   if(n_trk_bin_<0){
      std::cerr << "wrong dataset naem" << std::endl;
   }

   // variables
   double v2_PD0[ana::nMass][n_trk_bin_];
   double v2_PD0_err[ana::nMass][n_trk_bin_];

   double V2_PD0[ana::nMass][n_trk_bin_];
   double V2_PD0_err[ana::nMass][n_trk_bin_];
   double V2_REF[n_trk_bin_];
   double V2_REF_err[n_trk_bin_];

   double N_ass[ana::nMass][n_trk_bin_];
   double N_ass_err[ana::nMass][n_trk_bin_];
   double N_ass_low[ana::nMass];
   double N_ass_low_err[ana::nMass];

   double yields_jet[ana::nMass][n_trk_bin_];
   double yields_jet_err[ana::nMass][n_trk_bin_];
   double yields_jet_low[ana::nMass];

   double V2_PD0_low[ana::nMass];
   double V2_PD0_low_err[ana::nMass];

   double V2_Sub_PD0[ana::nMass][n_trk_bin_];
   double V2_Sub_PD0_err[ana::nMass][n_trk_bin_];

   double v2_Sub_PD0[ana::nMass][n_trk_bin_];
   double v2_Sub_PD0_err[ana::nMass][n_trk_bin_];

   // mass
   TH3D* hDcaVsMassAndMva[n_trk_bin_]; 

   // short range
   TH1D* hMult_ass[n_trk_bin_];
   TH1D* hMult_ass_low = nullptr;
   TH2D* h2DSignal_D0_low[ana::nMass];
   TH2D* h2DBackground_D0_low[ana::nMass];

   // long range
   TH2D* h2DSignal_D0[ana::nMass][n_trk_bin_];
   TH2D* h2DBackground_D0[ana::nMass][n_trk_bin_];
   TH1D* hMult_raw_D0[ana::nMass][n_trk_bin_]; // wrong normalized constant
   TH1D* hMass_D0[ana::nMass][n_trk_bin_];

   // read files
   TFile* f1 = new TFile(input_d0);
   TFile* f2 = new TFile(input_ref);
   TFile* f3 = nullptr;

   // read short range
   if(dataset != "PAMB" || dataset != "PPMB"){
      f3 = new TFile(input_d0_low_mult);
      if(!f3->IsOpen()) {
         cout << "not found " << input_d0_low_mult << endl;
         return;
      }
      f3->GetObject("hMult_ass_trk0", hMult_ass_low);
      if(!hMult_ass_low) {
      }
      for(int imass=0; imass<ana::nMass; imass++){
         f3->GetObject(Form("hSignal_mass%d_trk0", imass), h2DSignal_D0_low[imass]);
         f3->GetObject(Form("hBackground_mass%d_trk0", imass), h2DBackground_D0_low[imass]);
         if(!h2DSignal_D0_low[imass]) {
            cout << "not found hSignal_low" << endl;
            return;
         }
         if(!h2DBackground_D0_low[imass]) {
            cout << "not found hBackground_low" << endl;
            return;
         }
      }
   }else{
      f1->GetObject("hMult_ass_trk0", hMult_ass_low);
      if(!hMult_ass_low) {
      }
      for(int imass=0; imass<ana::nMass; imass++){
         f1->GetObject(Form("hSignal_mass%d_trk0", imass), h2DSignal_D0_low[imass]);
         f1->GetObject(Form("hBackground_mass%d_trk0", imass), h2DBackground_D0_low[imass]);
         if(!h2DSignal_D0_low[imass]) {
            cout << "not found" << endl;
            return;
         }
         if(!h2DBackground_D0_low[imass]) {
            cout << "not found" << endl;
            return;
         }
      }
   }

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      f1->GetObject(Form("hMult_ass_trk%d", i_trk_bin_), hMult_ass[i_trk_bin_]);
      if(!hMult_ass[i_trk_bin_]) {
      }
   }

   // read long range
   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         h2DSignal_D0[imass][i_trk_bin_] = (TH2D*) f1->Get(Form("hSignal_mass%d_trk%d", imass, i_trk_bin_));
         h2DBackground_D0[imass][i_trk_bin_] = (TH2D*) f1->Get(Form("hBackground_mass%d_trk%d", imass, i_trk_bin_));
         hMult_raw_D0[imass][i_trk_bin_] = (TH1D*) f1->Get(Form("hMult_raw_D0_mass%d_trk%d", imass, i_trk_bin_));
         hMass_D0[imass][i_trk_bin_] = (TH1D*) f1->Get(Form("hMassD0_mass%d_trk%d", imass, i_trk_bin_));
      }
   }

   // read mass
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      hDcaVsMassAndMva[i_trk_bin_] = (TH3D*) f1->Get(Form("hDcaVsMassAndMva_trk%d", i_trk_bin_));
   }

   // scaled by event number
   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         double nMult_D0= hMult_raw_D0[imass][i_trk_bin_]->GetEntries()-hMult_raw_D0[imass][i_trk_bin_]->GetBinContent(1);
         h2DSignal_D0[imass][i_trk_bin_]->Scale(1./nMult_D0);
      }
   }

   // long range
   TH1D* hNeg[ana::nMass][n_trk_bin_];
   TH1D* hPos[ana::nMass][n_trk_bin_];
   TH1D* hDeltaPhi[ana::nMass][n_trk_bin_];
   TH1D* hDeltaPhi_low[ana::nMass];

   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         hDeltaPhi[imass][i_trk_bin_] = proj1D_longrange(h2DSignal_D0[imass][i_trk_bin_], h2DBackground_D0[imass][i_trk_bin_], Form("deltaPhi_mass%d_trk%d", imass, i_trk_bin_));
      }

      hDeltaPhi_low[imass] = proj1D_longrange(h2DSignal_D0_low[imass], h2DBackground_D0_low[imass], Form("deltaPhi_mass%d_low", imass));
   }

   TH2D* h2DSignal_Ref[n_trk_bin_];
   TH2D* h2DBackground_Ref[n_trk_bin_];
   TH1D* hDeltaPhi_Ref[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      h2DSignal_Ref[i_trk_bin_] = (TH2D*) f2->Get(Form("hSignal_Ref_trk%d", i_trk_bin_));
      h2DBackground_Ref[i_trk_bin_] = (TH2D*) f2->Get(Form("hBackground_Ref_trk%d", i_trk_bin_));
      TH1D* hMult = (TH1D*) f2->Get(Form("hMult_trk%d", i_trk_bin_));
      long int nevents = hMult->GetEntries();

      hDeltaPhi_Ref[i_trk_bin_] = proj1D_longrange(h2DSignal_Ref[i_trk_bin_], h2DBackground_Ref[i_trk_bin_], Form("deltaPhi_Ref_trk%d", i_trk_bin_));
      hDeltaPhi_Ref[i_trk_bin_]->Scale(1./nevents);
   }
   

   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";
   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         std::string str = input_d0;
         auto found = str.find("/");
         while(found!=std::string::npos){
            str.replace(found, 1, "_");
            found = str.find("/");
         }
         while(str.at(0) == '.'){
            str.erase(0, 1);
         }
         while(str.at(0) == '_'){
            str.erase(0, 1);
         }
         std::string ntrk;
         auto trkRange = ana::get_Mult_Edges(dataset);
         if(trkRange.size() && trkRange[i_trk_bin_+1]!=std::numeric_limits<unsigned int>::max()) 
            ntrk = std::string(Form("%u#leqN_{trk}^{offline}<%u", trkRange[i_trk_bin_], trkRange[i_trk_bin_+1]));
         if(trkRange.size() && trkRange[i_trk_bin_+1]==std::numeric_limits<unsigned int>::max()) 
            ntrk = std::string(Form("N_{trk}^{offline}>%u", trkRange[i_trk_bin_]));
         auto func = draw1D_longrange(hDeltaPhi[imass][i_trk_bin_],
               Form("../plots/v2vsNtrk/y%.1f/%s/%s_deltaPhi_mass%d_trk%d.png", 
                  yMax, dataset.c_str(), str.c_str(), imass, i_trk_bin_),
               ntrk.c_str(),
               Form("%.1f<p_{T}<%.1fGeV, %.1f<y<%.1f", pTMin, pTMax, yMin, yMax), 
               Form("%.3f<mass<%.3fGeV", ana::massbin[imass], ana::massbin[imass+1])
               );
         V2_PD0[imass][i_trk_bin_] = func.GetParameter(2);
         V2_PD0_err[imass][i_trk_bin_] = func.GetParError(2);
         N_ass[imass][i_trk_bin_] = func.GetParameter(0);
         N_ass_err[imass][i_trk_bin_] = func.GetParError(0);
      }
   }

   for(int imass=0; imass<ana::nMass; imass++){
      std::string str = input_d0;
      auto found = str.find("/");
      while(found!=std::string::npos){
         str.replace(found, 1, "_");
         found = str.find("/");
      }
      while(str.at(0) == '.'){
         str.erase(0, 1);
      }
      while(str.at(0) == '_'){
         str.erase(0, 1);
      }
      std::string ntrk;
      auto func = draw1D_longrange(hDeltaPhi_low[imass],
            Form("../plots/v2vsNtrk/y%.1f/%s/%s_deltaPhi_mass%d_low.png", 
               yMax, dataset.c_str(), str.c_str(), imass),
            "N_{trk}^{offline} < 35",
            Form("%.1f<p_{T}<%.1fGeV, %.1f<y<%.1f", pTMin, pTMax, yMin, yMax), 
            Form("%.3f<mass<%.3fGeV", ana::massbin[imass], ana::massbin[imass+1])
            );
      V2_PD0_low[imass] = func.GetParameter(2);
      V2_PD0_low_err[imass] = func.GetParError(2);
      N_ass_low[imass] = func.GetParameter(0);
      N_ass_low_err[imass] = func.GetParError(0);
   }

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      std::string str = input_d0;
      auto found = str.find("/");
      while(found!=std::string::npos){
         str.replace(found, 1, "_");
         found = str.find("/");
      }
      while(str.at(0) == '.'){
         str.erase(0, 1);
      }
      while(str.at(0) == '_'){
         str.erase(0, 1);
      }
      std::string ntrk;
      auto trkRange = ana::get_Mult_Edges(dataset);
      if(trkRange.size() && trkRange[i_trk_bin_+1]!=std::numeric_limits<unsigned int>::max()) 
         ntrk = std::string(Form("%u#leqN_{trk}^{offline}<%u", trkRange[i_trk_bin_], trkRange[i_trk_bin_+1]));
      if(trkRange.size() && trkRange[i_trk_bin_+1]==std::numeric_limits<unsigned int>::max()) 
            ntrk = std::string(Form("N_{trk}^{offline}>%u", trkRange[i_trk_bin_]));
      auto func_ref = draw1D_longrange(hDeltaPhi_Ref[i_trk_bin_],
            Form("../plots/v2vsNtrk/%s_ref_deltaPhi_trk%d.png", 
            str.c_str(), i_trk_bin_),
            ntrk.c_str(),
            "0.3<p_{T}<3.0GeV |#eta|<2.4", 
            ""
            );
      V2_REF[i_trk_bin_] = func_ref.GetParameter(2);
      V2_REF_err[i_trk_bin_] = func_ref.GetParError(2);
   }

   TGraphErrors* g_v2_[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_v2_[i_trk_bin_] = new TGraphErrors(ana::nMass);
      g_v2_[i_trk_bin_]->SetName(Form("g_v2_trk%d", i_trk_bin_));
      for(int imass=0; imass<ana::nMass; imass++){
         double temp = V2_PD0[imass][i_trk_bin_]/ sqrt(V2_REF[i_trk_bin_]);
         double temp_err = temp* sqrt(pow(V2_PD0_err[imass][i_trk_bin_]/ V2_PD0[imass][i_trk_bin_], 2) 
               + pow(0.5*V2_REF_err[i_trk_bin_]/ V2_REF[i_trk_bin_], 2));

         v2_PD0[imass][i_trk_bin_] = temp;
         v2_PD0_err[imass][i_trk_bin_] = temp_err;

         g_v2_[i_trk_bin_]->SetPoint(imass, hMass_D0[imass][i_trk_bin_]->GetMean(), v2_PD0[imass][i_trk_bin_]);
         g_v2_[i_trk_bin_]->SetPointError(imass, 0, v2_PD0_err[imass][i_trk_bin_]);
      }
   }

   TH1D* hMass[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      hMass[i_trk_bin_] = hDcaVsMassAndMva[i_trk_bin_]->ProjectionX(Form("hmass_trk%d", i_trk_bin_));
   }

   TH1D* hNtrk[n_trk_bin_];
   TH1D* hKET[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      hNtrk[i_trk_bin_] = (TH1D*) f1->Get(Form("hNtrk_trk%d", i_trk_bin_));
      hKET[i_trk_bin_] = (TH1D*) f1->Get(Form("hKET_trk%d", i_trk_bin_));
   }

   TFile fout_v2(Form("%s_v2.root", input_d0), "recreate");

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_v2_[i_trk_bin_]->Write();
      hMass[i_trk_bin_]->Write();
      hNtrk[i_trk_bin_]->Write();
      hKET[i_trk_bin_]->Write();
   }
   fout_v2.Close();
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      delete g_v2_[i_trk_bin_];
      delete hDeltaPhi_Ref[i_trk_bin_];
      for(int imass=0; imass<ana::nMass; imass++){
         delete hDeltaPhi[imass][i_trk_bin_];
      }
   }

   // jet correlation
   int negBin_Short = h2DSignal_D0[0][0]->GetXaxis()->FindBin(-1.) ;
   int posBin_Short = h2DSignal_D0[0][0]->GetXaxis()->FindBin(1.) ;

   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         TH1D* hDeltaPhi_jet = h2DSignal_D0[imass][i_trk_bin_]->ProjectionY(
                  Form("deltaPhi_jet_mass%d_trk%d", imass, i_trk_bin_), negBin_Short, posBin_Short);

         TH1D* temp = h2DBackground_D0[imass][i_trk_bin_]->ProjectionY(
               Form("temp_mass%d_trk%d_temp", imass, i_trk_bin_), negBin_Short, posBin_Short);

         int center = h2DBackground_D0[imass][i_trk_bin_]->FindBin(0., 0.);
         hDeltaPhi_jet->Divide(temp);
         hDeltaPhi_jet->Scale(h2DBackground_D0[imass][i_trk_bin_]->GetBinContent(center) / temp->GetBinWidth(1));

         yields_jet[imass][i_trk_bin_] = hDeltaPhi_jet->IntegralAndError(0, 10000, yields_jet_err[imass][i_trk_bin_], "width");

         delete hDeltaPhi_jet;
         delete temp;
      }
      TH1D* hDeltaPhi_jet_low = h2DSignal_D0_low[imass]->ProjectionY(
               Form("deltaPhi_jet_low_mass%d_trk0", imass), negBin_Short, posBin_Short);
      TH1D* temp = h2DBackground_D0_low[imass]->ProjectionY(
               Form("temp_mass%d_trk0", imass), negBin_Short, posBin_Short);

      int center = h2DBackground_D0_low[imass]->FindBin(0., 0.);
      hDeltaPhi_jet_low->Divide(temp);
      hDeltaPhi_jet_low->Scale(h2DBackground_D0_low[imass]->GetBinContent(center) / temp->GetBinWidth(1));

      yields_jet_low[imass] = hDeltaPhi_jet_low->Integral("width");
      delete hDeltaPhi_jet_low;
      delete temp;
   }

   TGraphErrors* g_v2_sub[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_v2_sub[i_trk_bin_] = new TGraphErrors(ana::nMass);
      g_v2_sub[i_trk_bin_]->SetName(Form("g_v2_trk%d", i_trk_bin_));
      for(int imass=0; imass<ana::nMass; imass++){
         V2_Sub_PD0[imass][i_trk_bin_] = V2_PD0[imass][i_trk_bin_] - V2_PD0_low[imass]*N_ass_low[imass]/yields_jet_low[imass] / N_ass[imass][i_trk_bin_] * yields_jet[imass][i_trk_bin_];
         V2_Sub_PD0_err[imass][i_trk_bin_] = sqrt(pow(V2_PD0_err[imass][i_trk_bin_], 2) - 
               pow(V2_PD0_low_err[imass]*N_ass_low[imass]/yields_jet_low[imass] / N_ass[imass][i_trk_bin_] * yields_jet[imass][i_trk_bin_], 2));
         double temp = V2_Sub_PD0[imass][i_trk_bin_]/ sqrt(V2_REF[i_trk_bin_]);
         double temp_err = temp* sqrt(pow(V2_Sub_PD0_err[imass][i_trk_bin_]/ V2_Sub_PD0[imass][i_trk_bin_], 2) 
               + pow(0.5*V2_REF[i_trk_bin_]/ V2_REF[i_trk_bin_], 2));

         v2_Sub_PD0[imass][i_trk_bin_] = temp;
         v2_Sub_PD0_err[imass][i_trk_bin_] = temp_err;

         g_v2_sub[i_trk_bin_]->SetPoint(imass, hMass_D0[imass][i_trk_bin_]->GetMean(), v2_Sub_PD0[imass][i_trk_bin_]);
         g_v2_sub[i_trk_bin_]->SetPointError(imass, 0, v2_Sub_PD0_err[imass][i_trk_bin_]);
      }
   }
   TFile fout_v2_sub(Form("%s_v2_sub.root", input_d0), "recreate");

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_v2_sub[i_trk_bin_]->Write();
      hMass[i_trk_bin_]->Write();
      hNtrk[i_trk_bin_]->Write();
      hKET[i_trk_bin_]->Write();
   }
   fout_v2_sub.Close();

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      delete g_v2_sub[i_trk_bin_];
   }

   // ugly
   /*
   TGraphErrors* g_jets[3];
   for(int i=0; i<3; i++){
      cout << "test" << endl;
      g_jets[i] = new TGraphErrors(14);
      g_jets[i]->SetName(Form("g_v2_trk%d", i));
      for(int j=0; j<14; j++){
         g_jets[i]->SetPoint(i, hMass_D0[j][i]->GetMean(), yields_jet[j][i]);
         g_jets[i]->SetPointError(i, 0, yields_jet_err[j][i]);
      }
      g_jets[i]->Write();
      hMass[i]->Write();
      hNtrk[i]->Write();
      hKET[i]->Write();
      delete g_jets[i];
   }


   fout_consts.Close();
   */

   delete f1;
   delete f2;
   if(f3) delete f3;
   
   return;
}

TH1D* proj1D_longrange(TH2* h2DSignal, TH2* h2DBackground, const char* name)
{
   int negBinMin = 0;
   int negBinMax = h2DSignal->GetXaxis()->FindBin(-1.* deltaEtaBound)-1 ;
   int posBinMin = h2DSignal->GetXaxis()->FindBin(1.* deltaEtaBound)+1 ;
   int posBinMax = h2DSignal->GetXaxis()->GetNbins()+1;
   TH1D* hNeg = h2DSignal->ProjectionY("hneg", negBinMin, negBinMax);
   TH1D* hPos = h2DSignal->ProjectionY("hpos", posBinMin, posBinMax);
   hNeg->Add(hPos);

   TH1D* temp_neg = h2DBackground->ProjectionY("hneg_bkg", negBinMin, negBinMax);
   TH1D* temp_pos = h2DBackground->ProjectionY("hpos_bkg", posBinMin, posBinMax);
   temp_neg->Add(temp_pos);

   //int center = h2DBackground->FindBin(0., 0.);
   int center = temp_neg->FindBin(0.);
   hNeg->Divide(temp_neg);
   //hNeg->Scale(h2DBackground->GetBinContent(center) / temp_neg->GetBinWidth(1));
   hNeg->Scale(temp_neg->GetBinContent(center) / temp_neg->GetBinWidth(1));

   delete hPos;
   delete temp_neg;
   delete temp_pos;

   hNeg->SetName(name);

   return hNeg;
}

TF1 draw1D_longrange(TH1* hDeltaPhi, const char* name, 
      const char* cut1, const char* cut2, const char* cut3)
{
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";
   TCanvas* c_deltaPhi = new TCanvas("c_deltaPhi", "", 550, 450);
   c_deltaPhi->SetBottomMargin(0.14);
   TF1 func("deltaPhi", function.c_str(), -3.14159*0.5, 3.14159*1.5);
   func.SetParameter(0, hDeltaPhi->GetMaximum());
   func.SetParameter(1, 0.1);
   func.SetParameter(2, 0.1);
   func.SetParameter(3, 0.1);

   hDeltaPhi->SetMarkerStyle(20);

   hDeltaPhi->Fit(&func, "q");
   hDeltaPhi->Fit(&func, "q");
   hDeltaPhi->Fit(&func, "m q");
   hDeltaPhi->Fit(&func, "m q E");
   auto fitResult = hDeltaPhi->Fit(&func, "m S E q");

   hDeltaPhi->SetTitle(";#Delta#phi;");
   hDeltaPhi->GetXaxis()->CenterTitle();
   hDeltaPhi->GetXaxis()->SetTitleSize(0.05);

   hDeltaPhi->Draw();
   gPad->Update();
   TPaveStats* pave = (TPaveStats*) hDeltaPhi->FindObject("stats");
   pave->SetX1NDC(0.16);
   pave->SetX2NDC(0.56);
   gPad->Modified();
   gPad->Update();

   TLatex ltx;
   ltx.SetTextSize(0.035);
   ltx.SetTextFont(42);
   ltx.DrawLatexNDC(0.58, 0.42, cut1);
   ltx.DrawLatexNDC(0.58, 0.34, cut2);
   ltx.DrawLatexNDC(0.58, 0.26, cut3);

   c_deltaPhi->Print(name);

   delete c_deltaPhi;
   return func;
}

TH1D* proj1D_shortrange(TH2* h2DSignal, TH2* h2DBackground)
{
   return nullptr;
}

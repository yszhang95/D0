// read low multiplicity files, histograms


const double deltaEtaBound = 1;

TH1D* proj1D_longrange(TH2*, TH2*, const char*);
std::pair<double, double> proj1D_shortrange_yields(TH2*, TH2*, const char*, TCanvas* c, const int& ipad, TH1D*);
std::pair<double, double> proj1D_longrange_yields(TH2*, TH2*, const char*, TCanvas* c, const int& ipad, TH1D*);

std::pair<double, double> calYields(const pair<double, double>&, const pair<double,double>&);

TF1 draw1D_longrange(TH1*,
      const char*, const char* , const char*, TCanvas*c, const int& ipad);

string ouput_prefix(const char* name);

void proj1D_v2vsNtrk_PD0_Process(const char* input_d0= "",
      const char* input_ref = "", 
const string dataset="", 
const char* input_d0_low_mult="", 
const char* input_ref_low_mult="",
const float pTMin=0., const float pTMax=0., 
const float yMin =0., const float yMax =0.)
{

   //gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   string str = ouput_prefix(input_d0);

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

   double V2_PD0_low[ana::nMass];
   double V2_PD0_low_err[ana::nMass];
   double V2_REF_low;
   double V2_REF_low_err;

   double N_ass[ana::nMass][n_trk_bin_];
   double N_ass_err[ana::nMass][n_trk_bin_];
   double N_ass_ref[n_trk_bin_];
   double N_ass_ref_err[n_trk_bin_];

   double N_ass_low[ana::nMass];
   double N_ass_low_err[ana::nMass];
   double N_ass_ref_low;
   double N_ass_ref_low_err;

   double yields_jet[ana::nMass][n_trk_bin_];
   double yields_jet_err[ana::nMass][n_trk_bin_];
   double yields_jet_ref[n_trk_bin_];
   double yields_jet_ref_err[n_trk_bin_];

   double yields_jet_low[ana::nMass];
   double yields_jet_low_err[ana::nMass];
   double yields_jet_ref_low;
   double yields_jet_ref_low_err;


   double V2_Sub_PD0[ana::nMass][n_trk_bin_];
   double V2_Sub_PD0_err[ana::nMass][n_trk_bin_];
   double V2_Sub_REF[n_trk_bin_];
   double V2_Sub_REF_err[n_trk_bin_];

   double v2_Sub_PD0[ana::nMass][n_trk_bin_];
   double v2_Sub_PD0_err[ana::nMass][n_trk_bin_];

   // mass
   TH3D* hDcaVsMassAndMva[n_trk_bin_]; 

   // low multiplicity 
   TH2D* h2DSignal_D0_low[ana::nMass];
   TH2D* h2DBackground_D0_low[ana::nMass];
   TH1D* hMult_raw_D0_low[ana::nMass];
   TH1D* hMass_D0_low[ana::nMass];

   TH2D* h2DSignal_Ref_low;
   TH2D* h2DBackground_Ref_low;
   TH1D* hMult_Ref_low;

   // higher multiplicity
   TH2D* h2DSignal_D0[ana::nMass][n_trk_bin_];
   TH2D* h2DBackground_D0[ana::nMass][n_trk_bin_];
   TH1D* hMult_raw_D0[ana::nMass][n_trk_bin_]; 
   TH1D* hMass_D0[ana::nMass][n_trk_bin_];

   TH2D* h2DSignal_Ref[n_trk_bin_];
   TH2D* h2DBackground_Ref[n_trk_bin_];
   TH1D* hMult_Ref[n_trk_bin_];

   // open files
   TFile* f1 = new TFile(input_d0);
   TFile* f2 = new TFile(input_ref);
   TFile* f3 = nullptr;
   TFile* f4 = nullptr;


   // open low multiplicity
   if(dataset != "PAMB" || dataset != "PPMB"){
      f3 = new TFile(input_d0_low_mult);
      if(!f3->IsOpen()) {
         cout << "not found " << input_d0_low_mult << endl;
         return;
      }
      cout << "opened " << input_d0_low_mult << endl;
      f4 = new TFile(input_ref_low_mult);
      if(!f4->IsOpen()) {
         cout << "not found " << input_ref_low_mult << endl;
         return;
      }
      cout << "opened " << input_ref_low_mult << endl;
   }else{
      f3 = f1;
      f4 = f2;
   }

   // read low multiplicity histograms
   for(int imass=0; imass<ana::nMass; imass++){
      TH2D* htemp_s;
      TH2D* htemp_b;
      TH1D* htemp_mult;
      TH1D* htemp_mass;
      f3->GetObject(Form("hSignal_mass%d_trk0", imass), htemp_s);
      f3->GetObject(Form("hBackground_mass%d_trk0", imass), htemp_b);
      f3->GetObject(Form("hMult_raw_D0_mass%d_trk0", imass), htemp_mult);
      f3->GetObject(Form("hMassD0_mass%d_trk0", imass), htemp_mass);
      if(!htemp_s) {
         cout << "not found hSignal_low" << endl;
         return;
      }
      if(!htemp_b) {
         cout << "not found hBackground_low" << endl;
         return;
      }
      if(!htemp_mult) {
         cout << "not found hmult_low" << endl;
         return;
      }
      if(!htemp_mass) {
         cout << "not found hmass_low" << endl;
         return;
      }
      h2DSignal_D0_low[imass] = new TH2D(*htemp_s);
      h2DSignal_D0_low[imass]->SetName(Form("hSignal_low_mass%d", imass));
      h2DBackground_D0_low[imass] = new TH2D(*htemp_b);
      h2DBackground_D0_low[imass]->SetName(Form("hBackground_low_mass%d", imass));
      hMult_raw_D0_low[imass] = new TH1D(*htemp_mult);
      hMult_raw_D0_low[imass]->SetName(Form("hMult_raw_D0_low_mass%d", imass));
      hMass_D0_low[imass] = new TH1D(*htemp_mass);
      hMass_D0_low[imass]->SetName(Form("hMass_D0_low_mass%d", imass));
   }
   if(true){
      TH2D* htemp_s;
      TH2D* htemp_b;
      TH1D* htemp_mult;
      f4->GetObject("hSignal_Ref_trk0", htemp_s);
      f4->GetObject("hBackground_Ref_trk0", htemp_b);
      f4->GetObject("hMult_trk0", htemp_mult);
      if(!htemp_s) {
         cout << "not found hSignal_low" << endl;
         return;
      }
      if(!htemp_b) {
         cout << "not found hBackground_low" << endl;
         return;
      }
      if(!htemp_mult) {
         cout << "not found hmult_low" << endl;
         return;
      }
      h2DSignal_Ref_low = new TH2D(*htemp_s);
      h2DSignal_Ref_low->SetName("hSignal_Ref_low");
      h2DBackground_Ref_low = new TH2D(*htemp_b);
      h2DBackground_Ref_low->SetName("hBackground_Ref_low");
      hMult_Ref_low = new TH1D(*htemp_mult);
      hMult_Ref_low->SetName("hMult_Ref_low");
   }
   // scaled by event number
   for(int imass=0; imass<ana::nMass; imass++){
      double nMult_D0= hMult_raw_D0_low[imass]->Integral(2, 100000);
      h2DSignal_D0_low[imass]->Scale(1./nMult_D0);
   }
   if(true){
      auto nevents = hMult_Ref_low->Integral(3, 100000);
      h2DSignal_Ref_low->Scale(1./nevents);
   }

   // read mass, higher multiplicity
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      hDcaVsMassAndMva[i_trk_bin_] = (TH3D*) f1->Get(Form("hDcaVsMassAndMva_trk%d", i_trk_bin_));
   }

   // read higher multiplicity
   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         h2DSignal_D0[imass][i_trk_bin_] = (TH2D*) f1->Get(Form("hSignal_mass%d_trk%d", imass, i_trk_bin_));
         h2DBackground_D0[imass][i_trk_bin_] = (TH2D*) f1->Get(Form("hBackground_mass%d_trk%d", imass, i_trk_bin_));
         hMult_raw_D0[imass][i_trk_bin_] = (TH1D*) f1->Get(Form("hMult_raw_D0_mass%d_trk%d", imass, i_trk_bin_));
         hMass_D0[imass][i_trk_bin_] = (TH1D*) f1->Get(Form("hMassD0_mass%d_trk%d", imass, i_trk_bin_));
      }
   }
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      h2DSignal_Ref[i_trk_bin_] = (TH2D*) f2->Get(Form("hSignal_Ref_trk%d", i_trk_bin_));
      h2DBackground_Ref[i_trk_bin_] = (TH2D*) f2->Get(Form("hBackground_Ref_trk%d", i_trk_bin_));
      hMult_Ref[i_trk_bin_] = (TH1D*) f2->Get(Form("hMult_trk%d", i_trk_bin_));
   }

   // scaled by event number
   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         double nMult_D0= hMult_raw_D0[imass][i_trk_bin_]->Integral(2, 100000);
         h2DSignal_D0[imass][i_trk_bin_]->Scale(1./nMult_D0);
      }
   }
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      auto nevents = hMult_Ref[i_trk_bin_]->Integral(3, 100000);
      h2DSignal_Ref[i_trk_bin_]->Scale(1./nevents);
   }

   // long range
   TH1D* hDeltaPhi[ana::nMass][n_trk_bin_];
   TH1D* hDeltaPhi_low[ana::nMass];

   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         hDeltaPhi[imass][i_trk_bin_] = proj1D_longrange(h2DSignal_D0[imass][i_trk_bin_], h2DBackground_D0[imass][i_trk_bin_], Form("deltaPhi_mass%d_trk%d", imass, i_trk_bin_));
      }
      hDeltaPhi_low[imass] = proj1D_longrange(h2DSignal_D0_low[imass], h2DBackground_D0_low[imass], Form("deltaPhi_mass%d_low", imass));
   }

   TH1D* hDeltaPhi_Ref[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      hDeltaPhi_Ref[i_trk_bin_] = proj1D_longrange(h2DSignal_Ref[i_trk_bin_], h2DBackground_Ref[i_trk_bin_], Form("deltaPhi_Ref_trk%d", i_trk_bin_));
   }
   TH1D* hDeltaPhi_Ref_low;
   if(true){
      hDeltaPhi_Ref_low = proj1D_longrange(h2DSignal_Ref_low, h2DBackground_Ref_low, "deltaPhi_Ref_low");
   }
   
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";

   TCanvas* c_lr[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      c_lr[i_trk_bin_] = new TCanvas(Form("c_lr_%d", i_trk_bin_), "", 550*3, 500*5);
      c_lr[i_trk_bin_]->Divide(3, 5);
   }

   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
         std::string ntrk;
         auto trkRange = ana::get_Mult_Edges(dataset);
         if(trkRange.size() && trkRange[i_trk_bin_+1]!=std::numeric_limits<unsigned int>::max()) 
            ntrk = std::string(Form("%u#leqN_{trk}^{offline}<%u", trkRange[i_trk_bin_], trkRange[i_trk_bin_+1]));
         if(trkRange.size() && trkRange[i_trk_bin_+1]==std::numeric_limits<unsigned int>::max()) 
            ntrk = std::string(Form("N_{trk}^{offline}>%u", trkRange[i_trk_bin_]));
         auto func = draw1D_longrange(hDeltaPhi[imass][i_trk_bin_],
               ntrk.c_str(),
               Form("%.1f<p_{T}<%.1fGeV, %.1f<y<%.1f", pTMin, pTMax, yMin, yMax), 
               Form("%.3f<mass<%.3fGeV", ana::massbin[imass], ana::massbin[imass+1]),
               c_lr[i_trk_bin_], imass+1
               );
         V2_PD0[imass][i_trk_bin_] = func.GetParameter(2);
         V2_PD0_err[imass][i_trk_bin_] = func.GetParError(2);
         N_ass[imass][i_trk_bin_] = func.GetParameter(0);
         N_ass_err[imass][i_trk_bin_] = func.GetParError(0);
      }
   }

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      c_lr[i_trk_bin_]->Print(
               Form("../plots/v2vsNtrk/y%.1f/%s/%s_deltaPhi_trk%d.png", 
                  yMax, dataset.c_str(), str.c_str(), i_trk_bin_)
            );
      c_lr[i_trk_bin_]->Print(
               Form("../plots/v2vsNtrk/y%.1f/%s/%s_deltaPhi_trk%d.pdf", 
                  yMax, dataset.c_str(), str.c_str(), i_trk_bin_)
            );
   }

   TCanvas* c_lr_low = new TCanvas("c_lr_low", "", 550*3, 550*5);
   c_lr_low->Divide(3, 5);

   for(int imass=0; imass<ana::nMass; imass++){
      auto func = draw1D_longrange(hDeltaPhi_low[imass],
            "N_{trk}^{offline} < 35",
            Form("%.1f<p_{T}<%.1fGeV, %.1f<y<%.1f", pTMin, pTMax, yMin, yMax), 
            Form("%.3f<mass<%.3fGeV", ana::massbin[imass], ana::massbin[imass+1]),
            c_lr_low, imass+1
            );
      V2_PD0_low[imass] = func.GetParameter(2);
      V2_PD0_low_err[imass] = func.GetParError(2);
      N_ass_low[imass] = func.GetParameter(0);
      N_ass_low_err[imass] = func.GetParError(0);
   }
   c_lr_low->Print(
      Form("../plots/v2vsNtrk/y%.1f/%s/%s_deltaPhi_low.png", 
      yMax, dataset.c_str(), str.c_str())
         );
   c_lr_low->Print(
      Form("../plots/v2vsNtrk/y%.1f/%s/%s_deltaPhi_low.pdf", 
      yMax, dataset.c_str(), str.c_str())
         );

   TCanvas* c_lr_ref[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      c_lr_ref[i_trk_bin_] = new TCanvas(Form("c_lr_ref_%d", i_trk_bin_), "", 550*3, 500*5);
   }

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      std::string ntrk;
      auto trkRange = ana::get_Mult_Edges(dataset);
      if(trkRange.size() && trkRange[i_trk_bin_+1]!=std::numeric_limits<unsigned int>::max()) 
         ntrk = std::string(Form("%u#leqN_{trk}^{offline}<%u", trkRange[i_trk_bin_], trkRange[i_trk_bin_+1]));
      if(trkRange.size() && trkRange[i_trk_bin_+1]==std::numeric_limits<unsigned int>::max()) 
            ntrk = std::string(Form("N_{trk}^{offline}>%u", trkRange[i_trk_bin_]));
      auto func_ref = draw1D_longrange(hDeltaPhi_Ref[i_trk_bin_],
            ntrk.c_str(),
            "0.3<p_{T}<3.0GeV |#eta|<2.4", 
            "",
            c_lr_ref[i_trk_bin_], 1
            );
      V2_REF[i_trk_bin_] = func_ref.GetParameter(2);
      V2_REF_err[i_trk_bin_] = func_ref.GetParError(2);
      N_ass_ref[i_trk_bin_] = func_ref.GetParameter(0);
      N_ass_ref_err[i_trk_bin_] = func_ref.GetParError(0);
   }
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      c_lr_ref[i_trk_bin_]->Print(
            Form("../plots/v2vsNtrk/%s_ref_deltaPhi_trk%d.png", 
            str.c_str(), i_trk_bin_)
         );
      c_lr_ref[i_trk_bin_]->Print(
            Form("../plots/v2vsNtrk/%s_ref_deltaPhi_trk%d.pdf", 
            str.c_str(), i_trk_bin_)
         );
   }

   TCanvas* c_lr_ref_low = new TCanvas;
   if(true){
      auto func_ref = draw1D_longrange(hDeltaPhi_Ref_low,
            "N_{trk}^{offline} < 35",
            "0.3<p_{T}<3.0GeV |#eta|<2.4", 
            "",
            c_lr_ref_low, 1
            );
      V2_REF_low = func_ref.GetParameter(2);
      V2_REF_low_err= func_ref.GetParError(2);
      N_ass_ref_low = func_ref.GetParameter(0);
      N_ass_ref_low_err= func_ref.GetParError(0);
   }
   c_lr_ref_low->Print(
            Form("../plots/v2vsNtrk/%s_ref_deltaPhi_low.png", 
            str.c_str())
         );
   c_lr_ref_low->Print(
            Form("../plots/v2vsNtrk/%s_ref_deltaPhi_low.pdf", 
            str.c_str())
         );

   // started filling
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
      hMass[i_trk_bin_] = (TH1D*) f1->Get(Form("hMass_%d", i_trk_bin_));
      hMass[i_trk_bin_]->SetName(Form("hmass_trk%d", i_trk_bin_));
   }

   TH1D* hMass_low;
   if(true){
      TH1D* temp = (TH1D*)f3->Get("hMass_0");
      hMass_low = new TH1D(*temp);
      hMass_low->SetName("hMass_low");
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
   }

   // fix jet yields ratio and then subtract v2
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      V2_Sub_REF[i_trk_bin_] = V2_REF[i_trk_bin_] - V2_REF_low*N_ass_ref_low/ N_ass_ref[i_trk_bin_] * 1.5;
      double temp_2nd_err = V2_REF_low*N_ass_ref_low/ N_ass_ref[i_trk_bin_]*1.5 * sqrt(
               pow(V2_REF_low_err/V2_REF_low, 2) + pow(N_ass_ref_low_err/N_ass_ref_low, 2) + pow(N_ass_ref_err[i_trk_bin_]/N_ass_ref[i_trk_bin_], 2)
               );
      V2_Sub_REF_err[i_trk_bin_] = sqrt(
               pow(V2_REF_err[i_trk_bin_], 2) + pow(temp_2nd_err, 2)
            );
      for(int imass=0; imass<ana::nMass; imass++){
         V2_Sub_PD0[imass][i_trk_bin_] = V2_PD0[imass][i_trk_bin_] - V2_PD0_low[imass]*N_ass_low[imass]/ N_ass[imass][i_trk_bin_] * 1.5;

         double temp_2nd_err_pd0 = V2_PD0_low[imass]*N_ass_low[imass]/ N_ass[imass][i_trk_bin_] * 1.5 * sqrt(
               pow(V2_PD0_low_err[imass]/V2_PD0_low[imass] , 2) //+ 
               //pow(N_ass_low_err[imass]/N_ass_low[imass], 2) + 
               //pow(N_ass_err[imass][i_trk_bin_]/N_ass[imass][i_trk_bin_], 2)
               );

         V2_Sub_PD0_err[imass][i_trk_bin_] = sqrt(pow(V2_PD0_err[imass][i_trk_bin_], 2) + pow(temp_2nd_err_pd0, 2));

         double temp = V2_Sub_PD0[imass][i_trk_bin_]/ sqrt(V2_Sub_REF[i_trk_bin_]);
         double temp_err = temp* sqrt(pow(V2_Sub_PD0_err[imass][i_trk_bin_]/ V2_Sub_PD0[imass][i_trk_bin_], 2) 
               + pow(0.5*V2_Sub_REF_err[i_trk_bin_]/ V2_Sub_REF[i_trk_bin_], 2));

         v2_Sub_PD0[imass][i_trk_bin_] = temp;
         v2_Sub_PD0_err[imass][i_trk_bin_] = temp_err;
      }
   }
   TGraphErrors* g_v2_sub_fixed[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_v2_sub_fixed[i_trk_bin_] = new TGraphErrors(ana::nMass);
      g_v2_sub_fixed[i_trk_bin_]->SetName(Form("g_v2_trk%d", i_trk_bin_));
      for(int imass=0; imass<ana::nMass; imass++){
         g_v2_sub_fixed[i_trk_bin_]->SetPoint(imass, hMass_D0[imass][i_trk_bin_]->GetMean(), v2_Sub_PD0[imass][i_trk_bin_]);
         g_v2_sub_fixed[i_trk_bin_]->SetPointError(imass, 0, v2_Sub_PD0_err[imass][i_trk_bin_]);
      }
   }
   TFile fout_v2_sub_fixed(Form("%s_v2_sub_fixed.root", input_d0), "recreate");

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_v2_sub_fixed[i_trk_bin_]->Write();
      hMass[i_trk_bin_]->Write();
      hNtrk[i_trk_bin_]->Write();
      hKET[i_trk_bin_]->Write();
   }
   fout_v2_sub_fixed.Close();

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      delete g_v2_sub_fixed[i_trk_bin_];
   }


   // to fit the yields and then subtract non flow
   // low mult, ref
   TCanvas* c_sry_ref_low = new TCanvas("c_sry_ref_low", "", 550, 550);
   c_sry_ref_low->Divide(1,1);
   TCanvas* c_lry_ref_low = new TCanvas("c_lry_ref_low", "", 550, 550);
   c_lry_ref_low->Divide(1,1);
   TH1D* h_sry_ref_low;
   TH1D* h_lry_ref_low;
   if(true){
      auto sr_pair = proj1D_shortrange_yields(h2DSignal_Ref_low, h2DBackground_Ref_low, 
            Form("../plots/v2vsNtrk/%s_ref_low_sr.png", str.c_str()), c_sry_ref_low, 1, h_sry_ref_low);
      auto lr_pair = proj1D_longrange_yields(h2DSignal_Ref_low, h2DBackground_Ref_low, 
            Form("../plots/v2vsNtrk/%s_ref_low_lr.png", str.c_str()), c_lry_ref_low, 1, h_lry_ref_low);
      auto yields_pair = calYields(sr_pair, lr_pair);
      yields_jet_ref_low = yields_pair.first;
      yields_jet_ref_low_err = yields_pair.second;
   }
   c_sry_ref_low->Print(
            Form("../plots/v2vsNtrk/%s_ref_low_sr.png", str.c_str())
            );
   c_lry_ref_low->Print(
            Form("../plots/v2vsNtrk/%s_ref_low_lr.png", str.c_str())
            );

   TCanvas* c_sry_ref[n_trk_bin_];
   TCanvas* c_lry_ref[n_trk_bin_];
   TH1D* h_sry_ref[n_trk_bin_];
   TH1D* h_lry_ref[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      c_sry_ref[i_trk_bin_] = new TCanvas(Form("c_sry_ref_%d", i_trk_bin_), "", 550, 550);
      c_sry_ref[i_trk_bin_]->Divide(1, 1);
      c_lry_ref[i_trk_bin_] = new TCanvas(Form("c_lry_ref_%d", i_trk_bin_), "", 550, 550);
      c_lry_ref[i_trk_bin_]->Divide(1, 1);
   }
   // higher mult, ref
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      auto sr_pair = proj1D_shortrange_yields(h2DSignal_Ref[i_trk_bin_], h2DBackground_Ref[i_trk_bin_], 
            Form("../plots/v2vsNtrk/%s_ref_deltaPhi_trk%d_sr.png", str.c_str(), i_trk_bin_), 
            c_sry_ref[i_trk_bin_], 1, h_sry_ref[i_trk_bin_]);
      auto lr_pair = proj1D_longrange_yields(h2DSignal_Ref[i_trk_bin_], h2DBackground_Ref[i_trk_bin_], 
            Form("../plots/v2vsNtrk/%s_ref_deltaPhi_trk%d_lr.png", str.c_str(), i_trk_bin_),
            c_lry_ref[i_trk_bin_], 1, h_lry_ref[i_trk_bin_]);
      auto yields_pair = calYields(sr_pair, lr_pair);
      yields_jet_ref[i_trk_bin_]= yields_pair.first;
      yields_jet_ref_err[i_trk_bin_] = yields_pair.second;
   }
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      c_sry_ref[i_trk_bin_]->Print(
            Form("../plots/v2vsNtrk/%s_ref_deltaPhi_trk%d_sr.png", str.c_str(), i_trk_bin_)
            );
      c_lry_ref[i_trk_bin_]->Print(
            Form("../plots/v2vsNtrk/%s_ref_deltaPhi_trk%d_lr.png", str.c_str(), i_trk_bin_)
            );
   }

   // low mult, d0
   TCanvas* c_sry_low = new TCanvas("c_sry_low", "", 550*3, 550*5);
   c_sry_low->Divide(3,5);
   TCanvas* c_lry_low = new TCanvas("c_lry_low", "", 550*3, 550*5);
   c_lry_low->Divide(3,5);
   TH1D* h_sry_low;
   TH1D* h_lry_low;
   for(int imass=0; imass<ana::nMass; imass++){
      string str = ouput_prefix(input_d0);
      auto sr_pair = proj1D_shortrange_yields(h2DSignal_D0_low[imass], h2DBackground_D0_low[imass], 
            Form("../plots/v2vsNtrk/%s_deltaPhi_low_mass%d_sr.png", str.c_str(), imass),
            c_sry_low, imass+1, h_sry_low);
      auto lr_pair = proj1D_longrange_yields(h2DSignal_D0_low[imass], h2DBackground_D0_low[imass], 
            Form("../plots/v2vsNtrk/%s_deltaPhi_low_mass%d_lr.png", str.c_str(), imass),
            c_lry_low, imass+1, h_lry_low);
      auto yields_pair = calYields(sr_pair, lr_pair);
      yields_jet_low[imass]= yields_pair.first;
      yields_jet_low_err[imass] = yields_pair.second;
      cout << yields_jet_low[imass] << endl;
   }
   c_sry_low->Print(
            Form("../plots/v2vsNtrk/%s_deltaPhi_low_sr.png", str.c_str())
         );
   c_lry_low->Print(
            Form("../plots/v2vsNtrk/%s_deltaPhi_low_lr.png", str.c_str())
         );

   // higher mult, d0
   TCanvas* c_sry[n_trk_bin_];
   TCanvas* c_lry[n_trk_bin_];
   TH1D* h_sry[n_trk_bin_];
   TH1D* h_lry[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      c_sry[i_trk_bin_] = new TCanvas(Form("c_sry_low_%d", i_trk_bin_), "", 550*3, 550*5);
      c_sry[i_trk_bin_]->Divide(3, 5);
      c_lry[i_trk_bin_]= new TCanvas(Form("c_lry_low_%d", i_trk_bin_), "", 550*3, 550*5);
      c_lry[i_trk_bin_]->Divide(3, 5);
   }
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      for(int imass=0; imass<ana::nMass; imass++){
         string str = ouput_prefix(input_d0);
         auto sr_pair = proj1D_shortrange_yields(h2DSignal_D0[imass][i_trk_bin_], h2DBackground_D0[imass][i_trk_bin_], 
               Form("../plots/v2vsNtrk/%s_deltaPhi_mass%d_trk%d_sr.png", str.c_str(), imass, i_trk_bin_),
               c_sry[i_trk_bin_], imass+1, h_sry[i_trk_bin_]
               );
         auto lr_pair = proj1D_longrange_yields(h2DSignal_D0[imass][i_trk_bin_], h2DBackground_D0[imass][i_trk_bin_], 
               Form("../plots/v2vsNtrk/%s_deltaPhi_mass%d_trk%d_lr.png", str.c_str(), imass, i_trk_bin_),
               c_lry[i_trk_bin_], imass+1, h_lry[i_trk_bin_]
                  );
         auto yields_pair = calYields(sr_pair, lr_pair);
         yields_jet[imass][i_trk_bin_]= yields_pair.first;
         yields_jet_err[imass][i_trk_bin_] = yields_pair.second;
      }
   }
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      c_sry[i_trk_bin_]->Print(
               Form("../plots/v2vsNtrk/%s_deltaPhi_trk%d_sr.png", str.c_str(), i_trk_bin_)
            );
      c_lry[i_trk_bin_]->Print(
               Form("../plots/v2vsNtrk/%s_deltaPhi_trk%d_lr.png", str.c_str(), i_trk_bin_)
            );
   }

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      for(int imass=0; imass<ana::nMass; imass++){
         cout << yields_jet[imass][i_trk_bin_]/yields_jet_low[imass]<< endl;
      }
   }

   // calculate v_sub
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      V2_Sub_REF[i_trk_bin_] = V2_REF[i_trk_bin_] - V2_REF_low*N_ass_ref_low/ N_ass_ref[i_trk_bin_] * yields_jet_ref[i_trk_bin_]/yields_jet_ref_low;
      double temp_2nd_err = V2_REF_low*N_ass_ref_low/ N_ass_ref[i_trk_bin_]* yields_jet_ref[i_trk_bin_]/yields_jet_ref_low * sqrt(
               pow(V2_REF_low_err/V2_REF_low, 2)// + pow(N_ass_ref_low_err/N_ass_ref_low, 2) + pow(N_ass_ref_err[i_trk_bin_]/N_ass_ref[i_trk_bin_], 2)
               );
      V2_Sub_REF_err[i_trk_bin_] = sqrt(
               pow(V2_REF_err[i_trk_bin_], 2) + pow(temp_2nd_err, 2)
            );
      for(int imass=0; imass<ana::nMass; imass++){
         V2_Sub_PD0[imass][i_trk_bin_] = V2_PD0[imass][i_trk_bin_] - V2_PD0_low[imass]*N_ass_low[imass]/ N_ass[imass][i_trk_bin_] * yields_jet[imass][i_trk_bin_]/yields_jet_low[imass];

         double temp_2nd_err_pd0 = V2_PD0_low[imass]*N_ass_low[imass]/ N_ass[imass][i_trk_bin_] * yields_jet[imass][i_trk_bin_]/yields_jet_low[imass] * sqrt(
               pow(V2_PD0_low_err[imass]/V2_PD0_low[imass] , 2) //+ 
               //pow(N_ass_low_err[imass]/N_ass_low[imass], 2) + 
               //pow(N_ass_err[imass][i_trk_bin_]/N_ass[imass][i_trk_bin_], 2)
               );

         V2_Sub_PD0_err[imass][i_trk_bin_] = sqrt(pow(V2_PD0_err[imass][i_trk_bin_], 2) + pow(temp_2nd_err_pd0, 2));

         double temp = V2_Sub_PD0[imass][i_trk_bin_]/ sqrt(V2_Sub_REF[i_trk_bin_]);
         double temp_err = temp* sqrt(pow(V2_Sub_PD0_err[imass][i_trk_bin_]/ V2_Sub_PD0[imass][i_trk_bin_], 2) 
               + pow(0.5*V2_Sub_REF_err[i_trk_bin_]/ V2_Sub_REF[i_trk_bin_], 2));

         v2_Sub_PD0[imass][i_trk_bin_] = temp;
         v2_Sub_PD0_err[imass][i_trk_bin_] = temp_err;
      }
   }

   TGraphErrors* g_v2_sub[n_trk_bin_];
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_v2_sub[i_trk_bin_] = new TGraphErrors(ana::nMass);
      g_v2_sub[i_trk_bin_]->SetName(Form("g_v2_trk%d", i_trk_bin_));
      for(int imass=0; imass<ana::nMass; imass++){
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
   fout_v2.Close();

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      delete g_v2_sub[i_trk_bin_];
   }

   TFile fraw(Form("%s_raw.root", input_d0), "recreate");

   TGraphErrors* g_V2[n_trk_bin_]; // signal + background, nMass points, need to do fitting, filled
   TGraphErrors* g_V2_Ref[n_trk_bin_]; // signal, 1 point, done
   TGraphErrors* g_V2_low; // signal + background, nMass points, need to do fitting, filled, but to define hMass_D0_low
   TGraphErrors* g_V2_Ref_low; // signal, 1 point, done
   TGraphErrors* g_Jets[n_trk_bin_]; // signal + background, nMass points, need to do fitting, filled
   TGraphErrors* g_Jets_low; // signal + background, nMass points, need to do fitting, filled
   TGraphErrors* g_Jets_Ref[n_trk_bin_]; // signal, 1 point, done
   TGraphErrors* g_Jets_Ref_low; // signal, 1 point, done
   TGraphErrors* g_Nass[n_trk_bin_]; // signal, need to do average, nMass points, filled
   TGraphErrors* g_Nass_low; // signal, need to do average, nMass points, filled
   TGraphErrors* g_Nass_Ref[n_trk_bin_]; // signal, 1 point, done
   TGraphErrors* g_Nass_Ref_low; // signal, 1 point, done

   // high N
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_V2[i_trk_bin_] = new TGraphErrors(ana::nMass);
      g_V2_Ref[i_trk_bin_] = new TGraphErrors(1);
      g_Jets[i_trk_bin_] = new TGraphErrors(ana::nMass);
      g_Jets_Ref[i_trk_bin_] = new TGraphErrors(1);
      g_Nass[i_trk_bin_] = new TGraphErrors(ana::nMass);
      g_Nass_Ref[i_trk_bin_] = new TGraphErrors(1);
   }
   // low N
   g_V2_low = new TGraphErrors(ana::nMass);
   g_V2_Ref_low = new TGraphErrors(1);
   g_Jets_low = new TGraphErrors(ana::nMass);
   g_Jets_Ref_low = new TGraphErrors(1);
   g_Nass_low = new TGraphErrors(ana::nMass);
   g_Nass_Ref_low = new TGraphErrors(1);


   for(int imass=0; imass<ana::nMass; imass++){
      for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      // high N
         g_V2[i_trk_bin_]->SetPoint(imass, hMass_D0[imass][i_trk_bin_]->GetMean(), V2_PD0[imass][i_trk_bin_]);
         g_V2[i_trk_bin_]->SetPointError(imass, 0, V2_PD0_err[imass][i_trk_bin_]);

         g_Jets[i_trk_bin_]->SetPoint(imass, hMass_D0[imass][i_trk_bin_]->GetMean(), yields_jet[imass][i_trk_bin_]);
         g_Jets[i_trk_bin_]->SetPointError(imass, 0, yields_jet_err[imass][i_trk_bin_]);

         g_Nass[i_trk_bin_]->SetPoint(imass, hMass_D0[imass][i_trk_bin_]->GetMean(), N_ass[imass][i_trk_bin_]);
         g_Nass[i_trk_bin_]->SetPointError(imass, 0, N_ass_err[imass][i_trk_bin_]);
      }
      // low N
      g_V2_low->SetPoint(imass, hMass_D0_low[imass]->GetMean(), V2_PD0_low[imass]);
      g_V2_low->SetPointError(imass, 0, V2_PD0_low_err[imass]);

      g_Jets_low->SetPoint(imass, hMass_D0_low[imass]->GetMean(), yields_jet_low[imass]);
      g_Jets_low->SetPointError(imass, 0, yields_jet_low_err[imass]);

      g_Nass_low->SetPoint(imass, hMass_D0_low[imass]->GetMean(), N_ass_low[imass]);
      g_Nass_low->SetPointError(imass, 0, N_ass_low_err[imass]);
   }
   // high N
   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_V2_Ref[i_trk_bin_]->SetPoint(0, 1, V2_REF[i_trk_bin_]);
      g_V2_Ref[i_trk_bin_]->SetPointError(0, 1, V2_REF_err[i_trk_bin_]);

      g_Jets_Ref[i_trk_bin_]->SetPoint(0, 1, yields_jet_ref[i_trk_bin_]);
      g_Jets_Ref[i_trk_bin_]->SetPointError(0, 1, yields_jet_ref_err[i_trk_bin_]);

      g_Nass_Ref[i_trk_bin_]->SetPoint(0, 1, N_ass_ref[i_trk_bin_]);
      g_Nass_Ref[i_trk_bin_]->SetPointError(0, 1, N_ass_ref_err[i_trk_bin_]);
   }
   // low N
   g_V2_Ref_low->SetPoint(0, 1, V2_REF_low);
   g_V2_Ref_low->SetPointError(0, 1, V2_REF_low_err);

   g_Jets_Ref_low->SetPoint(0, 1, yields_jet_ref_low);
   g_Jets_Ref_low->SetPointError(0, 1, yields_jet_ref_low_err);

   g_Nass_Ref_low->SetPoint(0, 1, N_ass_ref_low);
   g_Nass_Ref_low->SetPointError(0, 1, N_ass_ref_low_err);

   cout <<"test" << endl;

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      g_V2[i_trk_bin_]->Write(Form("V2vsMass_trk%d", i_trk_bin_)); // signal + background, nMass points, need to do fitting, filled
      g_V2_Ref[i_trk_bin_]->Write(Form("V2_Ref_trk%d", i_trk_bin_)); // signal, 1 point, done
      g_Jets[i_trk_bin_]->Write(Form("JetsvsMass_trk%d", i_trk_bin_)); // signal + background, nMass points, need to do fitting, filled
      g_Jets_Ref[i_trk_bin_]->Write(Form("Jets_Ref_trk%d", i_trk_bin_)); // signal, 1 point, done
      g_Nass[i_trk_bin_]->Write(Form("NassvsMass_trk%d", i_trk_bin_)); // signal, need to do average, nMass points, filled
      g_Nass_Ref[i_trk_bin_]->Write(Form("Nass_Ref_trk%d", i_trk_bin_)); // signal, 1 point, done

      hMass[i_trk_bin_]->Write();
      hNtrk[i_trk_bin_]->Write();
   }
   g_V2_low->Write("V2_lowvsMass"); // signal + background, nMass points, need to do fitting, filled, but to define hMass_D0_low
   g_V2_Ref_low->Write("V2_Ref_low"); // signal, 1 point, done
   g_Jets_low->Write("Jets_lowvsMass"); // signal + background, nMass points, need to do fitting, filled
   g_Jets_Ref_low->Write("Jets_ref_low"); // signal, 1 point, done
   g_Nass_low->Write("Nass_lowvsMass"); // signal, need to do average, nMass points, filled
   g_Nass_Ref_low->Write("Nass_ref_low"); // signal, 1 point, done

   hMass_low->Write();

   fraw.Close();

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      delete g_V2[i_trk_bin_]; // signal + background, nMass points, need to do fitting, filled
      delete g_V2_Ref[i_trk_bin_]; // signal, 1 point, done
      delete g_Jets[i_trk_bin_]; // signal + background, nMass points, need to do fitting, filled
      delete g_Jets_Ref[i_trk_bin_]; // signal, 1 point, done
      delete g_Nass[i_trk_bin_]; // signal, need to do average, nMass points, filled
      delete g_Nass_Ref[i_trk_bin_]; // signal, 1 point, done
   }
   delete g_V2_low; // signal + background, nMass points, need to do fitting, filled, but to define hMass_D0_low
   delete g_V2_Ref_low; // signal, 1 point, done
   delete g_Jets_low; // signal + background, nMass points, need to do fitting, filled
   delete g_Jets_Ref_low; // signal, 1 point, done
   delete g_Nass_low; // signal, need to do average, nMass points, filled
   delete g_Nass_Ref_low; // signal, 1 point, done

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      delete hDeltaPhi_Ref[i_trk_bin_];
      for(int imass=0; imass<ana::nMass; imass++){
         delete hDeltaPhi[imass][i_trk_bin_];
      }
   }

   for(int imass=0; imass<ana::nMass; imass++){
      delete hDeltaPhi_low[imass];
      delete h2DSignal_D0_low[imass];
      delete h2DBackground_D0_low[imass];
      delete hMass_D0_low[imass];
      delete hMult_raw_D0_low[imass];
   }
   delete hDeltaPhi_Ref_low;
   delete h2DSignal_Ref_low;
   delete h2DBackground_Ref_low;
   delete hMult_Ref_low;

   delete hMass_low;

   for(int i_trk_bin_=0; i_trk_bin_<n_trk_bin_; i_trk_bin_++){
      delete c_lr[i_trk_bin_];
      delete c_lr_ref[i_trk_bin_];
      delete c_sry[i_trk_bin_];
      delete c_lry[i_trk_bin_];
      delete c_sry_ref[i_trk_bin_];
      delete c_lry_ref[i_trk_bin_];
   }
   delete c_lr_low;
   delete c_lr_ref_low;
   delete c_sry_low;
   delete c_lry_low;
   delete c_sry_ref_low;
   delete c_lry_ref_low;

   delete f1;
   delete f2;
   if(f3!=f1) delete f3;
   if(f4!=f2) delete f4;
   
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

   int center = h2DBackground->FindBin(0., 0.);
   hNeg->Divide(temp_neg);
   hNeg->Scale(h2DBackground->GetBinContent(center) / temp_neg->GetBinWidth(1)/h2DBackground->GetXaxis()->GetBinWidth(1));

   delete hPos;
   delete temp_neg;
   delete temp_pos;

   hNeg->SetName(name);

   return hNeg;
}

TF1 draw1D_longrange(TH1* hDeltaPhi,
      const char* cut1, const char* cut2, const char* cut3,
      TCanvas* c, const int& ipad)
{
   std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";
   c->cd(ipad);
   gPad->SetBottomMargin(0.14);
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

   return func;
}

std::pair<double, double> proj1D_shortrange_yields(TH2* h2DSignal, TH2* h2DBackground, const char* name, TCanvas* c, const int& ipad, TH1D* hsig)
{
   c->cd(ipad);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111111);
   gPad->SetLeftMargin(0.15);
   gPad->SetBottomMargin(0.15);
   int negBin = h2DSignal->GetXaxis()->FindBin(-1.);
   int posBin = h2DSignal->GetXaxis()->FindBin(1.);
   hsig = h2DSignal->ProjectionY("hsig", negBin, posBin);
   hsig->SetTitle(";#Delta#phi;dN/d(#Delta#phi)");

   TH1D* temp = h2DBackground->ProjectionY("hbkg", negBin, posBin);

   int center = h2DBackground->FindBin(0., 0.);
   hsig->Divide(temp);
   hsig->Scale(h2DBackground->GetBinContent(center) / temp->GetBinWidth(1)/h2DBackground->GetXaxis()->GetBinWidth(1));

   TF1 fitter("fitter", "[0]*x^2+[1]*x+[2]", 0.6, 2.2);
   fitter.SetParameters(1, 1, 1);

   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->SetMarkerStyle(20);
   hsig->Draw("E P");

   TLatex ltx;
   ltx.DrawLatexNDC(0.2, 0.8, "#font[42]{|#Delta#eta|<1}");

   double min = fitter.GetMinimum(0.6, 2.2);
   double minX = fitter.GetMinimumX(0.6, 2.2);

   TF1 minfunc("minfunc", "[0]", -0.5*TMath::Pi(), 1.5*TMath::Pi());
   minfunc.SetParameter(0, -min);
   hsig->Add(&minfunc);

   double yields_err;
   double yields = hsig->IntegralAndError(hsig->FindBin(0.), hsig->FindBin(minX), yields_err, "width");
   yields = 2*yields - hsig->GetBinContent(hsig->FindBin(0.)) * TMath::Pi()/16.;

   hsig->Add(&minfunc, -1);

   delete temp;

   return std::pair<double, double>(yields, yields_err);
}

std::pair<double, double> proj1D_longrange_yields(TH2* h2DSignal, TH2* h2DBackground, const char* name, TCanvas* c, const int& ipad , TH1D* hsig)
{
   c->cd(ipad);
   //double lw = 0.0;
   double lw = 0.1;
   double up = 2.0;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111111);
   c->SetLeftMargin(0.15);
   c->SetBottomMargin(0.15);

   int negBinMin = 1;
   int negBinMax = 10;
   int posBinMin = 24;
   int posBinMax = 33;

   TH1D* hNeg = h2DSignal->ProjectionY("hneg", negBinMin, negBinMax);
   TH1D* hPos = h2DSignal->ProjectionY("hpos", posBinMin, posBinMax);
   hNeg->Add(hPos);
   hsig = new TH1D(*hNeg);
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

   TF1 fitter("fitter", "[0]*x^2+[1]*x+[2]", lw, up);
   fitter.SetParameters(1, 1, 1);
   fitter.SetParLimits(0, 0, 10);

   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");
   fitter.ReleaseParameter(0);
   hsig->Fit(&fitter, "R q");
   hsig->Fit(&fitter, "R q");

   hsig->SetTitle(";#Delta#phi;dN/d(#Delta#phi)");
   hsig->SetMarkerStyle(20);
   hsig->Draw("E P");

   TLatex ltx;
   ltx.DrawLatexNDC(0.2, 0.8, "#font[42]{|#Delta#eta|>1}");

   double min = fitter.GetMinimum(lw, up);
   double minX = fitter.GetMinimumX(lw, up);

   TF1 minfunc("minfunc", "[0]", -0.5*TMath::Pi(), 1.5*TMath::Pi());
   minfunc.SetParameter(0, -min);
   hsig->Add(&minfunc);

   double yields_err;
   double yields = hsig->IntegralAndError(hsig->FindBin(0.), hsig->FindBin(minX), yields_err, "width");
   yields = 2*yields - hsig->GetBinContent(hsig->FindBin(0.)) * TMath::Pi()/16.;

   hsig->Add(&minfunc, -1);

   delete temp;

   return std::pair<double, double>(yields, yields_err);
}

std::pair<double, double> calYields(const pair<double, double>& sr, const pair<double,double>& lr)
{
   double yields = sr.first - lr.first;
   double yields_err = sqrt(pow(sr.second, 2) + pow(lr.second, 2));
   return pair<double, double>(yields, yields_err);
}

string ouput_prefix(const char* name)
{
   std::string str = name;
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
   return str;
}

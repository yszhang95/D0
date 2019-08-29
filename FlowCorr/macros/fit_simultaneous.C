#ifndef __MY_INCLUDE__
#define __MY_INCLUDE__
R__ADD_INCLUDE_PATH(../include);
R__ADD_INCLUDE_PATH(../src);
#include "myAnaConsts.h"
#include "functions.cxx"
#endif

#ifndef __FIT_SIMULTANEOUS__
#define __FIT_SIMULTANEOUS__
using namespace std;

int iparmassfit_poly3bkg_floatwidth[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
int iparvarfit_poly3bkg_floatwidth[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

struct GlobalChi2_poly3bkg_floatwidth {
    GlobalChi2_poly3bkg_floatwidth(  ROOT::Math::IMultiGenFunction & f1,
                                   ROOT::Math::IMultiGenFunction & f2) :
    fChi2_1(&f1), fChi2_2(&f2) {}

    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)
    double operator() (const double *par) const {
        double p1[13];
        for(int i = 0; i < 13; ++i) p1[i] = par[iparmassfit_poly3bkg_floatwidth[i]];

        double p2[16];
        for(int i = 0; i < 16; ++i) p2[i] = par[iparvarfit_poly3bkg_floatwidth[i]];

        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }

    const  ROOT::Math::IMultiGenFunction * fChi2_1;
    const  ROOT::Math::IMultiGenFunction * fChi2_2;
};

pair<double, double> fit_simultaneous_process(
      string dataset = "", 
      string input_d0_mass = "", string input_d0_mc = "", 
      string input_var = "", string var_name = "",
      string output = "",
      string type = "", int itrk=-1, int ipt=-1, 
      const vector<float>* ptbin=nullptr, const vector<unsigned int>* trkbin=nullptr
      //unordered_map<string, TGraphErrors*>* graph_map = nullptr
      )
{
   //gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   //gStyle->SetOptFit(1111);
   //
   double var = -1;
   double var_err = -1;
   
   double fit_range_low = 1.72;
   double fit_range_high = 2.0;
   double D0_mass = 1.8648;

   TString ntrk_str;
   if(trkbin->at(itrk+1) == numeric_limits<unsigned int>::max()) ntrk_str = Form("N_{trk}^{offline} #geq %u", trkbin->at(itrk));
   else ntrk_str = Form("%u #leq N_{trk}^{offline} < %u", trkbin->at(itrk), trkbin->at(itrk+1));

   TFile* fin_d0 = TFile::Open(input_d0_mass.c_str());
   if(isFailed(fin_d0, "d0 mass")) return pair<double, double>(var, var_err);
   
   TFile* fin_mc = TFile::Open(input_d0_mc.c_str());
   if(isFailed(fin_mc, "d0 mc mass")) return pair<double, double>(var, var_err);

   TFile* fin_var = TFile::Open(input_var.c_str());
   if(isFailed(fin_var, "var vs mass")) return pair<double, double>(var, var_err);

   TH1D* h_mc_match_signal = (TH1D*)fin_mc->Get("hMassPD0");
   TH1D* h_mc_match_all = (TH1D*)fin_mc->Get("hMassPD0_All");

   TH1D* h_data = (TH1D*)fin_d0->Get(Form("hMass _pt%d_trk%d%s", ipt, itrk, type.c_str()));

   TCanvas c("c", "", 800, 400);
   c.Divide(2, 1);
   c.cd(1)->SetTopMargin(0.06);
   c.cd(1)->SetLeftMargin(0.18);
   c.cd(1)->SetRightMargin(0.043);
   c.cd(1)->SetBottomMargin(0.145);
   c.cd(2)->SetTopMargin(0.06);
   c.cd(2)->SetLeftMargin(0.18);
   c.cd(2)->SetRightMargin(0.043);
   c.cd(2)->SetBottomMargin(0.145);

   TLatex tex;
   tex.SetNDC();
   tex.SetTextFont(42);
   tex.SetTextSize(0.045);
   tex.SetLineWidth(2);

   TLatex texCMS;
   texCMS.SetNDC();
   texCMS.SetTextFont(42);
   texCMS.SetTextSize(0.05);
   texCMS.SetTextAlign(12);

   c.cd(1);

   h_data->SetTitle("");
   h_data->SetMinimum(0);
   h_data->SetMarkerSize(0.8);
   h_data->SetMarkerStyle(20);
   h_data->SetLineWidth(1);
   h_data->SetOption("e");
   h_data->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
   h_data->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
   h_data->GetYaxis()->SetTitle("Entries / 5 MeV");
   h_data->GetXaxis()->CenterTitle();
   h_data->GetYaxis()->CenterTitle();
   h_data->GetXaxis()->SetTitleOffset(1.3);
   h_data->GetYaxis()->SetTitleOffset(2);
   h_data->GetXaxis()->SetLabelOffset(0.007);
   h_data->GetYaxis()->SetLabelOffset(0.007);
   h_data->GetXaxis()->SetTitleSize(0.045);
   h_data->GetYaxis()->SetTitleSize(0.045);
   h_data->GetXaxis()->SetTitleFont(42);
   h_data->GetYaxis()->SetTitleFont(42);
   h_data->GetXaxis()->SetLabelFont(42);
   h_data->GetYaxis()->SetLabelFont(42);
   h_data->GetXaxis()->SetLabelSize(0.04);
   h_data->GetYaxis()->SetLabelSize(0.04);

   h_data->GetXaxis()->SetNoExponent(true);
   ((TGaxis*)h_data->GetXaxis())->SetMaxDigits(7);

   h_data->SetMaximum(h_data->GetMaximum()*1.5);

   TF1* f = new TF1("fmass", "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x", fit_range_low, fit_range_high);
   f->SetLineColor(2);
   f->SetLineWidth(1);

   //first fit MC signal, swap and poly bkg set to 0

   f->SetParameter(0,100.);
   f->SetParameter(1,D0_mass);
   f->SetParameter(2,0.03);
   f->SetParameter(3,0.005);
   f->SetParameter(4,0.1);

   f->FixParameter(5,1);
   f->FixParameter(6,0); //always 0 in MC
   f->FixParameter(7,0.1); //does not really mater here as yield is fix to 0
   f->FixParameter(8,D0_mass); //does not really mater here as yield is fix to 0
   f->FixParameter(9,0);
   f->FixParameter(10,0);
   f->FixParameter(11,0);
   f->FixParameter(12,0);

   f->SetParLimits(2,0.01,0.1);
   f->SetParLimits(3,0.001,0.05);
   f->SetParLimits(4,0,1);
   f->SetParLimits(5,0,1);

   f->FixParameter(1,1.8648); //for first few attempt fix mean of gaussian to get reasonable estimation of other pars; later open it up
   h_mc_match_signal->Fit(f,"q","",fit_range_low,fit_range_high);
   h_mc_match_signal->Fit(f,"q","",fit_range_low,fit_range_high);
   f->ReleaseParameter(1); //now let gaussian mean float
   h_mc_match_signal->Fit(f,"L q","",fit_range_low,fit_range_high);
   h_mc_match_signal->Fit(f,"L q","",fit_range_low,fit_range_high);
   h_mc_match_signal->Fit(f,"L m","",fit_range_low,fit_range_high);

   //now fix signal double gaussian mean, sigma and gaus1,gaus2 yield ratio
   f->FixParameter(1,f->GetParameter(1));
   f->FixParameter(2,f->GetParameter(2));
   f->FixParameter(3,f->GetParameter(3));
   f->FixParameter(4,f->GetParameter(4));

   //now release swap bkg parameters to fit signal+swap MC
   f->ReleaseParameter(5);
   f->ReleaseParameter(7);
   f->ReleaseParameter(8);

   f->SetParameter(7,0.1);
   f->SetParameter(8,D0_mass);

   //fit signal+swap MC
   h_mc_match_all->Fit(f,"L q","",fit_range_low,fit_range_high);
   h_mc_match_all->Fit(f,"L q","",fit_range_low,fit_range_high);
   h_mc_match_all->Fit(f,"L q","",fit_range_low,fit_range_high);
   h_mc_match_all->Fit(f,"L q","",fit_range_low,fit_range_high);
   h_mc_match_all->Fit(f,"L m","",fit_range_low,fit_range_high);

   //now fix swap bkg parameters to fit data
   f->FixParameter(5,f->GetParameter(5));
   f->FixParameter(7,f->GetParameter(7));
   f->FixParameter(8,f->GetParameter(8));

   //now release poly bkg pars
   f->ReleaseParameter(9);
   f->ReleaseParameter(10);
   f->ReleaseParameter(11);
   f->ReleaseParameter(12);

   //now fit data
   h_data->Fit(f,"q","",fit_range_low,fit_range_high);
   h_data->Fit(f,"q","",fit_range_low,fit_range_high);
   f->ReleaseParameter(1); //allow data to have different mass peak mean than MC
   f->ReleaseParameter(6); //allow data to have different peak width than MC
   f->SetParameter(6,0);
   f->SetParLimits(6,-1.0,1.0);
   //f->FixParameter(5,1);
   h_data->Fit(f,"L q","",fit_range_low,fit_range_high);
   h_data->Fit(f,"L q","",fit_range_low,fit_range_high);
   h_data->Fit(f,"L E q","",fit_range_low,fit_range_high);
   auto result_mass = h_data->Fit(f,"L E m S","",fit_range_low,fit_range_high);

   //draw D0 signal separately
   TF1* f1 = new TF1("f_sig","[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
   f1->SetLineColor(kOrange-3);
   f1->SetLineWidth(1);
   f1->SetLineStyle(2);
   f1->SetFillColorAlpha(kOrange-3,0.3);
   f1->SetFillStyle(1001);
   f1->FixParameter(0,f->GetParameter(0));
   f1->FixParameter(1,f->GetParameter(1));
   f1->FixParameter(2,f->GetParameter(2));
   f1->FixParameter(3,f->GetParameter(3));
   f1->FixParameter(4,f->GetParameter(4));
   f1->FixParameter(5,f->GetParameter(5));
   f1->FixParameter(6,f->GetParameter(6));

   f1->Draw("LSAME");

   //draw swap bkg separately
   TF1* f2 = new TF1("f_swap",
         "[0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))+0*[1]*[2]*[3]*[4]", fit_range_low, fit_range_high);
   f2->SetLineColor(kGreen+4);
   f2->SetLineWidth(1);
   f2->SetLineStyle(1);
   f2->SetFillColorAlpha(kGreen+4,0.3);
   f2->SetFillStyle(1001);
   f2->FixParameter(0,f->GetParameter(0));
   f2->FixParameter(5,f->GetParameter(5));
   f2->FixParameter(6,f->GetParameter(6));
   f2->FixParameter(7,f->GetParameter(7));
   f2->FixParameter(8,f->GetParameter(8));

   f2->Draw("LSAME");

   TF1* f3 = new TF1("f_bkg","[9] + [10]*x + [11]*x*x + [12]*x*x*x+0*[0]*[1]*[2]*[3]*[4]*[5]*[6]*[7]*[8]", fit_range_low, fit_range_high);
   f3->SetLineColor(4);
   f3->SetLineWidth(1);
   f3->SetLineStyle(2);
   f3->FixParameter(9,f->GetParameter(9));
   f3->FixParameter(10,f->GetParameter(10));
   f3->FixParameter(11,f->GetParameter(11));
   f3->FixParameter(12,f->GetParameter(12));

   f3->Draw("LSAME");

   tex.DrawLatex(0.22,0.86,Form("%s", ntrk_str.Data()));
   tex.DrawLatex(0.22,0.80,Form("%.1f < p_{T} < %.1f GeV/c", ptbin->at(ipt), ptbin->at(ipt+1)));
   tex.DrawLatex(0.22,0.74,Form("|y| < 1.0"));
   tex.DrawLatex(0.22,0.68,Form("#chi_{2}/ndf = %.0f/%d",result_mass->Chi2(),result_mass->Ndf()));
   texCMS.DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
   texCMS.DrawLatex(0.62,0.97, "#scale[0.8]{pPb #sqrt{s_{NN}} = 8.16 TeV}");

   TLegend* leg = new TLegend(0.65,0.58,0.81,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.045);
   leg->SetTextFont(42);
   leg->SetFillStyle(0);
   leg->AddEntry(h_data,"data","p");
   leg->AddEntry(f,"Fit","L");
   leg->AddEntry(f1,"D^{0}+#bar{D^{#lower[0.2]{0}}} Signal","f");
   leg->AddEntry(f2,"K-#pi swap","f");
   leg->AddEntry(f3,"Combinatorial","l");
   leg->Draw("SAME");

   c.cd(2);

   TGraphErrors* var_data = (TGraphErrors*)fin_var->Get(Form("%s_pt%d_trk%d", var_name.c_str(), ipt, itrk));
   var_data->SetMarkerStyle(20);

   TF1* fmass_combinemassvarfit = new TF1("fmass_combinemassvarfit","[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x", fit_range_low, fit_range_high);
   TF1* fvar_combinemassvarfit = new TF1("fvar_combinemassvarfit", "( ( [0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) ) / ( [0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x ) ) * [13] + ( 1.0 - ( ( [0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) ) / ( [0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x ) ) ) * ( [14] + [15] * x)", fit_range_low, fit_range_high);

   fmass_combinemassvarfit->SetLineColor(2);
   fmass_combinemassvarfit->SetLineWidth(1);

   fvar_combinemassvarfit->SetLineColor(2);
   fvar_combinemassvarfit->SetLineWidth(1);

   ROOT::Math::WrappedMultiTF1 wfmass_combinemassvarfit(*fmass_combinemassvarfit,1);
   ROOT::Math::WrappedMultiTF1 wfvar_combinemassvarfit(*fvar_combinemassvarfit,1);

   ROOT::Fit::DataOptions opt;
   ROOT::Fit::DataRange range_massfit;

   range_massfit.SetRange(fit_range_low,fit_range_high);
   ROOT::Fit::BinData datamass(opt,range_massfit);
   ROOT::Fit::FillData(datamass, h_data);

   ROOT::Fit::DataRange range_varfit;
   range_varfit.SetRange(fit_range_low,fit_range_high);
   ROOT::Fit::BinData datavar(opt,range_varfit);
   ROOT::Fit::FillData(datavar, var_data);

   ROOT::Fit::Chi2Function chi2_B(datamass, wfmass_combinemassvarfit);
   ROOT::Fit::Chi2Function chi2_SB(datavar, wfvar_combinemassvarfit);

   GlobalChi2_poly3bkg_floatwidth globalChi2(chi2_B, chi2_SB);
   ROOT::Fit::Fitter fitter;

   const int Npar = 16;
   double par0[Npar];
   for( int ipar = 0; ipar < f->GetNpar(); ipar++ ) par0[ipar] = f->GetParameter(ipar);
   par0[13] = 0.01;
   par0[14] = 0.10;
   par0[15] = 0.05;

   fitter.Config().SetParamsSettings(Npar,par0);
   // fix parameter
   fitter.Config().ParSettings(2).Fix();
   fitter.Config().ParSettings(3).Fix();
   fitter.Config().ParSettings(4).Fix();
   fitter.Config().ParSettings(5).Fix();
   fitter.Config().ParSettings(7).Fix();
   fitter.Config().ParSettings(8).Fix();

   //fitter.Config().ParSettings(1).SetLimits(1.7, 2.0);

   fitter.Config().MinimizerOptions().SetPrintLevel(0);
   fitter.Config().SetMinimizer("Minuit2","Migrad");
   fitter.FitFCN(Npar,globalChi2,0,datamass.Size()+datavar.Size(),true);
   ROOT::Fit::FitResult result = fitter.Result();
   result.Print(std::cout);
   std::cout << "var, fit status: " << result.Status() << std::endl;
        
   fmass_combinemassvarfit->SetFitResult( result, iparmassfit_poly3bkg_floatwidth);
   fmass_combinemassvarfit->SetRange(range_massfit().first, range_massfit().second);
   fmass_combinemassvarfit->SetLineColor(kRed);
   h_data->GetListOfFunctions()->Add(fmass_combinemassvarfit);

   fvar_combinemassvarfit->SetFitResult( result, iparvarfit_poly3bkg_floatwidth);
   fvar_combinemassvarfit->SetRange(range_varfit().first, range_varfit().second);
   fvar_combinemassvarfit->SetLineColor(2);
   fvar_combinemassvarfit->SetLineStyle(2);
   var_data->GetListOfFunctions()->Add(fvar_combinemassvarfit);

   auto hist =  var_data->GetHistogram();
   auto y_elements = var_data->GetY();
   double max = *max_element(y_elements, y_elements+var_data->GetN());
   hist->SetLineWidth(0);
   hist->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
   hist->GetXaxis()->SetLimits(fit_range_low, fit_range_high);
   if(var_name == "Vn")hist->GetYaxis()->SetTitle("V_{2}");
   if(var_name == "jets")hist->GetYaxis()->SetTitle("Y_{jet}");
   hist->GetXaxis()->CenterTitle();
   hist->GetYaxis()->CenterTitle();
   hist->GetXaxis()->SetTitleOffset(1.3);
   hist->GetYaxis()->SetTitleOffset(2);
   hist->GetXaxis()->SetLabelOffset(0.007);
   hist->GetYaxis()->SetLabelOffset(0.007);
   hist->GetXaxis()->SetTitleSize(0.045);
   hist->GetYaxis()->SetTitleSize(0.045);
   hist->GetXaxis()->SetTitleFont(42);
   hist->GetYaxis()->SetTitleFont(42);
   hist->GetXaxis()->SetLabelFont(42);
   hist->GetYaxis()->SetLabelFont(42);
   hist->GetXaxis()->SetLabelSize(0.04);
   hist->GetYaxis()->SetLabelSize(0.04);
   hist->SetMinimum(0.00);
   hist->SetMaximum(2.0*max);
   var_data->SetTitle("");
   var_data->SetMarkerSize(0.8);
   var_data->SetLineWidth(1);
   hist->Draw();
   var_data->Draw("PE SAME");

   tex.DrawLatex(0.22,0.86,Form("%s", ntrk_str.Data()));
   tex.DrawLatex(0.22,0.80,Form("%.1f < p_{T} < %.1f GeV/c",  ptbin->at(ipt), ptbin->at(ipt+1)));
   tex.DrawLatex(0.22,0.74,Form("|y| < 1.0"));

   texCMS.DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
   texCMS.DrawLatex(0.62,0.97, "#scale[0.8]{pPb #sqrt{s_{NN}} = 8.16 TeV}");

   //tex.DrawLatex(0.22,0.68,Form("#chi_{2}/ndf = %.0f/%d",result_var->Chi2(), result_var->Ndf()));

   c.Print(Form("%s.png", output.c_str()));
   c.Print(Form("%s.pdf", output.c_str()));

   var = fvar_combinemassvarfit->GetParameter(13);
   var_err = fvar_combinemassvarfit->GetParError(13);

   delete leg;
   delete f;
   delete f1;
   delete f2;
   delete f3;

   return pair<double, double>(var, var_err);
}
#endif

void fit_simultaneous()
{
   vector<pair<double, double>> vec_Vn;
   vector<pair<double, double>> vec_jet;

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
   string cmdDirPlotsVn = 
   Form("if [ ! -d \"../plots/Vn\" ]; then\n"
         "    mkdir ../plots/Vn\n"
         "fi");
   string cmdDirPlotsJet = 
   Form("if [ ! -d \"../plots/Jet\" ]; then\n"
         "    mkdir ../plots/Jet\n"
         "fi");
   gSystem->Exec(cmdDirPlots.c_str());
   gSystem->Exec(cmdDirPlotsVn.c_str());
   gSystem->Exec(cmdDirPlotsJet.c_str());

   for(int iset=0; iset<1; iset++){
   //for(int iset=1; iset<2; iset++){
   //for(int iset=2; iset<3; iset++){

      cout << dataset[iset] << endl;

      auto trkbin = ana::get_Mult_Edges(dataset[iset]);
      int ntrk = trkbin.size() - 1;

      for(int itrk=0; itrk<ntrk; itrk++){
         for(int ipt=0; ipt<nPt; ipt++){
            TString output;
            if(trkbin.at(itrk+1) == numeric_limits<unsigned int>::max()) 
               output= Form("../plots/Vn/%s%u-inf_pT%.0f-%.0f%s", dataTrigger[iset].c_str(), 
                     trkbin[itrk], ptbin[ipt], ptbin[ipt+1], type[itype].c_str());
            else 
               output = Form("../plots/Vn/%s%u-%u_pT%.0f-%.0f%s", dataTrigger[iset].c_str(), 
                     trkbin[itrk], trkbin[itrk+1], ptbin[ipt], ptbin[ipt+1], type[itype].c_str());

            auto ret = fit_simultaneous_process(dataset[iset], Form("../data/corr2D_trg_pd0_%s.root", dataMult[iset].c_str()), 
                  Form("../MC/d0ana_hists_mass_pT%.1f-%.1f_y-1.0-1.0.root", ptbin[ipt], ptbin[ipt+1]),
                  Form("../data/%s_VarVsMass%s.root", dataMult[iset].c_str(),
                      type[itype].c_str()),
                  "Vn",
                  output.Data(), "", itrk, ipt, &ptbin, &trkbin
                  );
            vec_Vn.push_back(ret);
         }
      }

      for(int itrk=0; itrk<ntrk; itrk++){
         for(int ipt=0; ipt<nPt; ipt++){
            TString output;
            if(trkbin.at(itrk+1) == numeric_limits<unsigned int>::max()) 
               output= Form("../plots/Jet/%s%u-inf_pT%.0f-%.0f%s", dataTrigger[iset].c_str(), 
                     trkbin[itrk], ptbin[ipt], ptbin[ipt+1], type[itype].c_str());
            else 
               output = Form("../plots/Jet/%s%u-%u_pT%.0f-%.0f%s", dataTrigger[iset].c_str(), 
                     trkbin[itrk], trkbin[itrk+1], ptbin[ipt], ptbin[ipt+1], type[itype].c_str());

            auto ret = fit_simultaneous_process(dataset[iset], Form("../data/corr2D_trg_pd0_%s.root", dataMult[iset].c_str()), 
                  Form("../MC/d0ana_hists_mass_pT%.1f-%.1f_y-1.0-1.0.root", ptbin[ipt], ptbin[ipt+1]),
                  Form("../data/%s_VarVsMass%s.root", dataMult[iset].c_str(),
                      type[itype].c_str()),
                  "jets",
                  output.Data(), "", itrk, ipt, &ptbin, &trkbin
                  );
            vec_jet.push_back(ret);
         }
      }
      for(int itrk=0; itrk<ntrk; itrk++){
         for(int ipt=0; ipt<nPt; ipt++){
            printf("Vn = %.5f+/-%.5f\n", vec_Vn[itrk*nPt+ipt].first, vec_Vn[itrk*nPt+ipt].second);
            printf("jet = %.5f+/-%.5f\n", vec_jet[itrk*nPt+ipt].first, vec_jet[itrk*nPt+ipt].second);
         }
         printf("\n");
      }
   }
}

#include "../include/myAnaConsts.h"

int iparmassfit_poly3bkg_floatwidth[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
int iparvnfit_poly3bkg_floatwidth[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

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
        for(int i = 0; i < 16; ++i) p2[i] = par[iparvnfit_poly3bkg_floatwidth[i]];
        
        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }
    
    const  ROOT::Math::IMultiGenFunction * fChi2_1;
    const  ROOT::Math::IMultiGenFunction * fChi2_2;
};

void massfitVn_combine_pd0_ntrk()
{


   const char* input_mc= "../MC/d0ana_hists_mass_pT6.0-7.0_y-1.0-1.0.root";
   const char* input_data = "../data/corr2D_d0_d0ana_ntrk_pT6.0-7.0_y-1.0-1.0.root_v2.root";
   const char* output = "../data/v2vsNtrk_pd0_PAHM185-250_y1.0.root";
   const std::string dataset = "PAHM1-6";
   const float y=1.0;

   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L massfitVn_combine_pd0_ntrk_process.C");
   gInterpreter->ProcessLine(".L massfitVn_low_combine_pd0_ntrk_process.C");
   gInterpreter->ProcessLine(".L massfitJets_combine_pd0_ntrk_process.C");
   gInterpreter->ProcessLine(".L massfitJets_low_combine_pd0_ntrk_process.C");

   gInterpreter->ProcessLine(Form("massfitVn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT6.0-7.0_y-1.0-1.0.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2_sub_fixed.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2_sub.root", 
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_raw.root", 
            "../data/v2vsNtrk_pd0_PAHM185-250_pT6.0_7.0_y-1.0-1.0_new_raw.root", 
            "PAHM1-6", 
            1.0, 6.0, 7.0));

   gInterpreter->ProcessLine(Form("massfitVn_low_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT6.0-7.0_y-1.0-1.0.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2_sub_fixed.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2_sub.root", 
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_raw.root", 
            "../data/v2vsNtrk_pd0_PAHM185-250_pT6.0_7.0_y-1.0-1.0_new_raw.root", 
            "PAHM1-6", 
            1.0, 6.0, 7.0));

   gInterpreter->ProcessLine(Form("massfitJets_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT6.0-7.0_y-1.0-1.0.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2_sub_fixed.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2_sub.root", 
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_raw.root", 
            "../data/v2vsNtrk_pd0_PAHM185-250_pT6.0_7.0_y-1.0-1.0_new_raw.root", 
            "PAHM1-6", 
            1.0, 6.0, 7.0));
   gInterpreter->ProcessLine(Form("massfitJets_low_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT6.0-7.0_y-1.0-1.0.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2_sub_fixed.root", 
            //"../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_v2_sub.root", 
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk_new.root_raw.root", 
            "../data/v2vsNtrk_pd0_PAHM185-250_pT6.0_7.0_y-1.0-1.0_new_raw.root", 
            "PAHM1-6", 
            1.0, 6.0, 7.0));

   /*
   gInterpreter->ProcessLine(Form("massfitvn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT2.5-4.0_y-2.0-2.0.root", 
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root_v2.root", 
            "../data/v2vsNtrk_pd0_PAHM185-250_pT2.5_4.0_y-2.0-2.0.root", 
            "PAHM1-6", 
            2.0, 2.5, 4.0));

   gInterpreter->ProcessLine(Form("massfitvn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT2.5-4.0_y-2.0-2.0.root", 
            "../data/corr2D_trg_pd0_PAHM150-185_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root_v2.root", 
            "../data/v2vsNtrk_pd0_PAHM150-185_pT2.5_4.0_y-2.0-2.0.root", 
            "PAHM0", 
            2.0, 2.5, 4.0));

   gInterpreter->ProcessLine(Form("massfitvn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT2.5-4.0_y-2.0-2.0.root", 
            "../data/corr2D_trg_pd0_PAHM250-inf_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root_v2.root", 
            "../data/v2vsNtrk_pd0_PAHM250-inf_pT2.5_4.0_y-2.0-2.0.root", 
            "PAHM7", 
            2.0, 2.5, 4.0));

   gInterpreter->ProcessLine(Form("massfitvn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT2.5-4.0_y-2.0-2.0.root", 
            "../data/corr2D_trg_pd0_PAMB0-150_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root_v2.root", 
            "../data/v2vsNtrk_pd0_PAMB0-150_pT2.5_4.0_y-2.0-2.0.root", 
            "PAMB", 
            2.0, 2.5, 4.0));
            */
   /*
   gInterpreter->ProcessLine(Form("massfitvn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT3.0-5.0_y-2.0-2.0_pp.root",
            "../data/corr2D_trg_pd0_PPHM80-inf_d0ana_pT3.0-5.0_y-2.0-2.0_ntrk.root_v2.root",
            "../data/v2vsNtrk_pd0_PPHM80-inf_pT3.0_5.0_y-2.0-2.0.root", 
            "PPHM", 
            2.0, 3.0, 5.0));
            */
}

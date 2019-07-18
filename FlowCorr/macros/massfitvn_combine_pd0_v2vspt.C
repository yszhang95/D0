#include "../include/myAnaConsts.h"
void massfitvn_combine_pd0_v2vspt()
{
   const char* input_mc= "../MC/d0ana_hists_mass_pT6.0-7.0_y-1.0-1.0.root";
   const char* input_data = "../data/corr2D_d0_d0ana_ntrk_pT6.0-7.0_y-1.0-1.0.root_v2.root";
   const char* output = "../data/v2vsNtrk_pd0_PAHM185-250_y1.0.root";
   const std::string dataset = "PAHM1-6";
   const float y=1.0;

   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L massfitvn_combine_pd0_v2vspt_process.C");

   gInterpreter->ProcessLine(Form("massfitvn_combine_pd0_v2vspt_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT1.5-8.0_y-1.0-1.0_ptbin.root",
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT1.5-8.0_y-1.0-1.0_v2vspt.root_v2.root",
            "../data/v2vspt_pd0_PAHM185-250_pT1.5-8.0_y-1.0-1.0.root", 
            "PAHM1-6", 
            1.0, 1.5, 8.0));
}

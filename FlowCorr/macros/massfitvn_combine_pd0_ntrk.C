#include "../include/myAnaConsts.h"
void massfitvn_combine_pd0_ntrk()
{

   /*
   const char* input_mc= "../MC/d0ana_hists_mass_integrated.root";
   const char* input_data = "../data/corr2D_PAMB0-150_PD0_v2vsNtrk_2.0.root_v2.root";
   const char* output = "../data/v2vsNtrk_pd0_PAMB0-150_y2.0.root";
   const std::string dataset = "PAMB";
   const float y=2.0;
   */

   /*
   const char* input_mc= "../MC/d0ana_hists_mass_integrated.root";
   const char* input_data = "../data/corr2D_PAMB0-150_PD0_v2vsNtrk_2.0.root_v2_sub.root";
   const char* output = "../data/v2vsNtrk_pd0_PAMB0-150_v2_sub_y2.0.root";
   const std::string dataset = "PAMB";
   const float y=2.0;
   */

   /*
   const char* input_mc= "../MC/d0ana_hists_mass_integrated.root";
   const char* input_data = "consts.root";
   const char* output = "consts.root_jets.root";
   const std::string dataset = "PAMB";
   const float y=2.0;
   */

   const char* input_mc= "../MC/d0ana_hists_mass_pT6.0-7.0_y-1.0-1.0.root";
   const char* input_data = "../data/corr2D_d0_d0ana_ntrk_pT6.0-7.0_y-1.0-1.0.root_v2.root";
   const char* output = "../data/v2vsNtrk_pd0_PAHM185-250_y1.0.root";
   const std::string dataset = "PAHM1-6";
   const float y=1.0;

   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L massfitvn_combine_pd0_ntrk_process.C");

   /*
   gInterpreter->ProcessLine(Form("massfitvn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT6.0-7.0_y-1.0-1.0.root", 
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk.root_v2.root", 
            "../data/v2vsNtrk_pd0_PAHM185-250_pT6.0_7.0_y-1.0-1.0.root", 
            "PAHM1-6", 
            1.0, 6.0, 7.0));

   */
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
   gInterpreter->ProcessLine(Form("massfitvn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
            "../MC/d0ana_hists_mass_pT3.0-5.0_y-2.0-2.0_pp.root",
            "../data/corr2D_trg_pd0_PPHM80-inf_d0ana_pT3.0-5.0_y-2.0-2.0_ntrk.root_v2.root",
            "../data/v2vsNtrk_pd0_PPHM80-inf_pT3.0_5.0_y-2.0-2.0.root", 
            "PPHM", 
            2.0, 3.0, 5.0));
}

#include "../include/myAnaConsts.h"
 //void massfitvn_combine_pd0_ntrk_process(const char* input_mc = "../MC/d0ana_hists_mass_integrated.root",
              //const char* input_data = "../data/corr2D_PAMB0-150_PD0_v2vsNtrk_2.0.root_v2.root",
                     //const char* output = "../data/v2vsNtrk_pd0_PAMB0-150_y2.0.root",
                            //const std::string dataset = "PAMB",
                                   //const float y = 2.0
                                          //)

void massfitvn_combine_pd0_ntrk(const char* input_mc= "../MC/d0ana_hists_mass_integrated.root",
      const char* input_data = "../data/corr2D_PAMB0-150_PD0_v2vsNtrk_2.0.root_v2.root", 
      const char* output = "../data/v2vsNtrk_pd0_PAMB0-150_y2.0.root",
      const std::string dataset = "PAMB",
      const float y=2.0)
{
   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L massfitvn_combine_pd0_ntrk_process.C");
   gInterpreter->ProcessLine(Form("massfitvn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f)", input_mc, input_data, output, dataset.c_str(), y));
}
   

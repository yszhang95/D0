#include "../include/myAnaConsts.h"

void proj1D_v2vsNtrk_PD0(const char* input_d0= "../data/corr2D_PAMB0-150_PD0_v2vsNtrk_2.0.root",
      const char* input_ref = "../data/corr2D_PAMB0-150_REF_v2vsNtrk_2.0.root", 
const float y=2.0, const string dataset="PAMB")
{
   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L proj1D_v2vsNtrk_PD0_Process.C");
   //gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", input_d0, input_ref, dataset.c_str(), input_d0, 1.5, 8.0, -2., 2.));
   gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_d0_d0ana_ntrk_pT6.0-7.0_y-1.0-1.0.root", "../data/corr2D_PAHM185-250_REF_v2vsNtrk_2.0.root",
            "PAHM1-6", input_d0, 6.0, 7.0, -1., 1.));
   //gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", %f, \"%s\")", input_d0, input_ref, y, dataset.c_str()));
}
   

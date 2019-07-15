#include "../include/myAnaConsts.h"

void proj1D_v2vsNtrk_PD0(const char* input_d0= "../data/corr2D_PAMB0-150_PD0_v2vsNtrk_2.0.root",
      const char* input_ref = "../data/corr2D_PAMB0-150_REF_v2vsNtrk_2.0.root", 
const float y=2.0, const string dataset="PAMB")
{
   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L proj1D_v2vsNtrk_PD0_Process.C");
   gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", %f, \"%s\")", input_d0, input_ref, y, dataset.c_str()));
   //gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", %f, \"%s\")", input_d0, input_ref, y, dataset.c_str()));
}
   

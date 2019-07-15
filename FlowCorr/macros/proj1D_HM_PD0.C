#include "../include/myAnaConsts.h"

void proj1D_HM_PD0(const char* input_d0= "../data/corr2D_PAHM185-250_d0ana_HM_2.0.root",
      const char* input_ref = "../data/corr2D_PAHM185-250_REF_HM_2.0.root", 
const float y=2.0, const string dataset="PAHM1-6")
{
   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L proj1D_HM_PD0_Process.C");
   gInterpreter->ProcessLine(Form("proj1D_HM_PD0_Process(\"%s\", \"%s\", %f, \"%s\")", input_d0, input_ref, y, dataset.c_str()));
}
   
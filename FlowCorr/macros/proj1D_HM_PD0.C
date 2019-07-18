#include "../include/myAnaConsts.h"

void proj1D_HM_PD0(const char* input_d0= "",
      const char* input_ref = "",
      const string dataset="", const char* input_d0_low_mult="",
      const float pTMin=0., const float pTMax=0., 
      const float yMin =0., const float yMax =0.)
{
   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L proj1D_HM_PD0_Process.C");
   gInterpreter->ProcessLine(Form("proj1D_HM_PD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT1.5-8.0_y-1.0-1.0_v2vspt.root",
            "../data/corr2D_trg_ref_PAHM185-250_d0ana_v2vspt.root", "PAHM1-6", "", 
            1.5, 8.0, -1.0, 1.0));
}
   

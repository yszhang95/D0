#include "../include/myAnaConsts.h"

void proj1D_HM_NPD0(const char* input_d0= "",
      const char* input_ref = "",
      const string dataset="", const char* input_d0_low_mult="",
      const float pTMin=0., const float pTMax=0., 
      const float yMin =0., const float yMax =0.)
{
   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L proj1D_HM_NPD0_Process.C");

   /*
   gInterpreter->ProcessLine(Form("proj1D_HM_NPD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new.root",
            "../data/corr2D_trg_ref_PAHM185-250_d0ana_v2vspt.root", "PAHM1-6", "", 
            2.0, 8.0, -1.0, 1.0));
            */

   /*
   gInterpreter->ProcessLine(Form("proj1D_HM_NPD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-2.0-2.0.root",
            "../data/corr2D_trg_ref_PAHM185-250_d0ana_v2vspt.root", "PAHM1-6", "", 
            2.0, 8.0, -2.0, 2.0));
            */

   gInterpreter->ProcessLine(Form("proj1D_HM_NPD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_npd0_PAMB0-185_v2vspt_pT2.0-8.0_y-1.0-1.0.root",
            "../data/corr2D_trg_ref_PAMB0-185_d0ana_v2vspt.root", "PAMB", "", 
            2.0, 8.0, -1.0, 1.0));
}
   

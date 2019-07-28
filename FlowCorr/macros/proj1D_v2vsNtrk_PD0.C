#include "../include/myAnaConsts.h"

void proj1D_v2vsNtrk_PD0()
{
   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L proj1D_v2vsNtrk_PD0_Process.C");
   gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT4.0-6.0_y-1.0-1.0_ntrk_new.root", "../data/corr2D_trg_ref_PAHM185-250_d0ana_ntrk.root",
            "PAHM1-6", "../data/corr2D_trg_pd0_PAMB0-185_d0ana_pT4.0-6.0_y-1.0-1.0_ntrk_new.root",
            "../data/corr2D_trg_ref_PAMB0-185_d0ana_ntrk.root",
            4.0, 6.0, -1., 1.));
}

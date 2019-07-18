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
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT6.0-7.0_y-1.0-1.0_ntrk.root", "../data/corr2D_trg_ref_PAHM185-250_d0ana_ntrk.root",
            "PAHM1-6", "../data/corr2D_trg_pd0_PAMB0-150_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root", 6.0, 7.0, -1., 1.));
   //gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", %f, \"%s\")", input_d0, input_ref, y, dataset.c_str()));
   gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_pd0_PAHM185-250_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root", "../data/corr2D_trg_ref_PAHM185-250_d0ana_ntrk.root",
            "PAHM1-6", "../data/corr2D_trg_pd0_PAMB0-150_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root",
            2.5, 4.0, -2., 2.));
   gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_pd0_PAHM150-185_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root", "../data/corr2D_trg_ref_PAHM150-185_d0ana_ntrk.root",
            "PAHM1-6", "../data/corr2D_trg_pd0_PAMB0-150_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root",
            2.5, 4.0, -2., 2.));
   gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_pd0_PAHM250-inf_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root", "../data/corr2D_trg_ref_PAHM250-inf_d0ana_ntrk.root",
            "PAHM7", "../data/corr2D_trg_pd0_PAMB0-150_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root",
            2.5, 4.0, -2., 2.));
   gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_pd0_PAHM150-185_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root", "../data/corr2D_trg_ref_PAHM150-185_d0ana_ntrk.root",
            "PAHM0", "../data/corr2D_trg_pd0_PAMB0-150_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root",
            2.5, 4.0, -2., 2.));
   gInterpreter->ProcessLine(Form("proj1D_v2vsNtrk_PD0_Process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
            "../data/corr2D_trg_pd0_PAMB0-150_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root", "../data/corr2D_trg_ref_PAMB0-150_d0ana_ntrk.root",
            "PAMB", "../data/corr2D_trg_pd0_PAMB0-150_d0ana_pT2.5-4.0_y-2.0-2.0_ntrk.root",
            2.5, 4.0, -2., 2.));
}
   

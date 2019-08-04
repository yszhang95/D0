#include "../include/myAnaConsts.h"

void proj1D_v2vsNtrk_PD0_new()
{
   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L proj1D_v2vsNtrk_PD0_Process_new.C");

   const string dataset[] = {
      "PAMB", "PAHM1-6", "PAHM7"
   };
   const string mult[] = {
      "PAMB0-185", "PAHM185-250", "PAHM250-inf"
   };
   const string appendix[] = {
      "", "_loose", "_tight"
   };

   const int mode = 0;

   const float pTMin[] = {
      2., 4., 6.
   };
   const float pTMax[] = {
      4., 6., 8.
   };
   const float yMin = -1.;
   const float yMax = 1.;

   for(int iset=0; iset<1; iset++){
		for(int i=0; i<1; i++){
		   string input_d0 = string(
		         Form("../data/corr2D_trg_pd0_%s_d0ana_pT%.1f-%.1f_y%.1f-%.1f_ntrk_new.root%s", 
		            mult[iset].c_str(), pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
		      );
		   string input_ref = string(
		         Form("../data/corr2D_trg_ref_%s_d0ana_ntrk.root", mult[iset].c_str())
		      );
		   string input_d0_low = string(
		         Form("../data/corr2D_trg_pd0_%s_d0ana_pT%.1f-%.1f_y%.1f-%.1f_ntrk_new.root%s", 
		            mult[0].c_str(), pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
		      );
		   string input_ref_low = string(
		         Form("../data/corr2D_trg_ref_%s_d0ana_ntrk.root", mult[0].c_str())
		      );
		   gInterpreter->ProcessLine(
	            Form("proj1D_v2vsNtrk_PD0_Process_new(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %f)", 
		         input_d0.c_str(), input_ref.c_str(), dataset[iset].c_str(), input_d0_low.c_str(), input_ref_low.c_str(),
		         pTMin[i], pTMax[i], yMin, yMax)
	         );
		   }
   }
}

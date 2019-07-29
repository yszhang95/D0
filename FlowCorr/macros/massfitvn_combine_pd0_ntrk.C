#include "../include/myAnaConsts.h"

int iparmassfit_poly3bkg_floatwidth[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
int iparvnfit_poly3bkg_floatwidth[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

struct GlobalChi2_poly3bkg_floatwidth {
    GlobalChi2_poly3bkg_floatwidth(  ROOT::Math::IMultiGenFunction & f1,
                                   ROOT::Math::IMultiGenFunction & f2) :
    fChi2_1(&f1), fChi2_2(&f2) {}
    
    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)
    double operator() (const double *par) const {
        double p1[13];
        for(int i = 0; i < 13; ++i) p1[i] = par[iparmassfit_poly3bkg_floatwidth[i]];
        
        double p2[16];
        for(int i = 0; i < 16; ++i) p2[i] = par[iparvnfit_poly3bkg_floatwidth[i]];
        
        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }
    
    const  ROOT::Math::IMultiGenFunction * fChi2_1;
    const  ROOT::Math::IMultiGenFunction * fChi2_2;
};

void massfitvn_combine_pd0_ntrk()
{

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

   gInterpreter->ProcessLine(".include ../include");
   gInterpreter->ProcessLine(".L ../src/functions.cxx");
   gInterpreter->ProcessLine(".L massfitVn_combine_pd0_ntrk_process.C");
   gInterpreter->ProcessLine(".L massfitVn_low_combine_pd0_ntrk_process.C");
   gInterpreter->ProcessLine(".L massfitJets_combine_pd0_ntrk_process.C");
   gInterpreter->ProcessLine(".L massfitJets_combine_pd0_ntrk_process_test.C");
   gInterpreter->ProcessLine(".L massfitJets_low_combine_pd0_ntrk_process.C");
   gInterpreter->ProcessLine(".L massfitJets_low_combine_pd0_ntrk_process_test.C");
   gInterpreter->ProcessLine(".L fitNass.C");
   gInterpreter->ProcessLine(".L calPD0v2_sub.C");

   for(int iset=1; iset<2; iset++){
		for(int i=1; i<2; i++){
	      string input_mc = string(
	            Form("../MC/d0ana_hists_mass_pT%.1f-%.1f_y%.1f-%.1f.root%s", 
	               pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
	            );
		   string input_d0 = string(
		         Form("../data/corr2D_trg_pd0_%s_d0ana_pT%.1f-%.1f_y%.1f-%.1f_ntrk_new.root%s_raw.root", 
		            mult[iset].c_str(), pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
		         );
         string input_d0_low = string(
		         Form("../data/corr2D_trg_pd0_%s_d0ana_pT%.1f-%.1f_y%.1f-%.1f_ntrk_new.root%s_raw.root", 
		            mult[0].c_str(), pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
               );

         string output_Vn = string(
	            Form("../data/VnvsNtrk_pd0_%s_pT%.1f-%.1f_y%.1f-%.1f_new_raw.root%s", 
                  mult[iset].c_str(), pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
               );
         string output_Vn_low = string(
	            Form("../data/VnlowvsNtrk_pd0_%s_pT%.1f-%.1f_y%.1f-%.1f_new_raw.root%s", 
                  mult[0].c_str(), pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
               );
         string output_Jets = string(
	            Form("../data/JetsvsNtrk_pd0_%s_pT%.1f-%.1f_y%.1f-%.1f_new_raw.root%s", 
                  mult[iset].c_str(), pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
               );
         string output_Jets_low = string(
	            Form("../data/JetslowvsNtrk_pd0_%s_pT%.1f-%.1f_y%.1f-%.1f_new_raw.root%s", 
                  mult[0].c_str(), pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
               );
         string output_nass = string(
               Form("../data/nass_pd0_%s_pT%.1f-%.1f_y%.1f-%.1f_new_raw.root%s", 
                  mult[iset].c_str(), pTMin[i], pTMax[i], yMin, yMax, appendix[mode].c_str())
                  );

         /*
	      gInterpreter->ProcessLine(Form("massfitVn_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
               input_mc.c_str(), input_d0.c_str(), output_Vn.c_str(), dataset[iset].c_str(), yMax, pTMin[i], pTMax[i]));
	
	      gInterpreter->ProcessLine(Form("massfitVn_low_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
               input_mc.c_str(), input_d0.c_str(), output_Vn_low.c_str(), dataset[iset].c_str(), yMax, pTMin[i], pTMax[i]));
               */
	
	      gInterpreter->ProcessLine(Form("massfitJets_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
               input_mc.c_str(), input_d0.c_str(), output_Jets.c_str(), dataset[iset].c_str(), yMax, pTMin[i], pTMax[i]));
	      gInterpreter->ProcessLine(Form("massfitJets_combine_pd0_ntrk_process_test(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
               input_mc.c_str(), input_d0.c_str(), output_Jets.c_str(), dataset[iset].c_str(), yMax, pTMin[i], pTMax[i]));
	
         /*
	      gInterpreter->ProcessLine(Form("massfitJets_low_combine_pd0_ntrk_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
               input_mc.c_str(), input_d0.c_str(), output_Jets_low.c_str(), dataset[iset].c_str(), yMax, pTMin[i], pTMax[i]));
	      gInterpreter->ProcessLine(Form("massfitJets_low_combine_pd0_ntrk_process_test(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f)", 
               input_mc.c_str(), input_d0.c_str(), output_Jets_low.c_str(), dataset[iset].c_str(), yMax, pTMin[i], pTMax[i]));
               */
	
         /*
	      gInterpreter->ProcessLine(Form("fitNass(\"%s\", \"%s\", \"%s\", %f, %f, %f)", 
               input_d0.c_str(), output_nass.c_str(), dataset[iset].c_str(), yMax, pTMin[i], pTMax[i]));
	
	      gInterpreter->ProcessLine(Form("calPD0v2_sub(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\")",
               output_Vn.c_str(), output_Vn_low.c_str(), output_Jets.c_str(), output_Jets_low.c_str(), 
               output_nass.c_str(), input_d0.c_str(), dataset[iset].c_str()
	            ));
               */
	   }
   }
}

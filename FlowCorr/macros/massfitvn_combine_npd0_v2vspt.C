int iparmassfit_poly3bkg_floatwidth[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
int iparvnfit_poly3bkg_floatwidth[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

const int nDca= 2;
const float dcacut[2] = {0.012, 0.010};


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

void massfitvn_combine_npd0_v2vspt(
      )
{
      const char* input_mc = "";
      const char* input_data = "";
      const char* output = "";
      const string dataset = "";
      const float yMax = 0.0;
      const float pTMin = 0.0;
      const float pTMax = 0.0;
      const int   idca = 0;

      gInterpreter->ProcessLine(".include ../include");
      gInterpreter->ProcessLine(".L ../src/functions.cxx");
      gInterpreter->ProcessLine(".L massfitVn_combine_npd0_v2vspt_process_new.C");
      gInterpreter->ProcessLine(".L massfitVn_low_combine_npd0_v2vspt_process_new.C");
      gInterpreter->ProcessLine(".L massfitJets_combine_npd0_v2vspt_process.C");
      gInterpreter->ProcessLine(".L massfitJets_low_combine_npd0_v2vspt_process.C");
      gInterpreter->ProcessLine(".L fitNass_npd0.C");

      gInterpreter->ProcessLine(
            Form("massfitVn_combine_npd0_v2vspt_process_new(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin2_tight.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_tight_raw.root",
            "../data/V2vspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca0_new_binning_raw.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 0));


      gInterpreter->ProcessLine(
            Form("massfitVn_combine_npd0_v2vspt_process_new(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin2_tight.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_tight_raw.root",
            "../data/V2vspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca1_new_binning_raw.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 1));

      gInterpreter->ProcessLine(
            Form("massfitVn_low_combine_npd0_v2vspt_process_new(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin2_tight.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_tight_raw.root",
            "../data/V2_lowvspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca0_new_binning_raw.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 0));

      gInterpreter->ProcessLine(
            Form("massfitVn_low_combine_npd0_v2vspt_process_new(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin2_tight.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_tight_raw.root",
            "../data/V2_lowvspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca1_new_binning_raw.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 1));


      gInterpreter->ProcessLine(
            Form("massfitJets_combine_npd0_v2vspt_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin2_tight.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_tight_raw.root",
            "../data/Jetsvspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca0_new_binning_raw.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 0));

      gInterpreter->ProcessLine(
            Form("massfitJets_combine_npd0_v2vspt_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin2_tight.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_tight_raw.root",
            "../data/Jetsvspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca1_new_binning_raw.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 1));

      gInterpreter->ProcessLine(
            Form("massfitJets_low_combine_npd0_v2vspt_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin2_tight.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_tight_raw.root",
            "../data/Jets_lowvspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca0_new_binning_raw.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 0));
      gInterpreter->ProcessLine(
            Form("massfitJets_low_combine_npd0_v2vspt_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin2_tight.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_tight_raw.root",
            "../data/Jets_lowvspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca1_new_binning_raw.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 1));

      gInterpreter->ProcessLine(
            Form("fitNass_npd0(\"%s\", \"%s\")", "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new_binning.root_tight_raw.root", "../data/nass_npd0_pT2.0-8.0_y-1.0-1.0_new_binning.root")
            );

}

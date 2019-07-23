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
      gInterpreter->ProcessLine(".L massfitvn_combine_npd0_v2vspt_process.C");

      /*
      gInterpreter->ProcessLine(
            Form("massfitvn_combine_npd0_v2vspt_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new.root_v2.root",
            "../data/v2vspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca0_new.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 0));

      gInterpreter->ProcessLine(
            Form("massfitvn_combine_npd0_v2vspt_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin.root",
            "../data/corr2D_trg_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_new.root_v2.root",
            "../data/v2vspt_npd0_PAHM185-250_pT2.0-8.0_y-1.0-1.0_dca1_new.root",
            "PAHM1-6",
            1.0, 2.0, 8.0, 1));
            */

      gInterpreter->ProcessLine(
            Form("massfitvn_combine_npd0_v2vspt_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin.root",
            "../data/corr2D_trg_npd0_PAMB0-185_v2vspt_pT2.0-8.0_y-1.0-1.0.root_v2.root",
            "../data/v2vspt_npd0_PAMB0-185_pT2.0-8.0_y-1.0-1.0_dca0.root",
            "PAMB",
            1.0, 2.0, 8.0, 0));

      gInterpreter->ProcessLine(
            Form("massfitvn_combine_npd0_v2vspt_process(\"%s\", \"%s\", \"%s\", \"%s\", %f, %f, %f, %d)",
            "../MC/npd0ana1_hists_mass_pT2.0-8.0_y-1.0-1.0_ptbin.root",
            "../data/corr2D_trg_npd0_PAMB0-185_v2vspt_pT2.0-8.0_y-1.0-1.0.root_v2.root",
            "../data/v2vspt_npd0_PAMB0-185_pT2.0-8.0_y-1.0-1.0_dca1.root",
            "PAMB",
            1.0, 2.0, 8.0, 1));
}

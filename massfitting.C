//  this is a macro for fitting mass of d0 produced at cms                    
//  Here we have a short instruction for the macro                            
//  There are three parameters:                                               
//    The 1st paramter is the histogram represneting the data               
//    The 2nd paramter is the histogram represneting the signal of MC samples 
//    The 3rd paramter is the histogram represneting 
//        the signal and swap of MC samples  
//    The 4th paramter is the name of the fitting function            

#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TMath.h"

TF1 massfitting
(TH1* h_data, TH1* h_mc_match_signal, TH1* h_mc_match_all, const char* fname)
{
//   The full fitting function is constructed as follow
//   [0] is signal + swap yield;
//   [1] is common mean of double gaussian;
//   [2] is signal gaussian 1 sigma;
//   [3] is signal gaussian 2 sigma;
//   [4] is fractional signal gaussian 1 yield; 
//   1-[4] is fractional signal gaussian 2 yield;
//   [5] is fractional double gaussian signal yield, 
//   1-[5] is fractional swap yield;
//   [6] is a factor to let width of the gaussians to vary in data;
//   [7] is swap gaussian sigma;
//   [8] is swap gaussian mean;
//   [9-12] is 3rd order poly parameters
    const std::string double_gaussian = 
        "[4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))"
        "/(sqrt(2*3.14159)*[2]*(1.0 +[6]))"
        "+ (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))"
        "/(sqrt(2*3.14159)*[3]*(1.0 +[6]))";
    const std::string swap_gaussian =
        "TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))";
    const std::string bkg3rdpoly = "[9] + [10]*x + [11]*x*x + [12]*x*x*x";
    const std::string massfunc = "[0]*([5]*(" + double_gaussian +")"
        " + (1-[5])*(" + swap_gaussian + ")) + " + bkg3rdpoly;
    double fit_range_low = 1.7;
    double fit_range_high = 2.0;
    double D0_mass = 1.8648;
    TF1 f(fname, massfunc.c_str(), fit_range_low, fit_range_high);

    //first fit MC signal, swap and poly bkg set to 0
        
    f.SetParameter(0,100.);
    f.SetParameter(1,D0_mass);
    f.SetParameter(2,0.03);
    f.SetParameter(3,0.005);
    f.SetParameter(4,0.1);
    
    f.FixParameter(5,1);
    f.FixParameter(6,0); //always 0 in MC
    f.FixParameter(7,0.1); //does not really mater here as yield is fix to 0
    f.FixParameter(8,D0_mass); //does not really mater here 
                                // as yield is fix to 0
    f.FixParameter(9,0);
    f.FixParameter(10,0);
    f.FixParameter(11,0);
    f.FixParameter(12,0);
    
    f.SetParLimits(2,0.01,0.1);
    f.SetParLimits(3,0.001,0.05);
    f.SetParLimits(4,0,1);
    f.SetParLimits(5,0,1);
    
    f.FixParameter(1,1.8648); //for first few attempt fix mean of gaussian 
                               //to get reasonable estimation of other pars;
                               // later open it up
    h_mc_match_signal->Fit(&f,"q","",fit_range_low,fit_range_high);
    h_mc_match_signal->Fit(&f,"q","",fit_range_low,fit_range_high);
    f.ReleaseParameter(1); //now let gaussian mean float
    h_mc_match_signal->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_signal->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_signal->Fit(&f,"L m","",fit_range_low,fit_range_high);

    //now fix signal double gaussian mean, sigma and gaus1,gaus2 yield ratio
    f.FixParameter(1,f.GetParameter(1));
    f.FixParameter(2,f.GetParameter(2));
    f.FixParameter(3,f.GetParameter(3));
    f.FixParameter(4,f.GetParameter(4));
    
    //now release swap bkg parameters to fit signal+swap MC
    f.ReleaseParameter(5);
    f.ReleaseParameter(7);
    f.ReleaseParameter(8);
    
    f.SetParameter(7,0.1);
    f.SetParameter(8,D0_mass);
    
    //fit signal+swap MC
    h_mc_match_all->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_all->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_all->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_all->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_all->Fit(&f,"L m","",fit_range_low,fit_range_high);


    //now fix swap bkg parameters to fit data
    f.FixParameter(5,f.GetParameter(5));
    f.FixParameter(7,f.GetParameter(7));
    f.FixParameter(8,f.GetParameter(8));

    //now release poly bkg pars
    f.ReleaseParameter(9);
    f.ReleaseParameter(10);
    f.ReleaseParameter(11);
    f.ReleaseParameter(12);
    
    //now fit data
    f.SetParLimits(0, 0, 1e8);
    f.SetParLimits(1, 1.86, 1.87);

    //h_data->Fit(&f,"q","",fit_range_low,fit_range_high);
    //h_data->Fit(&f,"q","",fit_range_low,fit_range_high);
    h_data->Fit(&f,"q","",fit_range_low,fit_range_high);
    h_data->Fit(&f,"q","",fit_range_low,fit_range_high);

    f.ReleaseParameter(1); //allow data to have different 
                            //  mass peak mean than MC
    f.ReleaseParameter(6); //allow data to have different peak width than MC
    f.SetParameter(6,0);
    //f.SetParLimits(6,-1,1);
    f.SetParLimits(6,-0.5,0.5);
    //f.FixParameter(5,1);
    h_data->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_data->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_data->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_data->Fit(&f,"L m","",fit_range_low,fit_range_high);

    return f;
}

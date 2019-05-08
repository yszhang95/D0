#ifndef _MYANACONSTS_h
#define _MYANACONSTS_h
#include <string>
namespace ana{
   const int nuofY = 1;
   const float ybin[nuofY+1] = {-0.8, 0.8};
   const std::string whichtree[] = {"d0ana", "npd0ana", "npd0ana1", "npd0ana2"};

   const double mvaMin = 0.40;
   const double mvaMax = 0.80;
   const double mvaStep = 0.02;

   const int nuofpt = 8;
   const double ptbin[nuofpt+1]={1.5, 2.4, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0}; // GeV

   const int nuofDca = 12;
   const double dcaBin[nuofDca+1] = {0, 0.002, 0.004, 0.006, 0.008, 0.012, 0.016, 0.022, 0.028, 0.036, 0.044, 0.06, 0.08}; // cm
   const int nDca = 80;
   const double dcaMin = 0;
   const double dcaMax = 0.08;

   const double fit_range_low = 1.7;
   const double fit_range_high = 2.0;
   const double D0_mass = 1.8648;
}


int whichPt(const double pT){
   for(int ipt=0; ipt<ana::nuofpt; ipt++){
      if(ana::ptbin[ipt]<pT && pT<ana::ptbin[ipt+1]) return ipt;
   }
   return -1;
} 

int whichY(const double y){
   for(int iy=0; iy<ana::nuofY; iy++){
      if(ana::ybin[iy]<y && y>ana::ybin[iy]) return iy;
   }
   return -1;
} 
#endif

#include "myAnaConsts.h"

using namespace std;
string ana::findDCA(const double& DCA, const bool& isPromptD0)
{
   if(isPromptD0) {
      return "smallDCA";
   }else{
      bool isLarge = DCA > ana::dcaSep;
      if(isLarge) return "largeDCA";
      if(!isLarge) return "smallDCA";
   }
   return "";
}

int ana::findMassBin(const double& mass)
{
   for(int imass=0; imass<ana::nMass; imass++){
      if(ana::massbin[imass]<=mass && mass<ana::massbin[imass+1]) return imass;
   }
   return -1;
}

int ana::findPtBin(const double& pT)
{
   for(int ipt=0; ipt<ana::nPt; ipt++){
      if(ana::ptbin[ipt]<=pT && pT<ana::ptbin[ipt+1]) return ipt;
   }
   return -1;
}

int ana::findYBin(const double& y)
{
   for(int iy=0; iy<ana::nY; iy++){
      if(ana::ybin[iy]<=fabs(y) && fabs(y)<ybin[iy+1]) return iy;
   }
   return -1;
}

int ana::findZVtxBin(const float& z_vtx)
{
   for(unsigned int iz=0; iz<ana::nZ_Vtx_; iz++){
      if(ana::z_vtx_bin_[iz] <= z_vtx && 
            z_vtx < ana::z_vtx_bin_[iz+1])
         return iz;
   }
   return -1;
}

bool ana::pass_PD0_MVA(const float& pT, const float& mva){
   if(pT<4.) return mva > 0.58;
   else if(pT < 6. && pT >= 4.) return mva > 0.46;
   else if(pT < 8. && pT >= 6.) return mva > 0.28;
   else if(pT < 20. & pT >= 8.) return mva > 0.1;
   return false;
}

bool pass_NPD0_MVA(const float& pT, const float& mva){
   if(pT<4.) return mva > 0.52;
   else if(pT < 6. && pT >= 4.) return mva > 0.44;
   else if(pT < 8. && pT >= 6.) return mva > 0.32;
   else if(pT < 20. & pT >= 8.) return mva > 0.2;
   return false;
}

unsigned int findNtrkBin(const unsigned int& n_trk_bin)
{
}


#include "myAnaConsts.h"

#include <algorithm>

using namespace std;

int ana::findMassBin(const double& mass)
{
   for(int imass=0; imass<ana::nMass; imass++){
      if(ana::massbin[imass]<=mass && mass<ana::massbin[imass+1]) return imass;
   }
   return -1;
}

int ana::findPtBin(const double& pT, const vector<double>& ptbin)
{
   const unsigned int nPt = ptbin.size() - 1;
   for(unsigned int ipt=0; ipt<nPt; ipt++){
      if(ptbin[ipt]<=pT && pT<ptbin[ipt+1]) return ipt;
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

bool ana::pass_pPb2016_8TeV_PD0_MVA(const float& pT, const float& mva)
{
   if(pT<4.) return mva > 0.58;
   else if(pT < 6. && pT >= 4.) return mva > 0.46;
   else if(pT < 8. && pT >= 6.) return mva > 0.28;
   else if(pT < 20. && pT >= 8.) return mva > 0.1;
   return false;
}

bool ana::pass_pPb2016_8TeV_NPD0_MVA(const float& pT, const float& mva)
{
   if(pT<4.) return mva > 0.52;
   else if(pT < 6. && pT >= 4.) return mva > 0.44;
   else if(pT < 8. && pT >= 6.) return mva > 0.32;
   else if(pT < 20. && pT >= 8.) return mva > 0.2;
   return false;
}

bool ana::pass_pp2018_13TeV_PD0_MVA(const float& mva)
{
   return mva > 0.52;
}

int ana::Get_Trigger(const std::string& dataset)
{
   auto it = ana::dataset_trigger.find(dataset);
   if(it != ana::dataset_trigger.end()) return it->second;
   cerr << "type in the correct dataset!" << endl;
   return -1;
}

int ana::Get_N_nTrkBin(const std::string& dataset)
{
   auto it = ana::dataset_N_nTrkBin.find(dataset);
   if(it != ana::dataset_N_nTrkBin.end()) return it->second;
   cerr << "type in the correct dataset!" << endl;
   return -1;
}

int ana::findNtrkBin(const double& nTrkOffline, const int& trigger)
{
   switch(trigger) {
     case 0 :
        if(nTrkOffline<35 && nTrkOffline>=0) return 0;
        else if(nTrkOffline<90 && nTrkOffline>=35) return 1;
        else if(nTrkOffline<150 && nTrkOffline>=90) return 2;
        break;
     case 1 :
        if(nTrkOffline<185 && nTrkOffline>=150) return 0;
        break;
     case 2 :
        if(nTrkOffline<250 && nTrkOffline>=185) return 0;
        break;
     case 3 :
        if(nTrkOffline>=250) return 0;
        break;
     case 4:
        if(nTrkOffline<20) return 0;
        else if(nTrkOffline<40 && nTrkOffline>=20) return 1;
        else if(nTrkOffline<80 && nTrkOffline>=40) return 2;
        break;
     case 5 :
        if(nTrkOffline<100 && nTrkOffline>=80) return 0;
        break;
     case 6 :
        if(nTrkOffline>100) return 0;
         break;
   }
   return -1;
}

bool ana::isHM_PD0_DataSet(const string& dataset_name)
{
   auto it = find(ana::dataset_HM_PD0.begin(), ana::dataset_HM_PD0.end(), dataset_name);
   if(it != ana::dataset_HM_PD0.end()) return true;
   return false;
}

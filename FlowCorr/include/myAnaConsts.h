#ifndef __MYANACONST_h__
#define __MYANACONST_h__

#include <iostream>
#include <string>
#include <map>
#include <array>
#include <limits>

#include "TMath.h"

namespace ana{

   const double PI = TMath::Pi();

   const int nEtaBin = 33;
   const double etaBegin = -4.95;
   const double etaEnd = 4.95;

   const int nPhiBin = 32;
   const double phiBegin = -(0.5-1.0/32)*PI;
   const double phiEnd = (1.5-1.0/32)*PI;
   const double phiBegin_old = -0.5*PI;
   const double phiEnd_old = 1.5*PI;

   const double etaMin_ass_ = -2.4;
   const double etaMax_ass_ = 2.4;
   const double ptMin_ass_ = 0.3;
   const double ptMax_ass_ = 3.0;

   const int    multMin_PA_ = 185;
   const int    multMax_PA_ = 250;

   const int    multMin_low_PA_ = 0;
   const int    multMax_low_PA_ = 35;

   const int    multMin_PP_ = 100;
   const int    multMax_PP_ = 250;

   const int    multMin_low_PP_ = 0;
   const int    multMax_low_PP_ = 20;

   const bool   rejectDaughter_ = true;

   const double d0_dau_abs_eta_max_ = 2.4;
   const int    d0_dau_nhit_min_ = 11;
   const double d0_dau_pterr_max_ = 0.1;
   const double d0_dau_pt_min_ = 0.7;

   const double d0_pt_min_ = 1.5;
   const double d0_pt_max_ = 8.0;

   const double d0_y_min_  = -2;
   const double d0_y_max_  = 2;

   const int nMass = 14;
   const double massbin[nMass+1] = {1.72, 1.75, 1.78, 1.8, 1.82, 1.84, 1.85, 1.86, 1.865, 
      1.87, 1.88, 1.9, 1.92, 1.96, 2.0};
   const double massbin_HM[nMass+1] = {1.70, 1.74, 1.78, 1.8, 1.82, 1.84, 1.85, 1.86, 1.865, 
      1.87, 1.88, 1.9, 1.92, 1.96, 2.0};
   
   //const int nPt_NPD0_pPb = 6;
   //const double ptbin_NPD0_pPb[nPt_NPD0_pPb+1] = {2., 3., 4., 5., 6., 7., 8.};
   const int nPt_NPD0_pPb = 2;
   const double ptbin_NPD0_pPb[nPt_NPD0_pPb+1] = {2., 5., 8.};
   const int nPt_PD0_pPb = 8;
   const double ptbin_PD0_pPb[nPt_PD0_pPb+1] = {1.5, 2.4, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0};
   const int nPt_PD0_pp = 3;
   const double ptbin_PD0_pp[nPt_PD0_pp+1] = {2., 4., 6., 8.};
   
   const int nY = 1;
   const double ybin[nY+1] = {0., 2.0};

   const double dcaSep = 0.016;

   const unsigned int nMixedEvts = 10;

   const unsigned int nZ_Vtx_ = 15;
   const double z_vtx_bin_[nZ_Vtx_+1] = {
      -15, -13, -11, -9, -7, -5, -3, -1,
      1, 3, 5, 7, 9, 11, 13, 15
   };

   const std::string treeName[] = {
      "d0ana", "npd0ana", "npd0ana1"
   };

   const std::map<std::string, int> dataset_trigger = {
      {"PAMB", 0},
      {"PAHM1-6", 1},
      {"PAHM7", 2},
      {"PPMB", 3},
      {"PPHM", 4}, //80-100, >100
   };

   const std::array<std::string, 2> dataset_HM_PD0 = {"PAHM1-6", "PPHM"};
   const std::array<std::string, 2> dataset_Low_Mult_PD0 = {"PAMB", "PPMB"};
   const std::array<std::string, 1> dataset_HM_NPD0 = {"PAHM1-6"};
   const std::array<std::string, 1> dataset_Low_Mult_NPD0 = {"PAMB"};

   const std::map<std::string, int> dataset_N_nTrkBin = {
      {"PAMB", 4},
      {"PAHM1-6", 1},
      {"PAHM7", 1},
      {"PPMB", 3},
      {"PPHM", 2},
   };

   const std::map<std::string, int> dataset_PA_N_nTrkBin = {
      {"PAMB", 4},
      {"PAHM1-6", 1},
      {"PAHM7", 1},
   };

   const std::map<std::string, int> dataset_PP_N_nTrkBin = {
      {"PPMB", 3},
      {"PPHM", 2},
   };

   const unsigned int PA_Mult_Edges[] = {
      0, 35, 90, 120, 185, 250, std::numeric_limits<unsigned int>::max()
   };
   const unsigned int PP_Mult_Edges[] = {
      0, 20, 40, 80, 100, std::numeric_limits<unsigned int>::max()
   };

   const unsigned int PA_Mult_N = 3;
   const unsigned int PP_Mult_N = 2;
   const std::array<std::string, PA_Mult_N>PA_Mult_Order = {
      "PAMB", "PAHM1-6", "PAHM7"
   };
   const std::array<std::string, PP_Mult_N>PP_Mult_Order = {
      "PPMB", "PPHM" 
   };

   int findMassBin(const double&);
   int findMassBin_HM(const double&);
   int findPtBin(const double&, const std::vector<double>&);
   int findZVtxBin(const float&);

   bool pass_pPb2016_8TeV_PD0_MVA(const float&, const float&);
   bool pass_pPb2016_8TeV_PD0_tightMVA(const float&, const float&);
   bool pass_pPb2016_8TeV_PD0_looseMVA(const float&, const float&);

   bool pass_pPb2016_8TeV_NPD0_MVA(const float&, const float&);
   bool pass_pp2018_13TeV_PD0_MVA(const float&);

   int Get_Trigger(const std::string&);
   int Get_N_nTrkBin(const std::string&);
   int findNtrkBin(const double&, const int&);

   bool isHM_PD0_DataSet(const std::string&);
   bool isHM_NPD0_DataSet(const std::string&);

   bool isLow_Mult_PD0_DataSet(const std::string&);
   bool isLow_Mult_NPD0_DataSet(const std::string&);

   std::vector<unsigned int> get_Mult_Edges(const std::string&);
};

#endif

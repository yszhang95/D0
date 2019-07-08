#!/bin/bash

echo "setup cmssw"
cd /afs/cern.ch/user/y/yousen/public/pPb2016
source setup_CMSSW_8_0_31
#cd /home/yz144/pPb2016
#source setupfile8xy
cd CMSSW_8_0_31/src
eval `scramv1 runtime -sh`
cd D0/FlowCorr/batch
echo PWD: $PWD

./corr2D_trg_ref list/data_ref.list.$1 $2

#!/bin/bash

echo "setup cmssw"
cd /afs/cern.ch/user/y/yousen/public/pPb2016
source setup_CMSSW_8_0_31
cd CMSSW_8_0_31/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/user/y/yousen/work/pPb2016/D0/FlowCorr/batch
echo PWD: $PWD

../bin/corr2D_trg_ref_mult $1 $2 $3

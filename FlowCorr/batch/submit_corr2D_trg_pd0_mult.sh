#!/bin/bash

echo "setup cmssw"
cd /afs/cern.ch/user/y/yousen/public/pPb2016
source setup_CMSSW_8_0_31
cd CMSSW_8_0_31/src
eval `scramv1 runtime -sh`

# cd to work dir
cd /afs/cern.ch/user/y/yousen/work/pPb2016/D0/FlowCorr/batch
echo PWD: $PWD

#OUTDIR=/afs/cern.ch/user/y/yousen/work/pPb2016/HM185-250-PD0-v2vspt
#../bin/corr2D_trg_pd0 ../test/data_ref_HM.small.list PAHM1-6 ../eff/fEff.root $OUTDIR
# run
../bin/corr2D_trg_pd0_mult $1 $2 $3 $4 $5 $6 $7 $8

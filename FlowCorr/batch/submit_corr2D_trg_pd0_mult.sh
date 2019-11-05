#!/bin/bash

ROOTPATH=$9

echo "setup cmssw"
cd $ROOTPATH
source setup
export $SCRAM_ARCH

# cd to work dir
cd ${ROOTPATH}/D0/FlowCorr/batch
root -b -q
echo PWD: $PWD

#OUTDIR=/afs/cern.ch/user/y/yousen/work/pPb2016/HM185-250-PD0-v2vspt
#../bin/corr2D_trg_pd0 ../test/data_ref_HM.small.list PAHM1-6 ../eff/fEff.root $OUTDIR
# run
../bin/corr2D_trg_pd0_mult $1 $2 $3 $4 $5 $6 $7 $8 d0ana_newreduced

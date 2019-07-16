# Compiling
Makefile help one compile the source codes.

Before compiling, make sure dir, include, src, bin, lib, exist
# Executables 
Exes are in the dir bin, 
1. `corr2D_trg_pd0_mult` give the 2D correlation functions of prompt d0,
input parameters are the list of the ROOT files, dataset, effciency file, output dir, pTMin, pTMax, yMin, yMax

2. `corr2D_trg_ref_mult` give the 2D correlation functions  of ref particles, 
input parameters are the list of the ROOT files, dataset, output dir

3. `corr2D_trg_pd0` give the 2D correlation functions  of prompt-d0s, 
input parameters are the list of the ROOT files, dataset, effciency file, output dir, pTMin, pTMax, yMin, yMax

4. `corr2D_trg_npd0` give the 2D correlation functions  of non-prompt-d0s, 
input parameters are the list of the ROOT files, dataset, effciency file, output dir, pTMin, pTMax, yMin, yMax, csv file for dca cut

dcafile follows the format
```
pt0:pt1:pt2:
dcaCut0_0, dcaCut1_0, dcaCut2_0
dcaCut0_1, dcaCut1_1, dcaCut2_1
```
dcaCutX_X split the samples into different bin, (dcaCut0_0, dcaCut0_1)...(dcaCut0_X, inf), pt%d means different ptbin, defined in myAnaConsts.h

5. `corr2D_trg_ref` give the 2D correlation functions  of ref particles, 
input parameters are the list of the ROOT files, dataset, output dir, d0 tree index

datasets are listed, comments follows `//`:

```
PAMB    // 0-35, 35-90, 90-150
PAHM0   // 150-185
PAHM1-6 // 185-250
PAHM7   // 250-inf
PPMB    // MB 0-20, 20-40, 40-80
PPHM    // high mult 80-100 , 100 - inf
```

# Batch jobs
## out-dated

dir batch has the scrips for batch jobs, make sure `submit_corr2D_trg_d0/ref.sh` and `submit_jobs_trg_d0/ref.py`
when one want to submit jobs, edit `*.py`, and then type `python submit_XXXXX.py`, create a dir storing the list of files, another for output,
dir out, log, err, documenting the output from condor

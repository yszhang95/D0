# Compiling
Makefile help one compile the source codes.

Before compiling, make sure dir, include, src, bin, lib, exist
# Executables 
Exes are in the dir bin, `corr2D_trg_pd0_mult` give the 2D correlation functions seperately,
input parameter is the list of the ROOT files, dataset, effciency file, output dir
# Batch jobs
## out-dated

dir batch has the scrips for batch jobs, make sure `submit_corr2D_trg_d0/ref.sh` and `submit_jobs_trg_d0/ref.py`
when one want to submit jobs, edit `*.py`, and then type `python submit_XXXXX.py`, create a dir storing the list of files, another for output,
dir out, log, err, documenting the output from condor

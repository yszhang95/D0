import subprocess

tree = ["PromptD0", "NonPromptD0", "NonPromptD0"]
treeNumber = 0

#dataset_name = 'PAHM1-6'
#dataset_name = 'PAHM0'
#dataset_name = 'PAHM7'
dataset_name = 'PAMB'

dataset = {
      'PAMB'    : 'PAMB0-150.list',
      'PAHM0'   : 'PAHM150-185.list',
      'PAHM1-6' : 'PAHM185-250.list',
      'PAHM7'   : 'PAHM250-inf.list',
      'PPMB'    : 'PPMB0-80.list',
      'PPHM_1'  : 'PPHM80-100.list',
      'PPHM_2'  : 'PPHM100-inf.list'
      }
storage = {
      'PAMB'    : '/afs/cern.ch/user/y/yousen/work/pPb2016/MB0-150-PD0-v2vsNtrk',
      'PAHM0'   : '/afs/cern.ch/user/y/yousen/work/pPb2016/HM150-185-PD0-v2vsNtrk',
      'PAHM1-6' : '/afs/cern.ch/user/y/yousen/work/pPb2016/HM185-250-PD0-v2vsNtrk',
      'PAHM7'   : '/afs/cern.ch/user/y/yousen/work/pPb2016/HM250-inf-PD0-v2vsNtrk',
      'PPMB'    : '/afs/cern.ch/user/y/yousen/work/pp2018/MB0-80-PD0-v2vsNtrk',
      'PPHM_1'  : '/afs/cern.ch/user/y/yousen/work/pp2018/HM80-100-PD0-v2vsNtrk',
      'PPHM_2'  : '/afs/cern.ch/user/y/yousen/work/pp2018/HM100-inf-PD0-v2vsNtrk'
      }

sublist_number = {
      'PAMB'    : 458,
      'PAHM0'   : 32,
      'PAHM1-6' : 146,
      'PAHM7'   : 50,
      'PPMB'    : 0,
      'PPHM_1'  : 0,
      'PPHM_2'  : 0
      }

if sublist_number[dataset_name] != 0:
   name = "submit_corr2D_trg_pd0_%s_%s_v2vsNtrk.jdl" % (tree[treeNumber], dataset_name)
   f = open(name, "w")

   command_lines = '''universe   = vanilla
getenv     = True
executable = submit_corr2D_trg_pd0_mult.sh
arguments  = list/%s.000 %s ../eff/fEff_PD0_d0ana.root %s
log        = log/submit_corr2D_trg_pd0_%s_%s_v2vsNtrk.$(Process).log
output     = out/submit_corr2D_trg_pd0_%s_%s_v2vsNtrk.$(Process).out
error      = err/submit_corr2D_trg_pd0_%s_%s_v2vsNtrk.$(Process).err
requirements = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "workday"
queue
''' % (dataset[dataset_name], dataset_name, storage[dataset_name],
      tree[treeNumber], dataset_name,tree[treeNumber], dataset_name, tree[treeNumber], dataset_name)

   for i in range(1, sublist_number[dataset_name]):
      temp = '''
arguments  = list/%s.%03d %s ../eff/fEff_PD0_d0ana.root %s
queue
''' % (dataset[dataset_name], i, dataset_name, storage[dataset_name])
      command_lines += temp

   f.write(command_lines)
   f.close()
   subprocess.call(["condor_submit", name])
else:
   print "enter the correct dataset name!"

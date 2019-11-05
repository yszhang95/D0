import subprocess

tree = ["PromptD0", "NonPromptD0", "NonPromptD0"]
treeNumber = 0

#dataset_name = 'PAHM1-6'
dataset_name = 'PAHM7'
#dataset_name = 'PAMB'

dataset = {
      'PAMB'    : 'PAMB0-185_REF_Sep23.list',
      'PAHM1-6' : 'PAHM185-250_REF_Sep23.list',
      'PAHM7'   : 'PAHM250-inf_REF_Sep23.list'
      }
storage = {
      'PAMB'   : '/afs/cern.ch/user/y/yousen/work/pPb2016/PAMB0-185-PD0-v2vsNtrk_MC',
      'PAHM1-6'   : '/afs/cern.ch/user/y/yousen/work/pPb2016/PAHM185-250-PD0-v2vsNtrk_MC',
      'PAHM7'  : '/afs/cern.ch/user/y/yousen/work/pPb2016/PAHM250-inf-PD0-v2vsNtrk_MC'
      }

sublist_number = {
      'PAMB'      : 500,
      'PAHM1-6'   : 400,
      'PAHM7'   : 500,
      }

if sublist_number[dataset_name] != 0:
   name = "submit_corr2D_trg_ref_%s_%s_v2vsNtrk.jdl" % (tree[treeNumber], dataset_name)
   f = open(name, "w")

   command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = submit_corr2D_trg_ref_mult.sh
Arguments  = list/%s.000 %s %s
log        = log/submit_corr2D_trg_ref_%s_%s_v2vsNtrk.$(Process).log
output     = out/submit_corr2D_trg_ref_%s_%s_v2vsNtrk.$(Process).out
error      = err/submit_corr2D_trg_ref_%s_%s_v2vsNtrk.$(Process).err
requirements = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "tomorrow"
Queue
''' % (dataset[dataset_name], dataset_name, storage[dataset_name],
      tree[treeNumber], dataset_name,tree[treeNumber], dataset_name, tree[treeNumber], dataset_name)

   for i in range(1, sublist_number[dataset_name]):
      temp = '''
Arguments  = list/%s.%03d %s %s
Queue
''' % (dataset[dataset_name], i, dataset_name, storage[dataset_name])
      command_lines += temp

   f.write(command_lines)
   f.close()
   subprocess.call(["condor_submit", name]);
else:
   print "enter the correct dataset name!"

import subprocess

tree = ["PromptD0", "NonPromptD0", "NonPromptD0"]
treeNumber = 0

pTMin = 1.5
pTMax = 8.0
yMin  = -1.0
yMax  = 1.0

dataset_name = 'PAHM1-6'
#dataset_name = 'PAMB'

dataset = {
      'PAMB' : 'PAMB0-150.list', # too lazy to change it
      #'PAHM1-6' : 'PAHM185-250.list',
      'PAHM1-6' : 'PAHM185-250_new.list',
      'PPHM'    : 'PPHM80-inf.list'
      }
storage = {
      'PAMB' : '/afs/cern.ch/user/y/yousen/work/pPb2016/MB0-185-PD0-v2vspt',
      'PAHM1-6' : '/afs/cern.ch/user/y/yousen/work/pPb2016/HM185-250-PD0-v2vspt',
      'PPHM'  : '/afs/cern.ch/user/y/yousen/work/pp2018/HM80-inf-PD0-v2vspt'
      }

sublist_number = {
      'PAMB'    : 458,
      #'PAHM1-6' : 146,
      'PAHM1-6' : 240,
      'PPHM'  : 0,
      }

if sublist_number[dataset_name] != 0:
   name = "submit_corr2D_trg_pd0_%s_%s_HM_pT%.1f-%.1f_y%.1f-%.1f.jdl" % (tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax)
   f = open(name, "w")

   command_lines = '''universe   = vanilla
executable = submit_corr2D_trg_pd0.sh
arguments  = list/%s.000 %s ../eff/fEff.root %s %.1f %.1f %.1f %.1f
log        = log/submit_corr2D_trg_pd0_%s_%s_HM_pT%.1f-%.1f_y%.1f-%.1f.$(Process).log
output     = out/submit_corr2D_trg_pd0_%s_%s_HM_pT%.1f-%.1f_y%.1f-%.1f.$(Process).out
error      = err/submit_corr2D_trg_pd0_%s_%s_HM_pT%.1f-%.1f_y%.1f-%.1f.$(Process).err
requirements = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "workday"
queue
''' % (dataset[dataset_name], dataset_name, storage[dataset_name], pTMin, pTMax, yMin, yMax,
      tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax,
      tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax,
      tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax)

   for i in range(1, sublist_number[dataset_name]):
      temp = '''
arguments  = list/%s.%03d %s ../eff/fEff.root %s %.1f %.1f %.1f %.1f
queue
''' % (dataset[dataset_name], i, dataset_name, storage[dataset_name], pTMin, pTMax, yMin, yMax)
      command_lines += temp

   f.write(command_lines)
   f.close()
   subprocess.call(["condor_submit", name])
else:
   print "enter the correct dataset name!"

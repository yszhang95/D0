import subprocess

tree = ["PromptD0", "NonPromptD0", "NonPromptD0"]
treeNumber = 0

pTMin = 2.0
pTMax = 8.0
yMin  = -1.0
yMax  = 1.0

dataset_name = 'PAHM1-6'
#dataset_name = 'PAHM7'
#dataset_name = 'PAMB'

dataset = {
      'PAMB'    : 'PAMB0-185_Aug26.list',
      'PAHM1-6' : 'PAHM185-250_Aug26.list',
      'PAHM7'   : 'PAHM250-inf_Aug26.list',
      }
#storage = {
#      'PAMB'    : '/afs/cern.ch/user/y/yousen/work/pPb2016/PAMB0-185-PD0-v2vsNtrk',
#      'PAHM1-6' : '/afs/cern.ch/user/y/yousen/work/pPb2016/PAHM185-250-PD0-v2vsNtrk',
#      'PAHM7'   : '/afs/cern.ch/user/y/yousen/work/pPb2016/PAHM250-inf-PD0-v2vsNtrk',
#      }

storage = {
      'PAMB'    : '/afs/cern.ch/user/y/yousen/work/pPb2016/PAMB0-185-PD0-v2vsNtrk_MC',
      'PAHM1-6' : '/afs/cern.ch/user/y/yousen/work/pPb2016/PAHM185-250-PD0-v2vsNtrk_MC',
      'PAHM7'   : '/afs/cern.ch/user/y/yousen/work/pPb2016/PAHM250-inf-PD0-v2vsNtrk_MC',
      }

sublist_number = {
      'PAMB'    : 498,
      'PAHM1-6' : 191,
      'PAHM7'   : 37,
      }

path = {
      'PAMB' : '/afs/cern.ch/user/y/yousen/work/pPb2016',
      'PAHM1-6' : '/afs/cern.ch/user/y/yousen/work/pPb2016',
      'PAHM7' : '/afs/cern.ch/user/y/yousen/work/pPb2016',
      }

if sublist_number[dataset_name] != 0:
   name = "submit_corr2D_trg_pd0_mult_%s_%s_v2vsNtrk_pT%.1f-%.1f_y%.1f-%.1f_MC.jdl" % (tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax)
   f = open(name, "w")

   command_lines = '''universe   = vanilla
getenv     = True
executable = submit_corr2D_trg_pd0_mult.sh
arguments  = list/%s.000 %s ../eff/fEff.root %s %.1f %.1f %.1f %.1f %s
log        = log/submit_corr2D_trg_pd0_%s_%s_v2vsNtrk_pT%.1f-%.1f_y%.1f-%.1f.$(Process).log
output     = out/submit_corr2D_trg_pd0_%s_%s_v2vsNtrk_pT%.1f-%.1f_y%.1f-%.1f.$(Process).out
error      = err/submit_corr2D_trg_pd0_%s_%s_v2vsNtrk_pT%.1f-%.1f_y%.1f-%.1f.$(Process).err
requirements = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "workday"
queue
''' % (dataset[dataset_name], dataset_name, storage[dataset_name], pTMin, pTMax, yMin, yMax, path[dataset_name],
      tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax,
      tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax,
      tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax)

   for i in range(1, sublist_number[dataset_name]):
      temp = '''
arguments  = list/%s.%03d %s ../eff/fEff.root %s %.1f %.1f %.1f %.1f %s
queue
''' % (dataset[dataset_name], i, dataset_name, storage[dataset_name], pTMin, pTMax, yMin, yMax, path[dataset_name])
      command_lines += temp

   f.write(command_lines)
   f.close()
   subprocess.call(["condor_submit", name])
else:
   print "enter the correct dataset name!"

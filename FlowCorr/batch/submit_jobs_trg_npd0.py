import subprocess

tree = ["PromptD0", "NonPromptD0", "NonPromptD0"]
treeNumber = 2

#pTMin = 2.0
pTMin = 2.5
pTMax = 8.0
yMin  = -1.0
yMax  = 1.0

dcaCut_CSV = "dcacut.csv"

dataset_name = 'PAHM1-6'

dataset = {
      'PAHM1-6' : 'PAHM185-250.list',
      }
storage = {
      'PAHM1-6' : '/afs/cern.ch/user/y/yousen/work/pPb2016/HM185-250-NPD0-v2vspt',
      }

sublist_number = {
      'PAHM1-6' : 146,
      }

if sublist_number[dataset_name] != 0:
   name = "submit_corr2D_trg_npd0_%s_%s_HM_pT%.1f-%.1f_y%.1f-%.1f.jdl" % (tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax)
   f = open(name, "w")

   command_lines = '''universe   = vanilla
getenv     = True
executable = submit_corr2D_trg_npd0.sh
arguments  = list/%s.000 %s ../eff/fEff.root %s %.1f %.1f %.1f %.1f %s
log        = log/submit_corr2D_trg_npd0_%s_%s_HM_pT%.1f-%.1f_y%.1f-%.1f.$(Process).log
output     = out/submit_corr2D_trg_npd0_%s_%s_HM_pT%.1f-%.1f_y%.1f-%.1f.$(Process).out
error      = err/submit_corr2D_trg_npd0_%s_%s_HM_pT%.1f-%.1f_y%.1f-%.1f.$(Process).err
requirements = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "workday"
queue
''' % (dataset[dataset_name], dataset_name, storage[dataset_name], pTMin, pTMax, yMin, yMax, dcaCut_CSV,
      tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax,
      tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax,
      tree[treeNumber], dataset_name, pTMin, pTMax, yMin, yMax)

   for i in range(1, sublist_number[dataset_name]):
      temp = '''
arguments  = list/%s.%03d %s ../eff/fEff.root %s %.1f %.1f %.1f %.1f %s
queue
''' % (dataset[dataset_name], i, dataset_name, storage[dataset_name], pTMin, pTMax, yMin, yMax, dcaCut_CSV)
      command_lines += temp

   f.write(command_lines)
   f.close()
   subprocess.call(["condor_submit", name])
else:
   print "enter the correct dataset name!"

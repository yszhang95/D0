import subprocess

tree = ["PromptD0", "NonPromptD0", "NonPromptD0"]
treeNumber = 0

name = "submit_corr2D_trg_ref_%s_HM.jdl" % tree[treeNumber]
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = submit_corr2D_trg_ref.sh
Arguments  = list/PAHM185-250_REF.list.000 PAHM1-6 /afs/cern.ch/user/y/yousen/work/pPb2016/HM185-250-PD0-v2vspt %d
Log        = log/submit_corr2D_trg_ref_%s_HM1-6.$(Process).log
Output     = out/submit_corr2D_trg_ref_%s_HM1-6.$(Process).out
Error      = err/submit_corr2D_trg_ref_%s_HM1-6.$(Process).err
requirements = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "workday"
Queue
''' % (treeNumber, tree[treeNumber],tree[treeNumber],tree[treeNumber])

for i in range(1, 500):
   temp = '''
Arguments  = list/PAHM185-250_REF.list.%03d PAHM1-6 /afs/cern.ch/user/y/yousen/work/pPb2016/HM185-250-PD0-v2vspt %d
Queue
   ''' % (i, treeNumber)
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);

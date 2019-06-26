import subprocess

tree = ["PromptD0", "NonPromptD0", "NonPromptD0"]
treeNumber = 0
#treeNumber = 2

name = "submit_corr2D_trg_ref_%s.jdl" % tree[treeNumber]
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = submit_corr2D_trg_ref.sh
Arguments  = 00 %d
Log        = log/submit_corr2D_trg_ref_%s.$(Process).log
Output     = out/submit_corr2D_trg_ref_%s.$(Process).out
Error      = err/submit_corr2D_trg_ref_%s.$(Process).err
requirements = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "workday"
Queue
''' % (treeNumber, tree[treeNumber],tree[treeNumber],tree[treeNumber])

for i in range(1, 90):
   temp = '''
Arguments  = %02d %d
Queue
   ''' % (i, treeNumber)
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);

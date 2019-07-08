import subprocess

tree = ["PromptD0", "NonPromptD0", "NonPromptD0"]
#treeNumber = 0
treeNumber = 2

name = "submit_corr2D_trg_d0_%s.jdl" % tree[treeNumber]
f = open(name, "w")

command_lines = '''universe   = vanilla
getenv     = True
executable = submit_corr2D_trg_d0.sh
arguments  = 000 %d
log        = log/submit_corr2D_trg_d0_%s.$(Process).log
output     = out/submit_corr2D_trg_d0_%s.$(Process).out
error      = err/submit_corr2D_trg_d0_%s.$(Process).err
requirements = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "workday"
queue
''' % (treeNumber, tree[treeNumber],tree[treeNumber],tree[treeNumber])

for i in range(1, 168):
   temp = '''
arguments  = %03d %d
queue
   ''' % (i, treeNumber)
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name])

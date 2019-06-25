import subprocess
name = "submit_corr2D_trg_d0.jdl"
f = open(name, "w")

command_lines = '''universe   = vanilla
getenv     = True
executable = submit_corr2D_trg_d0.sh
arguments  = 00 0
log        = log/submit_corr2D_trg_d0.$(Process).log
output     = out/submit_corr2D_trg_d0.$(Process).out
error      = err/submit_corr2D_trg_d0.$(Process).err
requirements = (OpSysAndVer =?= "CentOS7")
+JobFlavour           = "workday"
queue
'''

for i in range(1, 89):
   temp = '''
arguments  = %02d 0
queue
   ''' % i
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);

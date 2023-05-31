#!/usr/bin/env python
import CompareFile as comf
import os
import sys
def red(txt):
    txt="\033[31m"+txt+"\033[0m"
    return txt
def green(txt):
    txt="\033[32m"+txt+"\033[0m"
    return txt
def yellow(txt):
    txt="\033[33m"+txt+"\033[0m"
    return txt
def check(filepair, threshold):
    testpass = True
    for files in filepair:
        testpass = testpass and comf.comparefile(files[0],files[1],threshold)
    return testpass
def runcandela(np):
    if np > 1:
        ierr=os.system("OMP_NUM_THREADS=2 ASAN_OPTIONS=protect_shadow_gap=0 mpirun -np "+str(2)+" ../../bin/candela > /dev/null")
    else:
        ierr=os.system("OMP_NUM_THREADS=2 ASAN_OPTIONS=protect_shadow_gap=0 ../../bin/candela > /dev/null")
    if ierr == 0:
        return True
    else:
        return False

if __name__ == "__main__":
    lines=[]
    comf.ReadFile("compare",lines)
    threshold = 0.01
    enable_mpi = False
    mpicompile = sys.argv[1]
    filepair=[]
    # Read compare
    for line in lines:
        linelist=line.replace('\n','').split()
        # print(linelist)
        if len(linelist) == 0:
            continue
        elif len(linelist) != 2:
            print("Wrong compare file.")
            exit(1)
        elif linelist[0] == "threshold":
            threshold = float(linelist[1])
        elif linelist[0] == "enable_mpi":
            if linelist[1] == "ON" or linelist[1] == "on":
                enable_mpi = True
        else:
            filepair.append(linelist)
    # Check
    testpass=True
    testpass = testpass and runcandela(1)
    testpass = testpass and check(filepair, threshold)
    if testpass:
        print("  One Processor      "+green( "PASS"))
    else:
        print("  One Processor      "+red("FAIL"))
    if enable_mpi and (mpicompile == "ON" or mpicompile == "on"):
        testpass = testpass and runcandela(2)
        testpass = testpass and check(filepair, threshold)
        testpass = testpass and runcandela(3)
        testpass = testpass and check(filepair, threshold)
        if testpass:
            print("  Multi Processor    "+green("PASS"))
        else:
            print("  Multi Processor    "+red("FAIL"))

#!/usr/bin/env python
import CompareFile as comf
import os
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
        ierr=os.system("mpirun -np "+str(2)+" ../../bin/candela > /dev/null")
    else:
        ierr=os.system("../../bin/candela > /dev/null")
    if ierr == 0:
        return True
    else:
        return False

if __name__ == "__main__":
    lines=[]
    comf.ReadFile("compare",lines)
    threshold = 0.01
    enable_mpi = False
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
        print("  Serial version   "+green( "PASS"))
    else:
        print("  Serial version   "+red("FAIL"))
    if enable_mpi:
        testpass = testpass and runcandela(2)
        testpass = testpass and check(filepair, threshold)
        testpass = testpass and runcandela(3)
        testpass = testpass and check(filepair, threshold)
        if testpass:
            print("  MPI version      "+green("PASS"))
        else:
            print("  MPI version      "+red("FAIL"))

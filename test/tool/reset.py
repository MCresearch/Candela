#!/usr/bin/env python
import CompareFile as comf
import os

if __name__ == "__main__":
    lines=[]
    comf.ReadFile("compare",lines)
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
        elif linelist[0] != "threshold" and linelist[0] != "enable_mpi":
            filepair.append(linelist)
            
    # Reset results
    for line in filepair:
        os.system('mv '+line[1]+' '+line[0])

#!/usr/bin/env python
import os,sys

# comparefile(file1,file2,epsilon)
# file1, file2: files to compare
# threshold: the threshold to check if two numbers are equal


def ReadFile(file1,lines):
    if os.path.isfile(file1):
        with open(file1) as f1: 
            for line in f1.readlines(): lines.append(line)
    else:
        print("Error: can not find file %s" % file1)
        sys.exit(1)

def IsFloat(x):
    try:
        float(x)
        return True
    except:
        return False

def IsComplex(x):
    x = str(x.strip())
    if x == '': return False
    if x[0] == '(' and x[-1] == ')' and len(x.split(',')) == 2:
        realp = x.split(',')[0][1:]
        imagp = x.split(',')[1][:-1]
        try:
            realp = float(realp)
            imagp = float(imagp)
        except:
            return False
    else:
        return False
    return [realp,imagp]

def ExitError(file1, file2, iline,line1,line2,jnumber=-1):
    if jnumber < 0:
        print('Error: line %d\n %s:%s %s:%s'%(iline,file1,line1,file2,line2))
    else:
        print('Error: line %d, column %d \n %s: %s\n %s: %s'%(iline,jnumber+1,file1,line1,file2,line2))
    return False


def comparefile(file1,file2,threshold):
    lines1 = []
    lines2 = []
    ReadFile(file1,lines1)
    ReadFile(file2,lines2)
   
    for i in range(min(len(lines1),len(lines2))):
        if lines1[i].strip() == lines2[i].strip():
            continue
        elif '' in [lines1[i].strip(),lines2[i].strip()]: 
            return ExitError(file1,file2,i,lines1[i],lines2[i])
        elif len(lines1[i].split()) != len(lines2[i].split()): 
            return ExitError(file1,file2,i,lines1[i],lines2[i])
        else:
            sline1 = lines1[i].split()
            sline2 = lines2[i].split()
            for j in range(len(sline1)):
                if IsComplex(sline1[j]) and IsComplex(sline2[j]):
                    dat1 = (x1[0]**2+x1[1]**2)**0.5
                    dat2 = (x2[0]**2+x2[1]**2)**0.5
                    if abs(dat1-dat2)/dat1 > threshold: 
                        return ExitError(file1,file2,i,sline1[j],sline2[j],j)
                elif IsFloat(sline1[j]) and IsFloat(sline2[j]):
                    dat1 = float(sline1[j])
                    dat2 = float(sline2[j])
                    if abs(dat1-dat2)/dat1 > threshold: 
                        return ExitError(file1,file2,i,sline1[j],sline2[j],j)
                else: 
                    return ExitError(file1,file2,i,sline1[j],sline2[j],j)
    return True

#!/usr/bin/env python
# -*-coding:utf-8 -*-

__author__ = "Meng Ye <yemeng77@gmail.com>"
__licence__ = "GNU General Public License"
__version__ = "0.0.1"

import numpy as np
from scipy.interpolate import UnivariateSpline

def readEline():
    k_st = []
    E_st = []    
    with open("E_line_W2K2", "r") as fEline:
        line = fEline.readline()
        line=line.strip()
        data = line.split()
        [nband,nkpt] = map(int, data)
        for i in range(nband):
            line = fEline.readline()
            line=line.strip()
            data = line.split()
            [iband,ikpt] = map(int, data)
            iline = ikpt/5
            if(ikpt%5 != 0):
                iline = iline + 1
            tmp = []
            for j in range(iline):
                line = fEline.readline()
                data = line.strip().split()
                tmp.extend(map(float, data))
            E_st.append(tmp)
            for j in range(iline):
                line = fEline.readline()
            tmp = []
            for j in range(iline):
                line = fEline.readline()
                data = line.strip().split()
                tmp.extend(map(float, data))
            k_st.append(tmp)
    return nband,k_st,E_st
    

def Elinespline(nband,k_st,E_st,Ew):
    iband = []
    kroot = []
    Eder = []
    for band_index in range(nband):
        if len(k_st[band_index]) <= 3:
            continue        
        k = k_st[band_index]
        E = map(lambda x: x - Ew, E_st[band_index])
        spl = UnivariateSpline(k,E,s=0)
        kroots = spl.roots()
        k_roots = [x for x in kroots if (x>min(k) and x < max(k))]        
        if len(k_roots) > 0:
            iband.extend([(band_index+1)]*len(k_roots))
            kroot.extend(k_roots)
            spld = spl.derivative()
            E_der = spld(k_roots)
            Eder.extend(E_der)
    return iband,kroot,Eder


def main():
    with open("etot.input","r") as f:
        line = f.readlines()[-1].strip()
    data = line.split()
    E1 = float(data[2].strip(","))
    E2 = float(data[3].strip(","))
    Ntrans = int(data[4].strip(","))
    dV = float(data[6].strip(","))
    nband,k_st,E_st = readEline()
    with open("dEdk_line_W2K2","w") as f:
        for Ew in np.linspace(E1,E2,Ntrans,endpoint=True):
            print >> f, "{:16.12f}".format(Ew)
            iband,kroot,Eder = Elinespline(nband,k_st,E_st,Ew+dV/2.0)
            print >> f, "{:16.12f}\t{:6d}".format(Ew+dV/2.0,len(iband))
            for i in range(len(iband)):
                print >> f, "{:6d}\t{:16.12f}\t{:.12e}".format(iband[i],kroot[i],Eder[i])
            iband,kroot,Eder = Elinespline(nband,k_st,E_st,Ew-dV/2.0)
            print >> f, "{:16.12f}\t{:6d}".format(Ew-dV/2.0,len(iband))
            for i in range(len(iband)):
                print >> f, "{:6d}\t{:16.12f}\t{:.12e}".format(iband[i],kroot[i],Eder[i])


if __name__ == '__main__':
    main()

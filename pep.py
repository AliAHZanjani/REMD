#!/usr/bin/env python

import os
import math
import sys
import numpy as np
import subprocess
import collections as coll
from   multiprocessing import Pool
global lambdas
global nodeList

n_replicas  = 64
rep_indices = range(n_replicas)

try:
    nodeList = subprocess.check_output(['cat $OAR_NODEFILE'],shell=True)
    nlist    = nodeList.decode()
    nlist    = nlist.split('\n')

    if len(nlist)-1 != n_replicas:
        print("warning! replicas doesn't match number of nodes.")  
    text_file = open("Out0put.txt", "w")
    for i in range(len(nlist)-1):
        text_file.write("nodelist[%d] = %s \n" %(i,nlist[i]))
    text_file.close()
       
except Exception as e:
    nodeList = []
    print("did not get nodelist: "+str(e))
lambdas = np.array([0.000000, 0.000008, 0.000032, 0.000073, 0.000129, 0.000202, 0.000290, 0.000395, 0.000516, 0.000653, 0.000806, 0.000976, 0.001161, 0.001363, 0.001580, 0.001814, 0.002064, 0.002330, 0.002612, 0.002911, 0.003225, 0.003556, 0.003902, 0.004265, 0.004644, 0.005039, 0.005450, 0.005878, 0.006321, 0.006781, 0.007256, 0.007748, 0.008256, 0.008780, 0.009320, 0.009877, 0.010449, 0.011038, 0.011642, 0.012263, 0.012900, 0.013553, 0.014222, 0.014908, 0.015609, 0.016327, 0.017060, 0.017810, 0.018576, 0.019358, 0.020156, 0.020971, 0.021801, 0.022648, 0.023510, 0.024389, 0.025284, 0.026195, 0.027122, 0.028066, 0.029025, 0.030001, 0.030992, 0.032000])
PathDir   = "/work/users/ahakamizanjani/Jobs/IL7kinematics=/"
for i in range(n_replicas):
    l=np.loadtxt(PathDir+"Replica%d/lseq.dat"%i)
    if (np.size(l) > 1):
        lambdas[i]=l[len(l)-1]
    else:
        lambdas[i]=l
#===============================================================================

segment_t   = "0.03"
Rgc         = 0.0019872036 # kcal/(K molecule)
beta	    = 1.0 / (Rgc * 300.0)

def rep(my_index) :

    if len(nlist) != 0:
        cmd  = "oarsh %s " % nlist[my_index]
    else:
        cmd  = ""

    exefile   = PathDir+"peptide1"
    initop    = PathDir+"Peptides.top"
    inifile   = PathDir+"Replica%d/Peptides.rst7" % my_index
    traj      = PathDir+"Replica%d/trajectory.traj" % my_index
    pot       = PathDir+"Replica%d/potentiel.dat" % my_index
    lseq      = PathDir+"Replica%d/lseq.dat" % my_index
    lseqAd    = PathDir+"Replica%d/lseqAd.dat" % my_index
    ghost     = PathDir+"Replica%d/ghost.pdb" % my_index
    outfile    = PathDir+"Replica%d/nabMD.dat" % my_index

    re_lambda = "%f" % lambdas[my_index]

    arg	      = subprocess.check_output (["oarsh", "%s " % nlist[my_index] ,exefile, "-iniConf", inifile, "-iniTop", initop, "-q", "1", "-Traj", traj, "-Pot", pot, "-t", segment_t, "-l", re_lambda, "-lseq", lseq, "-lseqAd", lseqAd, "-outf", outfile, "-gh", ghost, "|grep -v 'Using carbon'"])


    arg    = arg.split()
    retval = np.asarray([float(arg[0]), float(arg[1]), float(arg[2])])
    return retval
#########################################################################
def repLambda(my_index):
    l=np.loadtxt(PathDir+"Replica"+str(my_index)+"/lseqAd.dat")

    last100l=l[np.size(l)-5:np.size(l)]
    last500l=l[np.size(l)-25:np.size(l)]
    c=0
    for i in range (np.size(last100l)-1):
        if (last100l[i] != last100l[i+1]):
            c += 1
    ls=np.sort(last500l);
    un=np.unique(ls);
    cou=coll.Counter(ls)
    return [my_index, c, cou, l]
#########################################################################
for count in range(25) :
    pool = Pool(processes = 64)
    y_parallel = pool.map(rep, rep_indices)
    pool.close()
    pool.join()
    ePot       = np.asarray(y_parallel)  
    Delta      = np.zeros((n_replicas,n_replicas))
    for i in range(n_replicas) :
        for j in range(n_replicas) :
            Delta[i,j] =  ((lambdas[j] * ePot[i,1] - lambdas[i] * ePot[i,1]) +\
                           (lambdas[i] * ePot[j,1] - lambdas[j] * ePot[j,1]))
    n_exch       = 10 * n_replicas * n_replicas
    neglogflips  = -1 * np.log(np.random.random((n_exch)))
    pair_seq     = np.random.randint(0,n_replicas,(n_exch,2)) 

    for ii in range(n_exch) :
        if pair_seq[ii,0] == pair_seq[ii,1] :
            continue
        i = pair_seq[ii,0]
        j = pair_seq[ii,1]
        if Delta[i,j] <= 0. or neglogflips[ii] > beta * Delta[i,j] :
            Delta[i,j] *= -1
            Delta[j,i] *= -1
            temp       = lambdas[i]
            lambdas[i] = lambdas[j]
            lambdas[j] = temp

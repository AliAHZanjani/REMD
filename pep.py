#!/usr/bin/env python3

import os
import shutil
import math
import sys
import numpy as np
import subprocess
import collections as coll
from   multiprocessing import Pool

haveScoop=True
try:
    from scoop import futures
except:
    print("couldn't import scoop")
    haveScoop=False

global lambdas
global nodeList

"""
This is a master program to run the replica exchange, using the SLURM / OAR schedulers on a supercomputer,
or no scheduler (just fork processes).

"""

##parallelisation strategy requires some inputs
PathDir     = "./"
n_replicas  =  64
first_start = "./Peptides.rst7"

##nanoseconds
segment_t   = "0.020"

##1 femtosecond timestep, seems to be needed even with shake.
min_dt_ps   = 0.001

use_sched   = "None"
if len(sys.argv) > 1:
    for arg in sys.argv[1:]:
        #print("processing arg: ", arg)
        a1 = arg.split("=")[0]
        a2 = arg.split("=")[1]
        if   a1 == "PathDir":
                PathDir = a2
        elif a1 == "n_replicas":
                n_replicas = int(a2)
        elif a1 == "segment_t":
                segment_t  = a2
        elif a1 == "use_sched":
                use_sched  = a2
        else:
             print("unrecognised argument: '%s'" % arg)
else:
    print("received no command line arguments")

##set some (global) constants
Rgc         = 0.0019872036 # kcal/(K molecule)
beta	    = 1.0 / (Rgc * 300.0)

rep_indices = list(np.arange(n_replicas, dtype=int))

print("testing environment")
nodelist = []
if use_sched == "OAR":
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
       print("did not get OAR nodelist: "+str(e))
elif use_sched == "SLURM":
   ##no code needed here.
   #print("using SLURM scheduler")
   if haveScoop is False:
       print("Can't run on slurm without scoop task manager.")
       exit(1)
else:
   print("using no scheduler, just forking jobs.")
   

##(positive part of) this list of reaction coordinates was the product of some optimisation by Ali.
lambdas = np.array([ -0.000290, -0.000202, -0.000129, -0.000073, -0.000032, -0.000008,
0.000000, 0.000008, 0.000032, 0.000073, 0.000129, 0.000202, 0.000290, 0.000395, 0.000516, 0.000653, 0.000806, 0.000976, 0.001161, 0.001363, 0.001580, 0.001814, 0.002064, 0.002330, 0.002612, 0.002911, 0.003225, 0.003556, 0.003902, 0.004265, 0.004644, 0.005039, 0.005450, 0.005878, 0.006321, 0.006781, 0.007256, 0.007748, 0.008256, 0.008780, 0.009320, 0.009877, 0.010449, 0.011038, 0.011642, 0.012263, 0.012900, 0.013553, 0.014222, 0.014908, 0.015609, 0.016327, 0.017060, 0.017810, 0.018576, 0.019358, 0.020156, 0.020971, 0.021801, 0.022648, 0.023510, 0.024389, 0.025284, 0.026195, 0.027122, 0.028066, 0.029025, 0.030001, 0.030992, 0.032000])

lambdas = lambdas[:n_replicas]

##experimental: try some T-remd with H-remd.
Tscales = np.ones_like(lambdas)
Tscales[0] = 1.16
Tscales[1] = 1.15
Tscales[2] = 1.14
Tscales[3] = 1.12
Tscales[4] = 1.10
Tscales[5] = 1.05

print("testing environment")
nodelist = []
if use_sched == "OAR":
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
       print("did not get OAR nodelist: "+str(e))
elif use_sched == "SLURM":
   ##no code needed here.
   #print("using SLURM scheduler")
   if haveScoop is False:
       print("Can't run on slurm without scoop task manager.")
       exit(1)
else:
   print("using no scheduler, just forking jobs.")
   

##(positive part of) this list of reaction coordinates was the product of some optimisation by Ali.
lambdas = np.array([ -0.000290, -0.000202, -0.000129, -0.000073, -0.000032, -0.000008,
0.000000, 0.000008, 0.000032, 0.000073, 0.000129, 0.000202, 0.000290, 0.000395, 0.000516, 0.000653, 0.000806, 0.000976, 0.001161, 0.001363, 0.001580, 0.001814, 0.002064, 0.002330, 0.002612, 0.002911, 0.003225, 0.003556, 0.003902, 0.004265, 0.004644, 0.005039, 0.005450, 0.005878, 0.006321, 0.006781, 0.007256, 0.007748, 0.008256, 0.008780, 0.009320, 0.009877, 0.010449, 0.011038, 0.011642, 0.012263, 0.012900, 0.013553, 0.014222, 0.014908, 0.015609, 0.016327, 0.017060, 0.017810, 0.018576, 0.019358, 0.020156, 0.020971, 0.021801, 0.022648, 0.023510, 0.024389, 0.025284, 0.026195, 0.027122, 0.028066, 0.029025, 0.030001, 0.030992, 0.032000])
lambdas = lambdas[:n_replicas]
>>>>>>> cbdd19894f6c77f14ceaae3b043dd1b2a0af004f


#===============================================================================

def rep( in_tup ) :

    ##print("in rep: ", my_index)
    my_index   = in_tup[0]
    swap_index = in_tup[1]

    re_lambda   = "%f" % lambdas[my_index]
    re_tscale   = "%f" % Tscales[my_index]

    exe       = PathDir+"../peptide1"
    initop    = PathDir+"Peptides.top"
    ghost     = PathDir+"ghost.pdb"                 
    template  = PathDir+"template.pdb"                 
    ini_rst   = PathDir+"Replica%d/Peptides.rst7"   % swap_index
    out_rst   = PathDir+"Replica%d/Peptides.rst7"   % my_index
    traj      = PathDir+"Replica%d/trajectory.traj" % my_index
    pot       = PathDir+"Replica%d/potentiel.dat"   % my_index
    outfile   = PathDir+"Replica%d/nabMD.dat"       % my_index


    ##############################################3
    ##Main run: use OAR scheduler, or just fork.

    argSet = [exe, \
               "-iniConf", ini_rst, "-outConf", out_rst, \
               "-iniTop", initop, "-q", "1",\
               "-Traj", traj, "-Pot", pot, "-t", segment_t, "-l", re_lambda, \
               "-sy",  template, "-dt", "%.8f" % min_dt_ps, \
               "-outf", outfile, "-gh", ghost, "-T", re_tscale]


    try:
         arg = subprocess.check_output( [" ".join(argSet)], shell=True   )
         ##arg = subprocess.check_output (argSet + ["|grep -v 'Using carbon'"])
    except Exception as e:  
         print("failed subprocess.  Args were: ")
         for a in argSet:
             print("%s " % a, end="")
         print("\n")
         raise( e ) 
    ####################################       
    try:
        arg    = arg.split()
        retval = np.asarray([float(arg[0]), float(arg[1]), float(arg[2])])
    except Exception as e:
        print("couldn't process returned data from subprocess, this was:", arg)
        print("Expecting three floats, related to energies.")
        raise( e )

    return retval


##serial below this point: above code is thread-global initialisation
if __name__ == "__main__":

    for l in lambdas:
       est_dt = 0.002
       if l > 0.:
          est_w  = math.sqrt(  0.695 * l / 2e-26);
          est_dt = 1e12 * 0.1 * 2 * 3.14159 / est_w;
       elif l < 0.:
          est_w  = math.sqrt( -0.695 * l / 2e-26);
          est_dt = 1e12 * 0.1 * 2 * 3.14159 / est_w;
       print("estimated timesteps given restraint strengths %.6f = %.8e ps:" % (l, est_dt))


###set up replicas if not present.
    for my_index in range(n_replicas) :
        inifile   = PathDir+"Replica%d/Peptides.rst7"   % my_index
        if not os.path.isfile(inifile):
            print("initting replica dir: %i" % my_index)
            try:
               if not os.path.exists(PathDir+"/Replica%d" % my_index):
                   os.makedirs(PathDir+"/Replica%d" % my_index)
               shutil.copy2(first_start, inifile)
            except Exception as e:
               raise( e )

#########################################################################
## MAIN LOOP
#########################################################################
    try:
        with open( PathDir+"/swap_record.dat", "r" ) as swapLog:
            for line in swapLog: 
               pass
        L = line.split()
        swapSet = np.zeros((len(L)))
        for i in range(len(L)):
           swapSet[i] = int(L[i])
    except:
        swapSet = np.arange(n_replicas)
        swapLog = open( PathDir+"/swap_record.dat", "w" )
        swapLog.write("# ")
        for s in swapSet:
           swapLog.write("%.6f " % lambdas[s])
        swapLog.write("\n")
        swapLog.close()

    ##by default each run starts from its own restart.
    inp_indices  = np.arange((len(swapSet)))
    ##see $AMBERHOME/lib/constants.f90
    amber_kB  = 0.001987215873 
    amber_kBT = 300 * amber_kB

    for count in range(10) :

        y_parallel = list(futures.map( rep,  zip(rep_indices, inp_indices) ))

        #####do the replica exchange
        ePot       = np.asarray(y_parallel)  

        ##more exchange checks == faster equilibration.  Nothing wrong with doing 10N^2 exchanges.
        n_exch       = 10 * n_replicas * n_replicas

        ##RNG call for flip probabilities
        neglogflips  = -1 * np.log(np.random.random((n_exch)))

        ##random set of pairs to exchange.
        pair_seq     =  np.random.randint(0,n_replicas,(n_exch,2)) 

        ##by default each run starts from its own restart.
        inp_indices  = np.arange((len(swapSet)))

        for ii in range(n_exch) :
            if pair_seq[ii,0] == pair_seq[ii,1] :
                continue
            i = pair_seq[ii,0]
            j = pair_seq[ii,1]

            ibetai = amber_kBT * Tscales[i]
            ibetaj = amber_kBT * Tscales[j]

            betai = 1./ ibetai
            betaj = 1./ ibetaj

            ##if we do the swap, there is a change in the restraint energy
            delta = (lambdas[i]*ePot[j,1]+ePot[j,0]+ePot[j,2])*betai - (lambdas[i]*ePot[i,1]+ePot[i,0]+ePot[i,2])*betai +\
                    (lambdas[j]*ePot[i,1]+ePot[i,0]+ePot[i,2])*betaj - (lambdas[j]*ePot[j,1]+ePot[j,0]+ePot[j,2])*betaj

            if delta <= 0. or neglogflips[ii] > delta:
                swapSet[i], swapSet[j] = swapSet[j], swapSet[i]
                ePot[i,1],   ePot[j,1] = ePot[j,1],   ePot[i,1]
                ePot[i,0],   ePot[j,0] = ePot[j,0],   ePot[i,0]
                ePot[i,2],   ePot[j,2] = ePot[j,2],   ePot[i,2]
                inp_indices[i], inp_indices[j] = inp_indices[j], inp_indices[i]

        ###need something like an atomic write of the swap info: if it gets corrupted we 
        ###have some major hassle.
        swapLog = open( PathDir+"/swap_record.dat", "a" )
        for s in swapSet:
             swapLog.write("%i " % s)
        swapLog.write("\n") 
        swapLog.close()
###########################END of file.




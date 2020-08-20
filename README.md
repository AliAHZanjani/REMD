# REMD
Replica Exchange Molecular Dynamics for Peptide Aggregation

This code is associated to the PhD project of Ali Zanjani, thesis available here: https://orbilu.uni.lu/handle/10993/40678

If you use or are inspired by this code, please cite:

Ali Asghar Hakami Zanjani, Nicholas P. Reynolds, Afang Zhang, Tanja Schilling, Raffaele Mezzenga, Joshua T. Berryman,
Amyloid Evolution: Antiparallel Replaced by Parallel,
Biophysical Journal,
Volume 118, Issue 10,
2020,
Pages 2526-2536,
ISSN 0006-3495,
https://doi.org/10.1016/j.bpj.2020.03.023.

This code requires the free simulation package AmberTools to be installed, specifically it relies on the NAB molecular manipulation language (part of the AmberTools package).

The code runs a replica exchange calculation designed to accelerate aggregation of peptides into amyloid structures, without *in theory* biasing the observed structures.


## usage

To run a replica exchange calculation, the script pep.py is used.  This employs a parallel scheduler such as SLURM or OAR (only SLURM is recently tested), via the python "scoop" package, to make multiple short serial AMBER simulations, and then to redistribute start structures so as to accelerate convergence.

pep.py should also work without a scheduler, just forking processes using python multiprocessing: this is only for testing, in order to get converged results you will need to set the number of replicas to something like 70 (thus using minimum 70 cores) and run for a long time.  It should be possible to compile an MPI or even a GPU version of the serial part of the code, for further acceleration, but this has not yet been done.

Usage is to make a directory containing the files: Peptides.top, Peptides.rst7, ghost.pdb, template.pdb.  "ghost.pdb" is a pdb structure of a single peptide in an extended conformation.  "template.pdb" should be a pdb file of the whole system (N=20 to 60 peptides), the conformation doesn't matter.  Peptides.top and Peptides.rst7 are the system topology and starting configuration, in standard AMBER formats.  Having done this, you can run the replica exchange:

Without a scheduler, the command is: 

python $REMD_PATH/pep.py n_replicas=$YOUR_N_REPLICAS segment_t=0.010 

With SLURM via scoop, you need a launch script, something like:

```bash
##########################################################
#!/bin/bash -l
##get hostnames for this job:
HOSTFILE=$(pwd)/hostnames.dat
./getHosts.bash $SLURM_JOB_ID > $HOSTFILE

##build a wrapper for executing python with whatever environment needs to be set

SCOOP_WRAPPER=$(pwd)/scoop-python.sh

cat << EOF > $SCOOP_WRAPPER
#!/bin/bash -l
export SLURM_NTASKS=${SLURM_NTASKS}
EOF
echo 'python $@' >> $SCOOP_WRAPPER
chmod +x $SCOOP_WRAPPER

exe=../pep.py
python -m scoop --hostfile $HOSTFILE -n ${SLURM_NTASKS} --python-interpreter $SCOOP_WRAPPER \
           $exe n_replicas=${SLURM_NTASKS} segment_t=0.010 use_sched=SLURM 
##########################################################
```

To test the serial part of the code, you might find it useful to execute the nab binary directly:

peptide1  [-l lambda] [-q] [1]  [-boxR Radius] [-kWall wallConstant] [-t time(ns)] [-prStp printstep(ps)] -iniConf file.rst7 -iniTop file -Traj file.trj -Pot file.dat -lseq lambdaSeq.dat -lseqAd lambdaSeqAd.dat -gh ghost.pdb


## Inputs

Variables which might need to be changed for different systems are set to standard values in pep.py, then passed as arguments to peptide1.  Some effort was expended in finding a good set of lambda values; boundary conditions (size of the spherical confinement) however is a fairly arbitrary choice.  Work out the effective concentration that you want to simulate at, but beware that there is no acceleration of diffusion-driven sampling, so at low concentrations the code will spend a lot of time waiting for peptides to bump into eachother.

### Restraint strength parameter lambda

[-l lambda]
[-lseq lambdaseq.dat]

Either apply a single lambda, or load a set of lambdas from lambdaSeq.dat filename.

### Boundary conditions

[-boxR Radius]
[-kWall wallConstant]

Size and hardness of the spherical confinement of the simulation.

### Miscellaneous
 [ -q ] Partially surpress logging
 
### Run length
 
[-t time(ns)] [-prStp printstep(ps)]




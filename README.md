# REMD
Replica Exchange Molecular Dynamics for Peptide Aggregation

This code is associated to the PhD project of Ali Zanjani, thesis available here: https://orbilu.uni.lu/handle/10993/40678

It requires AMBER to be installed, specifically it relies on the NAB molecular manipulation language (provided in the ambertools package).

The code runs a replica exchange calculation designed to accelerate aggregation of peptides into amyloid structures, without *in theory* biasing the observed structures.


## usage

peptide1  [-l lambda] [-q] [1]  [-boxR Radius] [-kWall wallConstant] [-t time(ns)] [-prStp printstep(ps)] -iniConf file.rst7 -iniTop file -Traj file.trj -Pot file.dat -lseq lambdaSeq.dat -lseqAd lambdaSeqAd.dat -gh ghost.pdb


## inputs

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




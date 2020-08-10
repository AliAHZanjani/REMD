#!/bin/bash -l 

##tidy up previous tesst runs
rm -f out.log pot.dat test.trj
cp PTIYSLLL.rst.bak PTIYSLLL.rst

##do the run.
../peptide1 -t 0.01 -prStp 1 -iniConf PTIYSLLL.rst -iniTop PTIYSLLL_mb3.top \
	    -Traj test.trj -outf out.log -Pot pot.dat -lseq testLseq.dat -lseqAd testLseq_ad.dat \
	    -gh PTIYSLLL.pdb -sy PTIYSLLL.pdb

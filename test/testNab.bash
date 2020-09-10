#!/bin/bash -l 

##tidy up previous test runs
rm -f out.log pot.dat test.trj
cp FTIYSLLL.rst.bak FTIYSLLL.rst

##do the run.
../peptide1 -t 0.01 -prStp 1 -iniConf FTIYSLLL.rst -iniTop FTIYSLLL_mb3.top \
	    -Traj test.trj -outf out.log -Pot pot.dat -lseq testLseq.dat -lseqAd testLseq_ad.dat \
	    -gh FTIYSLLL.pdb -sy FTIYSLLL.pdb -S 88779876

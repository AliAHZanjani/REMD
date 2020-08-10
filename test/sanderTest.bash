#!/bin/bash -l 


sander=$AMBERHOME/bin/sander


cat > sander_test.in  <<EOF
 test of gbsa dynamics
 &cntrl
   nstlim=10, cut=99.0, igb=5,
   ntpr=1, ntwr=1000, 
   ntt=3, gamma_ln = 1000., 
   ntc=2, ntf=2, tol=0.000001,
   ntb=0, 
   offset=0.09, gbsa=1, ig=71277,
   ioutfm=1,
 /
EOF



$sander -O -c PTIYSLLL.crd -p PTIYSLLL_mb3.top -x sander_test.trj -o sander.out -i sander_test.in -r PTIYSLLL.rstnc

$AMBERHOME/bin/cpptraj PTIYSLLL_mb3.top << EOF
trajin PTIYSLLL.rstnc
trajout PTIYSLLL.rst restart 
go
EOF


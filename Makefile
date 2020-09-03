
AMBERHOME = /home/users/jberryman/amber19/amber
CC = icc 

##load environemnt: "module load numlib/imkl"

### linking command used to build nab:
#icc -DCC='"icc"' -DCPP='"ucpp -l"' -DFLIBS='"-lsff -lpbsa -lrism -lfftw3 -larpack -llapack -lblas  -lnetcdf  -lifport -lifcore -lsvml "' \
	-std=gnu99 -fPIC  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ -DHASGZ -DHASBZ2 -D__PLUMED_HAS_DLOPEN      \
	-shared-intel



all: NabFuncs.nab peptide1.nab CFuncs.c 
	$(CC) -c -O3 -I $(AMBERHOME)/include CFuncs.c -o CFuncs.o
	$(AMBERHOME)/bin/nab	-c NabFuncs.nab -o NabFuncs.o
	$(AMBERHOME)/bin/nab    -c peptide1.nab -o peptide1.o
	icc -static -L$(AMBERHOME)/lib peptide1.o CFuncs.o NabFuncs.o -lnab -lsff -lpbsa -lrism -lfftw3 -larpack -llapack -lblas  -lnetcdf -lcifparse  -lifport -lifcore -lsvml -std=gnu99 -fPIC -shared-intel -o peptide1 

	


##icc -static -L$(AMBERHOME)/lib peptide1.o CFuncs.o NabFuncs.o -lnab -lsff -lpbsa -lrism -lfftw3 -larpack -llapack -lblas  -lnetcdf -lcifparse  -lifport -lifcore -lsvml -std=gnu99 -fPIC -shared-intel -o clinked

#$(AMBERHOME)/bin/nab -static CFuncs.o NabFuncs.nab peptide1.nab -o peptide1

.PHONY: clean
clean:
	rm NabFuncs.c peptide1.c *.o

.PHONY: testdir
testdir: all
	mkdir -p test
	cd test
	ln -s ../peptide1


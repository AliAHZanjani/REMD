
all: NabFuncs.nab peptide1.nab CFuncs.c 
	nab CFuncs.c NabFuncs.nab peptide1.nab -o peptide1

.PHONY: clean
clean:
	rm NabFuncs.c peptide1.c *.o




AMBERHOME = /home/josh/amberCheckout/amber20_checkout/amber20

all: NabFuncs.nab peptide1.nab CFuncs.c
	gcc -c -O4 -I $(AMBERHOME)/include CFuncs.c -o CFuncs.o
	$(AMBERHOME)/bin/nab CFuncs.o NabFuncs.nab peptide1.nab -o peptide1

.PHONY: clean
clean:
	rm NabFuncs.c peptide1.c *.o



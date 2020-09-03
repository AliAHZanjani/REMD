//C functions

INT_T   init_mymme( REAL_T mo_xyz[], MOLECULE_T * *strandmol, INT_T *nPeps, INT_T *numatm, INT_T *atmsPerPep, REAL_T *lambd, REAL_T *boxle, REAL_T *kwa,\
	        	INT_T heavy[], INT_T *outsteps,STRING_T * *OutTraj, STRING_T * *Poten );

REAL_T  my_mme(REAL_T x[], REAL_T f[], INT_T *iter);

INT_T   fileclose(STRING_T * *ghname);

INT_T	printepot();

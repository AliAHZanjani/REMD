//NAB FUNCTIONS

int getNewConf(float  oldGh_xyz[0],     /* input coordinates of old ghost */
               float    newGh_xyz[0],   /* output of getNewConf, new ghost coordinates written into this */
               molecule strand, string strn, string atomGroup);

int gtNewGhost(float molc_xyz[0], int nPeps, int atmsPerPep, molecule strand, float strand_xyz[0], molecule temp, float temp_xyz[0], float singleGh_xyz[0], float ghost_xyz[0]);

int CreateStrands(molecule molc, int nPeps, molecule stran[0], point mgc[0]);

int MergeStrands(molecule stran[0], int nPeps, molecule molc);

int init_mymme(float mo_xyz[0],molecule strandmol, int nPeps, int numatm,int atmsPerPep, float lambd, float boxle,float kwa, int heavy[0], int outsteps, string OutTraj, string Poten);

float my_mme(float x[0], float f[0], int iter);

int fileclose(string ghname);

int printepot();



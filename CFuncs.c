//standard nab includes
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
#include "CFuncs.h"

//standard nab externs
extern char NAB_rsbuf[];
static int  mytaskid, numtasks;

//globals for peptide code
static STRING_T     *stri1  = NULL;
static STRING_T     *atmGrp = NULL;
static MOLECULE_T	*global_strand, *global_temp, *global_tempmol;
static POINT_T		*global_centre;
static REAL_T		*globalm_xyz, *globalgh_xyz, *globalPgh_xyz, *globalstrand_xyz, *globaltemp_xyz, *global_oldGhxyz, *global_newGhxyz, *globalfP;
static REAL_T		global_lambda, global_boxl, global_kwall;
static INT_T		*global_IsHeavy, global_nstr, global_natom, global_nsatom, counter,outStep;
static FILE_T		*potfout, *tfout, *tgout;
static REAL_T		gEPOT, gEWALL,gEGHOST;
static INT_T		outcrd_fieldCount, ccccc;

//function to save globals
INT_T init_mymme( REAL_T mo_xyz[], MOLECULE_T * *strandmol, INT_T *nPeps, INT_T *numatm, INT_T *atmsPerPep, REAL_T *lambd, REAL_T *boxle, REAL_T *kwa, INT_T heavy[], INT_T *outsteps, STRING_T * *OutTraj, STRING_T * *Poten ){
    INT_T i;
    INT_T nPeps_buf;
    INT_T numatm_buf;
    INT_T atmsPerPep_buf;
    nPeps_buf      = *nPeps;
    numatm_buf     = *numatm;
    atmsPerPep_buf = *atmsPerPep;
    outStep        = *outsteps;
//-------------------------------------------------------------------
    DA_ALLOC( global_IsHeavy = ( INT_T * )malloc( 3 * numatm_buf * (sizeof( INT_T))), "init_mymme", "global_IsHeavy" );
    DA_ALLOC( globalm_xyz = ( REAL_T * )malloc( 3 * numatm_buf * ( sizeof( REAL_T ) ) ), "init_mymme", "globalm_xyz" );
    DA_ALLOC( globalgh_xyz = ( REAL_T * )malloc( 3 * numatm_buf * ( sizeof( REAL_T ) ) ), "init_mymme", "globalgh_xyz" );
    DA_ALLOC( globalPgh_xyz = ( REAL_T * )malloc( 3 * numatm_buf * ( sizeof( REAL_T ) ) ), "init_mymme", "globalPgh_xyz" );
    DA_ALLOC( globalstrand_xyz = ( REAL_T * )malloc( 3 * atmsPerPep_buf * ( sizeof( REAL_T ) ) ), "init_mymme", "globalstrand_xyz" );
    DA_ALLOC( globaltemp_xyz = ( REAL_T * )malloc( 3 * atmsPerPep_buf * ( sizeof( REAL_T ) ) ), "init_mymme", "globaltemp_xyz" );
    DA_ALLOC( global_oldGhxyz = ( REAL_T * )malloc( 3 * atmsPerPep_buf * ( sizeof( REAL_T ) ) ), "init_mymme", "global_oldGhxyz" );
    DA_ALLOC( global_newGhxyz = ( REAL_T * )malloc( 3 * atmsPerPep_buf * ( sizeof( REAL_T ) ) ), "init_mymme", "global_newGhxyz" );
    DA_ALLOC( global_centre = ( POINT_T * )malloc( nPeps_buf * ( sizeof( POINT_T ) ) ), "init_mymme", "global_centre" );
    DA_ALLOC( globalfP = ( REAL_T * )malloc( 3 * numatm_buf * ( sizeof( REAL_T ) ) ), "init_mymme", "globalfP" );
    global_strand = newmolecule(  );
	global_strand = copymolecule(*strandmol);
    setxyz_from_mol(&global_strand, NULL, global_oldGhxyz);
    global_temp   = newmolecule(  );
    global_temp   = copymolecule(*strandmol);
    setxyz_from_mol(&global_temp, NULL, global_newGhxyz);
//-------------------------------------------------------------------
	for( i=1; i<=3 * (*numatm);i++)
		global_IsHeavy[i-1] = heavy[i-1];
	global_nstr   = *nPeps;
    global_natom  = *numatm;
    global_nsatom = *atmsPerPep;
    global_lambda = *lambd;
    global_boxl   = *boxle;
    global_kwall  = *kwa;
//global_previousEpGhost = INIT_EPGHOST;
	counter=0;ccccc=0;
    potfout = safe_fopen( *Poten, "a" );
	fprintf( potfout, "counter\t	eGhost/lambda\t    eGhost\t    eWall\t     eW+eGh\t    ePot\t      ePotTotal\n");
	tfout = safe_fopen( *OutTraj,"a" );
    outcrd_fieldCount = 0;
    return( i );
};
//-------------------------------------------------------------------
REAL_T my_mme(REAL_T x[], REAL_T f[], INT_T *iter){
    INT_T	i,j;
	REAL_T	epTot, epGhost, epGhostP, delE, epGhost1, epG0, df, dfP, r2, r, invr, neglog;
//-------------------------------------------------------------------
	counter++;
	epTot = mme(x, f, iter);
    for (i = 0; i < 3*global_natom; i++){
        globalfP[i] = f[i];
        }
    gtNewGhost( x,  &global_nstr,  &global_nsatom, &global_strand, globalstrand_xyz,  &global_temp, globaltemp_xyz, global_newGhxyz, globalgh_xyz);
    epG0	 = 0.0;
	epGhost  = 0.0;
	epGhostP = 0.0;
    epGhost1 = 0.0;
	for (i = 0; i < 3*global_natom; i++){
        if(global_IsHeavy[i]){
    		df           =  globalgh_xyz[i]-x[i];
            f[i]        -=  df * global_lambda;
            epGhost     +=  df * df;
        }
    }
	epG0	 = 0.5 * epGhost;
    epGhost *= 0.5 * global_lambda;
    epGhost1 = epGhost;  
//=====================================Spherical restraining potential
	for (i = 0; i < 3*global_natom; i += 3){
		r2 = x[i]*x[i] + x[i+1]*x[i+1] + x[i+2]*x[i+2];
        
        if (r2>(global_boxl*global_boxl)){
            r        =  sqrt(r2);
            invr     =  1./r;
            df       =  global_kwall*(global_boxl-r);
            f[i]    -=  df*x[i]*invr;
            f[i+1]  -=  df*x[i+1]*invr;
            f[i+2]  -=  df*x[i+2]*invr;
            epGhost +=  0.5*global_kwall*(r-global_boxl)*(r-global_boxl);
        }
    }
//=====================================
	if (counter % outStep == 1){
		fprintf( potfout, "%-7d\t%12.7f\t%12.7f\t%12.7f\t%12.7f\t%12.7f\t%12.7f\n",counter, epG0, epGhost1, epGhost-epGhost1, epGhost, epTot, epTot+epGhost );
        
        gEPOT   = epTot;
        gEGHOST = epG0;
        gEWALL  = epGhost-epGhost1;        
    }
    if (counter % (10*outStep) == 1){
        for( j = 1;j <= global_natom;j ++  ){
            fprintf( tfout, "%8.3lf", x[3 * j - 2 - 1] );
            outcrd_fieldCount += 1;
            if( outcrd_fieldCount % 10 == 0 ){
                fprintf( tfout, "\n" );
            }
            fprintf( tfout, "%8.3lf", x[3 * j - 1 - 1] );
            outcrd_fieldCount += 1;
            if( outcrd_fieldCount % 10 == 0 ){
                fprintf( tfout, "\n" );
            }
            fprintf( tfout, "%8.3lf", x[3 * j - 1] );
            outcrd_fieldCount += 1;
            if( outcrd_fieldCount % 10 == 0 ){
                fprintf( tfout, "\n" );
            }
        }
        if( outcrd_fieldCount % 10 != 0 ){
            fprintf( tfout, "\n" );
        }
        outcrd_fieldCount = 0;
    }
//=====================================
    return(epTot + epGhost);
};
INT_T   fileclose(STRING_T * *ghname){
	fclose( potfout );
	fclose( tfout );
	return 0;
};

INT_T   printepot(){
    printf("%f %f %f", gEPOT, gEGHOST, gEWALL);
    return 0;
};

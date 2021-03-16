/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdio.h>
#include<math.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"
#include "cell.h"
#include "vesicle.h"
#include "io.h"
#include "initial_distribution.h"
#include "frame.h"
#include "timestep.h"
#include "sh.h"

/** Entrance function to the program
  * @param argv is a number of parameters used in program call (including the program name
  * @param argc is a pointer to strings (character arrays) which holds the arguments
  * @returns returns 0 on success, any other number on fail.
*/
ts_bool saveAvgUlm2(ts_vesicle *vesicle);
int main(int argv, char *argc[]){
ts_uint i,j,k;
ts_vesicle *vesicle;
ts_double r0;
vesicle=initial_distribution_dipyramid(17,60,60,60,0.15);
//parsetape(vesicle,&i);

//similar to nmax in fortran code
ts_uint nmax;

//these four must come from parsetype!
vesicle->dmax=1.67*1.67;
vesicle->stepsize=0.15;
vesicle->clist->max_occupancy=8;
vesicle->bending_rigidity=25.0;
//fprintf(stderr,"xk=%f",vesicle->bending_rigidity);

	centermass(vesicle);
cell_occupation(vesicle);

//test if the structure is internally organized into cells correctly 
ts_uint cind;
for(i=0;i<vesicle->vlist->n;i++){
	cind=vertex_self_avoidance(vesicle, vesicle->vlist->vtx[i]);

	if(vesicle->clist->cell[cind]==vesicle->vlist->vtx[i]->cell){
		//fprintf(stdout,"(T) Idx match!\n");
	} else {
		fprintf(stderr,"(T) ***** Idx don't match!\n");

	}
}
//end test
vesicle->sphHarmonics=sph_init(vesicle->vlist, 21);

vesicle_volume(vesicle);
r0=getR0(vesicle);

preparationSh(vesicle,r0);
calculateYlmi(vesicle);
calculateUlm(vesicle);

//preloop:
ts_double vmsr, bfsr;
for(i=0;i<1000;i++){
	cell_occupation(vesicle);
	for(j=0;j<1000;j++){
		single_timestep(vesicle, &vmsr, &bfsr);
	}	
	centermass(vesicle);
	fprintf(stderr, "Preloop %d completed.\n",i+1);
}

nmax=1000;
for(i=0;i<nmax;i++){
	for(j=0;j<200;j++){
		cell_occupation(vesicle);
		for(k=0;k<5;k++){
		single_timestep(vesicle, &vmsr, &bfsr);
		}
		centermass(vesicle);
	}	
    vesicle_volume(vesicle);
    r0=getR0(vesicle);

    preparationSh(vesicle,r0);
    calculateYlmi(vesicle);
    calculateUlm(vesicle);

    storeUlm2(vesicle);
    saveAvgUlm2(vesicle);

	write_vertex_xml_file(vesicle,i);
	fprintf(stderr, "Loop %d out of %d completed.\n",i+1,nmax);

}

write_master_xml_file("test.pvd");
write_dout_fcompat_file(vesicle,"dout");
vesicle_free(vesicle);

return 0; //program finished perfectly ok. We return 0.
}




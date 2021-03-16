/* vim: set ts=4 sts=4 sw=4 noet : */
#include "general.h"
#include "sh.h"
#include <stdlib.h>
#include "vesicle.h"
#include "initial_distribution.h"
#include "frame.h"
#include <math.h>

int main(){
ts_vesicle *vesicle=initial_distribution_dipyramid(17,60,60,60,0.15);
//parsetape(vesicle,&i);

//these four must come from parsetype!
vesicle->dmax=1.67*1.67;
vesicle->stepsize=0.15;
vesicle->clist->max_occupancy=8;
vesicle->bending_rigidity=25.0;
fprintf(stderr,"xk=%f\n",vesicle->bending_rigidity);

centermass(vesicle);
vesicle->sphHarmonics=sph_init(vesicle->vlist, 21);
int i,j;
ts_double area;
for(i=1;i<=vesicle->sphHarmonics->l;i++){
    for(j=1;j<=2*i+1;j++){
    fprintf(stderr,"co(%d,%d)=%e\n",i,j,vesicle->sphHarmonics->co[i][j]);
    }
}
ts_double r0;
vesicle_volume(vesicle);
fprintf(stderr,"Volume=%e\n",vesicle->volume);
r0=getR0(vesicle);
fprintf(stderr,"r0=%e\n",r0);
area=0;
for(i=0;i<vesicle->tlist->n;i++){
	area+=vesicle->tlist->tria[i]->area;
}
fprintf(stderr,"area_dipyramid=%e\n",area);
fprintf(stderr,"Centroid=(%e,%e,%e)\n", vesicle->cm[0],vesicle->cm[1],vesicle->cm[2]);

preparationSh(vesicle,r0);
calculateYlmi(vesicle);
/*
ts_coord *coord=(ts_coord *)malloc(sizeof(ts_coord));
ts_double fi, theta;

  for(i=0;i<vesicle->vlist->n;i++){

    cart2sph(coord,vesicle->vlist->vtx[i]->x, vesicle->vlist->vtx[i]->y, vesicle->vlist->vtx[i]->z);
        fi=coord->e2;
        theta=coord->e3; 


	fprintf(stderr,"VTX(x,y,z,fi,theta)=%e,%e,%e,%e,%e ---> Ylmi(2,-2,%d)=%9.7e <---- data: omega=%e, r0=%e, plgndr(2,abs(-2),cos(theta))=%e, co(2,-2)=%e cos((m-l-1)*fi)=%e\n",vesicle->vlist->vtx[i]->x, vesicle->vlist->vtx[i]->y, vesicle->vlist->vtx[i]->z,fi,theta,i+1,vesicle->sphHarmonics->Ylmi[2][0][i], vesicle->vlist->vtx[i]->solAngle, vesicle->vlist->vtx[i]->relR, plgndr(2,abs(-2),cos(theta)), vesicle->sphHarmonics->co[2][1], cos(-2*fi));

	}

free(coord);
*/
calculateUlm(vesicle);
for(i=0;i<vesicle->sphHarmonics->l;i++){
    for(j=0;j<2*i+1;j++){
    fprintf(stderr,"ulm(%d,%d)=%e\n",i,j+1,vesicle->sphHarmonics->ulm[i][j]);
    }
}

storeUlm2(vesicle);

sph_free(vesicle->sphHarmonics);
vesicle_free(vesicle);
return 0;
}

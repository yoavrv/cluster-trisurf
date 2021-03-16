/* vim: set ts=4 sts=4 sw=4 noet : */
#include "general.h"
#include "vertex.h"
#include "initial_distribution.h"
#include "io.h"
#include "vesicle.h"
#include "sh.h"
#include "frame.c"
#include <math.h>
#include <stdlib.h>
int main(int argc, char *argv[]){

ts_fprintf(stdout,"SHdiscover was called with %d coefficients!\n",argc-1);
ts_uint n,i,j,l;
ts_int m;
ts_double fi,theta,r,Y;
ts_vesicle *vesicle=initial_distribution_dipyramid(17,60,60,60,0.15);
ts_vertex_list *vlist=vesicle->vlist;
centermass(vesicle);
ts_fprintf(stdout,"Vesicle has a CenterMass in %f,%f,%f\n",vesicle->cm[0],vesicle->cm[1], vesicle->cm[2]);

n=vlist->n;

ts_fprintf(stdout,"Tests\n");
ts_fprintf(stdout,"P(0,0,0.5)=%f (%f)\n",plgndr(0,0,0.5),1.0);
ts_fprintf(stdout,"P(1,0,0.5)=%f (%f)\n",plgndr(1,0,0.5),0.5);
ts_fprintf(stdout,"P(2,0,0.5)=%f (%f)\n",plgndr(2,0,0.5),0.5*(3*0.5*0.5-1));
ts_fprintf(stdout,"P(2,2,0.5)=%f (ni to:%f)\n",plgndr(2,2,0.5),0.5*(3*0.5*0.5-1));

ts_fprintf(stdout,"Y(0,0,pi/6,pi/4)=%f (%f)\n",shY(0,0,M_PI/6,M_PI/4),sqrt(1/(4*M_PI)));
ts_fprintf(stdout,"Y(1,0,pi/6,pi/4)=%f (%f)\n",shY(1,0,M_PI/6,M_PI/4),sqrt(3/(4*M_PI))*cos(M_PI/6));
ts_fprintf(stdout,"Y(1,0,4*pi/6,6*pi/4)=%f (%f)\n",shY(1,0,4*M_PI/6,6*M_PI/4),sqrt(3/(4*M_PI))*cos(4*M_PI/6));
ts_fprintf(stdout,"Y(1,1,pi/6,pi/4)=%f (%f)\n",shY(1,1,M_PI/6,M_PI/4),-sqrt(3/(8*M_PI))*sin(M_PI/6)*cos(M_PI/4));
ts_fprintf(stdout,"Y(2,0,pi/6,pi/4)=%f (%f)\n",shY(2,0,M_PI/6,M_PI/4),sqrt(5/(4*M_PI))*(3.0/2.0*cos(M_PI/6)*cos(M_PI/6)-1.0/2.0));
ts_fprintf(stdout,"Y(2,-2,pi/6,pi/4)=%f (0)\n",shY(2,-2,M_PI/6,M_PI/4));
ts_fprintf(stdout,"Y(2,2,pi/6,pi/3)=%f (%f)\n",shY(2,2,M_PI/6,M_PI/3), sqrt(15.0/(32.0*M_PI))*sin(M_PI/6)*sin(M_PI/6)*cos(2*M_PI/3));
	
	for(j=1;j<argc;j++){
		l=(int)sqrt(j-1); /* determine l from dataline */
		m=j-1-l*(l+1); /* determine m from dataline */
		ts_fprintf(stdout,"l=%d, m=%d, u=%s\n",l,m,argv[j]);
	}

/*we calculate new position of each vertex of vesicle */
for(i=0;i<n;i++){
	fi=atan2(vlist->vtx[i]->y, vlist->vtx[i]->x);
/*	theta=atan2(
	    sqrt(vlist->vtx[i]->data->x*vlist->vtx[i]->data->x + 
		vlist->vtx[i]->data->y*vlist->vtx[i]->data->y),
		vlist->vtx[i]->data->z 
	    ); */
	theta=acos(
		vlist->vtx[i]->z /
	    sqrt(vlist->vtx[i]->x*vlist->vtx[i]->x + 
		vlist->vtx[i]->y*vlist->vtx[i]->y+
		vlist->vtx[i]->z*vlist->vtx[i]->z)

		);



	r=0.0;
	for(j=1;j<argc;j++){
		l=(int)sqrt(j-1); /* determine l from dataline */
		m=j-1-l*(l+1); /* determine m from dataline */
		Y=shY(l,m,theta,fi);
		r+=fabs(atof(argv[j])*Y);
		/*ts_fprintf(stdout,"l=%d, m=%d, u=%s\n",l,m,argv[j]);*/
	}

	vlist->vtx[i]->z=fabs(r)*cos(theta);
	vlist->vtx[i]->x=fabs(r)*sin(theta)*cos(fi);
	vlist->vtx[i]->y=fabs(r)*sin(theta)*sin(fi);
}

write_vertex_xml_file(vesicle,0);
write_master_xml_file("test.pvd");


vesicle_free(vesicle);
return 0;
}

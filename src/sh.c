/* vim: set ts=4 sts=4 sw=4 noet : */
#include<math.h>
#include<stdlib.h>
#include "general.h"
#include "sh.h"
#include "io.h"
#include <string.h>


ts_spharm *sph_init(ts_vertex_list *vlist, ts_uint l){
    ts_uint j,i;
    ts_spharm *sph=(ts_spharm *)malloc(sizeof(ts_spharm));

    sph->N=0;
    /* lets initialize Ylm for each vertex. */
    sph->Ylmi=(ts_double ***)calloc(l,sizeof(ts_double **));
    for(i=0;i<l;i++){
            sph->Ylmi[i]=(ts_double **)calloc(2*i+1,sizeof(ts_double *));
            for(j=0;j<(2*i+1);j++){
                sph->Ylmi[i][j]=(ts_double *)calloc(vlist->n,sizeof(ts_double));
            }
    }
        
    /* lets initialize ulm */
    sph->ulm=(ts_double **)calloc(l,sizeof(ts_double *));
    for(j=0;j<l;j++){
        sph->ulm[j]=(ts_double *)calloc(2*j+1,sizeof(ts_double));
    }

    /* lets initialize sum of Ulm2 */
    sph->sumUlm2=(ts_double **)calloc(l,sizeof(ts_double *));
    for(j=0;j<l;j++){
        sph->sumUlm2[j]=(ts_double *)calloc(2*j+1,sizeof(ts_double));
    }

    /* lets initialize co */
//NOTE: C is has zero based indexing. Code is imported from fortran and to comply with original indexes we actually generate one index more. Also second dimension is 2*j+2 instead of 2*j+2. elements starting with 0 are useles and should be ignored!
    sph->co=(ts_double **)calloc(l+1,sizeof(ts_double *));
    for(j=0;j<=l;j++){
        sph->co[j]=(ts_double *)calloc(2*j+2,sizeof(ts_double));
    }

    sph->l=l;   

    /* Calculate coefficients that will remain constant during all the simulation */ 
   precomputeShCoeff(sph);
    
    return sph;
}


ts_bool sph_free(ts_spharm *sph){
    int i,j;
    if(sph==NULL) return TS_FAIL;
    for(i=0;i<sph->l;i++){
        if(sph->ulm[i]!=NULL) free(sph->ulm[i]);
        if(sph->sumUlm2[i]!=NULL) free(sph->sumUlm2[i]);
        if(sph->co[i]!=NULL) free(sph->co[i]);
    }
        if(sph->co[sph->l]!=NULL) free(sph->co[sph->l]);
    if(sph->co != NULL) free(sph->co);
    if(sph->ulm !=NULL) free(sph->ulm);

        if(sph->Ylmi!=NULL) {
            for(i=0;i<sph->l;i++){
                if(sph->Ylmi[i]!=NULL){
                    for(j=0;j<i*2+1;j++){
                        if(sph->Ylmi[i][j]!=NULL) free (sph->Ylmi[i][j]);
                    }
                    free(sph->Ylmi[i]);
                }
            }
            free(sph->Ylmi);
        }

    free(sph);
    return TS_SUCCESS;
}

/* Gives you legendre polynomials. Taken from NR, p. 254 */
ts_double plgndr(ts_int l, ts_int m, ts_double x){
	ts_double fact, pll, pmm, pmmp1, somx2;
	ts_int i,ll;

#ifdef TS_DOUBLE_DOUBLE
	if(m<0 || m>l || fabs(x)>1.0)
		fatal("Bad arguments in routine plgndr",1);
#endif
#ifdef TS_DOUBLE_FLOAT
	if(m<0 || m>l || fabsf(x)>1.0)
		fatal("Bad arguments in routine plgndr",1);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
	if(m<0 || m>l || fabsl(x)>1.0)
		fatal("Bad arguments in routine plgndr",1);
#endif
	pmm=1.0;
	if (m>0) {
#ifdef TS_DOUBLE_DOUBLE
		somx2=sqrt((1.0-x)*(1.0+x));
#endif
#ifdef TS_DOUBLE_FLOAT
		somx2=sqrtf((1.0-x)*(1.0+x));
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
		somx2=sqrtl((1.0-x)*(1.0+x));
#endif
		fact=1.0;
		for (i=1; i<=m;i++){
			pmm *= -fact*somx2;
			fact +=2.0;
		}
	}

	if (l == m) return pmm;
	else {
		pmmp1=x*(2*m+1)*pmm;
		if(l==(m+1)) return(pmmp1);
		else {
			pll=0; /* so it can not be uninitialized */
			for(ll=m+2;ll<=l;ll++){
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return(pll);
		}
	}
}


/** @brief: Precomputes coefficients that are required for spherical harmonics computations.

*/
ts_bool precomputeShCoeff(ts_spharm *sph){
    ts_int i,j,al,am;
    ts_double **co=sph->co;
    for(i=1;i<=sph->l;i++){
        al=i;
        sph->co[i][i+1]=sqrt((2.0*al+1.0)/2.0/M_PI);
        for(j=1;j<=i-1;j++){
            am=j;
            sph->co[i][i+1+j]=co[i][i+j]*sqrt(1.0/(al-am+1.0)/(al+am));
            sph->co[i][i+1-j]=co[i][i+1+j];
        }
        co[i][2*i+1]=co[i][2*i]*sqrt(1.0/(2.0*al));
        co[i][1]=co[i][2*i+1];
        co[i][i+1]=sqrt((2.0*al+1.0)/4.0/M_PI);
    }
    return TS_SUCCESS;

}


/** @brief: Computes Y(l,m,theta,fi) 
 *
 * Function calculates Y^l_m for vertex with given (\theta, \fi) coordinates in
 * spherical coordinate system.
 * @param l is an ts_int argument.
 * @param m is an ts_int argument.
 * @param theta is ts_double argument.
 * @param fi is a ts_double argument.
 *
 * (Miha's definition that is different from common definition for  factor srqt(1/(2*pi)) */
ts_double shY(ts_int l,ts_int m,ts_double theta,ts_double fi){
	ts_double fac1, fac2, K;
	int i;

	if(l<0 || m>l || m<-l)
		fatal("Error using shY function!",1);

	fac1=1.0;
	for(i=1; i<=l-abs(m);i++){
		fac1 *= i;
	}
	fac2=1.0;
	for(i=1; i<=l+abs(m);i++){
		fac2 *= i;
	}

	if(m==0){
		K=sqrt(1.0/(2.0*M_PI));
	}
	else if (m>0) {
		K=sqrt(1.0/(M_PI))*cos(m*fi);
	} 
	else {
		//K=pow(-1.0,abs(m))*sqrt(1.0/(2.0*M_PI))*cos(m*fi);
		if(abs(m)%2==0)
		K=sqrt(1.0/(M_PI))*cos(m*fi);
		else
		K=-sqrt(1.0/(M_PI))*cos(m*fi);
	}
	
	return K*sqrt((2.0*l+1.0)/2.0*(ts_double)(fac1/fac2))*plgndr(l,abs(m),cos(theta));	
}


/* Function transforms coordinates from cartesian to spherical coordinates
 * (r,phi, theta). */
ts_bool *cart2sph(ts_coord *coord, ts_double x, ts_double y, ts_double z){
    coord->coord_type=TS_COORD_SPHERICAL;
#ifdef TS_DOUBLE_DOUBLE
    coord->e1=sqrt(x*x+y*y+z*z);
    if(z==0) coord->e3=M_PI/2.0;
    else coord->e3=atan2(sqrt(x*x+y*y),z);
    coord->e2=atan2(y,x);
#endif
#ifdef TS_DOUBLE_FLOAT
    coord->e1=sqrtf(x*x+y*y+z*z);
    if(z==0) coord->e3=M_PI/2.0;
    else coord->e3=atanf(sqrtf(x*x+y*y)/z);
    coord->e2=atan2f(y,x);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    coord->e1=sqrtl(x*x+y*y+z*z);
    if(z==0) coord->e3=M_PI/2.0;
    else coord->e3=atanl(sqrtl(x*x+y*y)/z);
    coord->e2=atan2l(y,x);
#endif

    return TS_SUCCESS;
}


ts_bool sph2cart(ts_coord *coord){
    coord->coord_type=TS_COORD_CARTESIAN;
    ts_double x,y,z;

    x=coord->e1*cos(coord->e2)*sin(coord->e3);
    y=coord->e1*sin(coord->e2)*sin(coord->e3);
    z=coord->e1*cos(coord->e3);

    coord->e1=x;
    coord->e2=y;
    coord->e3=z;

    return TS_SUCCESS;
}


/* Function returns radius of the sphere with the same volume as vesicle (r0) */
ts_double getR0(ts_vesicle *vesicle){
    ts_double r0;
 #ifdef TS_DOUBLE_DOUBLE
   r0=pow(vesicle->volume*3.0/4.0/M_PI,1.0/3.0);
#endif
#ifdef TS_DOUBLE_FLOAT
   r0=powf(vesicle->volume*3.0/4.0/M_PI,1.0/3.0);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
   r0=powl(vesicle->volume*3.0/4.0/M_PI,1.0/3.0);
#endif
    return r0;
}


ts_bool preparationSh(ts_vesicle *vesicle, ts_double r0){
//TODO: before calling or during the call calculate area of each triangle! Can
//be also done after vertexmove and bondflip //
//DONE: in energy calculation! //
    ts_uint i,j;
    ts_vertex **vtx=vesicle->vlist->vtx;
    ts_vertex *cvtx;
    ts_triangle *ctri;
    ts_double centroid[3];
    ts_double r;
    for (i=0;  i<vesicle->vlist->n; i++){
        cvtx=vtx[i];
        //cvtx->projArea=4.0*M_PI/1447.0*(cvtx->x*cvtx->x+cvtx->y*cvtx->y+cvtx->z*cvtx->z)/r0/r0;
        cvtx->projArea=0.0;

        /* go over all triangles that have a common vertex i */
        for(j=0; j<cvtx->tristar_no; j++){
            ctri=cvtx->tristar[j];
            centroid[0]=(ctri->vertex[0]->x + ctri->vertex[1]->x + ctri->vertex[2]->x)/3.0;
            centroid[1]=(ctri->vertex[0]->y + ctri->vertex[1]->y + ctri->vertex[2]->y)/3.0;
            centroid[2]=(ctri->vertex[0]->z + ctri->vertex[1]->z + ctri->vertex[2]->z)/3.0;
        /* calculating projArea+= area(triangle)*cos(theta) */
#ifdef TS_DOUBLE_DOUBLE
            cvtx->projArea = cvtx->projArea + ctri->area*(-centroid[0]*ctri->xnorm - centroid[1]*ctri->ynorm - centroid[2]*ctri->znorm)/ sqrt(centroid[0]*centroid[0]+centroid[1]*centroid[1]+centroid[2]*centroid[2]);
#endif
#ifdef TS_DOUBLE_FLOAT
            cvtx->projArea = cvtx->projArea + ctri->area*(-centroid[0]*ctri->xnorm - centroid[1]*ctri->ynorm - centroid[2]*ctri->znorm)/ sqrtf(centroid[0]*centroid[0]+centroid[1]*centroid[1]+centroid[2]*centroid[2]);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
            cvtx->projArea = cvtx->projArea + ctri->area*(-centroid[0]*ctri->xnorm - centroid[1]*ctri->ynorm - centroid[2]*ctri->znorm)/ sqrtl(centroid[0]*centroid[0]+centroid[1]*centroid[1]+centroid[2]*centroid[2]);
#endif
        }

    cvtx->projArea=cvtx->projArea/3.0;
        //we dont store spherical coordinates of vertex, so we have to calculate
        //r(i) at this point.
#ifdef TS_DOUBLE_DOUBLE
    r=sqrt(cvtx->x*cvtx->x+cvtx->y*cvtx->y+cvtx->z*cvtx->z);
#endif
#ifdef TS_DOUBLE_FLOAT
    r=sqrtf(cvtx->x*cvtx->x+cvtx->y*cvtx->y+cvtx->z*cvtx->z);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    r=sqrtl(cvtx->x*cvtx->x+cvtx->y*cvtx->y+cvtx->z*cvtx->z);
#endif
    cvtx->relR=(r-r0)/r0;
    cvtx->solAngle=cvtx->projArea/r/r;
    }
    return TS_SUCCESS;
}



ts_bool calculateYlmi(ts_vesicle *vesicle){
    ts_int i,j,k;
    ts_spharm *sph=vesicle->sphHarmonics;
    ts_coord *coord=(ts_coord *)malloc(sizeof(ts_coord));
    ts_double fi, theta;
	ts_int m;
    ts_vertex *cvtx;
    for(k=0;k<vesicle->vlist->n;k++){
        cvtx=vesicle->vlist->vtx[k];
        sph->Ylmi[0][0][k]=sqrt(1.0/4.0/M_PI);
        cart2sph(coord,cvtx->x, cvtx->y, cvtx->z);
        fi=coord->e2;
        theta=coord->e3; 
        for(i=1; i<sph->l; i++){
            for(j=0;j<i;j++){
			m=j+1;
//Nastudiraj!!!!!
                sph->Ylmi[i][j][k]=sph->co[i][m]*cos((m-i-1)*fi)*pow(-1,m-i-1)*plgndr(i,abs(m-i-1),cos(theta));
		if(i==2 && j==0){
	/*	fprintf(stderr," **** vtx %d ****\n", k+1);
		fprintf(stderr,"m-i-1 =%d\n",m-i-1);
		fprintf(stderr,"fi =%e\n",fi);
		fprintf(stderr,"(m-i-1)*fi =%e\n",((ts_double)(m-i-1))*fi);
		fprintf(stderr,"-2*fi =%e\n",-2*fi);
		fprintf(stderr,"m =%d\n",m);
	
		fprintf(stderr," cos(m-i-1)=%e\n",cos((m-i-1)*fi));
		fprintf(stderr," cos(-2*fi)=%e\n",cos(-2*fi));
		fprintf(stderr," sph->co[i][m]=%e\n",sph->co[i][m]);
		fprintf(stderr," plgndr(i,abs(m-i-1),cos(theta))=%e\n",plgndr(i,abs(m-i-1),cos(theta)));
*/
		}
            }
//Nastudiraj!!!!!
		j=i;
		m=j+1;
                sph->Ylmi[i][j][k]=sph->co[i][m]*plgndr(i,0,cos(theta));
            for(j=i+1;j<2*i+1;j++){
			m=j+1;
//Nastudiraj!!!!!
                sph->Ylmi[i][j][k]=sph->co[i][m]*sin((m-i-1)*fi)*plgndr(i,m-i-1,cos(theta));
            }
        }

    }
    free(coord);
    return TS_SUCCESS;
}



ts_bool calculateUlm(ts_vesicle *vesicle){
    ts_uint i,j,k;
    ts_vertex *cvtx;
    for(i=0;i<vesicle->sphHarmonics->l;i++){
        for(j=0;j<2*i+1;j++) vesicle->sphHarmonics->ulm[i][j]=0.0;
    }

//TODO: call calculateYlmi !!!


    for(k=0;k<vesicle->vlist->n; k++){
        cvtx=vesicle->vlist->vtx[k];
        for(i=0;i<vesicle->sphHarmonics->l;i++){
            for(j=0;j<2*i+1;j++){
                vesicle->sphHarmonics->ulm[i][j]+= cvtx->solAngle*cvtx->relR*vesicle->sphHarmonics->Ylmi[i][j][k];
            }

        }
    }

    return TS_SUCCESS;
}





ts_bool storeUlm2(ts_vesicle *vesicle){

ts_spharm *sph=vesicle->sphHarmonics;
ts_int i,j;
for(i=0;i<sph->l;i++){
    for(j=0;j<2*i+1;j++){
	/* DEBUG fprintf(stderr,"sph->sumUlm2[%d][%d]=%e\n",i,j,sph->ulm[i][j]* sph->ulm[i][j]); */
        sph->sumUlm2[i][j]+=sph->ulm[i][j]* sph->ulm[i][j];
    }
}
	sph->N++;
return TS_SUCCESS;
}


ts_bool saveAvgUlm2(ts_vesicle *vesicle){

	FILE *fh;
    char filename[10000];
    strcpy(filename, command_line_args.path);
    strcat(filename, "sph2out.dat");
	fh=fopen(filename, "w");
	if(fh==NULL){
		err("Cannot open file %s for writing");
		return TS_FAIL;
	}

	ts_spharm *sph=vesicle->sphHarmonics;
	ts_int i,j;
	fprintf(fh,"l,\tm,\tulm^2avg\n");
	for(i=0;i<sph->l;i++){
    		for(j=0;j<2*i+1;j++){
		fprintf(fh,"%d,\t%d,\t%e\n", i, j-i, sph->sumUlm2[i][j]/(ts_double)sph->N);

    		}
    fprintf(fh,"\n");
	}
	fclose(fh);
	return TS_SUCCESS;
}

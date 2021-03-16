/* vim: set ts=4 sts=4 sw=4 noet : */
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_sf_legendre.h>

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_linalg.h>
#include "general.h"
#include "sh.h"
#include "shcomplex.h"


ts_spharm *complex_sph_init(ts_vertex_list *vlist, ts_uint l){
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
    sph->ulmComplex=(gsl_complex **)calloc(l,sizeof(gsl_complex *));
    for(j=0;j<l;j++){
        sph->ulm[j]=(ts_double *)calloc(2*j+1,sizeof(ts_double));
        sph->ulmComplex[j]=(gsl_complex *)calloc(2*j+1,sizeof(gsl_complex));
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

ts_bool complex_sph_free(ts_spharm *sph){
    int i,j;
    if(sph==NULL) return TS_FAIL;
    for(i=0;i<sph->l;i++){
        if(sph->ulm[i]!=NULL) free(sph->ulm[i]);
        if(sph->ulmComplex[i]!=NULL) free(sph->ulmComplex[i]);
        if(sph->sumUlm2[i]!=NULL) free(sph->sumUlm2[i]);
        if(sph->co[i]!=NULL) free(sph->co[i]);
    }
        if(sph->co[sph->l]!=NULL) free(sph->co[sph->l]);
    if(sph->co != NULL) free(sph->co);
    if(sph->ulm !=NULL) free(sph->ulm);
    if(sph->ulmComplex !=NULL) free(sph->ulmComplex);
    if(sph->sumUlm2 !=NULL) free(sph->sumUlm2);
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


ts_bool calculateUlmComplex(ts_vesicle *vesicle){
    ts_int i,j,k,m,l;
    ts_vertex *cvtx;
    ts_coord coord;
/* set all values to zero */
    for(i=0;i<vesicle->sphHarmonics->l;i++){
        for(j=0;j<2*i+1;j++) GSL_SET_COMPLEX(&(vesicle->sphHarmonics->ulmComplex[i][j]),0.0,0.0);
    }

    for(k=0;k<vesicle->vlist->n; k++){
        cvtx=vesicle->vlist->vtx[k];
	cart2sph(&coord,cvtx->x,cvtx->y,cvtx->z);
        for(i=0;i<vesicle->sphHarmonics->l;i++){
            for(j=0;j<2*i+1;j++){
		m=j-i;
		l=i;
		if(m>=0){	
	//	fprintf(stderr, "Racunam za l=%d, m=%d\n", l,m);
                vesicle->sphHarmonics->ulmComplex[i][j]=gsl_complex_add(vesicle->sphHarmonics->ulmComplex[i][j], gsl_complex_conjugate(gsl_complex_mul_real(gsl_complex_polar(1.0,(ts_double)m*coord.e2),cvtx->solAngle*cvtx->relR*gsl_sf_legendre_sphPlm(l,m,cos(coord.e3)))) );
		} else {
	//	fprintf(stderr, "Racunam za l=%d, abs(m=%d)\n", l,m);
                vesicle->sphHarmonics->ulmComplex[i][j]=gsl_complex_add(vesicle->sphHarmonics->ulmComplex[i][j], gsl_complex_conjugate(gsl_complex_mul_real(gsl_complex_polar(1.0,(ts_double)m*coord.e2),cvtx->solAngle*cvtx->relR*pow(-1,m)*gsl_sf_legendre_sphPlm(l,-m,cos(coord.e3)))) );

		}
            }
        }
    }
    return TS_SUCCESS;
}

char *Ulm2Complex2String(ts_vesicle *vesicle){
	ts_int i,j;
	char *strng=(char *)calloc(5000, sizeof(char));
	char tmpstrng[255];
	for(i=0;i<vesicle->sphHarmonics->l;i++){
    		for(j=i;j<2*i+1;j++){
			sprintf(tmpstrng,"%e ", gsl_complex_abs2(vesicle->sphHarmonics->ulmComplex[i][j]));
			strcat(strng,tmpstrng);
    		}
	}
	//strcat(strng,"\n");

	return strng;
}

ts_bool freeUlm2String(char *strng){
	free(strng);
	return TS_SUCCESS;
}


ts_bool storeUlmComplex2(ts_vesicle *vesicle){

	ts_spharm *sph=vesicle->sphHarmonics;
	ts_int i,j;
	for(i=0;i<sph->l;i++){
    		for(j=0;j<2*i+1;j++){
        		sph->sumUlm2[i][j]+=gsl_complex_abs2(sph->ulmComplex[i][j]);
    		}
	}
	sph->N++;
	return TS_SUCCESS;
}


ts_double calculateKc(ts_vesicle *vesicle, ts_int lmin, ts_int lmax){
    ts_int min=lmin;
    ts_int max=lmax; //vesicle->sphHarmonics->l-3;
    ts_long i,j;
    ts_double retval, bval;
    gsl_matrix *A=gsl_matrix_alloc(max-min,2);
    gsl_vector *tau=gsl_vector_alloc(2);
    gsl_vector *b=gsl_vector_alloc(max-min);
    gsl_vector *x=gsl_vector_alloc(2);
    gsl_vector *res=gsl_vector_alloc(max-min);

    //solving (A^T*A)*x=A^T*b
    //fill the data for matrix A and vector b
    for(i=min;i<max;i++){
            gsl_matrix_set(A, i-min,0,(ts_double)((i-1)*(i+2)));
            gsl_matrix_set(A, i-min,1,(ts_double)((i-1)*(i+2)*(i+1)*i));
//            fprintf(stderr,"%e %e\n", gsl_matrix_get(A,i-min,0), gsl_matrix_get(A,i-min,1));
            bval=0.0;
            //average for m from 0..l (only positive m's)
            for(j=0;j<=i;j++){
                bval+=vesicle->sphHarmonics->sumUlm2[i][(j+i)];
            }
                bval=bval/(ts_double)vesicle->sphHarmonics->N/(ts_double)(i+1);

            gsl_vector_set(b,i-min,1.0/bval);
//            fprintf(stderr,"%e\n", 1.0/gsl_vector_get(b,i-min));
    }
//    fprintf(stderr,"b[2]=%e\n",gsl_vector_get(b,1));
    gsl_linalg_QR_decomp(A,tau);
    gsl_linalg_QR_lssolve(A,tau,b,x,res);
//    fprintf(stderr,"kc=%e\n",gsl_vector_get(x,1));
    retval=gsl_vector_get(x,1);
    gsl_matrix_free(A);
    gsl_vector_free(tau);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_vector_free(res);
    
    return retval;
}

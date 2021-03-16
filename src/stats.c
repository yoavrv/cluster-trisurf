/* vim: set ts=4 sts=4 sw=4 noet : */
#include <stdio.h>
#include <time.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include "general.h"
 
void gyration_eigen(ts_vesicle *vesicle, ts_double *l1, ts_double *l2, ts_double *l3){
	gsl_matrix *matrix = gsl_matrix_alloc(3,3);
	gsl_vector *eig = gsl_vector_alloc(3);
	gsl_eigen_symm_workspace *workspace = gsl_eigen_symm_alloc(3);
	ts_uint i,j,k;
	ts_double mat[3][3];
	ts_double vec[3];

	for(i = 0; i < 3; i++)
        	for(j = 0; j < 3; j++)
            		mat[i][j]=0;
 
	for(k=0;k<vesicle->vlist->n;k++){
		vec[0]=vesicle->vlist->vtx[k]->x;
		vec[1]=vesicle->vlist->vtx[k]->y;
		vec[2]=vesicle->vlist->vtx[k]->z;
		for(i = 0; i < 3; i++)
        		for(j = 0; j <= i; j++)
				mat[i][j]+=vec[i]*vec[j];
	}

// Normalize gyration tensor:
	for(i = 0; i < 3; i++)
        	for(j = 0; j <= i; j++)
			mat[i][j]=mat[i][j]/(ts_double)vesicle->vlist->n;


// diagonal elements are copied twice!	
	for(i = 0; i < 3; i++)
        	for(j = 0; j <= i; j++){
            		gsl_matrix_set(matrix,i,j,mat[i][j]);
            		gsl_matrix_set(matrix,j,i,mat[i][j]);
	}
/*
for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
            printf("%g ", gsl_matrix_get(matrix, i, j));
        printf("\n");
    }
    printf("\n");
*/
	gsl_eigen_symm(matrix, eig, workspace);
	gsl_sort_vector(eig);
/* printf("** eigenvalues    \n");
    for(i = 0; i < 3; i++){
        printf("%3d %25.17e\n", i, gsl_vector_get(eig, i));
} */
	*l1=gsl_vector_get(eig,0);
	*l2=gsl_vector_get(eig,1);
	*l3=gsl_vector_get(eig,2);

	gsl_eigen_symm_free(workspace);
	gsl_matrix_free(matrix);
	gsl_vector_free(eig);
}

ts_ulong get_epoch(){
	time_t currTime;
	currTime = time(NULL);
	return (ts_ulong)currTime;
}

ts_bool get_area_volume(ts_vesicle *vesicle, ts_double *area, ts_double *volume){
	ts_uint i;
	*volume=0.0;
	*area=0.0;
	for(i=0;i<vesicle->tlist->n;i++){
		*volume+=vesicle->tlist->tria[i]->volume;
		*area+=vesicle->tlist->tria[i]->area;
	}
	return TS_SUCCESS;
} 

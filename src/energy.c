/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include "general.h"
#include "energy.h"
#include "vertex.h"
#include<math.h>
#include<stdio.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>


int cmpfunc(const void *x, const void *y)
{
	double diff=	fabs(*(double*)x) - fabs(*(double*)y);
	if(diff<0) return 1;
	else return -1;
}

/** @brief Wrapper that calculates energy of every vertex in vesicle
 *  
 *  Function calculated energy of every vertex in vesicle. It can be used in
 *  initialization procedure or in recalculation of the energy after non-MCsweep *  operations. However, when random move of vertex or flip of random bond occur *  call to this function is not necessary nor recommended. 
 *  @param *vesicle is a pointer to vesicle.
 *  @returns TS_SUCCESS on success.
*/
ts_bool mean_curvature_and_energy(ts_vesicle *vesicle){

    ts_uint i;
    
    ts_vertex_list *vlist=vesicle->vlist;
    ts_vertex **vtx=vlist->vtx;

    for(i=0;i<vlist->n;i++){
        if (!(vtx[i]->type==is_ghost_vtx)) {       
            energy_vertex(vesicle, vtx[i]);
        }
        
    }

    return TS_SUCCESS;
}

/** @brief Calculate energy of a bond (in models where energy is bond related)
 *
 *  This function is experimental and currently only used in polymeres calculation (PEGs or polymeres inside the vesicle).
 *
 *  @param *bond is a pointer to a bond between two vertices in polymere
 *  @param *poly is a pointer to polymere in which we calculate te energy of the bond
 *  @returns TS_SUCCESS on successful calculation
*/
inline ts_bool bond_energy(ts_bond *bond,ts_poly *poly){
//TODO: This value to be changed and implemented in data structure:
	ts_double d_relaxed=1.0;
	bond->energy=poly->k*pow(bond->bond_length-d_relaxed,2);
	return TS_SUCCESS;
};


inline ts_bool curvature_tensor_energy_vertex(ts_vesicle *vesicle, ts_vertex *vtx){
    // direct copy from Samo git repository
    // ...almost
    ts_uint jj, i, j;
    ts_double edge_vector_x[7]={0,0,0,0,0,0,0};
    ts_double edge_vector_y[7]={0,0,0,0,0,0,0};
    ts_double edge_vector_z[7]={0,0,0,0,0,0,0};
    ts_double edge_normal_x[7]={0,0,0,0,0,0,0};
    ts_double edge_normal_y[7]={0,0,0,0,0,0,0};
    ts_double edge_normal_z[7]={0,0,0,0,0,0,0};
    ts_double edge_binormal_x[7]={0,0,0,0,0,0,0};
    ts_double edge_binormal_y[7]={0,0,0,0,0,0,0};
    ts_double edge_binormal_z[7]={0,0,0,0,0,0,0};
    ts_double vertex_normal_x=0.0;
    ts_double vertex_normal_y=0.0;
    ts_double vertex_normal_z=0.0;
//    ts_triangle *triedge[2]={NULL,NULL};

    ts_uint nei,neip,neim;
    ts_vertex *it, *k, *kp,*km;
    ts_triangle *lm=NULL, *lp=NULL;
    ts_double sumnorm;
    ts_double temp_length;
    ts_double cross_x, cross_y, cross_z;

    ts_double Se11, Se21, Se22, Se31, Se32, Se33;
    ts_double Pv11, Pv21, Pv22, Pv31, Pv32, Pv33;
    ts_double Se12, Se13, Se23, Pv12, Pv13, Pv23;//test alias
    ts_double We;
    ts_double Av, We_Av;

	ts_double eigenval[3];

	gsl_matrix *gsl_Sv=gsl_matrix_alloc(3,3);
	gsl_vector *Sv_eigen=gsl_vector_alloc(3);
	gsl_eigen_symm_workspace *workspace=gsl_eigen_symm_alloc(3);

	ts_double mprod[7], phi[7], he[7];
	ts_double Sv[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    // Here edge vector is calculated
//    fprintf(stderr, "Vertex has neighbours=%d\n", vtx->neigh_no);




	Av=0;
	for(i=0; i<vtx->tristar_no; i++){
		vertex_normal_x=(vertex_normal_x - vtx->tristar[i]->xnorm*vtx->tristar[i]->area);
		vertex_normal_y=(vertex_normal_y - vtx->tristar[i]->ynorm*vtx->tristar[i]->area);
		vertex_normal_z=(vertex_normal_z - vtx->tristar[i]->znorm*vtx->tristar[i]->area);
		Av+=vtx->tristar[i]->area/3;
	}
	temp_length=sqrt(pow(vertex_normal_x,2)+pow(vertex_normal_y,2)+pow(vertex_normal_z,2));
	vertex_normal_x=vertex_normal_x/temp_length;
	vertex_normal_y=vertex_normal_y/temp_length;
	vertex_normal_z=vertex_normal_z/temp_length;

	Pv11=1-vertex_normal_x*vertex_normal_x;
	Pv22=1-vertex_normal_y*vertex_normal_y;
	Pv33=1-vertex_normal_z*vertex_normal_z;
	Pv21=vertex_normal_x*vertex_normal_y;
	Pv31=vertex_normal_x*vertex_normal_z;
	Pv32=vertex_normal_y*vertex_normal_z;
    Pv12=Pv21; Pv13=Pv31; Pv23=Pv32; //test

/*	if(vtx->idx==0){
		printf("Vertex normal for vertex %d: %f, %f, %f\n",vtx->idx,vertex_normal_x, vertex_normal_y, vertex_normal_z);
	}
*/
    for(jj=0;jj<vtx->neigh_no;jj++){
	edge_vector_x[jj]=vtx->neigh[jj]->x-vtx->x;
	edge_vector_y[jj]=vtx->neigh[jj]->y-vtx->y;
	edge_vector_z[jj]=vtx->neigh[jj]->z-vtx->z;

	//Here we calculate normalized edge vector

	temp_length=sqrt(edge_vector_x[jj]*edge_vector_x[jj]+edge_vector_y[jj]*edge_vector_y[jj]+edge_vector_z[jj]*edge_vector_z[jj]);
	edge_vector_x[jj]=edge_vector_x[jj]/temp_length;
	edge_vector_y[jj]=edge_vector_y[jj]/temp_length;
	edge_vector_z[jj]=edge_vector_z[jj]/temp_length;

	//end normalization
//	printf("(%f %f %f)\n", vertex_normal_x, vertex_normal_y, vertex_normal_z);
/*
	if(vtx->idx==0){
		printf("Edge vector for vertex %d (vector %d): %f, %f, %f\n",vtx->idx,jj,edge_vector_x[jj], edge_vector_y[jj], edge_vector_z[jj]);
	}
*/
	it=vtx;
	k=vtx->neigh[jj];
	nei=0;
    	for(i=0;i<it->neigh_no;i++){ // Finds the nn of it, that is k 
        	if(it->neigh[i]==k){
            	nei=i;
            	break;
        	}
    	}
    	neip=nei+1;  // I don't like it.. Smells like I must have it in correct order
    	neim=nei-1;
    	if(neip>=it->neigh_no) neip=0;
    	if((ts_int)neim<0) neim=it->neigh_no-1; /* casting is essential... If not
there the neim is never <0 !!! */
  //  fprintf(stderr,"The numbers are: %u %u\n",neip, neim);
    	km=it->neigh[neim];  // We located km and kp
    	kp=it->neigh[neip];

    	if(km==NULL || kp==NULL){
        	fatal("energy_vertex: cannot determine km and kp!",233);
    	}

   for(i=0;i<it->tristar_no;i++){
        for(j=0;j<k->tristar_no;j++){
            if((it->tristar[i] == k->tristar[j])){ //ce gre za skupen trikotnik
                if((it->tristar[i]->vertex[0] == km || it->tristar[i]->vertex[1]
== km || it->tristar[i]->vertex[2]== km )){
                lm=it->tristar[i];
         //       lmidx=i;
                }
                else
                {
                lp=it->tristar[i];
         //       lpidx=i;
                }

            }
        }
    }
if(lm==NULL || lp==NULL) fatal("energy_vertex: Cannot find triangles lm and lp!",233);

	//Triangle normals are NORMALIZED!

	sumnorm=sqrt( pow((lm->xnorm + lp->xnorm),2) + pow((lm->ynorm + lp->ynorm), 2) + pow((lm->znorm + lp->znorm), 2));

	edge_normal_x[jj]=(lm->xnorm+ lp->xnorm)/sumnorm;
	edge_normal_y[jj]=(lm->ynorm+ lp->ynorm)/sumnorm;
	edge_normal_z[jj]=(lm->znorm+ lp->znorm)/sumnorm;


	edge_binormal_x[jj]=(edge_normal_y[jj]*edge_vector_z[jj])-(edge_normal_z[jj]*edge_vector_y[jj]);
	edge_binormal_y[jj]=-(edge_normal_x[jj]*edge_vector_z[jj])+(edge_normal_z[jj]*edge_vector_x[jj]);
	edge_binormal_z[jj]=(edge_normal_x[jj]*edge_vector_y[jj])-(edge_normal_y[jj]*edge_vector_x[jj]);


	//mprod[jj]=it->x*(k->y*edge_vector_z[jj]-edge_vector_y[jj]*k->z)-it->y*(k->x*edge_vector_z[jj]-k->z*edge_vector_x[jj])+it->z*(k->x*edge_vector_y[jj]-k->y*edge_vector_x[jj]);
	mprod[jj]=lm->xnorm*(lp->ynorm*edge_vector_z[jj]-lp->znorm*edge_vector_y[jj]) - lm->ynorm*(lp->xnorm*edge_vector_z[jj]-lp->znorm*edge_vector_z[jj])+ lm->znorm*(lp->xnorm*edge_vector_y[jj]-lp->ynorm*edge_vector_x[jj]);

    cross_x = lm->ynorm*lp->znorm - lp->ynorm*lm->znorm;
    cross_y = lm->znorm*lp->xnorm - lp->znorm*lm->xnorm;
    cross_z = lm->xnorm*lp->ynorm - lp->xnorm*lm->ynorm;
	phi[jj]=copysign(atan2(sqrt(cross_x*cross_x+cross_y*cross_y+cross_z*cross_z),lm->xnorm*lp->xnorm+lm->ynorm*lp->ynorm+lm->znorm*lp->znorm-1e-10),-mprod[jj])+M_PI;
/*	if(vtx->idx==0){
		printf("Angle PHI vertex %d (angle %d): %f\n",vtx->idx,jj,phi[jj]);
	}
*/
//	printf("ACOS arg=%e\n", lm->xnorm*lp->xnorm+lm->ynorm*lp->ynorm+lm->znorm*lp->znorm);
	//he was multiplied with 2 before...
//	he[jj]=sqrt( pow((edge_vector_x[jj]),2) + pow((edge_vector_y[jj]), 2) + pow((edge_vector_z[jj]), 2))*cos(phi[jj]/2.0);
	he[jj]=temp_length*cos(phi[jj]/2.0);
//	printf("phi[%d]=%f\n", jj,phi[jj]);

    //lets try this identity: cos( (acos(b)+pi)/2) = -sqrt(1-b/2) 
    if (abs(cos(phi[jj]/2.0) - copysign(sqrt( (1-(lm->xnorm*lp->xnorm+lm->ynorm*lp->ynorm+lm->znorm*lp->znorm)+1e-15) /2), mprod[jj]))>1e-10){
        fprintf(stdout,"not equal: cos phi/2: %f, sqrt(1/2-cos(phi)/2) : %f, diff: %f\n mprod : %f\n",cos(phi[jj]/2),
         sqrt( (1 +1e-15 -( lm->xnorm*lp->xnorm+lm->ynorm*lp->ynorm+lm->znorm*lp->znorm))/2), 
         cos(phi[jj]/2.0)-sqrt( (1-(lm->xnorm*lp->xnorm+lm->ynorm*lp->ynorm+lm->znorm*lp->znorm)) /2), mprod[jj] );
        fatal("ouch\n",100);
    }

    
/*
	if(vtx->idx==0){
		printf("H operator of edge vertex %d (edge %d): %f\n",vtx->idx,jj,he[jj]);
	}
*/
	Se11=-edge_binormal_x[jj]*edge_binormal_x[jj]*he[jj];
	Se21=-edge_binormal_x[jj]*edge_binormal_y[jj]*he[jj];
	Se22=-edge_binormal_y[jj]*edge_binormal_y[jj]*he[jj];
	Se31=-edge_binormal_x[jj]*edge_binormal_z[jj]*he[jj];
	Se32=-edge_binormal_y[jj]*edge_binormal_z[jj]*he[jj];
	Se33=-edge_binormal_z[jj]*edge_binormal_z[jj]*he[jj];
    Se12=Se21; Se13=Se31; Se23=Se32; //test

	We=vertex_normal_x*edge_normal_x[jj]+vertex_normal_y*edge_normal_y[jj]+vertex_normal_z*edge_normal_z[jj];
	We_Av=We/Av;

/*
	Sv[0][0]+=We_Av* ( Pv11*(Pv11*Se11+Pv21*Se21+Pv31*Se31)+Pv21*(Pv11*Se21+Pv21*Se22+Pv31*Se32)+Pv31*(Pv11*Se31+Pv21*Se32+Pv31*Se33) );
	Sv[0][1]+=We_Av* (Pv21*(Pv11*Se11+Pv21*Se21+Pv31*Se31)+Pv22*(Pv11*Se21+Pv21*Se22+Pv31*Se32)+Pv32*(Pv11*Se31+Pv21*Se32+Pv31*Se33));
	Sv[0][2]+=We_Av* (Pv31*(Pv11*Se11+Pv21*Se21+Pv31*Se31)+Pv32*(Pv11*Se21+Pv21*Se22+Pv31*Se32)+Pv33*(Pv11*Se31+Pv21*Se32+Pv31*Se33));
	
	Sv[1][0]+=We_Av* (Pv11*(Pv21*Se11+Pv22*Se21+Pv32*Se31)+Pv21*(Pv21*Se21+Pv22*Se22+Pv32*Se32)+Pv31*(Pv21*Se31+Pv22*Se32+Pv32*Se33));
	Sv[1][1]+=We_Av* (Pv21*(Pv21*Se11+Pv22*Se21+Pv32*Se31)+Pv22*(Pv21*Se21+Pv22*Se22+Pv32*Se32)+Pv32*(Pv21*Se31+Pv22*Se32+Pv32*Se33));
	Sv[1][2]+=We_Av* (Pv31*(Pv21*Se11+Pv22*Se21+Pv32*Se31)+Pv32*(Pv21*Se21+Pv22*Se22+Pv32*Se32)+Pv33*(Pv21*Se31+Pv22*Se32+Pv32*Se33));

	Sv[2][0]+=We_Av* (Pv11*(Pv31*Se11+Pv32*Se21+Pv33*Se31)+Pv21*(Pv31*Se21+Pv32*Se22+Pv33*Se32)+Pv31*(Pv31*Se31+Pv32*Se32+Pv33*Se33));
	Sv[2][1]+=We_Av* (Pv21*(Pv31*Se11+Pv32*Se21+Pv33*Se31)+Pv22*(Pv31*Se21+Pv32*Se22+Pv33*Se32)+Pv32*(Pv31*Se31+Pv32*Se32+Pv33*Se33));
	Sv[2][2]+=We_Av* (Pv31*(Pv31*Se11+Pv32*Se21+Pv33*Se31)+Pv32*(Pv31*Se21+Pv32*Se22+Pv33*Se32)+Pv33*(Pv31*Se31+Pv32*Se32+Pv33*Se33));
    */
    Sv[0][0]+=We_Av* (Pv11*(Se11*Pv11+Se12*Pv21+Se13*Pv31)+Pv12*(Se21*Pv11+Se22*Pv21+Se23*Pv31)+Pv13*(Se31*Pv11+Se32*Pv21+Se33*Pv31));
	Sv[0][1]+=We_Av* (Pv11*(Se11*Pv12+Se12*Pv22+Se13*Pv32)+Pv12*(Se21*Pv12+Se22*Pv22+Se23*Pv32)+Pv13*(Se31*Pv12+Se32*Pv22+Se33*Pv32));
	Sv[0][2]+=We_Av* (Pv11*(Se11*Pv13+Se12*Pv23+Se13*Pv33)+Pv12*(Se21*Pv13+Se22*Pv23+Se23*Pv33)+Pv13*(Se31*Pv13+Se32*Pv23+Se33*Pv33));
	
	Sv[1][0]+=We_Av* (Pv21*(Se11*Pv11+Se12*Pv21+Se13*Pv31)+Pv22*(Se21*Pv11+Se22*Pv21+Se23*Pv31)+Pv23*(Se31*Pv11+Se32*Pv21+Se33*Pv31));
	Sv[1][1]+=We_Av* (Pv21*(Se11*Pv12+Se12*Pv22+Se13*Pv32)+Pv22*(Se21*Pv12+Se22*Pv22+Se23*Pv32)+Pv23*(Se31*Pv12+Se32*Pv22+Se33*Pv32));
	Sv[1][2]+=We_Av* (Pv21*(Se11*Pv13+Se12*Pv23+Se13*Pv33)+Pv22*(Se21*Pv13+Se22*Pv23+Se23*Pv33)+Pv23*(Se31*Pv13+Se32*Pv23+Se33*Pv33));

	Sv[2][0]+=We_Av* (Pv31*(Se11*Pv11+Se12*Pv21+Se13*Pv31)+Pv32*(Se21*Pv11+Se22*Pv21+Se23*Pv31)+Pv33*(Se31*Pv11+Se32*Pv21+Se33*Pv31));
	Sv[2][1]+=We_Av* (Pv31*(Se11*Pv12+Se12*Pv22+Se13*Pv32)+Pv32*(Se21*Pv12+Se22*Pv22+Se23*Pv32)+Pv33*(Se31*Pv12+Se32*Pv22+Se33*Pv32));
	Sv[2][2]+=We_Av* (Pv31*(Se11*Pv13+Se12*Pv23+Se13*Pv33)+Pv32*(Se21*Pv13+Se22*Pv23+Se23*Pv33)+Pv33*(Se31*Pv13+Se32*Pv23+Se33*Pv33));
//	printf("(%f %f %f); (%f %f %f); (%f %f %f)\n", edge_vector_x[jj], edge_vector_y[jj], edge_vector_z[jj], edge_normal_x[jj], edge_normal_y[jj], edge_normal_z[jj], edge_binormal_x[jj], edge_binormal_y[jj], edge_binormal_z[jj]);

    } // END FOR JJ
	gsl_matrix_set(gsl_Sv, 0,0, Sv[0][0]);
	gsl_matrix_set(gsl_Sv, 0,1, Sv[0][1]);
	gsl_matrix_set(gsl_Sv, 0,2, Sv[0][2]);
	gsl_matrix_set(gsl_Sv, 1,0, Sv[1][0]);
	gsl_matrix_set(gsl_Sv, 1,1, Sv[1][1]);
	gsl_matrix_set(gsl_Sv, 1,2, Sv[1][2]);
	gsl_matrix_set(gsl_Sv, 2,0, Sv[2][0]);
	gsl_matrix_set(gsl_Sv, 2,1, Sv[2][1]);
	gsl_matrix_set(gsl_Sv, 2,2, Sv[2][2]);

//	printf("Se= %f, %f, %f\n    %f, %f, %f\n    %f, %f, %f\n", Se11, Se21, Se31, Se21, Se22, Se32, Se31, Se32, Se33);
//	printf("Pv= %f, %f, %f\n    %f, %f, %f\n    %f, %f, %f\n", Pv11, Pv21, Pv31, Pv21, Pv22, Pv32, Pv31, Pv32, Pv33);
//	printf("Sv= %f, %f, %f\n    %f, %f, %f\n    %f, %f, %f\n", Sv[0][0], Sv[0][1], Sv[0][2], Sv[1][0], Sv[1][1], Sv[1][2], Sv[2][0], Sv[2][1], Sv[2][2]);


	gsl_eigen_symm(gsl_Sv, Sv_eigen, workspace);

//	printf("Eigenvalues: %f, %f, %f\n", gsl_vector_get(Sv_eigen, 0),gsl_vector_get(Sv_eigen, 1), gsl_vector_get(Sv_eigen, 2) );
//	printf("Eigenvalues: %f, %f, %f\n", gsl_matrix_get(evec, 0,0),gsl_matrix_get(evec, 0,1), gsl_matrix_get(evec, 0,2) );


	eigenval[0]= gsl_vector_get(Sv_eigen, 0);
	eigenval[1]= gsl_vector_get(Sv_eigen, 1);
	eigenval[2]= gsl_vector_get(Sv_eigen, 2);

	qsort(eigenval, 3, sizeof(ts_double), cmpfunc);
	if(vtx->idx==0){
	//printf("Eigenvalues: %f, %f, %f\n", eigenval[0], eigenval[1], eigenval[2] );
//	exit(0);
	}

    // Yoav : and the stuff I'm tracking
    vtx->nx = vertex_normal_x;
    vtx->ny = vertex_normal_y;
    vtx->nz = vertex_normal_z;
    //vtx->curvature = (eigenval[0] + eigenval[1])/2;
    //vtx->curvature2 = eigenval[0]*eigenval[1];
    vtx->curvature = (eigenval[0]+ eigenval[1]);
    vtx->curvature2 = eigenval[0]*eigenval[1];
	vtx->energy=4*vtx->xk*(pow(eigenval[0]+eigenval[1]-2*vtx->c,2))*Av;
    if (vtx->type&is_anisotropic_vtx && vtx->xk2!=0){
            vtx->energy += vtx->xk2 * Av * eigenval[0]*eigenval[1];
    }

	gsl_matrix_free(gsl_Sv);
	gsl_vector_free(Sv_eigen);
//	gsl_matrix_free(evec);
	gsl_eigen_symm_free(workspace);
	return TS_SUCCESS;
}


/** @brief Calculation of the bending energy of the vertex.
 *  
 *  Main function that calculates energy of the vertex \f$i\f$. Function returns \f$\frac{1}{2}(c_1+c_2-c)^2 s\f$, where \f$(c_1+c_2)/2\f$ is mean curvature,
 * \f$c/2\f$ is spontaneous curvature and \f$s\f$ is area per vertex \f$i\f$.
 *
 * Nearest neighbors (NN) must be ordered in counterclockwise direction for this function to work.
 *  Firstly NNs that form two neighboring triangles are found (\f$j_m\f$, \f$j_p\f$ and common \f$j\f$). Later, the scalar product of vectors \f$x_1=(\mathbf{i}-\mathbf{j_p})\cdot (\mathbf{i}-\mathbf{j_p})(\mathbf{i}-\mathbf{j_p})\f$, \f$x_2=(\mathbf{j}-\mathbf{j_p})\cdot  (\mathbf{j}-\mathbf{j_p})\f$  and \f$x_3=(\mathbf{j}-\mathbf{j_p})\cdot (\mathbf{i}-\mathbf{j_p})\f$  are calculated. From these three vectors the \f$c_{tp}=\frac{1}{\tan(\varphi_p)}\f$ is calculated, where \f$\varphi_p\f$ is the inner angle at vertex \f$j_p\f$. The procedure is repeated for \f$j_m\f$ instead of \f$j_p\f$ resulting in \f$c_{tn}\f$.
 *  
\begin{tikzpicture}{
\coordinate[label=below:$i$] (i) at (2,0);
\coordinate[label=left:$j_m$] (jm) at (0,3.7);
\coordinate[label=above:$j$] (j) at (2.5,6.4);
\coordinate[label=right:$j_p$] (jp) at (4,2.7);

\draw (i) -- (jm) -- (j) -- (jp) -- (i) -- (j);

\begin{scope}
\path[clip] (jm)--(i)--(j);
\draw (jm) circle (0.8);
\node[right] at (jm) {$\varphi_m$};
\end{scope}

\begin{scope}
\path[clip] (jp)--(i)--(j);
\draw (jp) circle (0.8);
\node[left] at (jp) {$\varphi_p$};
\end{scope}

%%vertices
\draw [fill=gray] (i) circle (0.1);
\draw [fill=white] (j) circle (0.1);
\draw [fill=white] (jp) circle (0.1);
\draw [fill=white] (jm) circle (0.1);
%\node[draw,circle,fill=white] at (i) {};
\end{tikzpicture}

 * The curvature is then calculated as \f$\mathbf{h}=\frac{1}{2}\Sigma_{k=0}^{\mathrm{neigh\_no}} c_{tp}^{(k)}+c_{tm}^{(k)} (\mathbf{j_k}-\mathbf{i})\f$, where \f$c_{tp}^{(k)}+c_{tm}^k=2\sigma^{(k)}\f$ (length in dual lattice?) and the previous equation can be written as \f$\mathbf{h}=\Sigma_{k=0}^{\mathrm{neigh\_no}}\sigma^{(k)}\cdot(\mathbf{j}-\mathbf{i})\f$ (See Kroll, p. 384, eq 70).
 *
 * From the curvature the enery is calculated by equation \f$E=\frac{1}{2}\mathbf{h}\cdot\mathbf{h}\f$.
 * @param *vtx is a pointer to vertex at which we want to calculate the energy
 * @returns TS_SUCCESS on successful calculation.
*/
inline ts_bool energy_vertex(ts_vesicle *vesicle, ts_vertex *vtx){
    if (vesicle->tape->type_of_curvature_model==10) {
        return curvature_tensor_energy_vertex(vesicle, vtx);
    }
    
    ts_uint jj;
    ts_uint jjp,jjm;
    ts_vertex *j,*jp, *jm;
    ts_triangle *jt;
    ts_double s=0.0,xh=0.0,yh=0.0,zh=0.0,txn=0.0,tyn=0.0,tzn=0.0;
    ts_double x1,x2,x3,ctp,ctm,tot,xlen;
    ts_double h,ht,norml;
    ts_double angle_sum=0;
    ts_double a_dot_b, a_cross_b_x, a_cross_b_y, a_cross_b_z, mag_a_cross_b;
    ts_bool model=vesicle->tape->type_of_curvature_model;
    
    for(jj=1; jj<=vtx->neigh_no;jj++){
        jjp=jj+1;
        if(jjp>vtx->neigh_no) jjp=1;
        jjm=jj-1;
        if(jjm<1) jjm=vtx->neigh_no;
        j=vtx->neigh[jj-1];
        jp=vtx->neigh[jjp-1];
        jm=vtx->neigh[jjm-1];
        jt=vtx->tristar[jj-1];
        x1=vtx_distance_sq(vtx,jp); //shouldn't be zero!
        x2=vtx_distance_sq(j,jp); // shouldn't be zero!
        //x1=pow(vtx->x-jp->x,2)+pow(vtx->y-jp->y,2)+pow(vtx->z-jp->z,2);
        //x2=pow(j->x-jp->x,2)+pow(j->y-jp->y,2)+pow(j->z-jp->z,2);
        x3=(j->x-jp->x)*(vtx->x-jp->x)+
           (j->y-jp->y)*(vtx->y-jp->y)+
           (j->z-jp->z)*(vtx->z-jp->z);
        
#ifdef TS_DOUBLE_DOUBLE
        ctp=x3/sqrt(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_FLOAT
        ctp=x3/sqrtf(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        ctp=x3/sqrtl(x1*x2-x3*x3);
#endif
        x1=vtx_distance_sq(vtx,jm);
        x2=vtx_distance_sq(j,jm);
        //x1=pow(vtx->x-jm->x,2)+pow(vtx->y-jm->y,2)+pow(vtx->z-jm->z,2);
        //x2=pow(j->x-jm->x,2)+pow(j->y-jm->y,2)+pow(j->z-jm->z,2);
        x3=(j->x-jm->x)*(vtx->x-jm->x)+
           (j->y-jm->y)*(vtx->y-jm->y)+
           (j->z-jm->z)*(vtx->z-jm->z);
#ifdef TS_DOUBLE_DOUBLE
        ctm=x3/sqrt(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_FLOAT
        ctm=x3/sqrtf(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        ctm=x3/sqrtl(x1*x2-x3*x3);
#endif
        tot=ctp+ctm;
        tot=0.5*tot;

        //testing
        xlen=vtx_distance_sq(j,vtx);
        //xlen=pow(vtx->x-j->x,2)+pow(vtx->y-j->y,2)+pow(vtx->z-j->z,2);
/*
#ifdef  TS_DOUBLE_DOUBLE 
        vtx->bond[jj-1]->bond_length=sqrt(xlen); 
#endif
#ifdef  TS_DOUBLE_FLOAT
        vtx->bond[jj-1]->bond_length=sqrtf(xlen); 
#endif
#ifdef  TS_DOUBLE_LONGDOUBLE 
        vtx->bond[jj-1]->bond_length=sqrtl(xlen); 
#endif

        vtx->bond[jj-1]->bond_length_dual=tot*vtx->bond[jj-1]->bond_length;
*/
        s+=tot*xlen;
        xh+=tot*(j->x - vtx->x);
        yh+=tot*(j->y - vtx->y);
        zh+=tot*(j->z - vtx->z);
        txn+=jt->xnorm;
        tyn+=jt->ynorm;
        tzn+=jt->znorm;

        if (model!=0){
            if (model==2 || (model==1 && vtx->type & is_anisotropic_vtx && vtx->xk2!=0)){
                // angle stuff
                //angle_sum += atan(ctp) + atan(ctm); // simple but slow!
                /// get the angle m-vtx-j
                // atan2(|axb|,a*b) was recommended at mathwork forum (cosin has small angle problems, and still need a sqrt)
                a_dot_b = (jm->x-vtx->x)*(j->x-vtx->x)+
                          (jm->y-vtx->y)*(j->y-vtx->y)+
                          (jm->z-vtx->z)*(j->z-vtx->z);
                a_cross_b_x = (jm->y-vtx->y)*(j->z-vtx->z)-(jm->z-vtx->z)*(j->y-vtx->y);
                a_cross_b_y = (jm->z-vtx->z)*(j->x-vtx->x)-(jm->x-vtx->x)*(j->z-vtx->z);
                a_cross_b_z = (jm->x-vtx->x)*(j->y-vtx->y)-(jm->y-vtx->y)*(j->x-vtx->x);
                mag_a_cross_b = sqrt(pow(a_cross_b_x,2)+pow(a_cross_b_y,2)+pow(a_cross_b_z,2));
                angle_sum += atan2(mag_a_cross_b, a_dot_b);
            }
        }
    }
    
    h=xh*xh+yh*yh+zh*zh;
    ht=txn*xh+tyn*yh + tzn*zh;
    s=s/4.0; 
#ifdef TS_DOUBLE_DOUBLE
    if(ht>=0.0) {
        vtx->curvature=sqrt(h);
    } else {
        vtx->curvature=-sqrt(h);
    }
#endif
#ifdef TS_DOUBLE_FLOAT
    if(ht>=0.0) {
        vtx->curvature=sqrtf(h);
    } else {
        vtx->curvature=-sqrtf(h);
    }
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    if(ht>=0.0) {
        vtx->curvature=sqrtl(h);
    } else {
        vtx->curvature=-sqrtl(h);
    }
#endif
    //also great point to update normal: note that the triangle normal is inwards
    norml=sqrt(txn*txn+tyn*tyn+tzn*tzn);
    vtx->nx=-txn/norml;
	vtx->ny=-tyn/norml;
	vtx->nz=-tzn/norml;
// c is forced curvature energy for each vertex. Should be set to zero for
// normal circumstances.
/* the following statement is an expression for $\frac{1}{2}\int(c_1+c_2-c_0^\prime)^2\mathrm{d}A$, where $c_0^\prime=2c_0$ (twice the spontaneous curvature)  */
    vtx->energy=vtx->xk* 0.5*s*(vtx->curvature/s-vtx->c)*(vtx->curvature/s-vtx->c);
    
    if (model!=0){
        if (model==2 || ( model==1 && vtx->type&is_anisotropic_vtx && vtx->xk2!=0)) {
            vtx->curvature2 = (2*M_PI- angle_sum)/s;
            if (vtx->type&is_anisotropic_vtx && vtx->xk2!=0){
                vtx->energy += vtx->xk2 * s * vtx->curvature2;
            }
        }
    }
    

    return TS_SUCCESS;
}



ts_bool sweep_attraction_bond_energy(ts_vesicle *vesicle){
	int i;
	for(i=0;i<vesicle->blist->n;i++){
		attraction_bond_energy(vesicle, vesicle->blist->bond[i]);
	}
	return TS_SUCCESS;
}


inline ts_bool attraction_bond_energy(ts_vesicle *vesicle, ts_bond *bond){

    if (vesicle->tape->type_of_bond_model==0){ // bonding type bond together
	    if((bond->vtx1->type&is_bonding_vtx && bond->vtx2->type&is_bonding_vtx)){
            // f(w1,w2)
		    bond->energy=-0.5*(bond->vtx1->w+bond->vtx2->w);
	    }
	    else {
		    bond->energy=0.0;
	    }
    }
    if (vesicle->tape->type_of_bond_model==1){ // bond by same type
	    if((bond->vtx1->type&is_bonding_vtx && bond->vtx2->type==bond->vtx1->type)){
            // f(w1,w2)
		    bond->energy=-0.5*(bond->vtx1->w+bond->vtx2->w);
	    }
	    else {
		    bond->energy=0.0;
	    }
    }
	return TS_SUCCESS;
}


ts_double direct_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
	// modified to include Vicsek Interaction and basic inhibition models
    // 0: force in normal direction
    // 1: inhibation F*=No_inactive/No_neigh
    // 2: inhibation F*=No_c>=0/No_neigh
    // 3: inhibition, no force if neighbor c0<0
    //16: Vicsek model, constant weight
    //17: Vicsek model, ~ 1/R

    
    // quit if there is no force
    if(fabs(vtx->f)<1e-15) return 0.0;

    ts_bool model=vesicle->tape->type_of_force_model;
    ts_double vicsek_strength=vesicle->tape->vicsek_strength;
    ts_double vicsek_radius= vesicle->tape->vicsek_radius;

    ts_double norml,ddp=0.0;
	//ts_double xnorm=0.0,ynorm=0.0,znorm=0.0;
	ts_double vixnorm=0.0,viynorm=0.0,viznorm=0.0;

    //stupid loop variables don't go in stupid for loop due to stupid C89 compiler
    ts_uint i, j, curr_dist;
    
    // estimated mallocation size for the cluster: roughly ~pi*r^2
    ts_uint max_vtx_seen=3*(( (int) vicsek_radius)+1)*(( (int) vicsek_radius)+1);
    // allocate the struct for the "seen vertex" (defined in general.h, functions in vertex.c)
    ts_seen_vertex *seen_vtx; //initalize in vicsek

    ts_uint No_neigh_activating; // number of activating neighbors
    ts_double inhibition_factor;



    // if vicsek type and vicsek model is relevant
    if ( vtx->type&is_vicsek_vtx && (model==16 || model==17) && fabs(vicsek_strength)>1e-15 && fabs(vicsek_radius)>1e-15) {
        
        //vicsek model
        //force directed by Vicsek sum-over-neighbors-normals

	
        //prime vertex normal
	    vixnorm=vtx->nx;
	    viynorm=vtx->ny;
	    viznorm=vtx->nz;

        //initialize seen_vtx
        seen_vtx = init_seen_vertex(max_vtx_seen);
        // we have now seen the prime vertex
        add_vtx_to_seen(seen_vtx, vtx);
        
    
        // Breadth first search using seen_vtx, layer by layer,
        // until reaching layer that is outside the maximum radius
        for (curr_dist=1 ; curr_dist<=vicsek_radius; curr_dist++){
            advance_seen_vertex_to_next_layer(seen_vtx);
            
            // The for loops are split by layers
            // The new "next" layer is being built from
            // the neighbors of the completed "current" layer
            // some of the checked neighbors may be
            // from the layer before that: the "previous" layer
            // ..............
            //  A1 --  A2 --  O1     O-not-yet added              seen_vtx[]:
            // /  \  /  |  \  / \    N-bext layer, not added yet          [...
            // C1--C2--N3--O1 - O1   A-next layer, already added            Q1
            // | / | \ | \ |  \ |    C-current, neighbor iterated layer     Q2
            // P1--P2--C3--N4 - O1   P-previous layer                       P1  <-previous
            // | \ | \ | \ |  \ |    Q-previous previous layer              P2
            // Q1--Q2--P3--C4 - N5                                          P3
            // ..............                                               C1  <-current
            //                                                              C2
            //                                                              C3
            //                                                              C4
            //                                                              A1  <-next
            //                                                              A2
            //                                                              _   <-n_top

            for (i=seen_vtx->n_curr; i<seen_vtx->n_next; i++) {
                // loop over the current layer: seen_vtx[i]
                

                for (j=0; j<seen_vtx->vtx[i]->neigh_no; j++) {
                    //loop over their neighbors: seen_vtx[i]->neigh[j]

                    // is this neighbor is not vicsek, skip to next neighbor
                    if (!(seen_vtx->vtx[i]->neigh[j]->type&is_vicsek_vtx)) continue;
                    //else{ rest of the j loop }

                    // has this neighbor been seen?
                    if (is_in_seen_vertex(seen_vtx,seen_vtx->vtx[i]->neigh[j]) && (curr_dist!=1)){
                        // if seen, skip to next neighbor in the j loop
                        // 1st layer is always new, so we can skip the check
                        continue;
                    }
                    else {
                        //new vertex! add to the next layer

                        add_vtx_to_seen(seen_vtx,seen_vtx->vtx[i]->neigh[j]);

                        //add normal to the sum
                        // Vicsek model 17: weight by 1/distance
                        if (vesicle->tape->type_of_force_model == 17) {
                            vixnorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nx / curr_dist;
                            viynorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->ny / curr_dist;
                            viznorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nz / curr_dist;
                        }
                        else {
                            vixnorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nx;
                            viynorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->ny;
                            viznorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nz;
                        }

                    }
                }
            }
            // finished this layer
            //if did not add any vertices: completed cluster
            if (seen_vtx->n_next==seen_vtx->n_top) break;
        }

        //having finished summing the normals, normalize the resulting vector
        norml=sqrt(vixnorm*vixnorm+viynorm*viynorm+viznorm*viznorm);
        vixnorm/=norml;
	    viynorm/=norml;
	    viznorm/=norml;

	    /*calculate ddp, Viscek force directed displacement*/
	    vtx->fx=vixnorm;
        vtx->fy=viynorm;
        vtx->fz=viznorm;

	    //don't forget to free! 
        seen_vertex_free(seen_vtx);

    //end if (!Vicsek)
    }
    else {
        //regular "force in normal direction"
        vtx->fx = vtx->nx;
        vtx->fy = vtx->ny;
        vtx->fz = vtx->nz;
    }
    
    // now we calculate the final force and work

    
    //we recalculate the force based on inhibition: we want the real force to be recorded!

    // HIV_Gag inhibition model
    if (model==1){
        // force ~ number of active neighbors
        No_neigh_activating = 0;
        for (i=0; i<vtx->neigh_no; i++){ 
            if ( !( vtx->neigh[i]->type & is_active_vtx ) ) {
                No_neigh_activating++;
            }
        }
        inhibition_factor = 1;
        inhibition_factor *= No_neigh_activating;  //take that, integer-integer division!
        inhibition_factor /= vtx->neigh_no;
        // change force by factor
        vtx->fx *= inhibition_factor;
        vtx->fy *= inhibition_factor;
        vtx->fz *= inhibition_factor;
    }

    if (model==2){
        // force ~ number of c>0 neighbors
        No_neigh_activating = 0;
        for (i=0; i<vtx->neigh_no; i++){
            if ( vtx->neigh[i]->c > -1e-15   ) {
                No_neigh_activating++;
            }
        }
        inhibition_factor = 1;
        inhibition_factor *= No_neigh_activating;  //take that, integer-integer division!
        inhibition_factor /= vtx->neigh_no;
        // change force by factor
        vtx->fx *= inhibition_factor;
        vtx->fy *= inhibition_factor;
        vtx->fz *= inhibition_factor;
    }

    if (model==3){
        // force disabled for any c<0 neighbors
        for (i=0; i<vtx->neigh_no; i++){
            if ( vtx->neigh[i]->c < -1e-15   ) {
                vtx->fx = 0;
                vtx->fy = 0;
                vtx->fz = 0;
                break;
            }
        }
    }

    /*calculate ddp, normal force directed displacement*/
	ddp=vtx->fx*(vtx->x-vtx_old->x)+vtx->fy*(vtx->y-vtx_old->y)+vtx->fz*(vtx->z-vtx_old->z);

    
    /*calculate dE*/
    return -vtx->f*ddp;
//end function
}

ts_double direct_force_from_Fz_balance(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
    // calculate the Fz force balance
    //     We could separat from direct_force_energy()
    //     since F_direct(old) = F_direct + Fz_balance \hat{z}
    //     W = F_direct * dX + Fz_balance * dz
    // Force is cached once and only once! : not recalculated!
    // Changing this behavior is weird: once, or always, work identically, 
    // once every 2 works identically, once every 64 works identically, 
    // but 65-128 don't!.
    ts_double f_diff=0;
    static ts_double Fz;
    static ts_char countdown=0;

    if (!countdown){
        Fz=total_force_on_vesicle(vesicle);
        countdown=64;//pow(2,20)-1; go back here in 64 steps
    }
    else{
        if (vtx->type&is_active_vtx) f_diff+=vtx->fz;
        if (vtx_old->type&is_active_vtx) f_diff-=vtx_old->fz;
        Fz+=f_diff;
        countdown-=1; //recalculate when reach 0
    }

    if (Fz>0){
	    return Fz*(vtx->z-vtx_old->z)/vesicle->vlist->n;
    }
    else{
        return 0;
    }	
	
}

inline ts_double total_force_on_vesicle(ts_vesicle *vesicle){
	ts_uint i;
	ts_double fz=0;
	/*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
    for (i=0;i<vesicle->vlist->n;i++){
		if(vesicle->vlist->vtx[i]->type&is_active_vtx){
			fz+=vesicle->vlist->vtx[i]->fz;
		}
	}
    return fz;
	
}

ts_double adhesion_energy_diff(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
    ts_double delta_energy=0;
    ts_double z=vtx->z, z_old=vtx_old->z;
    ts_double z0=vesicle->tape->z_adhesion, dz=vesicle->tape->adhesion_cuttoff;
    ts_double c0=vesicle->adhesion_center, r=vesicle->tape->adhesion_radius;
    //1 for step potential
	if(vesicle->tape->type_of_adhesion_model==1){

		if( (vtx->type&is_adhesive_vtx) && (abs(z-z0)<dz) ){
				delta_energy-=vtx->ad_w;
		}
		if( (vtx_old->type&is_adhesive_vtx) &&(abs(z_old-z0)<dz) ){
				delta_energy+=vtx_old->ad_w;
		}
	}

    //2 for parabolic potential
	else if(vesicle->tape->type_of_adhesion_model==2){

        // can't combine them well: each has (theoretically) different adhesion
		if( (vtx->type&is_adhesive_vtx) && ((z-z0)<=dz )){
				delta_energy-=(vtx->ad_w/pow(dz,2))*pow(z - dz,2);
		}
		if( (vtx_old->type&is_adhesive_vtx) && ((z_old-z0)>dz) ){
				delta_energy+=(vtx_old->ad_w/pow(dz,2))*pow(z_old - dz,2);
		}
	}

    //3 for sphrerical adhesion substrate with constant potential
	else if(vesicle->tape->type_of_adhesion_model==3){
		if( (vtx->type&is_adhesive_vtx) && (pow(pow(c0-z,2) + pow(vtx->x,2) + pow(vtx->y,2),0.5) - r < dz)){
			delta_energy-=vtx->ad_w;
		}
		if( (vtx_old->type&is_adhesive_vtx) && (pow(pow(c0-z_old,2) + pow(vtx_old->x,2) + pow(vtx_old->y,2),0.5) - r < dz)){
			delta_energy+=vtx_old->ad_w;
		}
	}
    
    //4 for cylindrical adhesive substrate with constant potential
	else if(vesicle->tape->type_of_adhesion_model==4){
		if( (vtx->type&is_adhesive_vtx) && (pow(pow(c0 -z,2) + pow(vtx->x,2),0.5) - r < dz)){
			delta_energy-=vtx->ad_w;
		}
		if( (vtx_old->type&is_adhesive_vtx) && (pow(pow(c0 - z_old,2) + pow(vtx_old->x,2),0.5) - r < dz)){
			delta_energy+=vtx_old->ad_w;
		}
	}

    return delta_energy;
}

void stretchenergy(ts_vesicle *vesicle, ts_triangle *triangle){
	triangle->energy=vesicle->tape->xkA0/2.0*pow((triangle->area/vesicle->tlist->a0-1.0),2);
}

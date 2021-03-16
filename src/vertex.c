/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include<stdio.h>

ts_bool vertex_list_assign_id(ts_vertex_list *vlist, ts_uint id){
	ts_uint i;	
	for(i=0;i<vlist->n;i++){
		vlist->vtx[i]->id = id;
	}
	return TS_SUCCESS;
}

ts_vertex_list *init_vertex_list(ts_uint N){	
	ts_int i;
    ts_vertex_list *vlist=(ts_vertex_list *)malloc(sizeof(ts_vertex_list));
    
	if(N==0){
		err("Initialized vertex list with zero elements. Pointer set to NULL");
        vlist->n=0;
		vlist->vtx=NULL;
		return vlist;
	}
	
    vlist->vtx=(ts_vertex **)calloc(N,sizeof(ts_vertex *));
    if(vlist->vtx==NULL)
        fatal("Fatal error reserving memory space for vertex list! Could number of requsted vertices be too large?", 100);
    for(i=0;i<N;i++) {
        vlist->vtx[i]=(ts_vertex *)calloc(1,sizeof(ts_vertex));
        vlist->vtx[i]->idx=i;

    /* initialize Ylm for spherical hamonics DONE in sh.c */
/*    for(i=0;i<l;i++){
        vlist->vtx[i]->Ylm[i]=(ts_double **)calloc(2*i+1,sizeof(ts_double *));
        for(j=0;j<(2*i+1);j++){
            clist->vtx[i]->Ylm[i][j]=(ts_double *)calloc(sizeof(ts_double));
        }
    }
*/


    }
    vlist->n=N;
	return vlist;
}

ts_bool vtx_add_neighbour(ts_vertex *vtx, ts_vertex *nvtx){
    ts_uint i;
    /* no neighbour can be null! */
    if(vtx==NULL || nvtx==NULL) return TS_FAIL;
    
    /*if it is already a neighbour don't add it to the list */
    for(i=0; i<vtx->neigh_no;i++){
        if(vtx->neigh[i]==nvtx) return TS_FAIL;
    }
    ts_uint nn=++vtx->neigh_no;
    vtx->neigh=(ts_vertex **)realloc(vtx->neigh, nn*sizeof(ts_vertex *));
    vtx->neigh[nn-1]=nvtx;
/* This was a bug in creating DIPYRAMID (the neighbours were not in right
 * order).
 */
    /* pa se sosedu dodamo vertex */
    /*if it is already a neighbour don't add it to the list */
/*
    for(i=0; i<nvtx->data->neigh_no;i++){
        if(nvtx->data->neigh[i]==vtx) return TS_FAIL;
    } 
    nn=++nvtx->data->neigh_no;
    nvtx->data->neigh=(ts_vertex **)realloc(nvtx->data->neigh, nn*sizeof(ts_vertex *));
    nvtx->data->neigh[nn-1]=vtx;
*/

    return TS_SUCCESS;
}

/* TODO: optimize this. test this. */
ts_bool vtx_remove_neighbour(ts_vertex *vtx, ts_vertex *nvtx){
/* find a neighbour */
/* remove it from the list while shifting remaining neighbours up */
    ts_uint i,j=0;
    for(i=0;i<vtx->neigh_no;i++){
//		fprintf(stderr,"neigh_addr=%ld\n", (long)vtx->neigh[i]);
        if(vtx->neigh[i]!=nvtx){
            vtx->neigh[j]=vtx->neigh[i];
            j++;
        }
    }
//	fprintf(stderr,"remove_neighbour: vtx1_addr=%ld, vtx2_addr=%ld\n",(long)vtx,(long)nvtx);
/* resize memory. potentionally time consuming */
    vtx->neigh_no--;
    vtx->neigh=(ts_vertex **)realloc(vtx->neigh,vtx->neigh_no*sizeof(ts_vertex *));
    if(vtx->neigh == NULL && vtx->neigh_no!=0)
        fatal("(1) Reallocation of memory failed during removal of vertex neighbour in vtx_remove_neighbour",100);
//fprintf(stderr,"first alloc");
/* repeat for the neighbour */
/* find a neighbour */
/* remove it from the list while shifting remaining neighbours up */
	j=0;
    for(i=0;i<nvtx->neigh_no;i++){
        if(nvtx->neigh[i]!=vtx){
            nvtx->neigh[j]=nvtx->neigh[i];
            j++;
        }
    }
/* resize memory. potentionally time consuming. */
//	fprintf(stderr,"Neigbours=%d\n",nvtx->neigh_no);
    nvtx->neigh_no--;
    		nvtx->neigh=(ts_vertex **)realloc(nvtx->neigh,nvtx->neigh_no*sizeof(ts_vertex *));
//	fprintf(stderr,"Neigbours=%d\n",nvtx->neigh_no);
    if(nvtx->neigh == NULL && nvtx->neigh_no!=0)
        fatal("(2) Reallocation of memory failed during removal of vertex neighbour in vtx_remove_neighbour",100);

    return TS_SUCCESS;
}



ts_bool vtx_add_bond(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2){
    ts_bond *bond;
    bond=bond_add(blist,vtx1,vtx2);
    if(bond==NULL) return TS_FAIL;
    vtx1->bond_no++;
    vtx2->bond_no++;
   // vtx2->data->bond_no++;

    vtx1->bond=(ts_bond **)realloc(vtx1->bond, vtx1->bond_no*sizeof(ts_bond *)); 
    vtx2->bond=(ts_bond **)realloc(vtx2->bond, vtx2->bond_no*sizeof(ts_bond *)); 
   // vtx2->data->bond=(ts_bond **)realloc(vtx2->data->bond, vtx2->data->bond_no*sizeof(ts_bond *)); 
    vtx1->bond[vtx1->bond_no-1]=bond;
    vtx2->bond[vtx2->bond_no-1]=bond;
   // vtx2->ata->bond[vtx2->data->bond_no-1]=bond;
    return TS_SUCCESS;
}

ts_bool vtx_add_cneighbour(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex *vtx2){
    ts_bool retval;
    retval=vtx_add_neighbour(vtx1,vtx2);
  //  retval=vtx_add_neighbour(vtx2,vtx1);
    if(retval==TS_SUCCESS)
    retval=vtx_add_bond(blist,vtx1,vtx2); 
    return retval;
}

/*TODO: write and optimize this urgently before use! */
ts_bool vtx_remove_cneighbour(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex
*vtx2){
//    ts_bool retval;
/* remove the bond */
//retval=vtx_remove_bond(blist,vtx1,vtx2);
/* remove the vertices */
    return TS_SUCCESS;
}


ts_bool vtx_free(ts_vertex  *vtx){
    if(vtx->neigh!=NULL)   free(vtx->neigh);
    if(vtx->tristar!=NULL) free(vtx->tristar);
    if(vtx->bond!=NULL)    free(vtx->bond);
    free(vtx);
    return TS_SUCCESS;
}

ts_bool vtx_list_free(ts_vertex_list *vlist){
    int i;
    for(i=0;i<vlist->n;i++){
		if(vlist->vtx[i]!=NULL) vtx_free(vlist->vtx[i]);
    }
    //free(*(vlist->vtx));
    free(vlist->vtx);
    free(vlist);
    return TS_SUCCESS;
}

inline ts_double vtx_distance_sq(ts_vertex *vtx1, ts_vertex *vtx2){
    ts_double dist;
#ifdef TS_DOUBLE_DOUBLE
    dist=pow(vtx1->x-vtx2->x,2) + pow(vtx1->y-vtx2->y,2) + pow(vtx1->z-vtx2->z,2);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    dist=powl(vtx1->x-vtx2->x,2) + powl(vtx1->y-vtx2->y,2) + powl(vtx1->z-vtx2->z,2);
#endif
#ifdef TS_DOUBLE_FLOAT
    dist=powf(vtx1->x-vtx2->x,2) + powf(vtx1->y-vtx2->y,2) + powf(vtx1->z-vtx2->z,2);
#endif
    return(dist);
}



ts_bool vtx_set_global_values(ts_vesicle *vesicle){ 
    ts_double xk=vesicle->bending_rigidity;
    ts_uint i; 
    for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]->xk=xk;
    }
    return TS_SUCCESS;
}

/** Calculates the triple product of vectors defined by vertices vtx1, vtx2 and vtx3, ($\mathrm{vtx}_1\cdot(\mathrm{vtx}_2\cross\mathrm{vtx}_3$):
 *  \begin{vmatrix}
 *  x_1 & y_1 & z_1 \\
 *  x_2-x_1 & y_2-y_1 & z_2-z_1\\
 *  x_3-x_1 & y_3-y_1 & z_3-z_1\\
 *  \end{vmatrix}
 *  where the vertices coordinates are denoted by corresponding vertex index number. Function is used to determine the orientation of area formed by triangle formed by the three given vertices.
 *
 *      @param vtx1 is first vertex, according to which the orientation is calculated
 *      @param vtx2 is the second vertex
 *      @param vtx3 is the third vertex
 *      @returns directionality of the area of the triangle formed by vertices vtx1, vtx2 and vtx3. It is positive if vtx1, vtx2 and vtx3 are oriented counter-clockwise.
*/
inline ts_double vtx_direct(ts_vertex *vtx1, ts_vertex *vtx2, ts_vertex *vtx3){
    ts_double dX2=vtx2->x-vtx1->x;
    ts_double dY2=vtx2->y-vtx1->y;
    ts_double dZ2=vtx2->z-vtx1->z;
    ts_double dX3=vtx3->x-vtx1->x;
    ts_double dY3=vtx3->y-vtx1->y;
    ts_double dZ3=vtx3->z-vtx1->z;
    ts_double direct=vtx1->x*(dY2*dZ3 -dZ2*dY3)+ 
        vtx1->y*(dZ2*dX3-dX2*dZ3)+
        vtx1->z*(dX2*dY3-dY2*dX3);
    return(direct);    
}


inline ts_bool vertex_add_tristar(ts_vertex *vtx, ts_triangle *tristarmem){
	vtx->tristar_no++;
	vtx->tristar=(ts_triangle **)realloc(vtx->tristar,vtx->tristar_no*sizeof(ts_triangle *));
	if(vtx->tristar==NULL){
			fatal("Reallocation of memory while adding tristar failed.",3);
	}
	vtx->tristar[vtx->tristar_no-1]=tristarmem;
	return TS_SUCCESS;
}


/* Insert neighbour is a function that is required in bondflip. It inserts a
 * neighbour exactly in the right place. */
inline ts_bool vtx_insert_neighbour(ts_vertex *vtx, ts_vertex *nvtx, ts_vertex *vtxm){
//nvtx is a vertex that is to be inserted after vtxm!
        ts_uint i,j,midx;
        vtx->neigh_no++;
        if(vtxm==NULL ||  nvtx==NULL || vtx==NULL)
                fatal("vertex_insert_neighbour: one of pointers has been zero.. Cannot proceed.",3);
        //We need to reallocate space! The pointer *neight must be zero if not having neighbours jey (if neigh_no was 0 at thime of calling
        vtx->neigh=realloc(vtx->neigh,vtx->neigh_no*sizeof(ts_vertex *));
        if(vtx->neigh == NULL){
            fatal("Reallocation of memory failed during insertion of vertex neighbour in vertex_insert_neighbour",3);
        }
        midx=0;
        for(i=0;i<vtx->neigh_no-1;i++) if(vtx->neigh[i]==vtxm) {midx=i; break;}
     //   fprintf(stderr,"midx=%d, vseh=%d\n",midx,vtx->neigh_no-2);
        if(midx==vtx->neigh_no-2) {
            vtx->neigh[vtx->neigh_no-1]=nvtx;
        } else {
            for(j=vtx->neigh_no-2;j>midx;j--) {
                vtx->neigh[j+1]=vtx->neigh[j];
//                vtx->bond_length[j+1]=vtx->bond_length[j];
//                vtx->bond_length_dual[j+1]=vtx->bond_length_dual[j];
            }
            vtx->neigh[midx+1]=nvtx;
        }
    return TS_SUCCESS;
}


/* vtx remove tristar is required in  bondflip. */
/* TODO: Check whether it is important to keep the numbering of tristar
 * elements in some order or not! */
inline ts_bool vtx_remove_tristar(ts_vertex *vtx, ts_triangle *tristar){
    ts_uint i,j=0;
    for(i=0;i<vtx->tristar_no;i++){
        if(vtx->tristar[i]!=tristar){
            vtx->tristar[j]=vtx->tristar[i];
            j++;
        }
    }
    vtx->tristar_no--;
    vtx->tristar=realloc(vtx->tristar,vtx->tristar_no*sizeof(ts_triangle *));
    if(vtx->neigh == NULL){
            fatal("Reallocation of memory failed during insertion of vertex neighbour in vertex_add_neighbour",3);
        }
    return TS_SUCCESS;
}



/* ****************************************************************** */
/* ***** New vertex copy operations. Inherently they are slow.  ***** */
/* ****************************************************************** */

ts_bool vtx_copy(ts_vertex *cvtx, ts_vertex *ovtx){
    memcpy((void *)cvtx,(void *)ovtx,sizeof(ts_vertex));
    cvtx->neigh=NULL;
    cvtx->neigh_no=0;
    cvtx->tristar_no=0;
    cvtx->bond_no=0;
    cvtx->tristar=NULL;
    cvtx->bond=NULL;
    cvtx->cell=NULL;
    return TS_SUCCESS;
}

ts_bool vtx_duplicate(ts_vertex *cvtx, ts_vertex *ovtx){
    memcpy((void *)cvtx,(void *)ovtx,sizeof(ts_vertex));
    return TS_SUCCESS;
}

//TODO: needs to be done
ts_vertex **vtx_neigh_copy(ts_vertex_list *vlist,ts_vertex *ovtx){
        return NULL;
}



ts_vertex_list *vertex_list_copy(ts_vertex_list *ovlist){
    ts_uint i;
    ts_vertex_list *vlist=(ts_vertex_list *)malloc(sizeof(ts_vertex_list));
    vlist=memcpy((void *)vlist, (void *)ovlist, sizeof(ts_vertex_list));
    ts_vertex **vtx=(ts_vertex **)malloc(vlist->n*sizeof(ts_vertex *));
    vlist->vtx=vtx;
    if(vlist->vtx==NULL)
        fatal("Fatal error reserving memory space for vertex list! Could number of requsted vertices be too large?", 100);
    for(i=0;i<vlist->n;i++) {
        vlist->vtx[i]=(ts_vertex *)calloc(1,sizeof(ts_vertex));
        vlist->vtx[i]->idx=i;
        vtx_copy(vlist->vtx[i],ovlist->vtx[i]);
    }

    return vlist;
}

/* vim: set ts=4 sts=4 sw=4 noet : */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"

#define p_diff(i, j) ((long int) (i) - (long int) (j))/ (long int) sizeof(*(i))

ts_bool vertex_list_assign_id(ts_vertex_list *vlist, ts_idx id){
	ts_idx i;	
	for(i=0;i<vlist->n;i++){
		vlist->vtx[i]->id = id;
	}
	return TS_SUCCESS;
}

ts_vertex_list *init_vertex_list(ts_idx N){	
	ts_idx i;
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
    /* 
    for(i=0;i<l;i++){
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

ts_seen_vertex *init_seen_vertex(ts_idx max_size){
    ts_seen_vertex *seen_vtx = (ts_seen_vertex *) malloc(sizeof(ts_seen_vertex));
    if (seen_vtx==NULL){
        fatal("Fatal error reserving memory space for seen_vertex!", 100);
        return NULL;
    }
    seen_vtx->size = max_size;
    seen_vtx->n_prev = 0;
    seen_vtx->n_curr = 0;
    seen_vtx->n_next = 0;
    seen_vtx->n_top = 0;
    seen_vtx->vtx=(ts_vertex **) malloc(max_size*sizeof(ts_vertex *));
    if (seen_vtx->vtx==NULL){
        fatal("Fatal error reserving memory space for vertex list in seen_vertex!", 100);
    }
    return seen_vtx;
}

ts_bool vtx_add_neighbour(ts_vertex *vtx, ts_vertex *nvtx){
    ts_small_idx i;
    /* no neighbour can be null! */
    if(vtx==NULL || nvtx==NULL) return TS_FAIL;
    
    /*if it is already a neighbour don't add it to the list */
    for(i=0; i<vtx->neigh_no;i++){
        if(vtx->neigh[i]==nvtx) return TS_FAIL;
    }
    ts_small_idx nn=++vtx->neigh_no;
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
    ts_small_idx i,j=0;
    for(i=0;i<vtx->neigh_no;i++){
        // fprintf(stderr,"neigh_addr=%ld\n", (long)vtx->neigh[i]);
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
    // fprintf(stderr,"Neigbours=%d\n",nvtx->neigh_no);
    nvtx->neigh_no--;
    nvtx->neigh=(ts_vertex **)realloc(nvtx->neigh,nvtx->neigh_no*sizeof(ts_vertex *));
    // fprintf(stderr,"Neigbours=%d\n",nvtx->neigh_no);
    if(nvtx->neigh == NULL && nvtx->neigh_no!=0)
        fatal("(2) Reallocation of memory failed during removal of vertex neighbour in vtx_remove_neighbour",100);

    return TS_SUCCESS;
}


// create and add bond between two vertices
ts_bool vtx_add_bond(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2){
    ts_small_idx i;
    ts_bond *bond;
    ts_bool bond_in_vtx1=0, bond_in_vtx2=0;
    bond=bond_add(blist,vtx1,vtx2); // registered bond, mallocate and register bond if none exists
    if(bond==NULL) return TS_FAIL;

    // check if bond is already in vtx1
    for(i=0; i<vtx1->bond_no; i++){
        if (bond==vtx1->bond[i]){
            bond_in_vtx1=1;
            break;
        };
    }
    if (!bond_in_vtx1){
        vtx1->bond_no++;
        vtx1->bond=(ts_bond **)realloc(vtx1->bond, vtx1->bond_no*sizeof(ts_bond *)); 
        vtx1->bond[vtx1->bond_no-1]=bond;
    }

    // check if bond is already in vtx2
    for(i=0; i<vtx2->bond_no; i++){
        if (bond==vtx2->bond[i]){
            bond_in_vtx2=1;
            break;
        }
    }
    if (!bond_in_vtx2){
        vtx2->bond_no++;
        vtx2->bond=(ts_bond **)realloc(vtx2->bond, vtx2->bond_no*sizeof(ts_bond *)); 
        vtx2->bond[vtx2->bond_no-1]=bond;
    }
    return TS_SUCCESS;
}

// add neighbors and connect with bond?
ts_bool vtx_add_cneighbour(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex *vtx2){
    ts_bool retval;
    retval=vtx_add_neighbour(vtx1,vtx2);
    // retval=vtx_add_neighbour(vtx2,vtx1);
    if(retval==TS_SUCCESS){
        retval=vtx_add_bond(blist,vtx1,vtx2);
    }
    return retval;
}


// add neighbors and connect with bond? doesn't work for some reason
ts_bool vtx_add_cneighbour2(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex *vtx2){
    ts_bool nei_in_v=0, bond_in_v=0, is_registered=0;
    ts_idx i,b=blist->n;
    ts_small_idx j;
    ts_bond *add_bond=NULL;
    if (vtx1==NULL || vtx2==NULL) return TS_FAIL;
    if (vtx1==vtx2) return TS_FAIL;
    for (i=0; i<blist->n; i++){
        if (blist->bond[i]!=NULL){
            if (in_bond(blist->bond[i], vtx1) && in_bond(blist->bond[i],vtx2)){
                is_registered += 1;
                b=i;
            }
        }
    }
    if (is_registered == 0){
        blist->n++;
        blist->bond = (ts_bond **) realloc(blist->bond,blist->n*sizeof(ts_bond *));
        add_bond = blist->bond[b];
        add_bond = (ts_bond*) malloc(sizeof(ts_bond));
        if(add_bond==NULL){
             fatal("Cannot allocate memory for additional ts_bond.",100);
        } else{
            add_bond->idx = b;
            add_bond->vtx1 = vtx1;
            add_bond->vtx2 = vtx2;
        }
    } else if (is_registered == 1){
        add_bond = blist->bond[b];
    } else {
        ts_fprintf(stdout,"bond is registered %d times", is_registered);
        fatal("bond registration failure",3);
    }
    
    for (j=0; j<vtx1->bond_no; j++){
        if (vtx1->bond[j]==add_bond){
            bond_in_v=1;
            break;
        }
    }
    if (!bond_in_v){
        if (vtx1->bond_no==0) {
            vtx_insert_bond_at(vtx1, add_bond, 0);
        } else {
            vtx_insert_bond_at(vtx1, add_bond, vtx1->bond_no-1);
        }
    }
    for (j=0; j<vtx1->neigh_no; j++){
        if (vtx1->neigh[j]==vtx2){
            nei_in_v=1;
            break;
        }
    }
    if (!nei_in_v){
        if (vtx1->neigh_no==0) {
            vtx1->neigh = (ts_vertex**) malloc(sizeof(ts_vertex *));
            vtx1->neigh[0] = vtx2;
        } else {
            vtx_insert_neigh_at(vtx1, vtx2, vtx1->neigh_no-1);
        }
    }
    nei_in_v=0;
    bond_in_v=0;
    for (j=0; j<vtx2->bond_no; j++){
        if (vtx2->bond[j]==add_bond){
            bond_in_v=1;
            break;
        }
    }
    if (!bond_in_v){
        if (vtx2->bond_no==0) {
            vtx_insert_bond_at(vtx2, add_bond, 0);
        } else {
            vtx_insert_bond_at(vtx2, add_bond, vtx2->bond_no-1);
        }
    }
    for (j=0; j<vtx2->neigh_no; j++){
        if (vtx2->neigh[j]==vtx1){
            nei_in_v=1;
            break;
        }
    }
    if (!nei_in_v){
        if (vtx2->neigh_no==0) {
            vtx2->neigh = (ts_vertex**) malloc(sizeof(ts_vertex*));
            vtx2->neigh[0] = vtx1;
        } else {
            vtx_insert_neigh_at(vtx2, vtx1, vtx2->neigh_no-1);
        }
    }
    ts_fprintf(stdout,"at %u, %u, creating bond %u\n", vtx1->idx, vtx2->idx, add_bond->idx);
    return TS_SUCCESS;
}

/*TODO: write and optimize this urgently before use! */
ts_bool vtx_remove_cneighbour(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex *vtx2){
    // ts_bool retval;
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
    ts_idx i;
    for(i=0;i<vlist->n;i++){
		if(vlist->vtx[i]!=NULL) vtx_free(vlist->vtx[i]);
    }
    //free(*(vlist->vtx));
    free(vlist->vtx);
    free(vlist);
    return TS_SUCCESS;
}

ts_bool seen_vertex_free(ts_seen_vertex *seen_vtx){
    free(seen_vtx->vtx);
    free(seen_vtx);
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
    // as it's set now: override any specific values! must use before initial distribution/ XML parsing!
    // (besides, if it's a real global value, it shouldn't be on the vertices anyway- in vesicle, depend on type)
    ts_idx i; 

    for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]->xk=vesicle->tape->xk0;
    }


    for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]->xk2=vesicle->tape->xk2;
    }

    // in case type and such lead to no curvature2, force, etc. being calculated
    for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]->mean_curvature=0;
        vesicle->vlist->vtx[i]->mean_curvature2=0;
        vesicle->vlist->vtx[i]->gaussian_curvature=0;
        vesicle->vlist->vtx[i]->gaussian_curvature2=0;
        vesicle->vlist->vtx[i]->energy=0;
        vesicle->vlist->vtx[i]->mean_energy=0;
        vesicle->vlist->vtx[i]->mean_energy2=0;
        vesicle->vlist->vtx[i]->gaussian_energy=0;
        vesicle->vlist->vtx[i]->gaussian_energy2=0;
        vesicle->vlist->vtx[i]->new_c1=0;
        vesicle->vlist->vtx[i]->new_c2=0;
        vesicle->vlist->vtx[i]->type=is_adhesive_vtx; // nonbonding, passive, adhesive, isotropic, nonedge,
		vesicle->vlist->vtx[i]->w=0;
		vesicle->vlist->vtx[i]->c=0;
		vesicle->vlist->vtx[i]->f=0;
		vesicle->vlist->vtx[i]->ad_w= vesicle->tape->adhesion_strength;
		vesicle->vlist->vtx[i]->d=0;  // curvature deviator
		vesicle->vlist->vtx[i]->xk = vesicle->tape->xk0;
		vesicle->vlist->vtx[i]->xk2 = 0; // Gauss-Bonet: we only need excess compare to the regular membrane
		vesicle->vlist->vtx[i]->nx=0; //normal
		vesicle->vlist->vtx[i]->ny=0;
		vesicle->vlist->vtx[i]->nz=0;
		vesicle->vlist->vtx[i]->fx=0; //force
		vesicle->vlist->vtx[i]->fy=0;
		vesicle->vlist->vtx[i]->fz=0;
		vesicle->vlist->vtx[i]->tx=0; //director
		vesicle->vlist->vtx[i]->ty=0;
		vesicle->vlist->vtx[i]->tz=0;
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
    ts_double direct=vtx1->x*(dY2*dZ3-dZ2*dY3)+ 
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
        ts_idx i,j,midx;
        vtx->neigh_no++;
        if(vtxm==NULL ||  nvtx==NULL || vtx==NULL)
            fatal("vertex_insert_neighbour: one of pointers has been zero.. Cannot proceed.",3);
        //We need to reallocate space! The pointer *neight must be zero if not having neighbours jey (if neigh_no was 0 at time of calling
        vtx->neigh=realloc(vtx->neigh,vtx->neigh_no*sizeof(ts_vertex *));
        if(vtx->neigh == NULL){
            fatal("Reallocation of memory failed during insertion of vertex neighbour in vertex_insert_neighbour",3);
        }
        midx=0;
        for(i=0;i<vtx->neigh_no-1;i++){
            if(vtx->neigh[i]==vtxm){
                midx=i;
                break;
            }
        }
        // fprintf(stderr,"midx=%d, vseh=%d\n",midx,vtx->neigh_no-2);
        if(midx==vtx->neigh_no-2) {
            vtx->neigh[vtx->neigh_no-1]=nvtx;
        } else {
            for(j=vtx->neigh_no-2;j>midx;j--) {
                vtx->neigh[j+1]=vtx->neigh[j];
                //  vtx->bond_length[j+1]=vtx->bond_length[j];
                //  vtx->bond_length_dual[j+1]=vtx->bond_length_dual[j];
            }
            vtx->neigh[midx+1]=nvtx;
        }
    return TS_SUCCESS;
}


/* vtx remove tristar is required in  bondflip. */
/* TODO: Check whether it is important to keep the numbering of tristar
 * elements in some order or not! */
inline ts_bool vtx_remove_tristar(ts_vertex *vtx, ts_triangle *tristar){
    ts_small_idx i,j=0;
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

// push vertex into the top of the current layer
ts_bool add_vtx_to_seen(ts_seen_vertex *seen_vtx, ts_vertex *vtx){
    if (seen_vtx->n_top == seen_vtx->size){
        // in case we run out of space: dynamically reallocate
        seen_vtx->size *= 2;
        seen_vtx->vtx = realloc(seen_vtx->vtx, sizeof(ts_vertex *) * seen_vtx->size);
        if (seen_vtx->vtx == NULL)
            fatal("Cannot reallocate memory to extend seen_vertex.", 100);
    }
    seen_vtx->vtx[seen_vtx->n_top] = vtx;
    seen_vtx->n_top++;
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
    ts_idx i;
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

// check if vertex is in first 3 layers
ts_bool is_in_seen_vertex(ts_seen_vertex *seen_vertex, ts_vertex *vtx){
    /* check if vertex was already accounted for in the "seen_vertex" list

    we don't need to go over all the vertices: as we iterate
    over the current layer, building the next one,
    we can only encounter neighbors that are new vertices, 
    vertices from the layer we're building, vertices from the currnet layer, 
    and vertice from the previous layers

    For example if we're checking neighbors of C2, we only encounter
    P1 and P2 from the previous layer, C1 and C3 from the current layer,
    A1 which we collected through C1 already in the next layer, and
    N2 and N3 in the next layer which which we haven't seen and want to add
    ..............
     A1 --  N2 --  O1     O-not added yet               seen_vtx[]:
    /  \  /  |  \  / \    N-next layer, not added yet          [...
    C1--C2--N3--O1 - O1   A-next layer, already added            Q1
    | / | \ | \ |  \ |    C-current, neighbor iterated layer     Q2
    P1--P2--C3--N4 - O1   P-previous layer                       P1  <-check_from
    | \ | \ | \ |  \ |    Q-previous previous layer              P2
    Q1--Q2--P3--C4 - N5                                          P3
    ..............                                               C1
                                                                 C2
                                                                 C3
                                                                 C4
                                                                 A1
                                                                  _ <-check_up_to
    */
    ts_idx i;
    for ( i=seen_vertex->n_prev ; i < seen_vertex->n_top ; i++){
        if (seen_vertex->vtx[i] == vtx) return 1;
    }
    return 0;
}

ts_bool advance_seen_vertex_to_next_layer(ts_seen_vertex *seen_vertex){
    seen_vertex->n_prev = seen_vertex->n_curr;
    seen_vertex->n_curr = seen_vertex->n_next;
    seen_vertex->n_next = seen_vertex->n_top;
    return TS_SUCCESS;
}

// swap triastar triangles at index i and j
ts_bool swap_triangles(ts_vertex* vtx, ts_small_idx i, ts_small_idx j){
    ts_triangle* temptri;
    if (i==j) return TS_SUCCESS;
    if (i >= vtx->tristar_no || j >= vtx->tristar_no){
        fatal("attempt to swap triangles outside of range tristar_no",3);
    }
    if (vtx->tristar[i] == NULL || vtx->tristar[j] == NULL){
        fatal("Attempt to swap tristars where one does not exist",3);
    }
    temptri = vtx->tristar[i];
    vtx->tristar[i] = vtx->tristar[j];
    vtx->tristar[j] = temptri;
    return TS_SUCCESS;
}

// swap bond location at index i and j
ts_bool swap_bonds(ts_vertex* vtx, ts_small_idx i, ts_small_idx j){
    ts_bond* tempbond;
    if (i==j) return TS_SUCCESS;
    if (i >= vtx->bond_no || j >= vtx->bond_no){
        fatal("attempt to swap bonds outside of range bond_no",3);
    }
    if (vtx->bond[i] == NULL || vtx->bond[j] == NULL){
        fatal("Attempt to swap bonds where one does not exist",3);
    }
    tempbond = vtx->bond[i];
    vtx->bond[i] = vtx->bond[j];
    vtx->bond[j] = tempbond;
    return TS_SUCCESS;
}

// check if triangle is ordered wrt to v1,v2 ordered
ts_bool tri_ordered(ts_triangle* t, ts_vertex* v1, ts_vertex* v2){
    return (    (t->vertex[0]==v1 && t->vertex[1]==v2) 
             || (t->vertex[1]==v1 && t->vertex[2]==v2)
             || (t->vertex[2]==v1 && t->vertex[0]==v2)) ;
}

// debug: print the order of triangles at a vertex, denoted by the two other neighbors
ts_bool print_tri_order(ts_vertex* vtx){
    ts_small_idx jj, jjm=2;
    for (jj=0; jj<vtx->tristar_no; jj++){
        if (vtx->tristar[jj]->vertex[0]==vtx) jjm=0;
        if (vtx->tristar[jj]->vertex[1]==vtx) jjm=1;
        // if (vtx->tristar[jj]->vertex[2]==vtx) jjm=2;
        fprintf(stdout,"(%ld, %ld), ", //(long int) vtx->tristar[jj]->vertex[jjm] - (long int) vtx, 
                                            p_diff(vtx->tristar[jj]->vertex[(jjm+1)%3], vtx),
                                            p_diff(vtx->tristar[jj]->vertex[(jjm+2)%3], vtx));
    }
    fprintf(stdout,"\n");
    return TS_SUCCESS;
}

ts_bool print_vertex_ordered(ts_vertex* vtx){
    ts_small_idx j;
    fprintf(stdout, "|%u|: (", vtx->idx);
    for (j=0; j<vtx->neigh_no; j++){
        fprintf(stdout, "%u,", vtx->neigh[j]->idx);
    }
    fprintf(stdout,")\n");
    for (j=0; j<vtx->tristar_no; j++){
        fprintf(stdout,"%u:{%u, %u, %u}, ", vtx->tristar[j]->idx, vtx->tristar[j]->vertex[0]->idx,
                                 vtx->tristar[j]->vertex[1]->idx, vtx->tristar[j]->vertex[2]->idx);
    }
    fprintf(stdout,"\n");
    for (j=0; j<vtx->bond_no; j++){
        fprintf(stdout,"%u:(%u, %u), ", vtx->bond[j]->idx, vtx->bond[j]->vtx1->idx,
                                 vtx->bond[j]->vtx2->idx);
    }
    fprintf(stdout,"\n");
    return TS_SUCCESS;
}

// order the triangles of the vertex according to the neighbors
// vtx->tristar[i] = {vtx, vtx->neigh[i], vtx->neigh[i+1]} (up to modulus)
ts_bool order_vertex_triangles(ts_vertex* vtx){

    ts_vertex* vl, *vr;
    ts_small_idx t, jj=0, li, ri, rri=0, lli=1;
    ts_triangle* jt;
    ts_bond* bi;
    if (vtx->tristar_no != vtx->neigh_no || vtx->bond_no != vtx->neigh_no){
        print_vertex_ordered(vtx);
        ts_fprintf(stdout, "vertex %u with neigh_no=%u, tristar_no=%u, bond_no=%u\n",vtx->idx,vtx->tristar_no, vtx->neigh_no, vtx->bond_no);
        fatal("Unable to order a vertex",3);
    }
    /* reorder the triangles: 
    - find first triangle with neighbors 0,1 , swap it to 0
    - keep a leftmost vertex 0 and rightmost vertex 1,
    - any triangle that has the leftmost, add to the left (No, No-1...), update leftmost, same with right (2,3,...)
    - repeat until rightmost and leftmost meet
    * fuse the first triangle search with the first left-right search
    */

    // find first, and also second and last triangles (0,1) (1,[2]),...([end],0) 
    vl = vtx->neigh[0];
    vr = vtx->neigh[1];
    for (t=0; t<vtx->tristar_no; t++){
        jt = vtx->tristar[t];
        if (in_tri(jt,vl)){
            if (in_tri(jt,vr)){
                jj = t;
            }
            else{
                lli = t;
            }
          
        }
        else if (in_tri(jt,vr)){
            rri = t;
        }  
    }
    swap_triangles(vtx, jj, 0); // move (0,1) to 0

    if (lli==0) lli=jj; // if 0 was occupied by (end,0), it is now at jj 
    swap_triangles(vtx, lli, vtx->tristar_no-1);

    if (rri==0) rri=jj; // if 0 was occupied by (1,2), it is now at jj
    if (rri==vtx->tristar_no-1) rri=lli; // if end was occupied by (1,2), it is now at lli
    swap_triangles(vtx, rri, 1);

    // now triangles can only be left of left or right of right
    ri = 2;
    li = vtx->neigh_no-1;
    while (ri<li){ 
        for (t=ri; t<li; t++){
            vl = vtx->neigh[li];
            vr = vtx->neigh[ri];
            jt = vtx->tristar[t];
            if (in_tri(jt, vl)){
                li-=1;
                swap_triangles(vtx, t, li);

                if (ri+1==li) break;
            }
            if (in_tri(jt, vr)){
                swap_triangles(vtx, t, ri);

                ri+=1;
                if (ri+1==li) break;
            }
        }
    }

    // and for bonds
    for (li=0; li<vtx->neigh_no; li++){
        vl = vtx->neigh[li];
        for (ri=li; ri<vtx->neigh_no; ri++){
            bi = vtx->bond[ri];
            if (in_bond(bi,vl)){
                swap_bonds(vtx, li, ri);
                break;
            }
        }
    }

    //for (i=0; i<vtx->neigh_no; i++){
    //    if (!in_tri(vtx->tristar[i], vtx->neigh[i]) && !in_tri(vtx->tristar[i], vtx->neigh[next_small(i,vtx->neigh_no)])){
    //        fatal("not ordered in triangles",3);
    //    }
    //    if (!in_bond(vtx->bond[i], vtx->neigh[i])){
    //        fatal("not ordered in bonds",3);
    //    }
    //}
    return TS_SUCCESS;
}


ts_bool assert_vtx_ordered(ts_vertex* vtx){
    ts_small_idx i;
    for (i=0; i<vtx->neigh_no; i++){
        if (!in_tri(vtx->tristar[i], vtx->neigh[i]) && !in_tri(vtx->tristar[i], vtx->neigh[next_small(i,vtx->neigh_no)])){
            print_vertex_ordered(vtx);
            fatal("not ordered in triangles",3);
        }
        if (!in_bond(vtx->bond[i], vtx->neigh[i])){
            print_vertex_ordered(vtx);
            fatal("not ordered in bonds",3);
        }
    }
    return TS_SUCCESS;
}

// add tristar to vertex at index- shift other vertex. use to maintain tristar order
ts_bool vtx_insert_tristar_at(ts_vertex *vtx, ts_triangle *tristarmem, ts_small_idx i){
    if (i > vtx->tristar_no){
        fatal("attempt to add tristar above tristar_no",3);
    }
	vtx->tristar_no++;
	vtx->tristar=(ts_triangle **)realloc(vtx->tristar,vtx->tristar_no*sizeof(ts_triangle *));
	if(vtx->tristar==NULL){
			fatal("Reallocation of memory while adding tristar failed.",3);
	}
    if (i+1 != vtx->tristar_no){ //no need to shift, don't want to tempt memmove(outside, edge, 0)
        memmove(vtx->tristar+i+1, vtx->tristar+i,(vtx->tristar_no-i-1)*sizeof(ts_triangle*));
    }
	vtx->tristar[i]=tristarmem;
	return TS_SUCCESS;
}

// add tristar from vertex at index- shift other vertex. use to maintain tristar order
ts_bool vtx_remove_tristar_at(ts_vertex *vtx, ts_small_idx i){

    if ( i >= vtx->tristar_no ){
        fatal("attempt to remove triangle above tristar_no",3);
    }
    if ( i+1 != vtx->tristar_no) {  //no need to shift, don't want to tempt memmove(edge, outside, 0)
        memmove(vtx->tristar+i, vtx->tristar+i+1, (vtx->tristar_no-i-1)*sizeof(ts_triangle*));
    }
    vtx->tristar_no--;
    vtx->tristar=(ts_triangle **)realloc(vtx->tristar,vtx->tristar_no*sizeof(ts_triangle *));
    if(vtx->tristar == NULL){
            fatal("Reallocation of memory failed during removal of tristar",3);
        }
    return TS_SUCCESS;
}

// It is the caller's responsibility to 1. add the vertex to the neighbor's list, 2. maintain order
ts_bool vtx_insert_neigh_at(ts_vertex *vtx, ts_vertex *vtxmem, ts_small_idx i){

    if ( i>vtx->neigh_no){
        fatal("attempt to add neighbor above neigh_no",3);
    }
	vtx->neigh_no++;
	vtx->neigh=(ts_vertex **)realloc(vtx->neigh,vtx->neigh_no*sizeof(ts_vertex *));
	if(vtx->neigh==NULL){
			fatal("Reallocation of memory while adding neighbor failed.",3);
	}
    if (i+1 != vtx->neigh_no) { //no need to shift, don't want to tempt memmove(outside, edge, 0)
        memmove(vtx->neigh+i+1, vtx->neigh+i,(vtx->neigh_no-i-1)*sizeof(ts_vertex*));
    }
	vtx->neigh[i]=vtxmem;
	return TS_SUCCESS;
}

// It is the caller's responsibility to remove the vertex from the neighbor's list
ts_bool vtx_remove_neigh_at(ts_vertex *vtx, ts_small_idx i){

    if ( i>=vtx->neigh_no){
        fatal("attempt to remove vertex above neigh_no",3);
    }
    if (i+1 != vtx->neigh_no) { //no need to shift, don't want to tempt memmove(edge, outside, 0)
        memmove(vtx->neigh+i, vtx->neigh+i+1,(vtx->neigh_no-i-1)*sizeof(ts_vertex*));
    }
    vtx->neigh_no--;
    vtx->neigh=(ts_vertex **)realloc(vtx->neigh,vtx->neigh_no*sizeof(ts_vertex *));
    if(vtx->neigh == NULL){
            fatal("Reallocation of memory failed during removal of neighbor",3);
        }
    return TS_SUCCESS;
}


// It is the caller's responsibility to 1. add the bond to the neighbor's list, 2. maintain order
ts_bool vtx_insert_bond_at(ts_vertex *vtx, ts_bond *bondmem, ts_small_idx i){

    if ( i>vtx->bond_no){
        fatal("attempt to add bond above bond_no",3);
    }
	vtx->bond_no++;
	vtx->bond=(ts_bond **)realloc(vtx->bond,vtx->bond_no*sizeof(ts_bond *));
	if(vtx->bond==NULL){
			fatal("Reallocation of memory while adding bond failed.",3);
	}
    if (i+1 != vtx->bond_no) { //no need to shift, don't want to tempt memmove(outside, edge, 0)
        memmove(vtx->bond+i+1, vtx->bond+i, (vtx->bond_no-i-1)*sizeof(ts_bond*));
    }
	vtx->bond[i]=bondmem;
	return TS_SUCCESS;
}

// It is the caller's responsibility to remove the bond from the neighbor's list
ts_bool vtx_remove_bond_at(ts_vertex *vtx, ts_small_idx i){

    if ( i>=vtx->bond_no){
        fatal("attempt to remove bond above bond_no",3);
    }
    if (i+1 != vtx->bond_no) { //no need to shift, don't want to tempt memmove(edge, outside, 0)
        memmove(vtx->bond+i, vtx->bond+i+1,(vtx->bond_no-i-1)*sizeof(ts_bond*));
    }
    vtx->bond_no--;
    vtx->bond=(ts_bond **)realloc(vtx->bond,vtx->bond_no*sizeof(ts_bond *));
    if(vtx->bond == NULL){
            fatal("Reallocation of memory failed during removal of bond",3);
        }
    return TS_SUCCESS;
}

// It is the caller's responsibility to 1. add to the neighbor's list, 2. maintain order
ts_bool vtx_insert_at(ts_vertex *vtx, ts_vertex *vtx_add, ts_bond* bond_add, ts_triangle* tri_add, ts_small_idx i){

    vtx_insert_neigh_at(vtx, vtx_add, i);
    vtx_insert_tristar_at(vtx, tri_add, i);
    vtx_insert_bond_at(vtx, bond_add, i);
	return TS_SUCCESS;
}

// It is the caller's responsibility to remove the vertex from the neighbor's list
ts_bool vtx_remove_at(ts_vertex *vtx, ts_small_idx i){

    vtx_remove_neigh_at(vtx, i);
    vtx_remove_tristar_at(vtx, i);
    vtx_remove_bond_at(vtx, i);
    return TS_SUCCESS;
}

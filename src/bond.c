/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include "general.h"
#include "vertex.h"

ts_bond_list *init_bond_list(){
    ts_bond_list *blist=(ts_bond_list *)malloc(sizeof(ts_bond_list));
	blist->n=0;
	blist->bond=NULL;
	return blist;
}

ts_bond  *bond_add(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex *vtx2){
    ts_uint i;
    /* no vertices must be null! */
    if(vtx1==NULL || vtx2==NULL) return NULL;
    /* Verify if the bond already exists... Don't do multiple bonds */
    for(i=0;i<blist->n;i++){
       if((blist->bond[i]->vtx1==vtx1 || blist->bond[i]->vtx1==vtx2)
        &&(blist->bond[i]->vtx2==vtx1 || blist->bond[i]->vtx2==vtx2))
        return blist->bond[i];
    }
	blist->n++;
	blist->bond=(ts_bond **)realloc(blist->bond,blist->n*sizeof(ts_bond *));
	if(blist->bond==NULL) fatal("Cannot reallocate memory for additional **ts_bond.",100);
    blist->bond[blist->n-1]=(ts_bond *)malloc(sizeof(ts_bond));
    if(blist->bond[blist->n-1]==NULL) fatal("Cannot allocate memory for additional *ts_bond.",100);
    
	//NOW insert vertices into data!	
	blist->bond[blist->n - 1]->vtx1=vtx1;	
	blist->bond[blist->n - 1]->vtx2=vtx2;
	blist->bond[blist->n - 1]->tainted=0;
	blist->bond[blist->n - 1]->idx=blist->n-1;
    //Should we calculate bond length NOW?
	
	return blist->bond[blist->n-1];
}


ts_bool bond_vector(ts_bond *bond){
	
	bond->x = bond->vtx1->x - bond->vtx2->x;
	bond->y = bond->vtx1->y - bond->vtx2->y;
	bond->z = bond->vtx1->z - bond->vtx2->z;

	return TS_SUCCESS;	
}


ts_bool bond_list_free(ts_bond_list *blist){
    ts_uint i;
    for(i=0;i<blist->n;i++){
    free(blist->bond[i]);
    }
    free(blist->bond);
    free(blist);
    return TS_SUCCESS;
}

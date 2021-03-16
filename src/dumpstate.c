/* vim: set ts=4 sts=4 sw=4 noet : */
#include <string.h>
#include "general.h"
#include <stdio.h>
#include "initial_distribution.h"
#include "vesicle.h"
#include "dumpstate.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "vertex.h"
#include "energy.h"

ts_vesicle *vtk2vesicle(char *filename, ts_tape *tape){

    ts_uint nshell=tape->nshell;
    ts_uint ncmax1=tape->ncxmax;
    ts_uint ncmax2=tape->ncymax;
    ts_uint ncmax3=tape->nczmax;
    ts_double stepsize=tape->stepsize;

    ts_uint no_vertices=5*nshell*nshell+2;
    ts_vesicle *vesicle=init_vesicle(no_vertices,ncmax1,ncmax2,ncmax3,stepsize);
    vesicle->nshell=nshell;
    parse_vtk(filename, vesicle);
    exit(1);
    return vesicle;
}


ts_bool parse_vtk(char *filename, ts_vesicle *vesicle){
    xmlDoc *doc;
    xmlNode *root_element=NULL;
    xmlNode *cur_node = NULL;
    doc = xmlReadFile(filename, NULL, 0);
    root_element=xmlDocGetRootElement(doc);
    cur_node=root_element->children;
    while(cur_node!=NULL){
//        fprintf(stderr,"Node name is: %s\n",cur_node->name);
        if(strcmp((char *)cur_node->name,"UnstructuredGrid")==0) break;
        cur_node=cur_node->next;
    }

    cur_node=cur_node->children;
    while(cur_node!=NULL){
//        fprintf(stderr,"Node name is: %s\n",cur_node->name);
        cur_node=cur_node->next;
        if(strcmp((char *)cur_node->name,"Piece")==0) break;
    }

    cur_node=cur_node->children;
    while(cur_node!=NULL){
        //fprintf(stderr,"Node name is: %s\n",cur_node->name);
        cur_node=cur_node->next;
        if(strcmp((char *)cur_node->name,"PointData")==0) vtk_index2vesicle(cur_node->children->next->children, vesicle);
        if(strcmp((char *)cur_node->name,"Points")==0) vtk_coordinates(cur_node->children->next->children,vesicle);
        if(strcmp((char *)cur_node->name,"Cells")==0) {
		vtk_neighbours(cur_node->children->next->children,vesicle);
        	break; //segfaults, because it finds another cells. Why?
	}	
    }
	//we have vertices and neighbour relations all set, but unordered. Let's do bonds, ordering, triangles...
	vesicle->vlist = vtk_sort_neighbours(vesicle->blist,vesicle->vlist);
	init_triangles(vesicle);
	init_triangle_neighbours(vesicle);
	init_common_vertex_triangle_neighbours(vesicle);
	init_normal_vectors(vesicle->tlist);
	mean_curvature_and_energy(vesicle);
	ts_fprintf(stdout,"restoring from vtk dump finished!\n");
    return TS_SUCCESS;
}


ts_bool vtk_index2vesicle(xmlNode *node, ts_vesicle *vesicle){
    //fprintf(stderr, "vsebina: %s\n",node->content);
    ts_uint i;
    char *token;
    token = strtok((char *)node->content, " ");
    for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]->idx=atoi(token);
        token=strtok(NULL," ");
    }
    //fprintf(stderr,"idx[11]=%d\n",vesicle->vlist->vtx[11]->idx);
    return TS_SUCCESS;
}

ts_bool vtk_coordinates(xmlNode *node, ts_vesicle *vesicle){
    //fprintf(stderr, "vsebina: %s\n",node->content);
    ts_uint i;
    char *token;
    token = strtok((char *)node->content, " ");
    for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]->x=atof(token);
//	fprintf(stderr,"%e, ",atof(token));
        token=strtok(NULL," ");
        vesicle->vlist->vtx[i]->y=atof(token);
        token=strtok(NULL,"\n");
        vesicle->vlist->vtx[i]->z=atof(token);
        if(i>0) token=strtok(NULL," ");
    }
    return TS_SUCCESS;

}

ts_bool vtk_neighbours(xmlNode *node, ts_vesicle *vesicle){
 //	fprintf(stderr, "vsebina: %s\n",node->content);
	ts_uint i;
	ts_uint idx1, idx2;
    char *token;
    token = strtok((char *)node->content, " ");
    for(i=0;i<3*(vesicle->vlist->n-2);i++){
	idx1=atoi(token);
        token=strtok(NULL,"\n");
	idx2=atoi(token);
        vtx_add_neighbour(vesicle->vlist->vtx[idx1], vesicle->vlist->vtx[idx2]);
        vtx_add_neighbour(vesicle->vlist->vtx[idx2], vesicle->vlist->vtx[idx1]);
	fprintf(stderr, "%d. povezujem %d in %d\n",i,idx1,idx2);
        if(i>0) token=strtok(NULL," ");
    }
	return TS_SUCCESS;
}

/* this is almost exact copy of init_sort_neighbours */
ts_vertex_list *vtk_sort_neighbours(ts_bond_list *blist,ts_vertex_list *vlist){
	ts_vertex **vtx=vlist->vtx -1; // take a look at dipyramid function for comment.
	ts_uint i,l,j,jj,jjj,k=0;   
    ts_double eps=0.001; // Take a look if EPS from math.h can be used

/*lets initialize memory for temporary vertex_list. Should we write a function instead */
    ts_vertex_list *tvlist=vertex_list_copy(vlist);
    ts_vertex **tvtx=tvlist->vtx -1;  /* again to compensate for 0-indexing */

	ts_double dist2; // Square of distance of neighbours
    ts_double direct; // Something, dont know what, but could be normal of some kind
	for(i=1;i<=vlist->n;i++){
		k++; // WHY i IS NOT GOOD??
       	vtx_add_cneighbour(blist,tvtx[k], tvtx[vtx[i]->neigh[0]->idx+1]); //always add 1st
       	jjj=1;
       	jj=1;
       	for(l=2;l<=vtx[i]->neigh_no;l++){
           	for(j=2;j<=vtx[i]->neigh_no;j++){
               	dist2=vtx_distance_sq(vtx[i]->neigh[j-1],vtx[i]->neigh[jj-1]);
               	direct=vtx_direct(vtx[i],vtx[i]->neigh[j-1],vtx[i]->neigh[jj-1]);
// TODO: check if fabs can be used with all floating point types!!
               	if( (direct>0.0) && (j!=jjj) ){
           			vtx_add_cneighbour(blist,tvtx[k],tvtx[vtx[i]->neigh[j-1]->idx+1]);
           			jjj=jj;
           			jj=j;
           		//	break;
           		}
       		}
       	}
	if(vtx[i]->neigh_no!=tvtx[i]->neigh_no){
		fprintf(stderr,"Hej! pri vtx=%d se razikuje stevilo sosedov (prej=%d, potem=%d)\n",i, vtx[i]->neigh_no, tvtx[i]->neigh_no);
		}
	}
	if(eps>dist2);
/* We use the temporary vertex for our main vertices and we abandon main
 * vertices, because their neighbours are not correctly ordered */
   // tvtx=vlist->vtx;
   // vlist->vtx=tvtx;
   // tvlist->vtx=vtx;
    vtx_list_free(vlist);
/* Let's make a check if the number of bonds is correct */
    if((blist->n)!=3*(tvlist->n-2)){
        ts_fprintf(stderr,"Number of bonds is %u should be %u!\n", blist->n, 3*(tvlist->n-2));
        fatal("Number of bonds is not 3*(no_vertex-2).",4);
    }

	return tvlist;
}



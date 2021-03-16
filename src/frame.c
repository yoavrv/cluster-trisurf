/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include "general.h"
#include "cell.h"
#include "frame.h"


#include "triangle.h"
ts_bool centermass(ts_vesicle *vesicle){
    ts_uint i,j, n=vesicle->vlist->n;
    ts_vertex **vtx=vesicle->vlist->vtx;
	ts_double temp_z_cm=0;
    vesicle->cm[0]=0;
    vesicle->cm[1]=0;
    vesicle->cm[2]=0;
    for(i=0;i<n;i++){
        vesicle->cm[0]+=vtx[i]->x;
        vesicle->cm[1]+=vtx[i]->y;
        vesicle->cm[2]+=vtx[i]->z; 
    } 
    vesicle->cm[0]/=(ts_float)n;
    vesicle->cm[1]/=(ts_float)n;
    vesicle->cm[2]/=(ts_float)n;

	//center mass for z component does not  change if we confine the vesicle with plates
	if(vesicle->tape->plane_confinement_switch){
		temp_z_cm=vesicle->cm[2];
		vesicle->cm[2]=0;
	}

	//center mass for z component does not  change when adhesion is switched on
	if(vesicle->tape->adhesion_switch){
		temp_z_cm=vesicle->cm[2];
		vesicle->cm[2]=0;
	}
	//center mass for x component does not  change for cylyndrical substrate
	if(vesicle->tape->type_of_adhesion_model==4){
		vesicle->cm[0]=0;
	}

	//center mass for x and y component does not  change for spherical substrate
	if(vesicle->tape->type_of_adhesion_model==3){
		vesicle->cm[0]=0;
		vesicle->cm[1]=0;
	}


    for(i=0;i<n;i++){
        vtx[i]->x-=vesicle->cm[0];
        vtx[i]->y-=vesicle->cm[1];
        vtx[i]->z-=vesicle->cm[2]; 
    } 
//move polymers for the same vector as we moved vesicle
	for(i=0;i<vesicle->poly_list->n;i++){
		for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
			vesicle->poly_list->poly[i]->vlist->vtx[j]->x-=vesicle->cm[0];
			vesicle->poly_list->poly[i]->vlist->vtx[j]->y-=vesicle->cm[1];
			vesicle->poly_list->poly[i]->vlist->vtx[j]->z-=vesicle->cm[2];
		}
    }
//move filaments for the same vector as we moved vesicle
	for(i=0;i<vesicle->filament_list->n;i++){
		for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
			vesicle->filament_list->poly[i]->vlist->vtx[j]->x-=vesicle->cm[0];
			vesicle->filament_list->poly[i]->vlist->vtx[j]->y-=vesicle->cm[1];
			vesicle->filament_list->poly[i]->vlist->vtx[j]->z-=vesicle->cm[2];
		}
    }
//move nucleus for the same vector as we moved vesicle
	vesicle->nucleus_center[0]-=vesicle->cm[0];
	vesicle->nucleus_center[1]-=vesicle->cm[1];
	vesicle->nucleus_center[2]-=vesicle->cm[2];

    vesicle->cm[0]=0.0;
    vesicle->cm[1]=0.0;
	if(vesicle->tape->plane_confinement_switch){
		vesicle->cm[2]=temp_z_cm;
	}

    for(i=0;i<vesicle->tlist->n;i++){
        triangle_normal_vector(vesicle->tlist->tria[i]);
    }


    return TS_SUCCESS;
}

ts_bool cell_occupation(ts_vesicle *vesicle){
    ts_uint i,j,cellidx, n=vesicle->vlist->n;

    cell_list_cell_occupation_clear(vesicle->clist);
    for(i=0;i<n;i++){
    cellidx=vertex_self_avoidance(vesicle, vesicle->vlist->vtx[i]);
//	already done in cell_add_vertex
//    vesicle->vlist->vtx[i]->cell=vesicle->clist->cell[cellidx];

    cell_add_vertex(vesicle->clist->cell[cellidx],vesicle->vlist->vtx[i]);
    }

//Add all polymers to cells
if(vesicle->poly_list!=NULL){
    for(i=0;i<vesicle->poly_list->n;i++){
	for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
    	cellidx=vertex_self_avoidance(vesicle, vesicle->poly_list->poly[i]->vlist->vtx[j]);
    	cell_add_vertex(vesicle->clist->cell[cellidx],vesicle->poly_list->poly[i]->vlist->vtx[j]);
	}
    }
}
//Add all filaments to cells
if(vesicle->filament_list!=NULL){
     for(i=0;i<vesicle->filament_list->n;i++){
	for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
    	cellidx=vertex_self_avoidance(vesicle, vesicle->filament_list->poly[i]->vlist->vtx[j]);
    	cell_add_vertex(vesicle->clist->cell[cellidx],vesicle->filament_list->poly[i]->vlist->vtx[j]);
	}
    }
}   

    return TS_SUCCESS;
}

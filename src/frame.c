/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include "general.h"
#include "cell.h"
#include "frame.h"
#include "triangle.h"

/** @brief restore center mass of vesicle to 0,0,?0
 *  
 *  remove centerm from all vertices
 *  in several cases, does not update cm from z changes
 *
 *  @returns TS_SUCCESS on success
*/
ts_bool centermass(ts_vesicle *vesicle){
    ts_idx i,j, n=vesicle->vlist->n;
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
	if(vesicle->tape->type_of_adhesion_model==model_plane_potential){
		temp_z_cm=vesicle->cm[2];
		vesicle->cm[2]=0;
	}
	//center mass for x and z component does not  change for cylyndrical substrate
	else if(vesicle->tape->type_of_adhesion_model==model_cylindrical_potential){
		temp_z_cm=vesicle->cm[2];
		vesicle->cm[2]=0;
		vesicle->cm[0]=0;
	} 
	//center mass for x y, and z component does not  change for spherical substrate
	else if(vesicle->tape->type_of_adhesion_model==model_spherical_potential){
		temp_z_cm=vesicle->cm[2];
		vesicle->cm[2]=0;
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

// update cell occupation (? which is not done each vertex move)
ts_bool cell_occupation(ts_vesicle *vesicle){
	ts_idx i, j, n=vesicle->vlist->n;
    ts_uint cellidx;
	ts_cell_list *clist=vesicle->clist; //aliases for less wide code
	ts_vertex_list *vlist;
	ts_poly_list *poly_list=vesicle->poly_list;
	ts_poly_list *filament_list=vesicle->filament_list;

    cell_list_cell_occupation_clear(clist);

	vlist = vesicle->vlist;
    for(i=0;i<n;i++){
    	cellidx=vertex_self_avoidance(vesicle, vlist->vtx[i]);
		//	already done in cell_add_vertex
		// vesicle->vlist->vtx[i]->cell=vesicle->clist->cell[cellidx];

    	cell_add_vertex(clist->cell[cellidx], vlist->vtx[i]);
    }

	//Add all polymers to cells
	if(poly_list!=NULL){
    	for(i=0;i<poly_list->n;i++){
			vlist = poly_list->poly[i]->vlist;
			for(j=0;j<vlist->n;j++){
    			cellidx=vertex_self_avoidance(vesicle, vlist->vtx[j]);
    			cell_add_vertex(clist->cell[cellidx],vlist->vtx[j]);
			}
    	}
	}

	//Add all filaments to cells
	if(filament_list!=NULL){
    	for(i=0;i<filament_list->n;i++){
			vlist = filament_list->poly[i]->vlist;
			for(j=0;j<vlist->n;j++){
    			cellidx=vertex_self_avoidance(vesicle, vlist->vtx[j]);
    			cell_add_vertex(clist->cell[cellidx], vlist->vtx[j]);
			}
    	}
	}   

    return TS_SUCCESS;
}

// same as cell_occupation but with special check to prevent
// overfilling a single cell with too many {x=0,y=0,z=0} vertices
ts_bool initialization_cell_occupation(ts_vesicle *vesicle){
	ts_bool is_initialized=0;
	ts_uint rolling = 0;
	ts_idx i, j, n=vesicle->vlist->n;
    ts_uint cellidx;
	ts_cell_list *clist=vesicle->clist; //aliases for less wide code
	ts_vertex_list *vlist;
	ts_poly_list *poly_list=vesicle->poly_list;
	ts_poly_list *filament_list=vesicle->filament_list;

    cell_list_cell_occupation_clear(clist);

	vlist = vesicle->vlist;
	for (i=0; i<n; i++){
		if (vlist->vtx[i]->x!=0 ||vlist->vtx[i]->y!=0 || vlist->vtx[i]->z!=0)
			is_initialized=1;
	}
    for(i=0;i<n;i++){
    	cellidx=vertex_self_avoidance(vesicle, vlist->vtx[i]);
		//	already done in cell_add_vertex
		// vesicle->vlist->vtx[i]->cell=vesicle->clist->cell[cellidx];
		if (is_initialized){
    		cell_add_vertex(clist->cell[cellidx], vlist->vtx[i]);
		}
		else{
			cell_add_vertex(clist->cell[rolling++%clist->cellno], vlist->vtx[i]);
		}
    }

	//Add all polymers to cells
	if(poly_list!=NULL){
    	for(i=0;i<poly_list->n;i++){
			vlist = poly_list->poly[i]->vlist;
			for (j=0; j<vlist->n; j++){
				if (vlist->vtx[j]->x!=0 || vlist->vtx[j]->y!=0 || vlist->vtx[j]->z!=0)
					is_initialized=1;
			}
			for(j=0;j<vlist->n;j++){
    			cellidx=vertex_self_avoidance(vesicle, vlist->vtx[j]);
				if (is_initialized){
    				cell_add_vertex(clist->cell[cellidx], vlist->vtx[i]);
				}
				else{
					cell_add_vertex(clist->cell[rolling++%clist->cellno], vlist->vtx[i]);
				}
			}
    	}
	}

	//Add all filaments to cells
	if(filament_list!=NULL){
    	for(i=0;i<filament_list->n;i++){
			vlist = filament_list->poly[i]->vlist;
			for (j=0; j<vlist->n; j++){
				if (vlist->vtx[j]->x!=0 || vlist->vtx[j]->y!=0 || vlist->vtx[j]->z!=0)
					is_initialized=1;
			}
			for(j=0;j<vlist->n;j++){
    			cellidx=vertex_self_avoidance(vesicle, vlist->vtx[j]);
				if (is_initialized){
    				cell_add_vertex(clist->cell[cellidx], vlist->vtx[i]);
				}
				else{
					cell_add_vertex(clist->cell[rolling++%clist->cellno], vlist->vtx[i]);
				}
			}
    	}
	}   

    return TS_SUCCESS;
}

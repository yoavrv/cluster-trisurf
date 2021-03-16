/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include "general.h"
#include "constvol.h"
#include "triangle.h"
#include "energy.h"
#include "vertex.h"
#include "cell.h"

ts_bool constvolume(ts_vesicle *vesicle, ts_vertex *vtx_avoid, ts_double Vol, ts_double *retEnergy, ts_vertex **vtx_moved_retval, ts_vertex **vtx_backup){
    ts_vertex *vtx_moved;
    ts_uint vtxind,i,j;
    ts_uint Ntries=100;
	ts_vertex *backupvtx;
    ts_double Rv, dh, dvol, volFirst, voldiff, oenergy,delta_energy;
    backupvtx=(ts_vertex *)calloc(sizeof(ts_vertex),10);
    ts_double l0 = (1.0 + sqrt(vesicle->dmax))/2.0; //make this a global constant if necessary
    for(i=0;i<Ntries;i++){
        vtxind=rand() % vesicle->vlist->n;
        vtx_moved=vesicle->vlist->vtx[vtxind];
        /* chosen vertex must not be a nearest neighbour. WASTODO: probably must.
         * DONE: solved below, when using different algorithm.
         * extend search in case of bondflip */
/*        if(vtx_moved==vtx_avoid) continue;
        for(j=0;j<vtx_moved->neigh_no;j++){
            if(vtx_moved->neigh[j]==vtx_avoid) continue;
        }
*/
/* different check of nearest neighbour (and second nearest neighbour) check.
 * Checking the distance between vertex and vertex to be moved to assure
 * constant volume. Solves upper todo problem. */

        if(vtx_distance_sq(vtx_moved,vtx_avoid)<16.0*vesicle->dmax){
            continue;
        }
        
	    memcpy((void *)&backupvtx[0],(void *)vtx_moved,sizeof(ts_vertex));

        //move vertex in specified direction. first try, test move!
        Rv=sqrt(pow(vtx_moved->x,2)+pow(vtx_moved->y,2)+pow(vtx_moved->z,2));
        dh=2.0*Vol/(sqrt(3.0)*l0*l0);
	    vtx_moved->x=vtx_moved->x*(1.0-dh/Rv);
	    vtx_moved->y=vtx_moved->y*(1.0-dh/Rv);
	    vtx_moved->z=vtx_moved->z*(1.0-dh/Rv);


//SKIPPING FIRST CHECK of CONSTRAINTS. This is not a final move.
        //check for constraints
/*          if(constvolConstraintCheck(vesicle, vtx_moved)==TS_FAIL){
		    vtx_moved=memcpy((void *)vtx_moved,(void *)&backupvtx[0],sizeof(ts_vertex));
//            continue;
                break;
        }  */
        // All checks OK!

        dvol=0.0;
        for(j=0;j<vtx_moved->tristar_no;j++){
            dvol-=vtx_moved->tristar[j]->volume;
        }
        volFirst=dvol;
        for(j=0;j<vtx_moved->tristar_no;j++){
            triangle_normal_vector(vtx_moved->tristar[j]);
            dvol+=vtx_moved->tristar[j]->volume;
        }
//TODO: here there is a bug. Don't know where, but preliminary success can
//happen sometimes. And when I do this checks, constant value is not achieved
//anymore
/*        voldiff=dvol-Vol;
        if(fabs(voldiff)/vesicle->volume < vesicle->tape->constvolprecision){
            //calculate energy, return change in energy...
             oenergy=vtx_moved->energy;
            energy_vertex(vtx_moved);
            delta_energy=vtx_moved->xk*(vtx_moved->energy - oenergy);
            //the same is done for neighbouring vertices
            for(j=0;j<vtx_moved->neigh_no;j++){
                oenergy=vtx_moved->neigh[j]->energy;
                energy_vertex(vtx_moved->neigh[j]);
                delta_energy+=vtx_moved->neigh[j]->xk*(vtx_moved->neigh[j]->energy-oenergy);
            }
            *retEnergy=delta_energy;
            *vtx_backup=backupvtx;
            *vtx_moved_retval=vtx_moved;
            fprintf(stderr, "Preliminary success\n");
            return TS_SUCCESS;
        }        */

//            fprintf(stderr, "Step 2 success\n");
        //do it again ;)
        dh=Vol*dh/dvol;
		vtx_moved=memcpy((void *)vtx_moved,(void *)&backupvtx[0],sizeof(ts_vertex));
        vtx_moved->x=vtx_moved->x*(1-dh/Rv);
	    vtx_moved->y=vtx_moved->y*(1-dh/Rv);
	    vtx_moved->z=vtx_moved->z*(1-dh/Rv);
        //check for constraints
        if(constvolConstraintCheck(vesicle, vtx_moved)==TS_FAIL){
	        vtx_moved=memcpy((void *)vtx_moved,(void *)&backupvtx[0],sizeof(ts_vertex));
            //also, restore normals
            for(j=0;j<vtx_moved->tristar_no;j++) triangle_normal_vector(vtx_moved->tristar[j]);
//            continue;
                break;
        }

        dvol=volFirst;
        for(j=0;j<vtx_moved->tristar_no;j++){
            triangle_normal_vector(vtx_moved->tristar[j]);
            dvol+=vtx_moved->tristar[j]->volume;
        }
        voldiff=dvol-Vol;
        if(fabs(voldiff)/vesicle->volume < vesicle->tape->constvolprecision){
            //calculate energy, return change in energy...
//            fprintf(stderr, "Constvol success! %e\n",voldiff);
            for(j=0;j<vtx_moved->neigh_no;j++){
        	memcpy((void *)&backupvtx[j+1],(void *)vtx_moved->neigh[j],sizeof(ts_vertex));
    	    }

            oenergy=vtx_moved->energy;
            energy_vertex(vtx_moved);
            delta_energy=vtx_moved->xk*(vtx_moved->energy - oenergy);
            //the same is done for neighbouring vertices
            for(j=0;j<vtx_moved->neigh_no;j++){
                oenergy=vtx_moved->neigh[j]->energy;
                energy_vertex(vtx_moved->neigh[j]);
                delta_energy+=vtx_moved->neigh[j]->xk*(vtx_moved->neigh[j]->energy-oenergy);
            }
            *retEnergy=delta_energy;
            *vtx_backup=backupvtx;
            *vtx_moved_retval=vtx_moved;
            return TS_SUCCESS;
        }        
    }
    free(backupvtx);
    return TS_FAIL;
}


ts_bool constvolConstraintCheck(ts_vesicle *vesicle, ts_vertex *vtx){ 
        ts_uint i;
        ts_double dist;
        ts_uint cellidx;
    	//distance with neighbours check
        for(i=0;i<vtx->neigh_no;i++){
            dist=vtx_distance_sq(vtx,vtx->neigh[i]);
            if(dist<1.0 || dist>vesicle->dmax) {
		    return TS_FAIL;
		    }
        }
        // Distance with grafted poly-vertex check:	
	    if(vtx->grafted_poly!=NULL){
		    dist=vtx_distance_sq(vtx,vtx->grafted_poly->vlist->vtx[0]);
            if(dist<1.0 || dist>vesicle->dmax) {
		    return TS_FAIL;
		    }
	    }

        // Nucleus penetration check:
	    if (vtx->x*vtx->x + vtx->y*vtx->y + vtx->z*vtx->z < vesicle->R_nucleus){
		    return TS_FAIL;
	    }

        //self avoidance check with distant vertices
	    cellidx=vertex_self_avoidance(vesicle, vtx);
    	//check occupation number
	    return cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx);
}



ts_bool constvolumerestore(ts_vertex *vtx_moved,ts_vertex *vtx_backup){
    ts_uint j;
	 memcpy((void *)vtx_moved,(void *)&vtx_backup[0],sizeof(ts_vertex));
    for(j=0;j<vtx_moved->tristar_no;j++) triangle_normal_vector(vtx_moved->tristar[j]);
     for(j=0;j<vtx_moved->neigh_no;j++){
        	   // memcpy((void *)vtx_moved->neigh[j],(void *)&vtx_backup[j+1],sizeof(ts_vertex));
                 energy_vertex(vtx_moved->neigh[j]);
	}

    free(vtx_backup);
    return TS_SUCCESS;
}

ts_bool constvolumeaccept(ts_vesicle *vesicle,ts_vertex *vtx_moved, ts_vertex *vtx_backup){
    ts_bool retval;
	ts_uint cellidx=vertex_self_avoidance(vesicle, vtx_moved);
    if(vtx_moved->cell!=vesicle->clist->cell[cellidx]){
		retval=cell_add_vertex(vesicle->clist->cell[cellidx],vtx_moved);
		if(retval==TS_SUCCESS) cell_remove_vertex(vtx_backup[0].cell,vtx_moved);
		
	}
    free(vtx_backup);

    return TS_SUCCESS;
}

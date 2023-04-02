/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include<math.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"
#include "vesicle.h"
#include "energy.h"
#include "timestep.h"
#include "cell.h"
//#include "io.h"
#include "io.h"
#include<stdio.h>
#include "vertexmove.h"
#include <string.h>
#include "constvol.h"

ts_bool single_verticle_timestep(ts_vesicle *vesicle,ts_vertex *vtx){
    ts_small_idx i,j;
    ts_double dist;
    ts_bool retval; 
    ts_cell_idx cellidx; 
    ts_double delta_energy, oenergy,dvol=0.0, darea=0.0, dstretchenergy=0.0;
    ts_double costheta,sintheta,phi,cosphi,sinphi,r, omega, cosomega, sinomega;
    ts_double tri_angle, tri_angle_old_min, tri_angle_new_min;
    //This will hold all the information of vtx and its neighbours
    ts_vertex backupvtx[20];
    // ts_vertex* *constvol_vtx_moved=NULL, *constvol_vtx_backup=NULL;
    ts_triangle *t1, *t2;
    memcpy((void *)&backupvtx[0],(void *)vtx,sizeof(ts_vertex));

    //Some stupid tests for debugging cell occupation!
    /*
    cellidx=vertex_self_avoidance(vesicle, vtx);
    if(vesicle->clist->cell[cellidx]==vtx->cell){
        fprintf(stderr,"Idx match!\n");
    } else {
        fprintf(stderr,"***** Idx don't match!\n");
        fatal("ENding.",1);
    }
    */

    // To understand the random movement, we need to take a look at the probability space:
    // We work in spherical coordinates: a probabity integral on a sphere of radius 1 is
    //
    //  3/(4pi) integral [0,1),[0,2pi),[0,pi) f(r,phi,theta) r^2 sin(theta) dr dphi dtheta
    //
    // We do a change of coordinates to variables in [0,1)
    //   q = r^3   |   v = phi/2pi   |    u = 0.5-cos(theta)/2
    // dq = 3r^2dr | dv = 1/2pi dphi | du = sin(theta)/2 dtheta
    // leading to jacobian transformation
    // 3/(4pi) r^2 sin(theta) dr dphi dtheta = (3r^2dr) (dphi/(2pi)) (sin(theta)dtheta/2) = dq dv du
    //    => integral [0,1),[0,1),[0, 1) f(q,v,u)   dq dv du 
    //
    // the coordinate in terms of equal variables in [0,1) i.e. rand48()
    //  r = q^(1/3)  |  phi = 2pi v   |    cos(theta) = 1 - 2u
    r=vesicle->stepsize*drand48(); // should be cubed root for equal volume probability, but maybe we prefer smaller moves
    phi=drand48()*2*M_PI;
    costheta=2*drand48()-1;
    sintheta=sqrt(1-pow(costheta,2)); // gcc doesn't seems to not know trigonometry
    cosphi=cos(phi);
    sinphi= (phi<M_PI)? sqrt(1-pow(cosphi,2)) : -sqrt(1-pow(cosphi,2)) ; // 0<phi<pi: sin(phi)>0, pi<phi<2pi: sin(phi)<0,
    vtx->x=vtx->x+r*sintheta*cosphi;
    vtx->y=vtx->y+r*sintheta*sinphi;
    vtx->z=vtx->z+r*costheta;

    
    //distance with neighbours check
    for(i=0;i<vtx->neigh_no;i++){
        dist=vtx_distance_sq(vtx,vtx->neigh[i]);
        if(dist<1.0 || dist>vesicle->dmax) {
            vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
            return TS_FAIL;
        }
    }

    // angle between triangles check must be done after update in energy!

    // Distance with grafted poly-vertex check:	
    if(vtx->grafted_poly!=NULL){
        dist=vtx_distance_sq(vtx,vtx->grafted_poly->vlist->vtx[0]);
        if(dist<1.0 || dist>vesicle->dmax) {
        vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
        return TS_FAIL;
        }
    }

    // TODO: Maybe faster if checks only nucleus-neighboring cells
    // Nucleus penetration check:
    //#define SQ(x) x*x
    if(vesicle->R_nucleus>0.0){
        if (   (vtx->x-vesicle->nucleus_center[0])*(vtx->x-vesicle->nucleus_center[0])
             + (vtx->y-vesicle->nucleus_center[1])*(vtx->y-vesicle->nucleus_center[1]) 
             + (vtx->z-vesicle->nucleus_center[2])*(vtx->z-vesicle->nucleus_center[2]) < vesicle->R_nucleus){
            vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
            return TS_FAIL;
        }
    } else if(vesicle->R_nucleusX>0.0){
        //	fprintf(stderr,"DEBUG, (Rx, Ry,Rz)^2=(%f,%f,%f)\n",vesicle->R_nucleusX, vesicle->R_nucleusY, vesicle->R_nucleusZ);
        //	if (SQ(vtx->x-vesicle->nucleus_center[0])/vesicle->R_nucleusX + SQ(vtx->y-vesicle->nucleus_center[1])/vesicle->R_nucleusY + SQ(vtx->z-vesicle->nucleus_center[2])/vesicle->R_nucleusZ < 1.0){
        if (   (vtx->x-vesicle->nucleus_center[0])*(vtx->x-vesicle->nucleus_center[0])/vesicle->R_nucleusX 
             + (vtx->y-vesicle->nucleus_center[1])*(vtx->y-vesicle->nucleus_center[1])/vesicle->R_nucleusY 
             + (vtx->z-vesicle->nucleus_center[2])*(vtx->z-vesicle->nucleus_center[2])/vesicle->R_nucleusZ < 1.0){
            //	if (SQ(vtx->x)/vesicle->R_nucleusX + SQ(vtx->y)/vesicle->R_nucleusY + SQ(vtx->z)/vesicle->R_nucleusZ < 1.0){
            vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
            return TS_FAIL;
        }

    }

    // plane confinement check whether the new position of vertex will be out of bounds
    if(vesicle->tape->plane_confinement_switch){
        if(vtx->z>vesicle->confinement_plane.z_max || vtx->z<vesicle->confinement_plane.z_min){
        vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
        return TS_FAIL;
        }

    }

    // adhesion check whether the new position of vertex will be out of bounds
    // parabolic potential will push things out instead
    if(vesicle->tape->adhesion_model==adhesion_step_potential){
        if( (adhesion_geometry_distance(vesicle, vtx) < 0) 
        && (adhesion_geometry_distance(vesicle, vtx)<adhesion_geometry_distance(vesicle, backupvtx)) ){
            vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
            return TS_FAIL;
        }

    }

    //#undef SQ
    //self avoidance check with distant vertices
    cellidx=vertex_self_avoidance(vesicle, vtx);
    //check occupation number
    retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx);

    if(retval==TS_FAIL){
        vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
        return TS_FAIL;
    } 
    
    //success
    //if all the tests are successful, then energy for vtx and neighbours is calculated
    for(i=0;i<vtx->neigh_no;i++){
        memcpy((void *)&backupvtx[i+1],(void *)vtx->neigh[i],sizeof(ts_vertex));
    }


    // remove current vtx values (for future update) (vesicle->prop -= vtx->prop; update(vtx), vesicle->prop += vtx->prop)
    if(vesicle->tape->pressure_switch == 1 || vesicle->tape->volume_switch>0){
        for(i=0;i<vtx->tristar_no;i++) dvol-=vtx->tristar[i]->volume;
    }

    if(vesicle->tape->area_switch==2 || vesicle->tape->volume_switch==4 ){
        for(i=0;i<vtx->tristar_no;i++) darea-=vtx->tristar[i]->area;
    
    }
    //stretching energy 1 of 3
    if(vesicle->tape->area_switch==1){
        for(i=0;i<vtx->tristar_no;i++) dstretchenergy-=vtx->tristar[i]->energy;
    }
    delta_energy=0;


    if(vesicle->tape->min_dihedral_angle_cosine>-1){
        // vesicle_volume(vesicle);
        // fprintf(stderr,"Volume in the beginning=%1.16e\n", vesicle->volume);

        //update the normals of triangles that share bead i.
        // combined with angle check!!!
        // for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]);
        // update the normals of the vertices is in the energy
        // angle between triangles must be large to prevent spikiness, 
        // which is small angle in the normals (provided they are oriented correctly)
        // but only if the angle was not already too small (fix existing bad angles)
        //
        tri_angle_old_min=1;
        tri_angle_new_min=1;
        // min old angles
        for(i=0; i<vtx->tristar_no;i++){
            t1 = vtx->tristar[i];   
            for(j=0; j<t1->neigh_no;j++){
                // technically we are double checking some angle, but hopefully it is easier and more parallel than checking
                t2 = t1->neigh[j];
                tri_angle = t1->xnorm*t2->xnorm + t1->ynorm*t2->ynorm + t1->znorm*t2->znorm; 
                tri_angle_old_min = fmin(tri_angle_old_min, tri_angle);
            }
        }
        // update normals
        for(i=0;i<vtx->tristar_no;i++){
            t1 = vtx->tristar[i];  
            triangle_normal_vector(t1);   
        }
        // min new angles
        for(i=0;i<vtx->tristar_no;i++){
            t1 = vtx->tristar[i];   
            for(j=0; j<t1->neigh_no;j++){
                t2 = t1->neigh[j];
                tri_angle = t1->xnorm*t2->xnorm + t1->ynorm*t2->ynorm + t1->znorm*t2->znorm; 
                tri_angle_new_min = fmin(tri_angle_new_min, tri_angle);
            }
        }

        //accept or reject
        if(  (tri_angle_new_min) < vesicle->tape->min_dihedral_angle_cosine && tri_angle_new_min < tri_angle_old_min) {
                // failure! too spiky (and step is not de-spiking)
                vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
                for(i=0;i<vtx->neigh_no;i++){
                    vtx->neigh[i]=memcpy((void *)vtx->neigh[i],(void *)&backupvtx[i+1],sizeof(ts_vertex));
                }
                for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]); // we need to un-do normal updates on the triangles, since they aren't saved

                return TS_FAIL;
            }
        }
    else {
       // update normals
       for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]); 
    }

    // rotate director
    if (vtx->type & is_anisotropic_vtx && !(vesicle->tape->type_of_curvature_model&to_not_rotate_directors)){
        // we do this before the parallel transport, because the new normal calculation is fused in the energy calculation
        // we may want to bias to the original direction
        omega=(2*drand48()-1)*M_PI;
        cosomega=cos(omega);
        // omega=2*drand48()-1
        // cosomega=1-2(omega*omega)
        // cosomega=1-4(omega*omega)+2(omega*omega*omega*omega) // even closer to a cosine, fudge 1-(2+a)x^2+ax^4 as long as a<=2 to keep cosoega>=-1
        sinomega= (omega<0)? sqrt(1-pow(cosomega,2)) : -sqrt(1-pow(cosomega,2)) ; // 0<phi<pi: sin(phi)>0, -pi<phi<0: sin(phi)<0,
        // rotation: d = cos()d + sin()dxn
        vtx->dx = cosomega*vtx->dx + sinomega*(vtx->dy*vtx->nz - vtx->dz*vtx->ny);
        vtx->dy = cosomega*vtx->dy + sinomega*(vtx->dz*vtx->nx - vtx->dx*vtx->nz);
        vtx->dz = cosomega*vtx->dz + sinomega*(vtx->dx*vtx->ny - vtx->dy*vtx->nx);
    }
    // bending energy of the vertex
    oenergy=vtx->energy;
    energy_vertex(vesicle, vtx);
    delta_energy=(vtx->energy - oenergy);
    //the same is done for neighbouring vertices
    for(i=0;i<vtx->neigh_no;i++){
        oenergy=vtx->neigh[i]->energy;
        energy_vertex(vesicle, vtx->neigh[i]);
        delta_energy+=(vtx->neigh[i]->energy-oenergy);
    }

    // bonding energy can change due to director updates
    for(i=0;i<vtx->neigh_no;i++){
        for(j=0;j<vtx->neigh[i]->bond_no;j++){
            if(vtx->neigh[i]->bond[j]->vtx1==vtx->neigh[prev_small(i,vtx->neigh_no)] || vtx->neigh[i]->bond[j]->vtx2==vtx->neigh[prev_small(i,vtx->neigh_no)] ){
                //skip bond with previous neighbor (which was calculated for it)
            }
            else{
                oenergy=vtx->neigh[i]->bond[j]->energy;
                attraction_bond_energy(vesicle, vtx->neigh[i]->bond[j]);
                delta_energy+=(vtx->neigh[i]->bond[j]->energy-oenergy);
            }
        }
    }


    // part 1 of 2 of volume and area update
    if(vesicle->tape->pressure_switch == 1 || vesicle->tape->volume_switch >0){
        for(i=0;i<vtx->tristar_no;i++) dvol+=vtx->tristar[i]->volume;
    }
    if(vesicle->tape->area_switch==2 || vesicle->tape->volume_switch == 4){
        for(i=0;i<vtx->tristar_no;i++) darea+=vtx->tristar[i]->area;
    }


    // volume area pressure energy and constraints
    retval = volume_pressure_area_energy_constraints(vesicle,&delta_energy,dvol,darea);
    if (retval == TS_FAIL){
        //restore old state.
        vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
        for(i=0;i<vtx->neigh_no;i++){
            vtx->neigh[i]=memcpy((void *)vtx->neigh[i],(void *)&backupvtx[i+1],sizeof(ts_vertex));
        }
        // unupdate triangles
        for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]);
        //unupdate bonds
        for(i=0;i<vtx->neigh_no;i++){
            for(j=0;j<vtx->neigh[i]->bond_no;j++){
                // no problem double updating, just double counting the energy
                attraction_bond_energy(vesicle, vtx->neigh[i]->bond[j]);
            }
        }
        return TS_FAIL;        
    }

    
    // vertices may be active, and have force applied on them
    // to do: think hard on stratonovich/ito approach
    if (vtx->type&is_active_vtx) delta_energy+=direct_force_energy(vesicle,vtx,backupvtx);
    
    // additionally, vesicle may have force balance, applying additional force
    // this is done separately here (since we can split the F)
    if (vesicle->tape->force_balance_along_z_axis==1){
        delta_energy+=direct_force_from_Fz_balance(vesicle,vtx,backupvtx);
    }
    

    //stretching energy 2 of 3
    if(vesicle->tape->area_switch==1){
        for(i=0;i<vtx->tristar_no;i++){ 
            stretchenergy(vesicle, vtx->tristar[i]);
            dstretchenergy+=vtx->tristar[i]->energy;
            }
    }

    delta_energy+=dstretchenergy;	
        
    /* No poly-bond energy for now!
    if(vtx->grafted_poly!=NULL){
        delta_energy+=
            (pow(sqrt(vtx_distance_sq(vtx, vtx->grafted_poly->vlist->vtx[0])-1),2)-
            pow(sqrt(vtx_distance_sq(&backupvtx[0], vtx->grafted_poly->vlist->vtx[0])-1),2)) *vtx->grafted_poly->k;
    }
    */

    // plane confinement energy due to compressing force
    if(vesicle->tape->plane_confinement_switch){
        if(vesicle->confinement_plane.force_switch){
            //substract old energy
            if(abs(vesicle->tape->plane_d-vesicle->confinement_plane.z_max)>1e-10) {
                delta_energy-=vesicle->tape->plane_F / pow(vesicle->confinement_plane.z_max-backupvtx[0].z,2);
                delta_energy+=vesicle->tape->plane_F / pow(vesicle->confinement_plane.z_max-vtx->z,2);
            }
        }
    }
    // change in energy due to adhesion
    delta_energy+=adhesion_energy_diff(vesicle, vtx, backupvtx);

    // fprintf(stderr, "DE=%f\n",delta_energy);
    //MONTE CARLOOOOOOOO
    // if(vtx->c!=0.0) printf("DE=%f\n",delta_energy);
    if(delta_energy>=0){
        if(exp(-delta_energy)< drand48()) { 
            //not accepted, reverting changes
            // fprintf(stderr,"MC failed\n");
            // Fz balance
            if (vesicle->tape->force_balance_along_z_axis==1) direct_force_from_Fz_balance(vesicle,backupvtx,vtx);

            //unupdate vertices
            vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
            for(i=0;i<vtx->neigh_no;i++){
                vtx->neigh[i]=memcpy((void *)vtx->neigh[i],(void *)&backupvtx[i+1],sizeof(ts_vertex));
            }
    
            //unupdate the normals of triangles that share bead i.
            for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]);
            //unupdate bond energy
            for(i=0;i<vtx->neigh_no;i++){
                for(j=0;j<vtx->neigh[i]->bond_no;j++){
                    // no problem double updating, just double counting the energy
                    attraction_bond_energy(vesicle, vtx->neigh[i]->bond[j]);
                }
            }

            //stretching energy 3 of 3
            if(vesicle->tape->area_switch==1){
                for(i=0;i<vtx->tristar_no;i++){ 
                    stretchenergy(vesicle,vtx->tristar[i]);
                }
            }

            // fprintf(stderr, "before vtx(x,y,z)=%e,%e,%e\n",constvol_vtx_moved->x, constvol_vtx_moved->y, constvol_vtx_moved->z);
            // if(vesicle->tape->volume_switch == 1){
            //     constvolumerestore(vesicle, constvol_vtx_moved,constvol_vtx_backup);
            // }


            return TS_FAIL; 
        }
    }
    //accepted	
    // fprintf(stderr,"MC accepted\n");
    // oldcellidx=vertex_self_avoidance(vesicle, &backupvtx[0]);
    if(vtx->cell!=vesicle->clist->cell[cellidx]){
        retval=cell_add_vertex(vesicle->clist->cell[cellidx],vtx);
        // if(retval==TS_SUCCESS) cell_remove_vertex(vesicle->clist->cell[oldcellidx],vtx);
        if(retval==TS_SUCCESS) cell_remove_vertex(backupvtx[0].cell,vtx);
        
    }


    // else if(vesicle->tape->volume_switch == 1){
    //     constvolumeaccept(vesicle,constvol_vtx_moved,constvol_vtx_backup);
    // }

    // part 2 of 2 of volume and area update
    if(vesicle->tape->pressure_switch == 1 || vesicle->tape->volume_switch > 0){
        vesicle->volume+=dvol;
    } 

    if(vesicle->tape->area_switch == 2 || vesicle->tape->volume_switch == 4){
        vesicle->area+=darea;
    }


    //END MONTE CARLOOOOOOO


    return TS_SUCCESS;
}


ts_bool single_poly_vertex_move(ts_vesicle *vesicle,ts_poly *poly,ts_vertex *vtx,ts_double *rn){
    ts_idx i;
    ts_bool retval; 
    ts_cell_idx cellidx; 
    // ts_double delta_energy;
    ts_double costheta,sintheta,phi,r;
    ts_double dist;
    //This will hold all the information of vtx and its neighbours
    ts_vertex backupvtx;
    // ts_bond backupbond[2];
    memcpy((void *)&backupvtx,(void *)vtx,sizeof(ts_vertex));

    //random move in a sphere with radius stepsize:
    r=vesicle->stepsize*rn[0];
    phi=rn[1]*2*M_PI;
    costheta=2*rn[2]-1;
    sintheta=sqrt(1-pow(costheta,2));
    vtx->x=vtx->x+r*sintheta*cos(phi);
    vtx->y=vtx->y+r*sintheta*sin(phi);
    vtx->z=vtx->z+r*costheta;


    //distance with neighbours check
    for(i=0;i<vtx->neigh_no;i++){
        dist=vtx_distance_sq(vtx,vtx->neigh[i]);
        if(dist<1.0 || dist>vesicle->dmax) {
            vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
            return TS_FAIL;
        }
    }

// Distance with grafted vesicle-vertex check:	
    if(vtx==poly->vlist->vtx[0]){
        dist=vtx_distance_sq(vtx,poly->grafted_vtx);
        if(dist<1.0 || dist>vesicle->dmax) {
        vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
        return TS_FAIL;
        }
    }


    //self avoidance check with distant vertices
    cellidx=vertex_self_avoidance(vesicle, vtx);
    //check occupation number
    retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx);
    
    if(retval==TS_FAIL){
        vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
        return TS_FAIL;
    } 


    //if all the tests are successful, then energy for vtx and neighbours is calculated
/* Energy ignored for now!
    delta_energy=0;
    for(i=0;i<vtx->bond_no;i++){
        memcpy((void *)&backupbond[i],(void *)vtx->bond[i],sizeof(ts_bond));

        vtx->bond[i]->bond_length=sqrt(vtx_distance_sq(vtx->bond[i]->vtx1,vtx->bond[i]->vtx2));
        bond_energy(vtx->bond[i],poly);
        delta_energy+= vtx->bond[i]->energy - backupbond[i].energy;
    }

    if(vtx==poly->vlist->vtx[0]){
        delta_energy+=
            (pow(sqrt(vtx_distance_sq(vtx, poly->grafted_vtx)-1),2)-
            pow(sqrt(vtx_distance_sq(&backupvtx, poly->grafted_vtx)-1),2)) *poly->k;
        
    }


    if(delta_energy>=0){
#ifdef TS_DOUBLE_DOUBLE
        if(exp(-delta_energy)< drand48() )
#endif
#ifdef TS_DOUBLE_FLOAT
        if(expf(-delta_energy)< (ts_float)drand48())
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        if(expl(-delta_energy)< (ts_ldouble)drand48())
#endif
        {
    //not accepted, reverting changes
    vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
    for(i=0;i<vtx->bond_no;i++){
    vtx->bond[i]=memcpy((void *)vtx->bond[i],(void *)&backupbond[i],sizeof(ts_bond));
    }

    return TS_FAIL; 
    }
    }
*/
        
//	oldcellidx=vertex_self_avoidance(vesicle, &backupvtx[0]);
    if(vtx->cell!=vesicle->clist->cell[cellidx]){
        retval=cell_add_vertex(vesicle->clist->cell[cellidx],vtx);
//		if(retval==TS_SUCCESS) cell_remove_vertex(vesicle->clist->cell[oldcellidx],vtx);
        if(retval==TS_SUCCESS) cell_remove_vertex(backupvtx.cell,vtx);	
    }
//	if(oldcellidx);
    //END MONTE CARLOOOOOOO
    return TS_SUCCESS;
}




ts_bool single_filament_vertex_move(ts_vesicle *vesicle,ts_poly *poly,ts_vertex *vtx,ts_double *rn){
    ts_idx i;
    ts_bool retval; 
    ts_cell_idx cellidx; 
    ts_double delta_energy;
    ts_double costheta,sintheta,phi,r;
    ts_double dist[2];
    //This will hold all the information of vtx and its neighbours
    ts_vertex backupvtx,backupneigh[2];
    ts_bond backupbond[2];

    //backup vertex:		
    memcpy((void *)&backupvtx,(void *)vtx,sizeof(ts_vertex));

    //random move in a sphere with radius stepsize:
    r=vesicle->stepsize*rn[0];
    phi=rn[1]*2*M_PI;
    costheta=2*rn[2]-1;
    sintheta=sqrt(1-pow(costheta,2));
    vtx->x=vtx->x+r*sintheta*cos(phi);
    vtx->y=vtx->y+r*sintheta*sin(phi);
    vtx->z=vtx->z+r*costheta;


    //distance with neighbours check
    for(i=0;i<vtx->bond_no;i++){
        dist[i]=vtx_distance_sq(vtx->bond[i]->vtx1,vtx->bond[i]->vtx2);
        if(dist[i]<1.0 || dist[i]>vesicle->dmax) {
            vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
            return TS_FAIL;
        }
    }

// TODO: Maybe faster if checks only nucleus-neighboring cells
// Nucleus penetration check:
    if (vtx->x*vtx->x + vtx->y*vtx->y + vtx->z*vtx->z < vesicle->R_nucleus){
        vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
        return TS_FAIL;
    }


    //self avoidance check with distant vertices
    cellidx=vertex_self_avoidance(vesicle, vtx);
    //check occupation number
    retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx);
    if(retval==TS_FAIL){
        vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
        return TS_FAIL;
    } 

    //backup bonds
    for(i=0;i<vtx->bond_no;i++){
        memcpy(&backupbond[i],vtx->bond[i], sizeof(ts_bond));
        vtx->bond[i]->bond_length=sqrt(dist[i]);
        bond_vector(vtx->bond[i]);
    }

    //backup neighboring vertices:
    for(i=0;i<vtx->neigh_no;i++){
        memcpy(&backupneigh[i],vtx->neigh[i], sizeof(ts_vertex));
    }
    
    //if all the tests are successful, then energy for vtx and neighbours is calculated
    delta_energy=0;
    
    if(vtx->bond_no == 2){
        vtx->energy = -(vtx->bond[0]->x*vtx->bond[1]->x + vtx->bond[0]->y*vtx->bond[1]->y + vtx->bond[0]->z*vtx->bond[1]->z)/vtx->bond[0]->bond_length/vtx->bond[1]->bond_length;
        delta_energy += vtx->energy - backupvtx.energy;
    }

    for(i=0;i<vtx->neigh_no;i++){
        if(vtx->neigh[i]->bond_no == 2){
            vtx->neigh[i]->energy = -(vtx->neigh[i]->bond[0]->x*vtx->neigh[i]->bond[1]->x + vtx->neigh[i]->bond[0]->y*vtx->neigh[i]->bond[1]->y + vtx->neigh[i]->bond[0]->z*vtx->neigh[i]->bond[1]->z)/vtx->neigh[i]->bond[0]->bond_length/vtx->neigh[i]->bond[1]->bond_length;
            delta_energy += vtx->neigh[i]->energy - backupneigh[i].energy;
        }
    }

    // poly->k is filament persistence length (in units l_min)
    delta_energy *= poly->k;

    if(delta_energy>=0){
#ifdef TS_DOUBLE_DOUBLE
        if(exp(-delta_energy)< drand48() )
#endif
#ifdef TS_DOUBLE_FLOAT
        if(expf(-delta_energy)< (ts_float)drand48())
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        if(expl(-delta_energy)< (ts_ldouble)drand48())
#endif
        {
    //not accepted, reverting changes
    vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
    for(i=0;i<vtx->neigh_no;i++){
        memcpy(vtx->neigh[i],&backupneigh[i],sizeof(ts_vertex));
    }
    for(i=0;i<vtx->bond_no;i++){
        vtx->bond[i]=memcpy((void *)vtx->bond[i],(void *)&backupbond[i],sizeof(ts_bond));
    }

    return TS_FAIL; 
    }
    }
    
    
//	oldcellidx=vertex_self_avoidance(vesicle, &backupvtx[0]);
    if(vtx->cell!=vesicle->clist->cell[cellidx]){
        retval=cell_add_vertex(vesicle->clist->cell[cellidx],vtx);
//		if(retval==TS_SUCCESS) cell_remove_vertex(vesicle->clist->cell[oldcellidx],vtx);
        if(retval==TS_SUCCESS) cell_remove_vertex(backupvtx.cell,vtx);	
    }
//	if(oldcellidx);
    //END MONTE CARLOOOOOOO
    return TS_SUCCESS;
}

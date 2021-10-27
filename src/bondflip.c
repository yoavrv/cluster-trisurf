/* vim: set ts=4 sts=4 sw=4 noet : */
#include <stdlib.h>
#include <math.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"
#include "vesicle.h"
#include "energy.h"
#include "timestep.h"
#include "cell.h"
#include "bondflip.h"
//#include "io.h"
#include<stdio.h>
#include<string.h>
#include "constvol.h"

ts_bool single_bondflip_timestep(ts_vesicle *vesicle, ts_bond *bond, ts_double *rn){
/*c  Vertex and triangle (lm and lp) indexing for bond flip:
c      +----- k-------+              +----- k ------+
c      |lm1 / | \ lp1 |              |lm1 /   \ lp1 |
c      |  /   |   \   |              |  /       \   |
c      |/     |     \ |     FLIP     |/    lm     \ |
c     km  lm  | lp   kp    --->      km ---------- kp  
c      |\     |     / |              |\    lp     / |  
c      |  \   |   /   |              |  \       /   |
c      |lm2 \ | / lp2 |              |lm2 \   / lp2 |
c      +------it------+              +----- it -----+
c
*/
    ts_vertex *it=bond->vtx1;
    ts_vertex *k=bond->vtx2;
    ts_small_idx nei,neip,neim;
    ts_small_idx i,j;
    ts_double oldenergy, delta_energy, dvol=0.0, darea=0.0;
    ts_triangle *lm=NULL,*lp=NULL, *lp1=NULL, *lm2=NULL;

    ts_vertex *kp,*km;

    ts_double delta_energy_cv;
    ts_vertex *constvol_vtx_moved, *constvol_vtx_backup;
    ts_bool retval;

    if (it->type==is_ghost_vtx && k->type==is_ghost_vtx) return TS_FAIL;

    if(it->neigh_no< 3) return TS_FAIL;
    if(k->neigh_no< 3) return TS_FAIL;
    if(k==NULL || it==NULL){
        fatal("In bondflip, number of neighbours of k or it is less than 3!",999);
    }

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
    if(nei==0) neim=it->neigh_no-1; /* casting is essential... If not there the neim is never <0 !!! */
    //  fprintf(stderr,"The numbers are: %u %u\n",neip, neim);
    km=it->neigh[neim];  // We located km and kp
    kp=it->neigh[neip];

    if(km==NULL || kp==NULL){
        fatal("In bondflip, cannot determine km and kp!",999);
    }

    //  fprintf(stderr,"I WAS HERE! after the 4 vertices are known!\n");

    /* test if the membrane is wrapped too much, so that kp is nearest neighbour of
    * km. If it is true, then don't flip! */

    for(i=0;i<km->neigh_no;i++){
        if(km->neigh[i] == kp) return TS_FAIL;
    }
    
    //   fprintf(stderr,"Membrane didn't wrap too much.. Continue.\n");
    /* if bond would be too long, return... */
    if(vtx_distance_sq(km,kp) > vesicle->dmax ) return TS_FAIL;
    //   fprintf(stderr,"Bond will not be too long.. Continue.\n");

    
    /* we make a bond flip. this is different than in original fortran */
    // find lm, lp
    // 1. step. We find lm and lp from k->tristar !
    for(i=0;i<it->tristar_no;i++){
        for(j=0;j<k->tristar_no;j++){
            if((it->tristar[i] == k->tristar[j])){ //ce gre za skupen trikotnik
                if((it->tristar[i]->vertex[0] == km || it->tristar[i]->vertex[1]== km || it->tristar[i]->vertex[2]== km )){
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
    if(lm==NULL || lp==NULL) fatal("ts_flip_bond: Cannot find triangles lm and lp!",999);

    //we look for important triangles lp1 and lm2.

    for(i=0;i<k->tristar_no;i++){
        for(j=0;j<kp->tristar_no;j++){
                if((k->tristar[i] == kp->tristar[j]) && k->tristar[i]!=lp){ // אם זה משולש משותף
                    lp1=k->tristar[i];
            }
        }
    }

    for(i=0;i<it->tristar_no;i++){
        for(j=0;j<km->tristar_no;j++){
            if((it->tristar[i] == km->tristar[j]) && it->tristar[i]!=lm){ //ce gre za skupen trikotnik
                    lm2=it->tristar[i];
            } 
        }
    }

    if(lm2==NULL || lp1==NULL) fatal("ts_flip_bond: Cannot find triangles lm2 and lp1!",999);


    // ####### Passed all base tests #######

    /* backup old structure */
    /* need to backup:
    * vertices k, kp, km, it
    * triangles lm, lp, lm2, lp1
    * bond
    */
    ts_vertex *bck_vtx[4];
    ts_triangle *bck_tria[4];
    ts_bond *bck_bond;
    ts_vertex *orig_vtx[]={k,it,kp,km};
    ts_triangle *orig_tria[]={lm,lp,lm2,lp1};

    //fprintf(stderr,"Backuping!!!\n");
	bck_bond=(ts_bond *)malloc(sizeof(ts_bond));
    for(i=0;i<4;i++){
    /*	fprintf(stderr,"vtx neigh[%d]=",i);
	    for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
	    fprintf(stderr,"\n");
    */
	    bck_vtx[i]=(ts_vertex *)malloc(sizeof(ts_vertex));
	    bck_tria[i]=(ts_triangle *)malloc(sizeof(ts_triangle));
    	memcpy((void *)bck_vtx[i],(void *)orig_vtx[i],sizeof(ts_vertex));
    	memcpy((void *)bck_tria[i],(void *)orig_tria[i],sizeof(ts_triangle));
	    /* level 2 pointers */

    	bck_vtx[i]->neigh=(ts_vertex **)malloc(orig_vtx[i]->neigh_no*sizeof(ts_vertex *));
    	bck_vtx[i]->tristar=(ts_triangle **)malloc(orig_vtx[i]->tristar_no*sizeof(ts_triangle *));
    	bck_vtx[i]->bond=(ts_bond **)malloc(orig_vtx[i]->bond_no*sizeof(ts_bond *));
    	bck_tria[i]->neigh=(ts_triangle **)malloc(orig_tria[i]->neigh_no*sizeof(ts_triangle *));

    	memcpy((void *)bck_vtx[i]->neigh,(void *)orig_vtx[i]->neigh,orig_vtx[i]->neigh_no*sizeof(ts_vertex *));
    	memcpy((void *)bck_vtx[i]->tristar,(void *)orig_vtx[i]->tristar,orig_vtx[i]->tristar_no*sizeof(ts_triangle *));
    	memcpy((void *)bck_vtx[i]->bond,(void *)orig_vtx[i]->bond,orig_vtx[i]->bond_no*sizeof(ts_bond *));
    
    	memcpy((void *)bck_tria[i]->neigh,(void *)orig_tria[i]->neigh,orig_tria[i]->neigh_no*sizeof(ts_triangle *));	
    }
	memcpy(bck_bond,bond,sizeof(ts_bond));
    //fprintf(stderr,"Backup complete!!!\n");
    /* end backup vertex */


    /* Save old energy */
    oldenergy=0;
    oldenergy+=k->energy;
    oldenergy+=kp->energy;
    oldenergy+=km->energy;
    oldenergy+=it->energy;
    oldenergy+=bond->energy; /* attraction with neighboring vertices, that have spontaneous curvature */
    //Neigbours of k, it, km, kp don't change its energy.

	if(vesicle->pswitch == 1 || vesicle->tape->constvolswitch>0){dvol = -lm->volume - lp->volume;}
    if(vesicle->tape->constareaswitch==2){darea=-lm->area-lp->area;} 
    /*    vesicle_volume(vesicle);
    fprintf(stderr,"Volume in the beginning=%1.16e\n", vesicle->volume);
    */
    /* fix data structure for flipped bond */
    ts_flip_bond(vesicle, k,it,km,kp, bond,lm, lp, lm2, lp1);


    /* Calculating the new energy */
    delta_energy=0;
    delta_energy+=k->energy;
    delta_energy+=kp->energy;
    delta_energy+=km->energy;
    delta_energy+=it->energy;
    delta_energy+=bond->energy; /* attraction with neighboring vertices, that have spontaneous curvature */
  //Neigbours of k, it, km, kp don't change its energy.
	if(vesicle->tape->stretchswitch==1){
		oldenergy+=lm->energy+lp->energy;
		stretchenergy(vesicle,lm);
		stretchenergy(vesicle,lp);
		delta_energy+=lm->energy+lp->energy;
	}

    delta_energy-=oldenergy;
	if(vesicle->pswitch == 1 || vesicle->tape->constvolswitch>0){
		dvol = dvol + lm->volume + lp->volume;
		if(vesicle->pswitch==1) delta_energy-= vesicle->pressure*dvol;
	}

    
 //   if(k->neigh_no<4 || kp->neigh_no<4 || km->neigh_no<4 || it->neigh_no<4){
//
//			//restore old state.
//			/* restoration procedure copied from few lines below */
//			    for(i=0;i<4;i++){
//			//			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
//				free(orig_vtx[i]->neigh);
//				free(orig_vtx[i]->tristar);
//				free(orig_vtx[i]->bond);
//				free(orig_tria[i]->neigh);
//				memcpy((void *)orig_vtx[i],(void *)bck_vtx[i],sizeof(ts_vertex));
//				memcpy((void *)orig_tria[i],(void *)bck_tria[i],sizeof(ts_triangle));
//			//			fprintf(stderr,"Restored vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
//				/* level 2 pointers are redirected*/
//			    }
//			    memcpy(bond,bck_bond,sizeof(ts_bond));
//			    for(i=0;i<4;i++){
//				free(bck_vtx[i]);
//				free(bck_tria[i]);
//			/*			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d =",i, orig_vtx[i]->neigh_no );
//				for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
//				fprintf(stderr,"\n"); */
//			    }
//			    free(bck_bond);
//			    return TS_FAIL;
//
//		
 //   }


    if(vesicle->tape->constareaswitch==2){
        darea=darea+lm->area+lp->area; 
/*check whether the dvol is gt than epsvol */
		if(fabs(vesicle->area+darea-A0)>epsarea){
			//restore old state.
			/* restoration procedure copied from few lines below */
			    for(i=0;i<4;i++){
			//			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
				free(orig_vtx[i]->neigh);
				free(orig_vtx[i]->tristar);
				free(orig_vtx[i]->bond);
				free(orig_tria[i]->neigh);
				memcpy((void *)orig_vtx[i],(void *)bck_vtx[i],sizeof(ts_vertex));
				memcpy((void *)orig_tria[i],(void *)bck_tria[i],sizeof(ts_triangle));
			//			fprintf(stderr,"Restored vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
				/* level 2 pointers are redirected*/
			    }
			    memcpy(bond,bck_bond,sizeof(ts_bond));
			    for(i=0;i<4;i++){
				free(bck_vtx[i]);
				free(bck_tria[i]);
			/*			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d =",i, orig_vtx[i]->neigh_no );
				for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
				fprintf(stderr,"\n"); */
			    }
			    free(bck_bond);
			    return TS_FAIL;

		}
    }




	if(vesicle->tape->constvolswitch == 2){
		/*check whether the dvol is gt than epsvol */
		if(fabs(vesicle->volume+dvol-V0)>epsvol){
			//restore old state.
			/* restoration procedure copied from few lines below */
			    for(i=0;i<4;i++){
			//			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
				free(orig_vtx[i]->neigh);
				free(orig_vtx[i]->tristar);
				free(orig_vtx[i]->bond);
				free(orig_tria[i]->neigh);
				memcpy((void *)orig_vtx[i],(void *)bck_vtx[i],sizeof(ts_vertex));
				memcpy((void *)orig_tria[i],(void *)bck_tria[i],sizeof(ts_triangle));
			//			fprintf(stderr,"Restored vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
				/* level 2 pointers are redirected*/
			    }
			    memcpy(bond,bck_bond,sizeof(ts_bond));
			    for(i=0;i<4;i++){
				free(bck_vtx[i]);
				free(bck_tria[i]);
			/*			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d =",i, orig_vtx[i]->neigh_no );
				for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
				fprintf(stderr,"\n"); */
			    }
			    free(bck_bond);
			    return TS_FAIL;

		}

	} else
    if(vesicle->tape->constvolswitch == 1){
        retval=constvolume(vesicle, it, -dvol, &delta_energy_cv, &constvol_vtx_moved,&constvol_vtx_backup);
        if(retval==TS_FAIL){
/* restoration procedure copied from few lines below */
            for(i=0;i<4;i++){
    //			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
                free(orig_vtx[i]->neigh);
                free(orig_vtx[i]->tristar);
                free(orig_vtx[i]->bond);
                free(orig_tria[i]->neigh);
                memcpy((void *)orig_vtx[i],(void *)bck_vtx[i],sizeof(ts_vertex));
                memcpy((void *)orig_tria[i],(void *)bck_tria[i],sizeof(ts_triangle));
    //			fprintf(stderr,"Restored vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
                /* level 2 pointers are redirected*/
            }
            memcpy(bond,bck_bond,sizeof(ts_bond));
            for(i=0;i<4;i++){
                free(bck_vtx[i]);
                free(bck_tria[i]);
    /*			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d =",i, orig_vtx[i]->neigh_no );
                for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
                fprintf(stderr,"\n"); */
            }
            free(bck_bond);
            return TS_FAIL;
        }
        delta_energy+=delta_energy_cv;
    }


/* MONTE CARLO */
    if(delta_energy>=0){
#ifdef TS_DOUBLE_DOUBLE
        if(exp(-delta_energy)< drand48())
#endif
#ifdef TS_DOUBLE_FLOAT
        if(expf(-delta_energy)< (ts_float)drand48())
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        if(expl(-delta_energy)< (ts_ldouble)drand48())
#endif
        {
            //not accepted, reverting changes
	    //restore all backups
//		fprintf(stderr,"Restoring!!!\n");
        if(vesicle->tape->constvolswitch == 1){
            constvolumerestore(vesicle, constvol_vtx_moved,constvol_vtx_backup);
        }

		for(i=0;i<4;i++){
//			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
			free(orig_vtx[i]->neigh);
			free(orig_vtx[i]->tristar);
			free(orig_vtx[i]->bond);
			free(orig_tria[i]->neigh);
			memcpy((void *)orig_vtx[i],(void *)bck_vtx[i],sizeof(ts_vertex));
			memcpy((void *)orig_tria[i],(void *)bck_tria[i],sizeof(ts_triangle));
//			fprintf(stderr,"Restored vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
			/* level 2 pointers are redirected*/
		}
		memcpy(bond,bck_bond,sizeof(ts_bond));

		for(i=0;i<4;i++){
			free(bck_vtx[i]);
			free(bck_tria[i]);
/*			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d =",i, orig_vtx[i]->neigh_no );
			for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
			fprintf(stderr,"\n"); */
		}

		free(bck_bond);

//		fprintf(stderr,"Restoration complete!!!\n");
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after fail=%1.16e\n", vesicle->volume);
	if(vesicle->tape->stretchswitch==1){
		stretchenergy(vesicle,lm);
		stretchenergy(vesicle,lp);
	}

		return TS_FAIL;
        }
    }
     /* IF BONDFLIP ACCEPTED, THEN RETURN SUCCESS! */
//            fprintf(stderr,"SUCCESS!!!\n");

    if(vesicle->tape->constvolswitch == 2){
	    vesicle->volume+=dvol;
    } else if(vesicle->tape->constvolswitch == 1){
        constvolumeaccept(vesicle,constvol_vtx_moved,constvol_vtx_backup);
    }
    if(vesicle->tape->constareaswitch==2){
        vesicle->area+=darea;
    }
	// delete all backups
	for(i=0;i<4;i++){
	    free(bck_vtx[i]->neigh);
	    free(bck_vtx[i]->bond);
	    free(bck_vtx[i]->tristar);
	    free(bck_vtx[i]);
 	    free(bck_tria[i]->neigh);
        free(bck_tria[i]);
	}
	free(bck_bond);

//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after success=%1.16e\n", vesicle->volume);
    return TS_SUCCESS;
}


ts_bool ts_flip_bond(ts_vesicle *vesicle, ts_vertex *k,ts_vertex *it,ts_vertex *km, ts_vertex *kp, ts_bond *bond,
                     ts_triangle *lm, ts_triangle *lp, ts_triangle *lm2, ts_triangle *lp1){

    ts_small_idx i; //lmidx, lpidx;
    if(k==NULL || it==NULL || km==NULL || kp==NULL){
        fatal("ts_flip_bond: You called me with invalid pointers to vertices",999);
    }
    
    // 2. step. We change the triangle vertices... (actual bond flip)
    for(i=0;i<3;i++) if(lm->vertex[i]== it) lm->vertex[i]= kp;
    for(i=0;i<3;i++) if(lp->vertex[i]== k) lp->vertex[i]= km;

    // 2a. step. If any changes in triangle calculations must be done, do it here!
    //   * normals are recalculated here
    triangle_normal_vector(lp);
    triangle_normal_vector(lm);

    // 3. step. Correct neighbours in vertex_list

    vtx_remove_neighbour(k,it);

    
    //Tukaj pa nastopi tezava... Kam dodati soseda?
    vtx_insert_neighbour(km,kp,k);
    vtx_insert_neighbour(kp,km,it);
    //pazi na vrstni red.


    // 3a. step. If any changes to ts_vertex, do it here!
    //   bond_length calculatons not required for it is done in energy.c

    // 4. step. Correct bond_list (don't know why I still have it!)
    bond->vtx1=km;
    bond->vtx2=kp;
    //fprintf(stderr,"4. step: bondlist corrected\n");


    // 5. step. Correct neighbouring triangles 
   
    triangle_remove_neighbour(lp,lp1);
    triangle_remove_neighbour(lp1,lp);
    triangle_remove_neighbour(lm,lm2);
    triangle_remove_neighbour(lm2,lm);
   
    triangle_add_neighbour(lm,lp1);    
    triangle_add_neighbour(lp1,lm);
    triangle_add_neighbour(lp,lm2);  //Vrstni red?!
    triangle_add_neighbour(lm2,lp);



    // 6. step. Correct tristar for vertices km, kp, k and it
    vertex_add_tristar(km,lp);  // Preveri vrstni red!
    vertex_add_tristar(kp,lm);
    vtx_remove_tristar(it,lm);
    vtx_remove_tristar(k,lp);

    // END modifications to data structure!



    // 7. step. Update energy
    energy_vertex(vesicle, k);
    energy_vertex(vesicle, kp);
    energy_vertex(vesicle, km);
    energy_vertex(vesicle, it);
    // Yoav: this also updates the normals

    // Yoav: Do we need to update force?

    // Yoav; No more sense in giving w to the energy: it's in the vertices now
    attraction_bond_energy(vesicle, bond);
    return TS_SUCCESS;
}



// helper functions for ordered versions

// find index i vtx->neigh[i]==nei. if nei is not in btx->neigh, return neigh_no
ts_small_idx find_neigh_idx(ts_vertex* vtx, ts_vertex* nei){
    ts_small_idx i;
    for(i=0;i<vtx->neigh_no;i++){
        if(vtx->neigh[i]==nei){
            return i;
            break;
        }
    }
    return vtx->neigh_no; // return neigh_no, a faulty index that should give error (rather than valid but wrong neigh_no-1)
}

// swap a triangle's vertex without reallocating
ts_bool swap_tri_vertex(ts_triangle* tri, ts_vertex* vtx_remove, ts_vertex* vtx_insert){
    ts_small_idx i;
    for(i=0;i<3;i++) {
        if(tri->vertex[i] == vtx_remove){
            tri->vertex[i] = vtx_insert;
            return TS_SUCCESS;
            break;
        }
    }
    return TS_FAIL;
}

// swap a triangle's neighbor without reallocating
ts_bool swap_tri_neigh(ts_triangle* tri, ts_triangle* tri_remove, ts_triangle* tri_insert){
    ts_small_idx i;
    for(i=0;i<3;i++) {
        if(tri->neigh[i] == tri_remove){
            tri->neigh[i] = tri_insert;
            return TS_SUCCESS;
            break;
        }
    }
    return TS_FAIL;
}

// same as swap_tri_neigh but reverse the neighbors cyclic order
// {a,b,c}->{b,d,c}, {a,b,c}->{a,c,d}, {a,b,c}->{d,b,a}
ts_bool swap_tri_neigh_and_reverse(ts_triangle* tri, ts_triangle* tri_remove, ts_triangle* tri_insert){
    if (tri->neigh[0] == tri_remove){
        tri->neigh[0] = tri->neigh[1];
        tri->neigh[1] = tri_insert;
        return TS_SUCCESS;
    }
    else if (tri->neigh[1] == tri_remove){
        tri->neigh[1] = tri->neigh[2];
        tri->neigh[2] = tri_insert;
        return TS_SUCCESS;
    }
    else if (tri->neigh[2] == tri_remove){
        tri->neigh[2] = tri->neigh[0];
        tri->neigh[0] = tri_insert;
        return TS_SUCCESS;
    }
    return TS_FAIL;
}

ts_bool single_bondflip_timestep_ordered(ts_vesicle *vesicle, ts_bond *bond, ts_double *rn){
/*c  Vertex and triangle (lm and lp) indexing for bond flip:
c      +----- k-------+              +----- k ------+
c      |lm1 / | \ lp1 |              |lm1 /   \ lp1 |
c      |  /   |   \   |              |  /       \   |
c      |/     |     \ |     FLIP     |/    lm     \ |
c     km  lm  | lp   kp    --->      km ---------- kp  
c      |\     |     / |              |\    lp     / |  
c      |  \   |   /   |              |  \       /   |
c      |lm2 \ | / lp2 |              |lm2 \   / lp2 |
c      +------it------+              +----- it -----+
c
*/
// same as bondflip but ensure vertices remained ordered in their tristar too
    ts_vertex *it=bond->vtx1;
    ts_vertex *k=bond->vtx2;
    ts_small_idx nei,neip,neim;
    ts_small_idx nei_k_at_km, nei_kp_at_k, nei_it_at_kp; // nei_k_at_it == neim
    ts_small_idx i;
    ts_double oldenergy, delta_energy, dvol=0.0, darea=0.0;
    ts_triangle *lm=NULL,*lp=NULL, *lp1=NULL, *lm2=NULL;

    ts_vertex *kp,*km;

    ts_double delta_energy_cv;
    ts_vertex *constvol_vtx_moved, *constvol_vtx_backup;
    ts_bool retval;

    if (it->type==is_ghost_vtx && k->type==is_ghost_vtx) return TS_FAIL;

    if(it->neigh_no< 3) return TS_FAIL;
    if(k->neigh_no< 3) return TS_FAIL;
    if(k==NULL || it==NULL){
        fatal("In bondflip, number of neighbours of k or it is less than 3!",999);
    }

    nei=find_neigh_idx(it, k);
    neip=next_small(nei, it->neigh_no);  //  must have it in correct order
    neim=prev_small(nei, it->neigh_no);

    km=it->neigh[neim];  // We located km and kp
    kp=it->neigh[neip];

    if(km==NULL || kp==NULL){
        fatal("In bondflip, cannot determine km and kp!",999);
    }
    if (km->type==is_ghost_vtx && kp->type==is_ghost_vtx) return TS_FAIL;

    if(it->type & is_edge_vtx || km->type & is_edge_vtx || k->type & is_edge_vtx || kp->type & is_edge_vtx){
        // not even bothering for the mixed case
        // hoping compiler can throw out repeated steps when composing
        retval = single_bondflip_timestep(vesicle, bond, rn);
        if (!(it->type & is_edge_vtx)) order_vertex_triangles(it);
        if (!(km->type & is_edge_vtx)) order_vertex_triangles(km);
        if (!( k->type & is_edge_vtx)) order_vertex_triangles(k);
        if (!(kp->type & is_edge_vtx)) order_vertex_triangles(kp);
        return retval;
    }

    // We now assume all vertices are ordered in neighbors and tristars 


    /* test if the membrane is wrapped too much, so that kp is nearest neighbour of
    * km. If it is true, then don't flip! */
    for(i=0;i<km->neigh_no;i++){
        if(km->neigh[i] == kp) return TS_FAIL;
    }
    //   fprintf(stderr,"Membrane didn't wrap too much.. Continue.\n");
    /* if bond would be too long, return... */
    if(vtx_distance_sq(km,kp) > vesicle->dmax ) return TS_FAIL;
    //   fprintf(stderr,"Bond will not be too long.. Continue.\n");

    
    /* we make a bond flip. this is different than in original fortran */
    // find lm, lp
    // 1. step. We find lm and lp from it->tristar !
    //much simpler, since in it->tristar lm(km,k) is associated with km and lp(k,kp) with k
    lm = it->tristar[neim];
    lp = it->tristar[nei];


    // we want to find the order of things around km, k, and kp:
    /*c  
    c      +----- k-------+              it = {...neim:km, nei:k, neip:kp}
    c      |lm1 / | \ lp1 |          {neim-1:lm2, neim:lm, nei:lp, nei+1:lp2}
    c      |  /   |   \   |             
    c      |/     |     \ |              km = {...x:k, x+1:it}
    c     km  lm  | lp   kp            {x-1:lm1, x:lm, x+1:lm2}
    c      |\     |     / |              k = {...y:kp, y+1:it, y+2:km}
    c      |  \   |   /   |            {y-1:lm1, y:lp, y+1:lm, y+2:lm1}
    c      |lm2 \ | / lp2 |              kp = {...z:it, z+1:k}
    c      +------it------+            {z-1:lp2, z:lp, z+1:lm1}
    c
    c we can find all vertices, triangles, and their index relation using only nei,x,y,z
    */
    //we look for important triangles lp1 and lm2.
    // lm2 can be found straight from it->tristar[neim-1]
    lm2 = it->tristar[prev_small(neim, it->tristar_no)];
    nei_k_at_km  = find_neigh_idx(km, k); //x
    nei_kp_at_k  = find_neigh_idx(k, kp); //y
    nei_it_at_kp = find_neigh_idx(kp, it);//z
    //lm2 = km->tristar[next_small(nei_k_at_km, km->tristar_no)]; // already done
    lp1 = k->tristar[prev_small(nei_kp_at_k, k->tristar_no)];
    
    // ####### Passed all base tests #######

    /* backup old structure */
    /* need to backup:
    * vertices k, kp, km, it
    * triangles lm, lp, lm2, lp1
    * bond
    */
    ts_vertex *bck_vtx[4];
    ts_triangle *bck_tria[4];
    ts_bond *bck_bond;
    ts_vertex *orig_vtx[]={k,it,kp,km};
    ts_triangle *orig_tria[]={lm,lp,lm2,lp1};

    //fprintf(stderr,"Backuping!!!\n");
	bck_bond=(ts_bond *)malloc(sizeof(ts_bond));
    for(i=0;i<4;i++){
    /*	fprintf(stderr,"vtx neigh[%d]=",i);
	    for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
	    fprintf(stderr,"\n");
    */
	    bck_vtx[i]=(ts_vertex *)malloc(sizeof(ts_vertex));
	    bck_tria[i]=(ts_triangle *)malloc(sizeof(ts_triangle));
    	memcpy((void *)bck_vtx[i],(void *)orig_vtx[i],sizeof(ts_vertex));
    	memcpy((void *)bck_tria[i],(void *)orig_tria[i],sizeof(ts_triangle));
	    /* level 2 pointers */

    	bck_vtx[i]->neigh=(ts_vertex **)malloc(orig_vtx[i]->neigh_no*sizeof(ts_vertex *));
    	bck_vtx[i]->tristar=(ts_triangle **)malloc(orig_vtx[i]->tristar_no*sizeof(ts_triangle *));
    	bck_vtx[i]->bond=(ts_bond **)malloc(orig_vtx[i]->bond_no*sizeof(ts_bond *));
    	bck_tria[i]->neigh=(ts_triangle **)malloc(orig_tria[i]->neigh_no*sizeof(ts_triangle *));

    	memcpy((void *)bck_vtx[i]->neigh,(void *)orig_vtx[i]->neigh,orig_vtx[i]->neigh_no*sizeof(ts_vertex *));
    	memcpy((void *)bck_vtx[i]->tristar,(void *)orig_vtx[i]->tristar,orig_vtx[i]->tristar_no*sizeof(ts_triangle *));
    	memcpy((void *)bck_vtx[i]->bond,(void *)orig_vtx[i]->bond,orig_vtx[i]->bond_no*sizeof(ts_bond *));
    
    	memcpy((void *)bck_tria[i]->neigh,(void *)orig_tria[i]->neigh,orig_tria[i]->neigh_no*sizeof(ts_triangle *));	
    }
	memcpy(bck_bond,bond,sizeof(ts_bond));
    //fprintf(stderr,"Backup complete!!!\n");
    /* end backup vertex */


    /* Save old energy */
    oldenergy=0;
    oldenergy+=k->energy;
    oldenergy+=kp->energy;
    oldenergy+=km->energy;
    oldenergy+=it->energy;
    oldenergy+=bond->energy; /* attraction with neighboring vertices, that have spontaneous curvature */
    //Neigbours of k, it, km, kp don't change its energy.

	if(vesicle->pswitch == 1 || vesicle->tape->constvolswitch>0){dvol = -lm->volume - lp->volume;}
    if(vesicle->tape->constareaswitch==2){darea=-lm->area-lp->area;} 
    /*    vesicle_volume(vesicle);
    fprintf(stderr,"Volume in the beginning=%1.16e\n", vesicle->volume);
    */
    
    // ####################################//
    /* fix data structure for flipped bond */
    // ####################################//
    ts_flip_bond_ordered(vesicle, bond, it, neim, km, nei_k_at_km, 
                         k, nei_kp_at_k, kp, nei_it_at_kp, lm, lp, lm2, lp1);


    /* Calculating the new energy */
    delta_energy=0;
    delta_energy+=k->energy;
    delta_energy+=kp->energy;
    delta_energy+=km->energy;
    delta_energy+=it->energy;
    delta_energy+=bond->energy; /* attraction with neighboring vertices, that have spontaneous curvature */
    //Neigbours of k, it, km, kp don't change its energy.
	if(vesicle->tape->stretchswitch==1){
		oldenergy+=lm->energy+lp->energy;
		stretchenergy(vesicle,lm);
		stretchenergy(vesicle,lp);
		delta_energy+=lm->energy+lp->energy;
	}

    delta_energy-=oldenergy;
	if(vesicle->pswitch == 1 || vesicle->tape->constvolswitch>0){
		dvol = dvol + lm->volume + lp->volume;
		if(vesicle->pswitch==1) delta_energy-= vesicle->pressure*dvol;
	}


    if(vesicle->tape->constareaswitch==2){
        darea=darea+lm->area+lp->area; 
/*check whether the dvol is gt than epsvol */
		if(fabs(vesicle->area+darea-A0)>epsarea){
			//restore old state.
			/* restoration procedure copied from few lines below */
			    for(i=0;i<4;i++){
			//			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
				free(orig_vtx[i]->neigh);
				free(orig_vtx[i]->tristar);
				free(orig_vtx[i]->bond);
				free(orig_tria[i]->neigh);
				memcpy((void *)orig_vtx[i],(void *)bck_vtx[i],sizeof(ts_vertex));
				memcpy((void *)orig_tria[i],(void *)bck_tria[i],sizeof(ts_triangle));
			//			fprintf(stderr,"Restored vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
				/* level 2 pointers are redirected*/
			    }
			    memcpy(bond,bck_bond,sizeof(ts_bond));
			    for(i=0;i<4;i++){
				free(bck_vtx[i]);
				free(bck_tria[i]);
			/*			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d =",i, orig_vtx[i]->neigh_no );
				for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
				fprintf(stderr,"\n"); */
			    }
			    free(bck_bond);
			    return TS_FAIL;

		}
    }




	if(vesicle->tape->constvolswitch == 2){
		/*check whether the dvol is gt than epsvol */
		if(fabs(vesicle->volume+dvol-V0)>epsvol){
			//restore old state.
			/* restoration procedure copied from few lines below */
			    for(i=0;i<4;i++){
			//			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
				free(orig_vtx[i]->neigh);
				free(orig_vtx[i]->tristar);
				free(orig_vtx[i]->bond);
				free(orig_tria[i]->neigh);
				memcpy((void *)orig_vtx[i],(void *)bck_vtx[i],sizeof(ts_vertex));
				memcpy((void *)orig_tria[i],(void *)bck_tria[i],sizeof(ts_triangle));
			//			fprintf(stderr,"Restored vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
				/* level 2 pointers are redirected*/
			    }
			    memcpy(bond,bck_bond,sizeof(ts_bond));
			    for(i=0;i<4;i++){
				free(bck_vtx[i]);
				free(bck_tria[i]);
			/*			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d =",i, orig_vtx[i]->neigh_no );
				for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
				fprintf(stderr,"\n"); */
			    }
			    free(bck_bond);
			    return TS_FAIL;

		}

	} else
    if(vesicle->tape->constvolswitch == 1){
        retval=constvolume(vesicle, it, -dvol, &delta_energy_cv, &constvol_vtx_moved,&constvol_vtx_backup);
        if(retval==TS_FAIL){
/* restoration procedure copied from few lines below */
            for(i=0;i<4;i++){
    //			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
                free(orig_vtx[i]->neigh);
                free(orig_vtx[i]->tristar);
                free(orig_vtx[i]->bond);
                free(orig_tria[i]->neigh);
                memcpy((void *)orig_vtx[i],(void *)bck_vtx[i],sizeof(ts_vertex));
                memcpy((void *)orig_tria[i],(void *)bck_tria[i],sizeof(ts_triangle));
    //			fprintf(stderr,"Restored vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
                /* level 2 pointers are redirected*/
            }
            memcpy(bond,bck_bond,sizeof(ts_bond));
            for(i=0;i<4;i++){
                free(bck_vtx[i]);
                free(bck_tria[i]);
    /*			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d =",i, orig_vtx[i]->neigh_no );
                for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
                fprintf(stderr,"\n"); */
            }
            free(bck_bond);
            return TS_FAIL;
        }
        delta_energy+=delta_energy_cv;
    }


/* MONTE CARLO */
    if(delta_energy>=0){
#ifdef TS_DOUBLE_DOUBLE
        if(exp(-delta_energy)< drand48())
#endif
#ifdef TS_DOUBLE_FLOAT
        if(expf(-delta_energy)< (ts_float)drand48())
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        if(expl(-delta_energy)< (ts_ldouble)drand48())
#endif
        {
            //not accepted, reverting changes
	        //restore all backups
            // fprintf(stderr,"Restoring!!!\n");
            if(vesicle->tape->constvolswitch == 1){
                constvolumerestore(vesicle, constvol_vtx_moved,constvol_vtx_backup);
            }

		    for(i=0;i<4;i++){
            //	fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
			free(orig_vtx[i]->neigh);
			free(orig_vtx[i]->tristar);
			free(orig_vtx[i]->bond);
			free(orig_tria[i]->neigh);
			memcpy((void *)orig_vtx[i],(void *)bck_vtx[i],sizeof(ts_vertex));
			memcpy((void *)orig_tria[i],(void *)bck_tria[i],sizeof(ts_triangle));
//			fprintf(stderr,"Restored vtx neigh[%d] with neighbours %d\n",i, orig_vtx[i]->neigh_no );
			/* level 2 pointers are redirected*/
		}
		memcpy(bond,bck_bond,sizeof(ts_bond));

		for(i=0;i<4;i++){
			free(bck_vtx[i]);
			free(bck_tria[i]);
/*			fprintf(stderr,"Restoring vtx neigh[%d] with neighbours %d =",i, orig_vtx[i]->neigh_no );
			for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
			fprintf(stderr,"\n"); */
		}

		free(bck_bond);

//		fprintf(stderr,"Restoration complete!!!\n");
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after fail=%1.16e\n", vesicle->volume);
	if(vesicle->tape->stretchswitch==1){
		stretchenergy(vesicle,lm);
		stretchenergy(vesicle,lp);
	}

		return TS_FAIL;
        }
    }
     /* IF BONDFLIP ACCEPTED, THEN RETURN SUCCESS! */
//            fprintf(stderr,"SUCCESS!!!\n");

    if(vesicle->tape->constvolswitch == 2){
	    vesicle->volume+=dvol;
    } else if(vesicle->tape->constvolswitch == 1){
        constvolumeaccept(vesicle,constvol_vtx_moved,constvol_vtx_backup);
    }
    if(vesicle->tape->constareaswitch==2){
        vesicle->area+=darea;
    }
	// delete all backups
	for(i=0;i<4;i++){
	    free(bck_vtx[i]->neigh);
	    free(bck_vtx[i]->bond);
	    free(bck_vtx[i]->tristar);
	    free(bck_vtx[i]);
 	    free(bck_tria[i]->neigh);
        free(bck_tria[i]);
	}
	free(bck_bond);

//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after success=%1.16e\n", vesicle->volume);
    return TS_SUCCESS;
}

// flip bond while preserving all orders.
ts_bool ts_flip_bond_ordered(ts_vesicle *vesicle, ts_bond *bond,
                             ts_vertex *it, ts_small_idx neim,
                             ts_vertex *km, ts_small_idx nei_k_at_km,
                             ts_vertex *k, ts_small_idx nei_kp_at_k,
                             ts_vertex *kp, ts_small_idx nei_it_at_kp,
                             ts_triangle *lm, ts_triangle *lp,
                             ts_triangle *lm2, ts_triangle *lp1){
    // we know order of things around km, k, and kp:
    // nei,x,y,z give all the relations between the components
    /*c  
    c      +----- k-------+              it = {...neim:km, nei:k, neip:kp}
    c      |lm1 / | \ lp1 |          {neim-1:lm2, neim:lm, nei:lp, nei+1:lp2}
    c      |  /   |   \   |             
    c      |/     |     \ |              km = {...x:k, x+1:it}
    c     km  lm  | lp   kp            {x-1:lm1, x:lm, x+1:lm2}
    c      |\     |     / |              k = {...y:kp, y+1:it, y+2:km}
    c      |  \   |   /   |            {y-1:lm1, y:lp, y+1:lm, y+2:lm1}
    c      |lm2 \ | / lp2 |              kp = {...z:it, z+1:k}
    c      +------it------+            {z-1:lp2, z:lp, z+1:lm1}
    c
    c     Flipping:
    c
    c           +----- k ------+        it = {...neim:km, nei:kp} -> remove nei
    c           |lm1 /   \ lp1 |    {neim-1:lm2, neim:lp, nei:lp2} -> remove neim, shift back
    c           |  /       \   |       
    c  FLIP     |/    lm     \ |        km = {...x:k, x+1:kp, x+2:it} -> insert kp at x+1
    c --->     km ----------  kp      {x-1:lm1, x:lm, x+1:lp, x+2:lm2} -> insert lp at x+1
    c           |\    lp     / |        k = {...y:kp, y+1::km} -> remove y+1
    c           |  \       /   |      {y-1:lm1, y:lm, y+1:lm1} -> remove , shift back
    c           |lm2 \   / lp2 |        kp = {...z:it, z+1:km. z+2:k} -> insert km at z+1
    c           +----- it -----+      {z-1:lp2, z:lp, z+1:lm1} -> insert lm at z+1
    c
    c  and insert/remove triangles as usual
    c  lm = {it,k,km}->{kp,k,km}   lp = {k, it, kp}->{km, it,kp}
    c  switching lm(it->kp) and lp(k->km) maintains whatever order they had
    c
    c  changing triangle triangle neighbors
    c
    c  lm={lm1.lp,lm2}->{lm1,lp, lp1}, lp={lp1,lp2,lm}->{lm2,lp2,lm}
    c  lm2={?1,?2,lm}->{?1,?2,lp}, lp1={?1,?2,lp}->{?1,?2,lm}
    c
    c  switching lm(lm2->lp1), lp(lp1->lm2) reverses the order!
    c  switching lm2(lm->lp), lp(lp->lm) maintains the order
    */
   
    ts_small_idx it_rem, km_add, k_rem, kp_add;
    if(k==NULL || it==NULL || km==NULL || kp==NULL){
        fatal("ts_flip_bond: You called me with invalid pointers to vertices",999);
    }
    
    // 1. step. Correct bond_list (Samo doesn't know why he still has it)
    bond->vtx1=km;
    bond->vtx2=kp;

    // 2. step. We change the triangle vertices... (actual bond flip)
    swap_tri_vertex(lm,it,kp);
    swap_tri_vertex(lp,k,km);


    // 3. step. remove and insert neighbors of vertices and the bond
    it_rem=next_small(neim, it->neigh_no); // next(idx,no) will change after inserting/removing neighbors
    km_add=next_small(nei_k_at_km, km->neigh_no);
    k_rem =next_small(nei_kp_at_k, k->neigh_no);
    kp_add=next_small(nei_it_at_kp, kp->neigh_no);

    vtx_remove_neigh_at(it, it_rem);
    vtx_remove_bond_at(it, it_rem);

    vtx_insert_neigh_at(km, kp, km_add);
    vtx_insert_bond_at(km, bond, km_add);

    vtx_remove_neigh_at(k, k_rem);
    vtx_remove_bond_at(k, k_rem);
    
    vtx_insert_neigh_at(kp, km, kp_add);
    vtx_insert_bond_at(kp, bond, kp_add);

    // 4. step. Correct tristar for vertices km, kp, k and it
    // we need to shift back the ones to remove:
    it->tristar[neim]=lp; // copy the remaining one on top of the removed
    vtx_remove_tristar_at(it, it_rem); // remove the original one, shifting the rest correctly
    vtx_insert_tristar_at(km, lp, km_add);
    k->tristar[nei_kp_at_k]=lm;
    vtx_remove_tristar_at(k, k_rem);
    vtx_insert_tristar_at(kp, lm, kp_add);


    // 5. step. Correct neighbouring triangles 
    swap_tri_neigh(lm2, lm, lp);
    swap_tri_neigh(lp1, lp, lm);
    swap_tri_neigh_and_reverse(lm, lm2, lp1);
    swap_tri_neigh_and_reverse(lp, lp1, lm2);

    // END modifications to data structure!

    // 7. step. Update properties of the structures: normals and energy
    triangle_normal_vector(lp);
    triangle_normal_vector(lm);
    energy_vertex(vesicle, k);
    energy_vertex(vesicle, kp);
    energy_vertex(vesicle, km);
    energy_vertex(vesicle, it);
    // energy also updates the normals
    // Yoav: Do we need to update force?

    // Yoav; No more sense in giving w to the energy: it's in the vertices now
    attraction_bond_energy(vesicle, bond);
    return TS_SUCCESS;
}

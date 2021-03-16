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
    ts_uint nei,neip,neim;
    ts_uint i,j;
    ts_double oldenergy, delta_energy, dvol=0.0, darea=0.0;
    ts_triangle *lm=NULL,*lp=NULL, *lp1=NULL, *lm2=NULL;

    ts_vertex *kp,*km;

    ts_double delta_energy_cv;
    ts_vertex *constvol_vtx_moved, *constvol_vtx_backup;
    ts_bool retval;

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
    if((ts_int)neim<0) neim=it->neigh_no-1; /* casting is essential... If not
there the neim is never <0 !!! */
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
if(lm==NULL || lp==NULL) fatal("ts_flip_bond: Cannot find triangles lm and lp!",999);

//we look for important triangles lp1 and lm2.

 for(i=0;i<k->tristar_no;i++){
        for(j=0;j<kp->tristar_no;j++){
                if((k->tristar[i] == kp->tristar[j]) && k->tristar[i]!=lp){ //ce gre za skupen trikotnik
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
  oldenergy+=k->xk* k->energy;
  oldenergy+=kp->xk* kp->energy;
  oldenergy+=km->xk* km->energy;
  oldenergy+=it->xk* it->energy;
  oldenergy+=bond->energy; /* attraction with neighboring vertices, that have spontaneous curvature */
  //Neigbours of k, it, km, kp don't change its energy.

	if(vesicle->pswitch == 1 || vesicle->tape->constvolswitch>0){dvol = -lm->volume - lp->volume;}
    if(vesicle->tape->constareaswitch==2){darea=-lm->area-lp->area;} 
/*    vesicle_volume(vesicle);
    fprintf(stderr,"Volume in the beginning=%1.16e\n", vesicle->volume);
*/
/* fix data structure for flipped bond */
    ts_flip_bond(k,it,km,kp, bond,lm, lp, lm2, lp1, vesicle->tape->w);


/* Calculating the new energy */
  delta_energy=0;
  delta_energy+=k->xk* k->energy;
  delta_energy+=kp->xk* kp->energy;
  delta_energy+=km->xk* km->energy;
  delta_energy+=it->xk* it->energy;
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
//		fprintf(stderr,"Restoring!!!\n");
        if(vesicle->tape->constvolswitch == 1){
            constvolumerestore(constvol_vtx_moved,constvol_vtx_backup);
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
/*	fprintf(stderr,"Afret backup deletion vtx neigh[%d]=",i);
	for(j=0;j<orig_vtx[i]->neigh_no;j++) fprintf(stderr," %d", orig_vtx[i]->neigh[j]->idx);
	fprintf(stderr,"\n");
*/	
	}
	free(bck_bond);

//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after success=%1.16e\n", vesicle->volume);
    return TS_SUCCESS;
}


ts_bool ts_flip_bond(ts_vertex *k,ts_vertex *it,ts_vertex *km, ts_vertex *kp,
ts_bond *bond, ts_triangle *lm, ts_triangle *lp, ts_triangle *lm2, ts_triangle *lp1, ts_double w_energy){

    ts_uint i; //lmidx, lpidx;
if(k==NULL || it==NULL || km==NULL || kp==NULL){
    fatal("ts_flip_bond: You called me with invalid pointers to vertices",999);
}
// 2. step. We change the triangle vertices... (actual bond flip)
    for(i=0;i<3;i++) if(lm->vertex[i]== it) lm->vertex[i]= kp;
    for(i=0;i<3;i++) if(lp->vertex[i]== k) lp->vertex[i]= km;
//fprintf(stderr,"2. step: actual bondflip made\n");
// 2a. step. If any changes in triangle calculations must be done, do it here!
//   * normals are recalculated here
    triangle_normal_vector(lp);
    triangle_normal_vector(lm);
//fprintf(stderr,"2a. step: triangle normals recalculated\n");
// 3. step. Correct neighbours in vertex_list


            vtx_remove_neighbour(k,it);
//            vtx_remove_neighbour(it,k);
//fprintf(stderr,"3. step (PROGRESS): removed k and it neighbours\n");
    
            //Tukaj pa nastopi tezava... Kam dodati soseda?
            vtx_insert_neighbour(km,kp,k);
            vtx_insert_neighbour(kp,km,it);
//            vertex_add_neighbour(km,kp); //pazi na vrstni red.
//            vertex_add_neighbour(kp,km);
//fprintf(stderr,"3. step: vertex neighbours corrected\n");

// 3a. step. If any changes to ts_vertex, do it here!
//   bond_length calculatons not required for it is done in energy.c

// 4. step. Correct bond_list (don't know why I still have it!)
            bond->vtx1=km;
            bond->vtx2=kp;
//fprintf(stderr,"4. step: bondlist corrected\n");


// 5. step. Correct neighbouring triangles 
   
    triangle_remove_neighbour(lp,lp1);
  //  fprintf(stderr,".\n");
    triangle_remove_neighbour(lp1,lp);
  //  fprintf(stderr,".\n");
    triangle_remove_neighbour(lm,lm2);
  //  fprintf(stderr,".\n");
    triangle_remove_neighbour(lm2,lm);
   
    triangle_add_neighbour(lm,lp1);    
    triangle_add_neighbour(lp1,lm);
    triangle_add_neighbour(lp,lm2);  //Vrstni red?!
    triangle_add_neighbour(lm2,lp);

//fprintf(stderr,"5. step: triangle neigbours corrected\n");


// 6. step. Correct tristar for vertices km, kp, k and it
            vertex_add_tristar(km,lp);  // Preveri vrstni red!
            vertex_add_tristar(kp,lm);
            vtx_remove_tristar(it,lm);
            vtx_remove_tristar(k,lp);
//fprintf(stderr,"6. step: tristar corrected\n");
  energy_vertex(k);
  energy_vertex(kp);
  energy_vertex(km);
  energy_vertex(it);
  attraction_bond_energy(bond, w_energy);
// END modifications to data structure!
    return TS_SUCCESS;
}

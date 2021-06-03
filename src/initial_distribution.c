/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "vesicle.h"
#include "vertex.h"
#include "triangle.h"
#include "initial_distribution.h"
#include "energy.h"
#include "poly.h"
#include "io.h"
#include "sh.h"
#include "shcomplex.h"

ts_vesicle *initial_distribution_dipyramid(ts_uint nshell, ts_uint ncmax1, ts_uint ncmax2, ts_uint ncmax3, ts_double stepsize){
	ts_fprintf(stdout,"Starting initial_distribution on vesicle with %u shells!...\n",nshell);
	ts_bool retval;
	ts_uint no_vertices=5*nshell*nshell+2;	
	ts_vesicle *vesicle=init_vesicle(no_vertices,ncmax1,ncmax2,ncmax3,stepsize);
	vesicle->nshell=nshell;
	//retval = vtx_set_global_values(vesicle);
	retval = pentagonal_dipyramid_vertex_distribution(vesicle->vlist);
	retval = init_vertex_neighbours(vesicle->vlist);
	vesicle->vlist = init_sort_neighbours(vesicle->blist,vesicle->vlist);
   // retval = init_vesicle_bonds(vesicle); // bonds are created in sort_neigh
	retval = init_triangles(vesicle);
	retval = init_triangle_neighbours(vesicle);
	retval = init_common_vertex_triangle_neighbours(vesicle);
	retval = init_normal_vectors(vesicle->tlist);
	retval = mean_curvature_and_energy(vesicle);
	ts_fprintf(stdout,"initial_distribution finished!\n");
	if(retval);
	return vesicle;
} 



ts_vesicle *create_vesicle_from_tape(ts_tape *tape){
	ts_vesicle *vesicle;

	vesicle=initial_distribution_dipyramid(tape->nshell,tape->ncxmax,tape->ncymax,tape->nczmax,tape->stepsize);
    	vesicle->tape=tape;
	vesicle->clist->dmin_interspecies = tape->dmin_interspecies*tape->dmin_interspecies;
	vesicle->poly_list=init_poly_list(tape->npoly,tape->nmono, vesicle->vlist, vesicle);
	set_vesicle_values_from_tape(vesicle);
	initial_population(vesicle,tape);
	return vesicle;
}

ts_bool set_vesicle_values_from_tape(ts_vesicle *vesicle){
	// Nucleus:
	ts_vertex *vtx;
	ts_tape *tape=vesicle->tape;
	vesicle->R_nucleus=tape->R_nucleus*tape->R_nucleus;
	vesicle->R_nucleusX=tape->R_nucleusX*tape->R_nucleusX;
	vesicle->R_nucleusY=tape->R_nucleusY*tape->R_nucleusY;
	vesicle->R_nucleusZ=tape->R_nucleusZ*tape->R_nucleusZ;
	vesicle->clist->dmin_interspecies = tape->dmin_interspecies*tape->dmin_interspecies;

	//Initialize grafted polymers (brush):
	//vesicle->poly_list=init_poly_list(tape->npoly,tape->nmono, vesicle->vlist, vesicle);
	vesicle->spring_constant=tape->kspring;
	poly_assign_spring_const(vesicle);

	//Initialize filaments (polymers inside the vesicle):
	vesicle->filament_list=init_poly_list(tape->nfil,tape->nfono, NULL, vesicle);
	poly_assign_filament_xi(vesicle,tape);

	ts_uint i,j;
	for(i=0;i<vesicle->filament_list->n;i++){
		for(j=0;j<vesicle->filament_list->poly[i]->blist->n;j++){
			bond_vector(vesicle->filament_list->poly[i]->blist->bond[j]);
			vesicle->filament_list->poly[i]->blist->bond[j]->bond_length = sqrt(vtx_distance_sq(vesicle->filament_list->poly[i]->blist->bond[j]->vtx1,vesicle->filament_list->poly[i]->blist->bond[j]->vtx2));
		}
	}

	for(i=0;i<vesicle->filament_list->n;i++){
		for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
			vtx = vesicle->filament_list->poly[i]->vlist->vtx[j];
			if(vtx->bond_no == 2){
			vtx->energy = -(vtx->bond[0]->x*vtx->bond[1]->x + vtx->bond[0]->y*vtx->bond[1]->y + vtx->bond[0]->z*vtx->bond[1]->z)/vtx->bond[0]->bond_length/vtx->bond[1]->bond_length;
			}
		}
	}

	for(i=0;i<vesicle->filament_list->n;i++){
		vertex_list_assign_id(vesicle->filament_list->poly[i]->vlist,TS_ID_FILAMENT);
	}


//Initialize the position of the adhesive surface
//	ts_double z_adh=1e10;
//	for(k=0;k<vesicle->vlist->n;k++){
//		if(vesicle->vlist->vtx[k]->z < z_adh) z_adh=vesicle->vlist->vtx[k]->z;
//	}	
//	vesicle->adhesion.z_adhesion = z_adh;


//	vesicle->spring_constant=tape->kspring;
//	poly_assign_spring_const(vesicle);

	
	vesicle->nshell=tape->nshell;
	vesicle->dmax=tape->dmax*tape->dmax; /* dmax^2 in the vesicle dmax variable */
	vesicle->bending_rigidity=tape->xk0;
	vtx_set_global_values(vesicle); /* make xk0 default value for every vertex */ 
//	ts_fprintf(stdout, "Tape setting: xk0=%e\n",tape->xk0);
	vesicle->stepsize=tape->stepsize;
	vesicle->clist->ncmax[0]=tape->ncxmax;
	vesicle->clist->ncmax[1]=tape->ncymax;
	vesicle->clist->ncmax[2]=tape->nczmax;
	vesicle->clist->max_occupancy=16; /* hard coded max occupancy? */

	vesicle->pressure= tape->pressure;
	vesicle->pswitch=tape->pswitch;
    if(tape->shc>0){
	    vesicle->sphHarmonics=complex_sph_init(vesicle->vlist,tape->shc);
    }
    else {
        vesicle->sphHarmonics=NULL;
    }

	vesicle->tlist->a0=sqrt(3)/4*pow((vesicle->tape->dmax+1.0)/2.0,2);  
    return TS_SUCCESS;

}


ts_bool initial_population(ts_vesicle *vesicle, ts_tape *tape){
	// modified by Yoav: initialize possibly heterogeneous vertices population
	ts_uint i,j, rndvtx, temp, n;
	ts_uint *indices;
	ts_vertex* vtx;
	n = vesicle->vlist->n;
	indices = (ts_uint*) malloc(n * sizeof(ts_uint));
	// create an array fo indices; copied a bunch from wikipedia "FIsher-Yates shuffle"
	for (i=0; i<n; i++){
		indices[i] = i;
	} // shuffle
	for(i=n-1; i>0 ; i--){
		rndvtx = rand() % (i+1);
		temp = indices[i];
		indices[i] = indices[rndvtx];
		indices[rndvtx] = temp;
	}
	// now we have an array of random indices: we can use it to populate the vertices however we like
	j=0;
	
	if (tape->number_of_vertices_with_c0>n) fatal("number of vertices with c0 larger than number of vertices!",100);
	for(i=j;i<j+tape->number_of_vertices_with_c0;i++){
		vtx = vesicle->vlist->vtx[indices[i]];
		vtx->type=is_bonding_vtx + is_active_vtx + is_adhesive_vtx; // bonds, active, adhesive, (not an)isotropic
		vtx->w=tape->w;
		vtx->c=tape->c0;
		vtx->f=tape->F;
		vtx->ad_w=tape->adhesion_strength;
		vtx->d=0;  // curvature deviator
		vtx->nx=0; //normal
		vtx->ny=0;
		vtx->nz=0;
		vtx->fx=0; //force
		vtx->fy=0;
		vtx->fz=0;
		vtx->tx=0; //director
		vtx->ty=0;
		vtx->tz=0;
	}
	j=i;
	// add however many types with for(i=j, i<j+n_next_type; i++){} j=i;

	//finally, regular non-protein vertices
	for(i=j; i<n; i++){
		vtx = vesicle->vlist->vtx[indices[i]];
		vtx->type=is_adhesive_vtx; // nonbonding, passive, adhesive, isotropic
		vtx->w=0;
		vtx->c=0;
		vtx->f=0;
		vtx->ad_w=tape->adhesion_strength;
		vtx->d=0;  // curvature deviator
		vtx->nx=0; //normal
		vtx->ny=0;
		vtx->nz=0;
		vtx->fx=0; //force
		vtx->fy=0;
		vtx->fz=0;
		vtx->tx=0; //director
		vtx->ty=0;
		vtx->tz=0;
	}
	free(indices);

	// This updates the normals using the energy_vertex
	mean_curvature_and_energy(vesicle);
	if(fabs(tape->w)>1e-16){ //if nonzero energy
		//	ts_fprintf(stderr,"Setting attraction between vertices with spontaneous curvature\n");
		sweep_attraction_bond_energy(vesicle);
	}
	return TS_SUCCESS;
}


ts_bool pentagonal_dipyramid_vertex_distribution(ts_vertex_list *vlist){
	/* Some often used relations */
	const ts_double s1= sin(2.0*M_PI/5.0);
	const ts_double s2= sin(4.0*M_PI/5.0);
	const ts_double c1= cos(2.0*M_PI/5.0);
	const ts_double c2= cos(4.0*M_PI/5.0);

	/* Calculates projection lenght of an edge bond to pentagram plane */
	const ts_double xl0=DEF_A0/(2.0*sin(M_PI/5.0));
#ifdef TS_DOUBLE_DOUBLE
	const ts_double z0=sqrt(pow(DEF_A0,2)-pow(xl0,2));
#endif
#ifdef TS_DOUBLE_FLOAT
	const ts_double z0=sqrtf(powf(DEF_A0,2)-powf(xl0,2));
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
	const ts_double z0=sqrtl(powl(DEF_A0,2)-powl(xl0,2));
#endif
//	const z0=sqrt(A0*A0 -xl0*xl0); /* I could use pow function but if pow is used make a check on the float type. If float then powf, if long double use powl */

/*placeholder for the pointer to vertex datastructure list... DIRTY: actual pointer points towards invalid address, one position before actual beginning of the list... This is to solve the difference between 1 based indexing in original program in fortran and 0 based indexing in C. All algorithms remain unchanged because of this!*/
	ts_vertex **vtx=vlist->vtx -1 ; 


	ts_uint nshell=(ts_uint)( sqrt((ts_double)(vlist->n-2)/5));
//	printf("nshell=%u\n",nshell);
	ts_uint i,n0; // some for loop prereq
	ts_int j,k;
	ts_double dx,dy; // end loop prereq

	/* topmost vertex */
	vtx[1]->x=0.0;
	vtx[1]->y=0.0;
	vtx[1]->z=z0*(ts_double)nshell;
	
	/* starting from to in circular order on pentagrams */	
	for(i=1;i<=nshell;i++){
		n0=2+5*i*(i-1)/2; //-1 would be for the reason that C index starts from 0 
		vtx[n0]->x=0.0;
		vtx[n0]->y=(ts_double)i*xl0;
		vtx[n0+i]->x=vtx[n0]->y*s1;
		vtx[n0+i]->y=vtx[n0]->y*c1;
		vtx[n0+2*i]->x=vtx[n0]->y*s2;
		vtx[n0+2*i]->y=vtx[n0]->y*c2;
		vtx[n0+3*i]->x=-vtx[n0+2*i]->x;
		vtx[n0+3*i]->y=vtx[n0+2*i]->y;
		vtx[n0+4*i]->x=-vtx[n0+i]->x;
		vtx[n0+4*i]->y=vtx[n0+i]->y;
	}

	/* vertexes on the faces of the dipyramid */
	for(i=1;i<=nshell;i++){
		n0=2+5*i*(i-1)/2; // -1 would be because of C!
		for(j=1;j<=i-1;j++){
			dx=(vtx[n0]->x-vtx[n0+4*i]->x)/(ts_double)i;
			dy=(vtx[n0]->y-vtx[n0+4*i]->y)/(ts_double)i;
			vtx[n0+4*i+j]->x=(ts_double)j*dx+vtx[n0+4*i]->x;
			vtx[n0+4*i+j]->y=(ts_double)j*dy+vtx[n0+4*i]->y;
		}
		for(k=0;k<=3;k++){ // I would be worried about zero starting of for
			dx=(vtx[n0+(k+1)*i]->x - vtx[n0+k*i]->x)/(ts_double) i;
			dy=(vtx[n0+(k+1)*i]->y - vtx[n0+k*i]->y)/(ts_double) i;
			for(j=1; j<=i-1;j++){
				vtx[n0+k*i+j]->x= (ts_double)j*dx+vtx[n0+k*i]->x;
				vtx[n0+k*i+j]->y= (ts_double)j*dy+vtx[n0+k*i]->y;
			} 
		} 
	}

	for(i=1;i<=nshell;i++){
		n0= 2+ 5*i*(i-1)/2;
		for(j=0;j<=5*i-1;j++){
		vtx[n0+j]->z= z0*(ts_double)(nshell-i);   // I would be worried about zero starting of for
		}
	}

/* for botom part of dipyramide we calculate the positions of vertices */
	for(i=2+5*nshell*(nshell+1)/2;i<=vlist->n;i++){
		vtx[i]->x=vtx[vlist->n - i +1]->x;
		vtx[i]->y=vtx[vlist->n - i +1]->y;
		vtx[i]->z=-vtx[vlist->n - i +1]->z;
	}

	for(i=1;i<=vlist->n;i++){
		for(j=1;j<=vlist->n;j++){
			if(i!=j && vtx_distance_sq(vtx[i],vtx[j])<0.001){
				printf("Vertices %u and %u are the same!\n",i,j);
			}
		}
	}
	return TS_SUCCESS;
}



ts_bool init_vertex_neighbours(ts_vertex_list *vlist){
	ts_vertex **vtx=vlist->vtx -1; // take a look at dipyramid function for comment.
	const ts_double eps=0.001; //TODO: find out if you can use EPS from math.h
	ts_uint i,j;
	ts_double dist2; // Square of distance of neighbours
	/*this is not required if we zero all data in vertex structure at initialization */
	/*if we force zeroing at initialization this for loop can safely be deleted */
	//for(i=1;i<=vlist->n;i++){
	//	vtx[i].neigh_no=0;
	//}
	for(i=1;i<=vlist->n;i++){
		for(j=1;j<=vlist->n;j++){
			dist2=vtx_distance_sq(vtx[i],vtx[j]);
			if( (dist2>eps) && (dist2<(DEF_A0*DEF_A0+eps))){ 
	//if it is close enough, but not too much close (solves problem of comparing when i==j)
				vtx_add_neighbour(vtx[i],vtx[j]);
			}
		}
	//		printf ("vertex %u ima %u sosedov!\n",i,vtx[i]->data->neigh_no);
	}

	return TS_SUCCESS;
}

// TODO: with new datastructure can be rewritten. Partially it is done, but it is complicated.
ts_vertex_list *init_sort_neighbours(ts_bond_list *blist,ts_vertex_list *vlist){
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
               	if( (fabs(dist2-DEF_A0*DEF_A0)<=eps) && (direct>0.0) && (j!=jjj) ){
           			vtx_add_cneighbour(blist,tvtx[k],tvtx[vtx[i]->neigh[j-1]->idx+1]);
           			jjj=jj;
           			jj=j;
           			break;
           		}
       		}
       	}	
	}
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


ts_bool init_vesicle_bonds(ts_vesicle *vesicle){
	ts_vertex_list *vlist=vesicle->vlist;
	ts_bond_list *blist=vesicle->blist;
	ts_vertex **vtx=vesicle->vlist->vtx - 1; // Because of 0 indexing
/* lets make correct clockwise ordering of in nearest neighbour list */
	ts_uint i,j,k;
	for(i=1;i<=vlist->n;i++){
		for(j=i+1;j<=vlist->n;j++){
			for(k=0;k<vtx[i]->neigh_no;k++){ // has changed 0 to < instead of 1 and <=
				if(vtx[i]->neigh[k]==vtx[j]){  //if addresses matches it is the same
					bond_add(blist,vtx[i],vtx[j]);
					break;
				}
			}
		}
	} 
/* Let's make a check if the number of bonds is correct */
    if((blist->n)!=3*(vlist->n-2)){
        ts_fprintf(stderr,"Number of bonds is %u should be %u!\n", blist->n, 3*(vlist->n-2));
        fatal("Number of bonds is not 3*(no_vertex-2).",4);
    }
	return TS_SUCCESS;
}



ts_bool init_triangles(ts_vesicle *vesicle){
	ts_uint i,j,jj,k;
	ts_vertex **vtx=vesicle->vlist->vtx -1; // difference between 0 indexing and 1 indexing
	ts_triangle_list *tlist=vesicle->tlist;
	ts_double dist, direct;
	ts_double eps=0.001; // can we use EPS from math.h?
	k=0;
	for(i=1;i<=vesicle->vlist->n;i++){
		for(j=1;j<=vtx[i]->neigh_no;j++){
			for(jj=1;jj<=vtx[i]->neigh_no;jj++){
		//		ts_fprintf(stderr,"%u: (%u,%u) neigh_no=%u ",i,j,jj,vtx[i].neigh_no);
        //      ts_fprintf(stderr,"%e, %e",vtx[i].neigh[j-1]->x,vtx[i].neigh[jj-1]->x);
				dist=vtx_distance_sq(vtx[i]->neigh[j-1],vtx[i]->neigh[jj-1]);
				direct=vtx_direct(vtx[i],vtx[i]->neigh[j-1],vtx[i]->neigh[jj-1]);				
// TODO: same as above				
				if(fabs(dist-DEF_A0*DEF_A0)<=eps && direct < 0.0 && vtx[i]->neigh[j-1]->idx+1 > i && vtx[i]->neigh[jj-1]->idx+1 >i){
					triangle_add(tlist,vtx[i],vtx[i]->neigh[j-1],vtx[i]->neigh[jj-1]);
				}	
			}	
		}
	}
/* We check if all triangles have 3 vertices and if the number of triangles
 * matches the theoretical value.
 */
	for(i=0;i<tlist->n;i++){
        k=0;
		for(j=0;j<3;j++){
			if(tlist->tria[i]->vertex[j]!=NULL)
            k++;
		}
            if(k!=3){
                fatal("Some triangles have less than 3 vertices..",4);
            }   
	} 
    if(tlist->n!=2*(vesicle->vlist->n -2)){
        ts_fprintf(stderr,"The number of triangles is %u but should be %u!\n",tlist->n,2*(vesicle->vlist->n -2));
        fatal("The number of triangles doesn't match 2*(no_vertex -2).",4);
    }
	return TS_SUCCESS;
}



ts_bool init_triangle_neighbours(ts_vesicle *vesicle){
	ts_uint i,j,nobo;
    ts_vertex *i1,*i2,*i3,*j1,*j2,*j3;
//	ts_vertex **vtx=vesicle->vlist->vtx -1; // difference between 0 indexing and 1 indexing
	ts_triangle_list *tlist=vesicle->tlist;
    ts_triangle **tria=tlist->tria -1;
    nobo=0;
    for(i=1;i<=tlist->n;i++){
        i1=tria[i]->vertex[0]; 
        i2=tria[i]->vertex[1]; 
        i3=tria[i]->vertex[2]; 
        for(j=1;j<=tlist->n;j++){
            if(j==i) continue;
            j1=tria[j]->vertex[0]; 
            j2=tria[j]->vertex[1]; 
            j3=tria[j]->vertex[2]; 
            if((i1==j1 && i3==j2) || (i1==j2 && i3==j3) || (i1==j3 && i3==j1)){
                    triangle_add_neighbour(tria[i],tria[j]);
                    nobo++;
            }
        }
    }
    for(i=1;i<=tlist->n;i++){
        i1=tria[i]->vertex[0]; 
        i2=tria[i]->vertex[1]; 
        i3=tria[i]->vertex[2]; 
        for(j=1;j<=tlist->n;j++){
            if(j==i) continue;
            j1=tria[j]->vertex[0]; 
            j2=tria[j]->vertex[1]; 
            j3=tria[j]->vertex[2]; 
            if((i1==j1 && i2==j3) || (i1==j3 && i2==j2) || (i1==j2 && i2==j1)){
                triangle_add_neighbour(tria[i],tria[j]);
                nobo++;
            }
        }
    }
    for(i=1;i<=tlist->n;i++){
        i1=tria[i]->vertex[0]; 
        i2=tria[i]->vertex[1]; 
        i3=tria[i]->vertex[2]; 
        for(j=1;j<=tlist->n;j++){
            if(j==i) continue;
            j1=tria[j]->vertex[0]; 
            j2=tria[j]->vertex[1]; 
            j3=tria[j]->vertex[2]; 
            if((i2==j1 && i3==j3) || (i2==j3 && i3==j2) || (i2==j2 && i3==j1)){
                triangle_add_neighbour(tria[i],tria[j]);
                nobo++;
            }
        }
    }
    if(nobo != vesicle->blist->n*2) {
            ts_fprintf(stderr,"Number of triangles= %u, number of bonds= %u\n",nobo/2, vesicle->blist->n);
            fatal("Number of triangle neighbour pairs differs from double the number of bonds!",4);
    }
    return TS_SUCCESS;
}


ts_bool init_common_vertex_triangle_neighbours(ts_vesicle *vesicle){
	ts_uint i,j,jp,k;
    ts_vertex *k1,*k2,*k3,*k4,*k5;
	ts_vertex **vtx=vesicle->vlist->vtx -1; // difference between 0 indexing and 1 indexing
	ts_triangle_list *tlist=vesicle->tlist;
    ts_triangle **tria=tlist->tria -1;

    for(i=1;i<=vesicle->vlist->n;i++){
        for(j=1;j<=vtx[i]->neigh_no;j++){
            k1=vtx[i]->neigh[j-1];
            jp=j+1;
            if(j == vtx[i]->neigh_no) jp=1;
            k2=vtx[i]->neigh[jp-1];
            for(k=1;k<=tlist->n;k++){		// VERY NON-OPTIMAL!!! too many loops (vlist.n * vtx.neigh * tlist.n )!
                k3=tria[k]->vertex[0];
                k4=tria[k]->vertex[1];
                k5=tria[k]->vertex[2];
//                ts_fprintf(stderr,"%u %u: k=(%u %u %u)\n",k1,k2,k3,k4,k5);
                if((vtx[i]==k3 && k1==k4 && k2==k5) ||
                (vtx[i]==k4 && k1==k5 && k2==k3) ||
                (vtx[i]==k5 && k1==k3 && k2==k4)){

//TODO: probably something wrong with neighbour distribution.
//                if(vtx[i]==k3 || vtx[i]==k4 || vtx[i]==k5){
    //                    if(i==6) ts_fprintf(stdout, "Vtx[%u] > Added to tristar!\n",i);
                    vertex_add_tristar(vtx[i],tria[k]);
                }
            }
        }
/*        ts_fprintf(stderr,"TRISTAR for %u (%u):",i-1,vtx[i].tristar_no);
        for(j=0;j<vtx[i].tristar_no;j++){
            ts_fprintf(stderr," %u,",vtx[i].tristar[j]->idx);
        }
        ts_fprintf(stderr,"\n"); */
    }
    return TS_SUCCESS;
}


ts_bool init_normal_vectors(ts_triangle_list *tlist){
	/* Normals point INSIDE vesicle */
	ts_uint k;
	ts_triangle **tria=tlist->tria -1; //for 0 indexing
	for(k=1;k<=tlist->n;k++){
		triangle_normal_vector(tria[k]);	
	}
	return TS_SUCCESS;
}



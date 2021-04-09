/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include "general.h"
#include "energy.h"
#include "vertex.h"
#include<math.h>
#include<stdio.h>


/** @brief Wrapper that calculates energy of every vertex in vesicle
 *  
 *  Function calculated energy of every vertex in vesicle. It can be used in
 *  initialization procedure or in recalculation of the energy after non-MCsweep *  operations. However, when random move of vertex or flip of random bond occur *  call to this function is not necessary nor recommended. 
 *  @param *vesicle is a pointer to vesicle.
 *  @returns TS_SUCCESS on success.
*/
ts_bool mean_curvature_and_energy(ts_vesicle *vesicle){

    ts_uint i;
    
    ts_vertex_list *vlist=vesicle->vlist;
    ts_vertex **vtx=vlist->vtx;

    for(i=0;i<vlist->n;i++){
        energy_vertex(vtx[i]);
        
    }

    return TS_SUCCESS;
}

/** @brief Calculate energy of a bond (in models where energy is bond related)
 *
 *  This function is experimental and currently only used in polymeres calculation (PEGs or polymeres inside the vesicle).
 *
 *  @param *bond is a pointer to a bond between two vertices in polymere
 *  @param *poly is a pointer to polymere in which we calculate te energy of the bond
 *  @returns TS_SUCCESS on successful calculation
*/
inline ts_bool bond_energy(ts_bond *bond,ts_poly *poly){
//TODO: This value to be changed and implemented in data structure:
	ts_double d_relaxed=1.0;
	bond->energy=poly->k*pow(bond->bond_length-d_relaxed,2);
	return TS_SUCCESS;
};

/** @brief Calculation of the bending energy of the vertex.
 *  
 *  Main function that calculates energy of the vertex \f$i\f$. Function returns \f$\frac{1}{2}(c_1+c_2-c)^2 s\f$, where \f$(c_1+c_2)/2\f$ is mean curvature,
 * \f$c/2\f$ is spontaneous curvature and \f$s\f$ is area per vertex \f$i\f$.
 *
 * Nearest neighbors (NN) must be ordered in counterclockwise direction for this function to work.
 *  Firstly NNs that form two neighboring triangles are found (\f$j_m\f$, \f$j_p\f$ and common \f$j\f$). Later, the scalar product of vectors \f$x_1=(\mathbf{i}-\mathbf{j_p})\cdot (\mathbf{i}-\mathbf{j_p})(\mathbf{i}-\mathbf{j_p})\f$, \f$x_2=(\mathbf{j}-\mathbf{j_p})\cdot  (\mathbf{j}-\mathbf{j_p})\f$  and \f$x_3=(\mathbf{j}-\mathbf{j_p})\cdot (\mathbf{i}-\mathbf{j_p})\f$  are calculated. From these three vectors the \f$c_{tp}=\frac{1}{\tan(\varphi_p)}\f$ is calculated, where \f$\varphi_p\f$ is the inner angle at vertex \f$j_p\f$. The procedure is repeated for \f$j_m\f$ instead of \f$j_p\f$ resulting in \f$c_{tn}\f$.
 *  
\begin{tikzpicture}{
\coordinate[label=below:$i$] (i) at (2,0);
\coordinate[label=left:$j_m$] (jm) at (0,3.7);
\coordinate[label=above:$j$] (j) at (2.5,6.4);
\coordinate[label=right:$j_p$] (jp) at (4,2.7);

\draw (i) -- (jm) -- (j) -- (jp) -- (i) -- (j);

\begin{scope}
\path[clip] (jm)--(i)--(j);
\draw (jm) circle (0.8);
\node[right] at (jm) {$\varphi_m$};
\end{scope}

\begin{scope}
\path[clip] (jp)--(i)--(j);
\draw (jp) circle (0.8);
\node[left] at (jp) {$\varphi_p$};
\end{scope}

%%vertices
\draw [fill=gray] (i) circle (0.1);
\draw [fill=white] (j) circle (0.1);
\draw [fill=white] (jp) circle (0.1);
\draw [fill=white] (jm) circle (0.1);
%\node[draw,circle,fill=white] at (i) {};
\end{tikzpicture}

 * The curvature is then calculated as \f$\mathbf{h}=\frac{1}{2}\Sigma_{k=0}^{\mathrm{neigh\_no}} c_{tp}^{(k)}+c_{tm}^{(k)} (\mathbf{j_k}-\mathbf{i})\f$, where \f$c_{tp}^{(k)}+c_{tm}^k=2\sigma^{(k)}\f$ (length in dual lattice?) and the previous equation can be written as \f$\mathbf{h}=\Sigma_{k=0}^{\mathrm{neigh\_no}}\sigma^{(k)}\cdot(\mathbf{j}-\mathbf{i})\f$ (See Kroll, p. 384, eq 70).
 *
 * From the curvature the enery is calculated by equation \f$E=\frac{1}{2}\mathbf{h}\cdot\mathbf{h}\f$.
 * @param *vtx is a pointer to vertex at which we want to calculate the energy
 * @returns TS_SUCCESS on successful calculation.
*/
inline ts_bool energy_vertex(ts_vertex *vtx){
    ts_uint jj;
    ts_uint jjp,jjm;
    ts_vertex *j,*jp, *jm;
    ts_triangle *jt;
    ts_double s=0.0,xh=0.0,yh=0.0,zh=0.0,txn=0.0,tyn=0.0,tzn=0.0;
    ts_double x1,x2,x3,ctp,ctm,tot,xlen;
    ts_double h,ht;
    for(jj=1; jj<=vtx->neigh_no;jj++){
        jjp=jj+1;
        if(jjp>vtx->neigh_no) jjp=1;
        jjm=jj-1;
        if(jjm<1) jjm=vtx->neigh_no;
        j=vtx->neigh[jj-1];
        jp=vtx->neigh[jjp-1];
        jm=vtx->neigh[jjm-1];
        jt=vtx->tristar[jj-1];
        x1=vtx_distance_sq(vtx,jp); //shouldn't be zero!
        x2=vtx_distance_sq(j,jp); // shouldn't be zero!
        x3=(j->x-jp->x)*(vtx->x-jp->x)+
           (j->y-jp->y)*(vtx->y-jp->y)+
           (j->z-jp->z)*(vtx->z-jp->z);
        
#ifdef TS_DOUBLE_DOUBLE
        ctp=x3/sqrt(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_FLOAT
        ctp=x3/sqrtf(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        ctp=x3/sqrtl(x1*x2-x3*x3);
#endif
        x1=vtx_distance_sq(vtx,jm);
        x2=vtx_distance_sq(j,jm);
        x3=(j->x-jm->x)*(vtx->x-jm->x)+
           (j->y-jm->y)*(vtx->y-jm->y)+
           (j->z-jm->z)*(vtx->z-jm->z);
#ifdef TS_DOUBLE_DOUBLE
        ctm=x3/sqrt(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_FLOAT
        ctm=x3/sqrtf(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        ctm=x3/sqrtl(x1*x2-x3*x3);
#endif
        tot=ctp+ctm;
        tot=0.5*tot;

        xlen=vtx_distance_sq(j,vtx);
/*
#ifdef  TS_DOUBLE_DOUBLE 
        vtx->bond[jj-1]->bond_length=sqrt(xlen); 
#endif
#ifdef  TS_DOUBLE_FLOAT
        vtx->bond[jj-1]->bond_length=sqrtf(xlen); 
#endif
#ifdef  TS_DOUBLE_LONGDOUBLE 
        vtx->bond[jj-1]->bond_length=sqrtl(xlen); 
#endif

        vtx->bond[jj-1]->bond_length_dual=tot*vtx->bond[jj-1]->bond_length;
*/
        s+=tot*xlen;
        xh+=tot*(j->x - vtx->x);
        yh+=tot*(j->y - vtx->y);
        zh+=tot*(j->z - vtx->z);
        txn+=jt->xnorm;
        tyn+=jt->ynorm;
        tzn+=jt->znorm;
    }
    
    h=xh*xh+yh*yh+zh*zh;
    ht=txn*xh+tyn*yh + tzn*zh;
    s=s/4.0; 
#ifdef TS_DOUBLE_DOUBLE
    if(ht>=0.0) {
        vtx->curvature=sqrt(h);
    } else {
        vtx->curvature=-sqrt(h);
    }
#endif
#ifdef TS_DOUBLE_FLOAT
    if(ht>=0.0) {
        vtx->curvature=sqrtf(h);
    } else {
        vtx->curvature=-sqrtf(h);
    }
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    if(ht>=0.0) {
        vtx->curvature=sqrtl(h);
    } else {
        vtx->curvature=-sqrtl(h);
    }
#endif
// c is forced curvature energy for each vertex. Should be set to zero for
// normal circumstances.
/* the following statement is an expression for $\frac{1}{2}\int(c_1+c_2-c_0^\prime)^2\mathrm{d}A$, where $c_0^\prime=2c_0$ (twice the spontaneous curvature)  */
    vtx->energy=0.5*s*(vtx->curvature/s-vtx->c)*(vtx->curvature/s-vtx->c);

    return TS_SUCCESS;
}



ts_bool sweep_attraction_bond_energy(ts_vesicle *vesicle){
	int i;
	for(i=0;i<vesicle->blist->n;i++){
		attraction_bond_energy(vesicle->blist->bond[i], vesicle->tape->w);
	}
	return TS_SUCCESS;
}


inline ts_bool attraction_bond_energy(ts_bond *bond, ts_double w){

	if(fabs(bond->vtx1->c)>1e-16 && fabs(bond->vtx2->c)>1e-16){
		bond->energy=-w;
	}
	else {
		bond->energy=0.0;
	}
	return TS_SUCCESS;
}


ts_double direct_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
	//Yoav- modified to include Vicsek Interaction
    
    //quit if this is a non-active vertex (or no vertex is active)
    if(fabs(vesicle->tape->F)<1e-15) return 0.0;
    if(vtx->c<1e-15) return 0.0;


    ts_double norml,ddp=0.0;
	ts_double xnorm=0.0,ynorm=0.0,znorm=0.0;
	ts_double vixnorm=0.0,viynorm=0.0,viznorm=0.0;

    // breadth-first search "seen" vertex list
    // should probably be a struct, but for now is kinda volatile
    // vertex list, with 4 locations: 
    // n_top, the top of the list, where new vertices are added
    // n_next_layer, first vertex of the next, under construction layer
    // n_curr_layer, first vertex of the last completed layer, used to construct the next
    // n_prev_layer, first vertex of the layer before the current one
    //
    // !!! mallocation size should be subject to change: roughly how many
    // vertex are expected to be counted, which highly depends on the
    // maximum distance and how crumpled the vesicle is, ~pi*r^2
    ts_uint max_vtx_seen=3*(( (int) vesicle->tape->vicsek_radius)+1)*(( (int) vesicle->tape->vicsek_radius)+1); //future: vesicle->tape->max_dist
    ts_vertex **seen_vtx=(ts_vertex**) malloc(sizeof(ts_vertex*)*max_vtx_seen);
    ts_uint n_prev_layer=0, n_curr_layer=0, n_next_layer=0, n_top=0;


    if (!vesicle->tape->vicsek_model || fabs(vesicle->tape->vicsek_strength)<1e-15 || fabs(vesicle->tape->vicsek_radius)<1e-15) {//no vicsek
        //regular "force in normal direction"

	
        //vertex normal
	    /*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
	    for(ts_uint i=0;i<vtx->tristar_no;i++){
	    	xnorm+=vtx->tristar[i]->xnorm;
	    	ynorm+=vtx->tristar[i]->ynorm;
	    	znorm+=vtx->tristar[i]->znorm;
    	}
	    /*normalize*/
        norml=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
        xnorm/=norml;
	    ynorm/=norml;
	    znorm/=norml;

	    /*calculate ddp, Viscek force directed displacement*/
	    ddp=xnorm*(vtx->x-vtx_old->x)+ynorm*(vtx->y-vtx_old->y)+znorm*(vtx->z-vtx_old->z);

    }
    else {//vicsek
        //force directed by Vicsek sum-over-neighbors-normals
        //(provided they are active too (c>0)

	
        //prime vertex normal
	    /*find normal of the vertex as sum of all the
        normals of the triangles surrounding it. */
	    for(ts_uint i=0;i<vtx->tristar_no;i++){
	    	xnorm+=vtx->tristar[i]->xnorm;
	    	ynorm+=vtx->tristar[i]->ynorm;
	    	znorm+=vtx->tristar[i]->znorm;
	    }
	    /*normalize*/
	    norml=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
	    vixnorm=xnorm/norml;
	    viynorm=ynorm/norml;
	    viznorm=znorm/norml;

        //seen the prime vertex
        seen_vtx[n_top]=vtx;
        n_top++;
        n_next_layer=n_top;
    
        // Breadth first search using seen_vtx and neighbors
        // while current distance from the prime vertex
        // is less than or equal to the maximum
        for (ts_uint curr_dist=1 ; curr_dist<=vesicle->tape->vicsek_radius; curr_dist++){
            
            // The for loops are split by layers
            // The new "next" layer is being built from
            // the neighbors of the completed "current" layer
            // some of the checked neighbors may be
            // from the layer before that "previous" layer
            // ..............
            //  A1 --  A2 --  O1     O-not-yet added              seen_vtx[]:
            // /  \  /  |  \  / \    N-bext layer, not added yet          [...
            // C1--C2--N3--O1 - O1   A-next layer, already added            Q1
            // | / | \ | \ |  \ |    C-current, neighbor iterated layer     Q2
            // P1--P2--C3--N4 - O1   P-previous layer                       P1  <-previous
            // | \ | \ | \ |  \ |    Q-previous previous layer              P2
            // Q1--Q2--P3--C4 - N5                                          P3
            // ..............                                               C1  <-current
            //                                                              C2
            //                                                              C3
            //                                                              C4
            //                                                              A1  <-next
            //                                                              A2
            //                                                              _   <-n_top

            for (ts_uint i=n_curr_layer; i<n_next_layer; i++) {
                // loop over the current layer: seen_vtx[i]
                

                for (ts_uint j=0; j<seen_vtx[i]->neigh_no; j++) {
                    //loop over their neighbors: seen_vtx[i]->neigh[j]

                    // is this neighbor is a bare vertex, skip to next neighbor
                    if (seen_vtx[i]->neigh[j]->c<1e-15) continue;
                    //else{ rest of the j loop }

                    // has this neighbor been seen?
                    if (was_vertex_seen(seen_vtx,seen_vtx[i]->neigh[j],n_prev_layer,n_top)){
                        continue; //if seen, skip to next neighbor in the loop
                    }
                    else {
                        //new vertex! add to the next layer

                        //add to seen_vtx (pending structure)
                        if (n_top >= max_vtx_seen) {
                            max_vtx_seen *= 2;
                            seen_vtx = realloc(seen_vtx, sizeof(ts_vertex *) * max_vtx_seen);
                            if (seen_vtx==NULL) fatal("Cannot reallocate memory to extend seen_vtx.",100);
                        }
                        seen_vtx[n_top] = seen_vtx[i]->neigh[j];

                        //calculate normal and add to the sum
                        xnorm = 0.0;
                        ynorm = 0.0;
                        znorm = 0.0;
                        for (ts_uint l = 0; l < seen_vtx[n_top]->tristar_no; l++){
                            xnorm += seen_vtx[n_top]->tristar[l]->xnorm;
                            ynorm += seen_vtx[n_top]->tristar[l]->ynorm;
                            znorm += seen_vtx[n_top]->tristar[l]->znorm;
                        }
                        /*normalize, and add to normal sum with weight by Vicsek strength*/
                        norml = sqrt(xnorm * xnorm + ynorm * ynorm + znorm * znorm);
                        vixnorm += vesicle->tape->vicsek_strength * xnorm / norml;
                        viynorm += vesicle->tape->vicsek_strength * ynorm / norml;
                        viznorm += vesicle->tape->vicsek_strength * znorm / norml;

                        n_top++; //next available spot
                    }
                }
            }
            // finished this layer
            if (n_next_layer==n_top) break; //did not add any vertices: completed cluster
            n_prev_layer=n_curr_layer;
            n_curr_layer=n_next_layer;
            n_next_layer=n_top;
        }

        //having finished summing the normals, normalize the resulting vector
        norml=sqrt(vixnorm*vixnorm+viynorm*viynorm+viznorm*viznorm);
        vixnorm/=norml;
	    viynorm/=norml;
	    viznorm/=norml;

	    /*calculate ddp, Viscek force directed displacement*/
	    ddp=vixnorm*(vtx->x-vtx_old->x)+viynorm*(vtx->y-vtx_old->y)+viznorm*(vtx->z-vtx_old->z);
	    /*calculate dE*/
    //	printf("ddp=%e",ddp);
	    
    //end else from if Vicsek
    }
    
    //don't forget to free! pending struct free_seen(seen_vtx)
    free(seen_vtx);
    return vesicle->tape->F*ddp;
//end function
}

ts_double direct_force_from_Fz_balance(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
	if(fabs(vesicle->tape->F)<1e-15) return 0.0;
    // calculate the Fz force balance separately
    // F_direct(old) = F_direct + Fz_balance \hat{z}
    // W = F_direct * dX + Fz_balance * dz
	ts_double ddp=0.0;
	ts_double Fz;
	force_per_vertex(vesicle, &Fz);

	ddp=-Fz*(vtx->z-vtx_old->z);

	return vesicle->tape->F*ddp;		
	
}

void force_per_vertex(ts_vesicle *vesicle, ts_double *Fz){
    ts_double norml;
	ts_uint i,j;
	ts_double fz=0.0;
	ts_double xnorm,ynorm,znorm;
	/*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
    for (j=0;j<vesicle->vlist->n;j++){
		if(vesicle->vlist->vtx[j]->c > 1e-15){
			xnorm=0.0;ynorm=0.0;znorm=0.0;
			for(i=0;i<vesicle->vlist->vtx[j]->tristar_no;i++){
				xnorm+=vesicle->vlist->vtx[j]->tristar[i]->xnorm;
				ynorm+=vesicle->vlist->vtx[j]->tristar[i]->ynorm;
				znorm+=vesicle->vlist->vtx[j]->tristar[i]->znorm;
			}
			/*normalize*/
			norml=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
			znorm/=norml;

			fz+=znorm;
		}
	}
	fz/=vesicle->vlist->n;
	if(fz<0.0000){
		*Fz=fz;
	}
	else if (fz>=0){
		*Fz=0.0;
	}
	
}

void stretchenergy(ts_vesicle *vesicle, ts_triangle *triangle){
	triangle->energy=vesicle->tape->xkA0/2.0*pow((triangle->area/vesicle->tlist->a0-1.0),2);
}

inline ts_bool was_vertex_seen(ts_vertex **seen_vtx, ts_vertex *vtx, ts_uint check_from, ts_uint check_up_to){
    // check if vertex was already accounted for in the "seen_vertex" list
    // 
    // we don't typically need to go over all the vertices: when we iterate
    // over neighbors of the current layer, building the next one, we only
    // need to check from the previous layer, the current layer, and the
    // layer we're building
    // if we're checking neighbors of C, we can only be at either P node,
    // C node, A node, or N node, and we've seen the P,C, and A nodes
    // ..............
    //  A1 --  A2 --  O1     O-not-yet added              seen_vtx[]:
    // /  \  /  |  \  / \    N-bext layer, not added yet          [...
    // C1--C2--N3--O1 - O1   A-next layer, already added            Q1
    // | / | \ | \ |  \ |    C-current, neighbor iterated layer     Q2
    // P1--P2--C3--N4 - O1   P-previous layer                       P1  <-check_from
    // | \ | \ | \ |  \ |    Q-previous previous layer              P2
    // Q1--Q2--P3--C4 - N5                                          P3
    // ..............                                               C1
    //                                                              C2
    //                                                              C3
    //                                                              C4
    //                                                              A1
    //                                                              A2
    //                                                              _ <-check_up_to


    for (ts_uint k=check_from; k<check_up_to; k++){
        if (seen_vtx[k]==vtx) return 1;
    }
    return 0;
}

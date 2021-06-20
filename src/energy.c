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
    ts_double h,ht,norml;
    ts_double angle_sum=0;
    ts_double a_dot_b, a_cross_b_x, a_cross_b_y, a_cross_b_z, mag_a_cross_b;
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
        //x1=pow(vtx->x-jp->x,2)+pow(vtx->y-jp->y,2)+pow(vtx->z-jp->z,2);
        //x2=pow(j->x-jp->x,2)+pow(j->y-jp->y,2)+pow(j->z-jp->z,2);
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
        //x1=pow(vtx->x-jm->x,2)+pow(vtx->y-jm->y,2)+pow(vtx->z-jm->z,2);
        //x2=pow(j->x-jm->x,2)+pow(j->y-jm->y,2)+pow(j->z-jm->z,2);
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

        //testing
        xlen=vtx_distance_sq(j,vtx);
        //xlen=pow(vtx->x-j->x,2)+pow(vtx->y-j->y,2)+pow(vtx->z-j->z,2);
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

        // angle stuff
        //angle_sum += atan(ctp) + atan(ctm); // simple but slow!
        /// get the angle m-vtx-j
        // atan2(|axb|,a*b) was recommended at mathwork forum (cosin has small angle problems, and still need a sqrt)
        a_dot_b = (jm->x-vtx->x)*(j->x-vtx->x)+
                  (jm->y-vtx->y)*(j->y-vtx->y)+
                  (jm->z-vtx->z)*(j->z-vtx->z);
        a_cross_b_x = (jm->y-vtx->y)*(j->z-vtx->z)-(jm->z-vtx->z)*(j->y-vtx->y);
        a_cross_b_y = (jm->z-vtx->z)*(j->x-vtx->x)-(jm->x-vtx->x)*(j->z-vtx->z);
        a_cross_b_z = (jm->x-vtx->x)*(j->y-vtx->y)-(jm->y-vtx->y)*(j->x-vtx->x);
        mag_a_cross_b = sqrt(pow(a_cross_b_x,2)+pow(a_cross_b_y,2)+pow(a_cross_b_z,2));
        angle_sum += atan2(mag_a_cross_b, a_dot_b);

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
    //also great point to update normal: note that the triangle normal is inwards
    norml=sqrt(txn*txn+tyn*tyn+tzn*tzn);
    vtx->nx=-txn/norml;
	vtx->ny=-tyn/norml;
	vtx->nz=-tzn/norml;
// c is forced curvature energy for each vertex. Should be set to zero for
// normal circumstances.
/* the following statement is an expression for $\frac{1}{2}\int(c_1+c_2-c_0^\prime)^2\mathrm{d}A$, where $c_0^\prime=2c_0$ (twice the spontaneous curvature)  */
    vtx->energy=vtx->xk* 0.5*s*(vtx->curvature/s-vtx->c)*(vtx->curvature/s-vtx->c);
    vtx->curvature2 = (2*M_PI- angle_sum)/s;
    if (vtx->type&is_anisotropic_vtx){
        vtx->energy += vtx->xk2 * s * vtx->curvature2;
    }

    return TS_SUCCESS;
}



ts_bool sweep_attraction_bond_energy(ts_vesicle *vesicle){
	int i;
	for(i=0;i<vesicle->blist->n;i++){
		attraction_bond_energy(vesicle->blist->bond[i]);
	}
	return TS_SUCCESS;
}


inline ts_bool attraction_bond_energy(ts_bond *bond){

	if((bond->vtx1->type&is_bonding_vtx && bond->vtx2->type&is_bonding_vtx)){
        // f(w1,w2)
		bond->energy=-0.5*(bond->vtx1->w+bond->vtx2->w);
	}
	else {
		bond->energy=0.0;
	}
	return TS_SUCCESS;
}


ts_double direct_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
	// modified to include Vicsek Interaction
    
    // quit if there is no point
    if(fabs(vtx->f)<1e-15) return 0.0;


    ts_double norml,ddp=0.0;
	//ts_double xnorm=0.0,ynorm=0.0,znorm=0.0;
	ts_double vixnorm=0.0,viynorm=0.0,viznorm=0.0;

    //stupid loop variables don't go in stupid for loop due to stupid C89 compiler
    ts_uint i, j, curr_dist;
    
    // estimated mallocation size for the cluster: roughly ~pi*r^2
    ts_uint max_vtx_seen=3*(( (int) vesicle->tape->vicsek_radius)+1)*(( (int) vesicle->tape->vicsek_radius)+1);
    // allocate the struct for the "seen vertex" (defined in general.h, functions in vertex.c)
    ts_seen_vertex *seen_vtx; //initalize in vicsek



    // if not vicsek type, or vicsek model is not relevant
    if ( !(vtx->type&is_vicsek_vtx) || !vesicle->tape->vicsek_model || fabs(vesicle->tape->vicsek_strength)<1e-15 || fabs(vesicle->tape->vicsek_radius)<1e-15) {//no vicsek
        //regular "force in normal direction"
        vtx->fx = vtx->nx;
        vtx->fy = vtx->ny;
        vtx->fz = vtx->nz;

    }
    else {
        //vicsek model
        //force directed by Vicsek sum-over-neighbors-normals

	
        //prime vertex normal
	    vixnorm=vtx->nx;
	    viynorm=vtx->ny;
	    viznorm=vtx->nz;

        //initialize seen_vtx
        seen_vtx = init_seen_vertex(max_vtx_seen);
        // we have now seen the prime vertex
        add_vtx_to_seen(seen_vtx, vtx);
        
    
        // Breadth first search using seen_vtx, layer by layer,
        // until reaching layer that is outside the maximum radius
        for (curr_dist=1 ; curr_dist<=vesicle->tape->vicsek_radius; curr_dist++){
            advance_seen_vertex_to_next_layer(seen_vtx);
            
            // The for loops are split by layers
            // The new "next" layer is being built from
            // the neighbors of the completed "current" layer
            // some of the checked neighbors may be
            // from the layer before that: the "previous" layer
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

            for (i=seen_vtx->n_curr; i<seen_vtx->n_next; i++) {
                // loop over the current layer: seen_vtx[i]
                

                for (j=0; j<seen_vtx->vtx[i]->neigh_no; j++) {
                    //loop over their neighbors: seen_vtx[i]->neigh[j]

                    // is this neighbor is not vicsek, skip to next neighbor
                    if (!(seen_vtx->vtx[i]->neigh[j]->type&is_vicsek_vtx)) continue;
                    //else{ rest of the j loop }

                    // has this neighbor been seen?
                    if (is_in_seen_vertex(seen_vtx,seen_vtx->vtx[i]->neigh[j]) && (curr_dist!=1)){
                        // if seen, skip to next neighbor in the j loop
                        // 1st layer is always new, so we can skip the check
                        continue;
                    }
                    else {
                        //new vertex! add to the next layer

                        add_vtx_to_seen(seen_vtx,seen_vtx->vtx[i]->neigh[j]);

                        //add normal to the sum
                        // Vicsek model 2: weight by 1/distance
                        if (vesicle->tape->vicsek_model == 2) {
                            vixnorm += vesicle->tape->vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nx / curr_dist;
                            viynorm += vesicle->tape->vicsek_strength * seen_vtx->vtx[i]->neigh[j]->ny / curr_dist;
                            viznorm += vesicle->tape->vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nz / curr_dist;
                        }
                        else {
                            vixnorm += vesicle->tape->vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nx;
                            viynorm += vesicle->tape->vicsek_strength * seen_vtx->vtx[i]->neigh[j]->ny;
                            viznorm += vesicle->tape->vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nz;
                        }

                    }
                }
            }
            // finished this layer
            //if did not add any vertices: completed cluster
            if (seen_vtx->n_next==seen_vtx->n_top) break;
        }

        //having finished summing the normals, normalize the resulting vector
        norml=sqrt(vixnorm*vixnorm+viynorm*viynorm+viznorm*viznorm);
        vixnorm/=norml;
	    viynorm/=norml;
	    viznorm/=norml;

	    /*calculate ddp, Viscek force directed displacement*/
	    vtx->fx=vixnorm;
        vtx->fy=viynorm;
        vtx->fz=viznorm;

	    //don't forget to free! 
        seen_vertex_free(seen_vtx);

    //end else from if (!Vicsek)
    }
    
    /*calculate ddp, normal force directed displacement*/
	ddp=vtx->fx*(vtx->x-vtx_old->x)+vtx->fy*(vtx->y-vtx_old->y)+vtx->fz*(vtx->z-vtx_old->z);
    
    /*calculate dE*/
    return -vtx->f*ddp;
//end function
}

ts_double direct_force_from_Fz_balance(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
    // calculate the Fz force balance
    //     We could separat from direct_force_energy()
    //     since F_direct(old) = F_direct + Fz_balance \hat{z}
    //     W = F_direct * dX + Fz_balance * dz
    // Force is cached once and only once! : not recalculated!
    // Changing this behavior is weird: once, or always, work identically, 
    // once every 2 works identically, once every 64 works identically, 
    // but 65-128 don't!.
    ts_double f_diff=0;
    static ts_double Fz;
    static ts_char countdown=0;

    if (!countdown){
        Fz=total_force_on_vesicle(vesicle);
        countdown=64;//pow(2,20)-1; go back here in 64 steps
    }
    else{
        if (vtx->type&is_active_vtx) f_diff+=vtx->fz;
        if (vtx_old->type&is_active_vtx) f_diff-=vtx_old->fz;
        Fz+=f_diff;
        countdown-=1; //recalculate when reach 0
    }

    if (Fz>0){
	    return Fz*(vtx->z-vtx_old->z)/vesicle->vlist->n;
    }
    else{
        return 0;
    }	
	
}

inline ts_double total_force_on_vesicle(ts_vesicle *vesicle){
	ts_uint i;
	ts_double fz=0;
	/*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
    for (i=0;i<vesicle->vlist->n;i++){
		if(vesicle->vlist->vtx[i]->type&is_active_vtx){
			fz+=vesicle->vlist->vtx[i]->fz;
		}
	}
    return fz;
	
}

ts_double adhesion_energy_diff(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
    ts_double delta_energy=0;
    ts_double z=vtx->z, z_old=vtx_old->z;
    ts_double z0=vesicle->tape->z_adhesion, dz=vesicle->tape->adhesion_cuttoff;
    ts_double c0=vesicle->adhesion_center, r=vesicle->tape->adhesion_radius;
    //1 for step potential
	if(vesicle->tape->type_of_adhesion_model==1){

		if( (vtx->type&is_adhesive_vtx) && (abs(z-z0)<dz) ){
				delta_energy-=vtx->ad_w;
		}
		if( (vtx_old->type&is_adhesive_vtx) &&(abs(z_old-z0)<dz) ){
				delta_energy+=vtx_old->ad_w;
		}
	}

    //2 for parabolic potential
	else if(vesicle->tape->type_of_adhesion_model==2){

        // can't combine them well: each has (theoretically) different adhesion
		if( (vtx->type&is_adhesive_vtx) && ((z-z0)<=dz )){
				delta_energy-=(vtx->ad_w/pow(dz,2))*pow(z - dz,2);
		}
		if( (vtx_old->type&is_adhesive_vtx) && ((z_old-z0)>dz) ){
				delta_energy+=(vtx_old->ad_w/pow(dz,2))*pow(z_old - dz,2);
		}
	}

    //3 for sphrerical adhesion substrate with constant potential
	else if(vesicle->tape->type_of_adhesion_model==3){
		if( (vtx->type&is_adhesive_vtx) && (pow(pow(c0-z,2) + pow(vtx->x,2) + pow(vtx->y,2),0.5) - r < dz)){
			delta_energy-=vtx->ad_w;
		}
		if( (vtx_old->type&is_adhesive_vtx) && (pow(pow(c0-z_old,2) + pow(vtx_old->x,2) + pow(vtx_old->y,2),0.5) - r < dz)){
			delta_energy+=vtx_old->ad_w;
		}
	}
    
    //4 for cylindrical adhesive substrate with constant potential
	else if(vesicle->tape->type_of_adhesion_model==4){
		if( (vtx->type&is_adhesive_vtx) && (pow(pow(c0 -z,2) + pow(vtx->x,2),0.5) - r < dz)){
			delta_energy-=vtx->ad_w;
		}
		if( (vtx_old->type&is_adhesive_vtx) && (pow(pow(c0 - z_old,2) + pow(vtx_old->x,2),0.5) - r < dz)){
			delta_energy+=vtx_old->ad_w;
		}
	}

    return delta_energy;
}

void stretchenergy(ts_vesicle *vesicle, ts_triangle *triangle){
	triangle->energy=vesicle->tape->xkA0/2.0*pow((triangle->area/vesicle->tlist->a0-1.0),2);
}

/* vim: set ts=4 sts=4 sw=4 noet : */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "general.h"
#include "energy.h"
#include "vertex.h"
#include "triangle.h"

// DEBUG
#include "io.h"
#include <string.h>


/** @brief Wrapper that calculates and saves energy of every vertex in vesicle
 *  
 *  Function calculated energy of every vertex in vesicle. It can be used in
 *  initialization procedure or in recalculation of the energy after non-MCsweep *  operations. However, when random move of vertex or flip of random bond occur *  call to this function is not necessary nor recommended. 
 *  @param *vesicle is a pointer to vesicle.
 *  @returns TS_SUCCESS on success.
*/
ts_bool sweep_vertex_curvature_energy(ts_vesicle *vesicle){
    ts_idx i;  
    ts_vertex_list *vlist=vesicle->vlist;
    ts_vertex **vtx=vlist->vtx;

    for(i=0;i<vlist->n;i++){
        if (!(vtx[i]->type==is_ghost_vtx)) {       
            vertex_curvature_energy(vesicle, vtx[i]);
        }
        
    }

    return TS_SUCCESS;
}

/** @brief Wrapper that calculate and saves force of every vertex in vesicle
 * 
 *  @param *vesicle is a pointer to vesicle.
 *  @returns TS_SUCCESS on success.
*/
ts_bool sweep_vertex_forces(ts_vesicle *vesicle){
    ts_idx i;
    ts_vertex_list *vlist=vesicle->vlist;
    ts_vertex **vtx=vlist->vtx;
    for(i=0;i<vlist->n;i++){
        if (!(vtx[i]->type==is_ghost_vtx)) {       
            direct_force_energy(vesicle,vtx[i],vtx[i]);
        }    
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


/** @brief Isotropic calculation of the bending energy of the vertex (subfunction of vertex_curvature_energy).
 *  
 *  Calculates energy of the vertex \f$i\f$. Function returns \f$\frac{1}{2}(c_1+c_2-c)^2 s\f$, where \f$(c_1+c_2)/2\f$ is mean curvature,
 * \f$c/2\f$ is spontaneous curvature and \f$s\f$ is area per vertex \f$i\f$.
 *
 * Nearest neighbors (NN) must be ordered in clockwise direction for this function to work.
 * First, NNs that form two neighboring triangles are found (\f$j_m\f$, \f$j_p\f$ and common \f$j\f$). Later, the scalar product of vectors \f$x_1=(\mathbf{i}-\mathbf{j_p})\cdot (\mathbf{i}-\mathbf{j_p})(\mathbf{i}-\mathbf{j_p})\f$, \f$x_2=(\mathbf{j}-\mathbf{j_p})\cdot  (\mathbf{j}-\mathbf{j_p})\f$  and \f$x_3=(\mathbf{j}-\mathbf{j_p})\cdot (\mathbf{i}-\mathbf{j_p})\f$  are calculated. From these three vectors the \f$c_{tp}=\frac{1}{\tan(\varphi_p)}\f$ is calculated, where \f$\varphi_p\f$ is the inner angle at vertex \f$j_p\f$. The procedure is repeated for \f$j_m\f$ instead of \f$j_p\f$ resulting in \f$c_{tn}\f$.
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
inline ts_bool laplace_beltrami_curvature_energy(ts_vesicle *vesicle, ts_vertex *vtx){
    ts_small_idx jj,jjm; // current, next (p) and prev (m) neighbor idx
    ts_vertex *j, *jm; // current and previous neighbor
    ts_triangle *tm, *tp; // triangle (vtx,jm,j) and (vtx,j,jp)
    ts_double s=0.0; // area
    ts_double xh=0.0,yh=0.0,zh=0.0,txn=0.0,tyn=0.0,tzn=0.0;
    ts_double lx, ly, lz, sigx, sigy, sigz, sigl, l_sqr;
    ts_double h,ht,norml;
    ts_double angle_sum=0;
    ts_double a_dot_b, a_cross_b_x, a_cross_b_y, a_cross_b_z, mag_a_cross_b;
    ts_flag model=vesicle->tape->curvature_model; // control how and what model we use to calculate energy: see enum curvature_model_type in general.h
    ts_bool do_angle_sum=0;
    do_angle_sum=(model&to_calculate_sum_angle 
                  && (!(model&to_use_sum_angle_for_kx2_only) || (model&to_use_sum_angle_for_kx2_only && vtx->xk2!=0))
                 );

    // we have 4 steps:
    // 1: iterate neighbors: calculate things on the dual lattice (+ normal + angle sum + area)
    // 2: get mean curvature using the dual lattice laplace-beltrami formula
    // 3: calculate gaussian curvature
    // 4: calcualte the bending energy

    // step 1. iterate over the neighbors
    // - calculate the normal
    // - calculate the dual lattice edge contribution to the curvature xh voodoo magic vector
    // - calculate the angle sum for the gaussian curvature
    for(jj=0; jj<vtx->neigh_no;jj++){
        jjm=prev_small(jj, vtx->neigh_no);
        j=vtx->neigh[jj];
        jm=vtx->neigh[jjm];

        // step 1.1 calculate the contribution to the vertex normal
        // txn normal points inwards!!
        tp=vtx->tristar[jj]; 
        tm=vtx->tristar[jjm]; 
        txn+=tp->xnorm;
        tyn+=tp->ynorm;
        tzn+=tp->znorm;
    
        // step 1.2 calculate cotangent of the edge, dual lattice edge triangles (vertex, edge middle, circumcenter-m), (vertex, edge middle, circumcenter-p)
        // These have the same angle as the opposing angle (half of a central angle = inscribed angle)
        lx = (j->x-vtx->x);
        ly = (j->y-vtx->y); // edge vector
        lz = (j->z-vtx->z);
        l_sqr = lx*lx + ly*ly + lz*lz; // half edge square
        sigx = tm->xcirc - (j->x+vtx->x)/2;
        sigy = tm->ycirc - (j->y+vtx->y)/2;
        sigz = tm->zcirc - (j->z+vtx->z)/2;
        sigl = lx*(sigy*tm->znorm - sigz*tm->ynorm) + ly*(sigz*tm->xnorm - sigx*tm->znorm) + lz*(sigx*tm->ynorm - sigy*tm->xnorm); // l*(sigxN)
        sigl *= -1; // the jm section is left handed
        s += 0.25*sigl; // A = 1/2 (sigma * l/1) = 1/4 sigl
        xh += sigl*(lx)/l_sqr;
        yh += sigl*(ly)/l_sqr;
        zh += sigl*(lz)/l_sqr;
        sigx = tp->xcirc - (j->x+vtx->x)/2;
        sigy = tp->ycirc - (j->y+vtx->y)/2;
        sigz = tp->zcirc - (j->z+vtx->z)/2;
        sigl = lx*(sigy*tp->znorm - sigz*tp->ynorm) + ly*(sigz*tp->xnorm - sigx*tp->znorm) + lz*(sigx*tp->ynorm - sigy*tp->xnorm); // sigma * l
        s += 0.25*sigl; // A = 1/2 (sigma * l/1) = 1/4 sigl
        xh += sigl*(lx)/l_sqr;
        yh += sigl*(ly)/l_sqr;
        zh += sigl*(lz)/l_sqr;

        // step 1.4 angle calculation for gaussian curvature
        // angle_sum += atan(ctp) + atan(ctm); // simple but slow!
        // instead: get the angle m-vtx-j
        // atan2(|axb|,a*b) was recommended at mathwork forum (cosin has small angle problems, and still need a sqrt)
        // possibly more complicated but better one from linked pdf (kahan)
        // TODO: sheck if there is a good formula using sigl/l^2  (maybe tan(a+b) = tan(a)+tan(b)/1-tan(a)tan(b), and we get tans from sigl1/l^2 sigl2*l2)
        if (do_angle_sum){
            a_dot_b = (jm->x-vtx->x)*(j->x-vtx->x)+
                        (jm->y-vtx->y)*(j->y-vtx->y)+
                        (jm->z-vtx->z)*(j->z-vtx->z);
            a_cross_b_x = (jm->y-vtx->y)*(j->z-vtx->z)-(jm->z-vtx->z)*(j->y-vtx->y);
            a_cross_b_y = (jm->z-vtx->z)*(j->x-vtx->x)-(jm->x-vtx->x)*(j->z-vtx->z);
            a_cross_b_z = (jm->x-vtx->x)*(j->y-vtx->y)-(jm->y-vtx->y)*(j->x-vtx->x);
            mag_a_cross_b = sqrt(pow(a_cross_b_x,2)+pow(a_cross_b_y,2)+pow(a_cross_b_z,2));
            angle_sum += atan2(mag_a_cross_b, a_dot_b);
        }

    } // end for jj neighbors

    // step 2 calculate the curvatures from the xh voodoo formula
    // xh voodoo magic vector has the mean curvature times area as the magnitude and a direction roughly(?) towards the center of curvature
    h=xh*xh + yh*yh + zh*zh; 
    ht=txn*xh + tyn*yh + tzn*zh; // direction of center of curvature with the normal i.e. convex or concave
    if(ht>=0.0) {
        vtx->mean_curvature= (sqrt(h)/s)/2;
    } else {
        vtx->mean_curvature=-(sqrt(h)/s)/2;
    }
    norml=sqrt(txn*txn+tyn*tyn+tzn*tzn); //also great point to update the vertex normal
    vtx->nx=-txn/norml;
    vtx->ny=-tyn/norml; // !!the triangle normal points inwards!!
    vtx->nz=-tzn/norml;
    
    // step 3 gaussian curvature and update for anisotropic vertex

    // step 3.1 calculate gaussian curvature using sum angle formula
    if (do_angle_sum){
        vtx->gaussian_curvature = (2*M_PI- angle_sum)/s;
    }

    //step 4: calculate the bending energy
    // step 4.1 isotropic curvature energies
    /* the following statement is an expression for $\frac{1}{2}\int(c_1+c_2-c_0^\prime)^2\mathrm{d}A$, where $c_0^\prime=2c_0$ (twice the spontaneous curvature)  */
    vtx->mean_energy=vtx->xk* 0.5*s*(2*vtx->mean_curvature-vtx->c)*(2*vtx->mean_curvature-vtx->c);
    vtx->gaussian_energy = do_angle_sum ? vtx->xk2 * s * vtx->gaussian_curvature : 0;
    vtx->energy = vtx->mean_energy + vtx->gaussian_energy;
    

    return TS_SUCCESS;
}

// Function that update the vertex with the tensor
// make sure that the curvature is on the right side!
ts_bool update_vertex_from_curvature_tensor(ts_vertex* vtx, ts_double Av,
                             ts_double s00,ts_double s01,ts_double s10,ts_double s11,
                             ts_double nx,ts_double ny,ts_double nz,
                             ts_double dx,ts_double dy,ts_double dz,
                             ts_double tx,ts_double ty,ts_double tz){

    //get curvature information from the 2x2 shape operator
    // 1. get eigenvectors in the director-tangent directions
    // 2.: get eigenvectors in real space
    // 3.: get energy information from the 2x2 shape difference
    ts_double dSd,tSt,dSt,tSd;
    ts_double tr,det,discrim_sqrt, lambda1,lambda2, length;
    ts_double eigen_vec1d,eigen_vec2d,eigen_vec1t,eigen_vec2t;
    // first, the primary matrix invariants, the trace and determinant
    dSd=s00;
    dSt=s01;
    tSd=s10;
    tSt=s11;
    tr  = dSd + tSt;
    det = dSd*tSt - dSt*tSd;
    vtx->mean_curvature = tr/2; // curvatures up to signs and factors of 2
    vtx->gaussian_curvature = det;
    // eigenvalues: trace determinant formula. We have real symmetric matrix, so positive discriminant
    discrim_sqrt = sqrt(tr*tr - 4*det);
    lambda1 = (tr + discrim_sqrt)/2;
    lambda2 = (tr - discrim_sqrt)/2;
    // ----------------------------------------------------
    // get eigenvectors in the director-tangent directions
    // ----------------------------------------------------
    // construct eigenvectors in the director-tangent plane: make sure we don't have any 0 vectors!
    // based on https://math.stackexchange.com/questions/4103294/is-there-a-closed-form-expression-for-the-eigenvectors-of-a-2x2-matrix
    //      for [ a b ], we use  v+ = [l1 - d]  v- =  [  -b  ]  vectors if a>d 
    //          [ b d ]               [   b  ],       [a - l2]    
    //                    and    v+ = [   b  ]   v- = [l2 - d]  if d>a
    //                                [l1 - a],       [   b  ]
    //      with the v+, v- making a left handed xy plane 
    //      playing in Mathematica shows switching these combinations lead to larger eigenvectors
    //      in the [ a 0 ]  b=0 case, this also assures we have ~[1,0] and ~[0,1] vectors
    //             [ 0 d ]
    //      ( since lambda1 = max[a,d] and lambda2 = min[a,d])
    //      for dSd==tSt and tSd==dSt!=0 i.e. [ [a,b],[b,a] ], mathematica shows both work
    //      the only problematic case is dSd==tSt and tSd==dSt==0 [[a,0],[0,a]]=aI, which is the degenerate case
    //
    // To order the eigenvalues by magnitude and keep the orientation, we can always take v+ and construct v- out of it
    // v- = [ -v+[1] , 
    //         v+[0]  ]
    // and only later 

    // a. degenerate case
    if(lambda1==lambda2){
        eigen_vec1d = 1; // we pick the director and tangent vectors as the eigenvectors
        eigen_vec1t = 0;
    }
    else{
        // b. nondegenerate case
        if(dSd>=tSt){  // a>d
            eigen_vec1d = lambda1-tSt;  //tSd==0 -> lambda1=dSd, eigen_vec1d!=0
            eigen_vec1t = tSd;
        }
        else { // dSd<tSt, d<a
            eigen_vec1d = dSt;
            eigen_vec1t = lambda1-dSd;; //tSd==0 -> lambda1=tSt, eigen_vec1t!=0
        }
        // normalize the eigenvectors
        length = sqrt(eigen_vec1d*eigen_vec1d + eigen_vec1t*eigen_vec1t);
        eigen_vec1d/=length;
        eigen_vec1t/=length;
    }

    // -----------------------------------------------
    // get eigenvalues and eigenvectors in real space
    // -----------------------------------------------

    //sort eigenvalues by absolute value
    if(abs(lambda1)>=abs(lambda2)){
        vtx->eig_v0 = lambda1;
        vtx->eig_v1 = lambda2;

        // lambda2's eigenvector is second
        eigen_vec2d = -eigen_vec1t;
        eigen_vec2t = eigen_vec1d;

    } else {
        vtx->eig_v0 = lambda2;
        vtx->eig_v1 = lambda1;

        // lambda2's eigenvector is first
        eigen_vec2d = eigen_vec1d;
        eigen_vec2t = eigen_vec1t;
        eigen_vec1d = eigen_vec2t;
        eigen_vec1t = -eigen_vec2d;

    }
    vtx->eig0[0] = eigen_vec1d*dx + eigen_vec1t*tx;
    vtx->eig0[1] = eigen_vec1d*dy + eigen_vec1t*ty;
    vtx->eig0[2] = eigen_vec1d*dz + eigen_vec1t*tz;

    vtx->eig1[0] = eigen_vec2d*dx + eigen_vec2t*tx;
    vtx->eig1[1] = eigen_vec2d*dy + eigen_vec2t*ty;
    vtx->eig1[2] = eigen_vec2d*dz + eigen_vec2t*tz;

    vtx->eig_v2 = 0; // The shape operator is annihilated in the normal direction
    
    vtx->eig2[0] = nx;
    vtx->eig2[1] = ny;
    vtx->eig2[2] = nz;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // get energy information from the 2x2 shape difference
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // spontaneous curvature and isotropy go here!
    dSd -= 0.5*(vtx->c + vtx->d);
    tSt -= 0.5*(vtx->c - vtx->d);
    tr = dSd + tSt;
    det = dSd*tSt - tSd * dSt;
    vtx->mean_energy = vtx->xk*Av* pow(tr,2)/2;
    vtx->gaussian_energy = vtx->xk2 * Av * det;
    vtx->energy = vtx->mean_energy + vtx->gaussian_energy;

    return TS_SUCCESS;
}


ts_bool error_correction_scheme(ts_double *s00, ts_double *s01, ts_double *s11, ts_double H, ts_double Kg){
    ts_double delta=0,a=0,beta=0;
    ts_double a_d=0;
    ts_double h,d, kg;
    h = (*s00+*s11)/2;
    d = (*s00-*s11)/2;
    kg = *s00*(*s11)-(*s01)*(*s01);
    delta = H-h;
    a_d = d*d*(1-(Kg-kg+h*h-H*H)/(d*d+pow(*s01,2)));
    // a_d = a_d>0? sqrt(a_d) : 0;
    a = (abs(a_d-d)<=abs(-a_d-d))? a_d-d : -a_d-d;
    beta = a*(*s01)/d;
    *s00 += delta+a;
    *s01 += beta;
    *s11 += delta-a;
    return TS_SUCCESS;
}

/** @brief anisotropic subfunction of vertex_curvature_energy.
 *
 *  This function is experimental, branch from vertex_energy for anisotropic proteins 
 *  to calculate the tensor-based curvature values c1,c2 and principle directions,
 *  calculate the energy, and save it all on the vertex
 * 
 *  vertex must remain ordered through initial_dist and through bondflips
 *
 *  @returns TS_SUCCESS on successful calculation
*/
inline ts_bool tensor_curvature_energy2(ts_vesicle *vesicle, ts_vertex *vtx){
    //  step 1. calculate the area assigned to the vertex and the vertex normal
    //      step 1.1: update vertex normal and director, create normal-director-tangent frame
    //  step 2. calculate and accumulate the shape operator triangle
    //  step 3. project to a 2x2 matrix in the surface tangent plane, obtain energy and curvature.

    // we hardcoded 10 neighbor limit!
    ts_small_idx i, ip;

    ts_double edge_e_x[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_e_y[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_e_z[10]={0,0,0,0,0,0,0,0,0,0};

    ts_double edge_middle_x[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_middle_y[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_middle_z[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_len[10]={0,0,0,0,0,0,0,0,0,0};

    ts_double edge_normal_x[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_normal_y[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_normal_z[10]={0,0,0,0,0,0,0,0,0,0};



    ts_double sigma_L_x[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double sigma_L_y[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double sigma_L_z[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double sigma_L_len[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double area_L[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double sigma_R_x[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double sigma_R_y[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double sigma_R_z[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double sigma_R_len[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double area_R[10]={0,0,0,0,0,0,0,0,0,0};

    ts_double vertex_normal_x=0.0;
    ts_double vertex_normal_y=0.0;
    ts_double vertex_normal_z=0.0;
    ts_double director_x=0.0; // vertex director
    ts_double director_y=0.0;
    ts_double director_z=0.0;
    ts_double tangent_x=0.0; // tangent vector in director-tangent-normal axis
    ts_double tangent_y=0.0;
    ts_double tangent_z=0.0;

    ts_double temp_length; // length of a vector for normalization
    ts_double dot_product; 

    ts_double Ne_sigma, Ni_sigma, Ni_e;

    ts_double Se11=0, Se21=0, Se22=0; // Edgewise shape operator
    ts_double Se12; //alias for clarity of symmetric matrices: hopefully the compiler removes them
    ts_double dSd,dSt,tSt,tSd; // vertex shape operator in the director-tangent plane

    ts_triangle *t, *tp;
    ts_double s;

    ts_double Sv11=0, Sv12=0, Sv13=0, Sv21=0, Sv22=0, Sv23=0, Sv31=0, Sv32=0, Sv33=0;

    // ################################################
    // step 1. prepare all vectors and arrays we need #
    // ################################################

    /* We have a vertex, and the edges are ordered clockwise
    
        i---------------------ip
        \ \                 / /  \
         \   \     t      /  /     \
          \    \        /   /         \
           \     \   /     /            \
            \     _O_     /               \
             \ _-^ | ^-_ /        tp        \
              \    |    /                     \
               \ L | S /                        \
          e[i]  \  |  /  e[ip]                    \
                 \ | /                      ____-ipp
                  \|/          ____-----^^^^
                  vtx -----^^^^
    
    
    */

    ts_double Av=0;

    // edge vectors
    for(i=0; i<vtx->neigh_no; i++){
        edge_e_x[i] = (vtx->x - vtx->neigh[i]->x)/2;
        edge_e_y[i] = (vtx->y - vtx->neigh[i]->y)/2;
        edge_e_z[i] = (vtx->z - vtx->neigh[i]->z)/2;
        edge_middle_x[i] = (vtx->x + vtx->neigh[i]->x)/2;
        edge_middle_y[i] = (vtx->y + vtx->neigh[i]->y)/2;
        edge_middle_z[i] = (vtx->z + vtx->neigh[i]->z)/2;
        edge_len[i] = sqrt(pow(edge_e_x[i],2) + pow(edge_e_y[i],2) + pow(edge_e_z[i],2));
        edge_e_x[i] /= edge_len[i];
        edge_e_y[i] /= edge_len[i];
        edge_e_z[i] /= edge_len[i];
    }

    // sigma and vertex normal
    // notice everything is ordered clockwise from outside perspective
    for (i=0; i<vtx->tristar_no; i++) {
        ip = next_small(i,vtx->neigh_no);
        t = vtx->tristar[i];
        sigma_L_x[i] = t->xcirc - edge_middle_x[i];
        sigma_L_y[i] = t->ycirc - edge_middle_y[i];
        sigma_L_z[i] = t->zcirc - edge_middle_z[i];
        sigma_L_len[i] = sqrt(pow(sigma_L_x[i],2)+pow(sigma_L_y[i],2)+pow(sigma_L_z[i],2));
        sigma_L_x[i] /= sigma_L_len[i];
        sigma_L_y[i] /= sigma_L_len[i];
        sigma_L_z[i] /= sigma_L_len[i];
        sigma_R_x[i] = t->xcirc - edge_middle_x[ip];
        sigma_R_y[i] = t->ycirc - edge_middle_y[ip];
        sigma_R_z[i] = t->zcirc - edge_middle_z[ip];
        sigma_R_len[i] = sqrt(pow(sigma_R_x[i],2)+pow(sigma_R_y[i],2)+pow(sigma_R_z[i],2));
        sigma_R_x[i] /= sigma_R_len[i];
        sigma_R_y[i] /= sigma_R_len[i];
        sigma_R_z[i] /= sigma_R_len[i];

        if (sigma_L_x[i]*edge_e_x[i]+sigma_L_y[i]*edge_e_y[i]+sigma_L_z[i]*edge_e_z[i]>0.00001){
            ts_fprintf(stdout,"At vertex %d neighbor %d",vtx->idx, i);
            ts_fprintf(stdout,"Sigma L: %f,%f,%f e: %f,%f,%f\n",sigma_L_x[i],sigma_L_y[i],sigma_L_z[i],
                                                              edge_e_x[i],edge_e_y[i],edge_e_z[i]);
            strcat(command_line_args.dump_fullfilename,"ary_emergency.bin");
            dump_state(vesicle,0);
            fatal("bad sigma calculation\n",231);
        }
        if (sigma_R_x[i]*edge_e_x[ip]+sigma_R_y[i]*edge_e_y[ip]+sigma_R_z[i]*edge_e_z[ip]>0.00001){
            ts_fprintf(stdout,"At vertex %d neighbor %d",vtx->idx, i);
            ts_fprintf(stdout,"Sigma R: %f,%f,%f e: %f,%f,%f\n",sigma_R_x[i],sigma_R_y[i],sigma_R_z[i],
                                                              edge_e_x[ip],edge_e_y[ip],edge_e_z[ip]);
            strcat(command_line_args.dump_fullfilename,"ary_emergency.bin");
            dump_state(vesicle,0);
            fatal("bad sigma calculation\n",231);
        }

        area_L[i] = 0.5 * sigma_L_len[i]*edge_len[i];
        area_R[i] = 0.5 * sigma_R_len[i]*edge_len[ip]; 
        s = area_R[i] + area_L[i];
        vertex_normal_x -= t->xnorm*s; // t->norm points inwards!
        vertex_normal_y -= t->ynorm*s;
        vertex_normal_z -= t->znorm*s;
        Av += s;
    }
    //DEBUG
    vtx->area = Av;
    // edge normals
    if (! (vtx->type & is_edge_vtx) || vtx->neigh_no==vtx->tristar_no) {
        for (i=0; i<vtx->neigh_no; i++) {
            ip = next_small(i,vtx->neigh_no);
            t = vtx->tristar[i];
            tp = vtx->tristar[ip];
            edge_normal_x[ip] = t->xnorm + tp->xnorm;
            edge_normal_y[ip] = t->ynorm + tp->ynorm;
            edge_normal_z[ip] = t->znorm + tp->znorm;
            temp_length = sqrt(pow(edge_normal_x[ip],2) + pow(edge_normal_y[ip],2) + pow(edge_normal_z[ip],2));
            edge_normal_x[ip] /= temp_length;
            edge_normal_y[ip] /= temp_length;
            edge_normal_z[ip] /= temp_length;
        }
    } else if ((vtx->type & is_edge_vtx) && vtx->neigh_no==vtx->tristar_no+1) {
        // if we have an edge: pretend the first and final edges go right out
        edge_normal_x[0] = -sigma_L_x[0];
        edge_normal_y[0] = -sigma_L_y[0];
        edge_normal_z[0] = -sigma_L_z[0];
        temp_length = sqrt(pow(edge_normal_x[0],2) + pow(edge_normal_y[0],2) + pow(edge_normal_z[0],2));
        edge_normal_x[0] /= temp_length;
        edge_normal_y[0] /= temp_length;
        edge_normal_z[0] /= temp_length;
        for (i=0; i<vtx->neigh_no-1; i++) {
            ip = next_small(i,vtx->neigh_no);
            t = vtx->tristar[i];
            tp = vtx->tristar[ip];
            edge_normal_x[ip] = t->xnorm + tp->xnorm;
            edge_normal_y[ip] = t->ynorm + tp->ynorm;
            edge_normal_z[ip] = t->znorm + tp->znorm;
            temp_length = sqrt(pow(edge_normal_x[ip],2) + pow(edge_normal_y[ip],2) + pow(edge_normal_z[ip],2));
            edge_normal_x[ip] /= temp_length;
            edge_normal_y[ip] /= temp_length;
            edge_normal_z[ip] /= temp_length;
            }
        i = vtx->tristar_no;
        ip = vtx->neigh_no;
        edge_normal_x[ip] = -sigma_R_x[i];
        edge_normal_y[ip] = -sigma_R_y[i];
        edge_normal_z[ip] = -sigma_R_z[i];
        temp_length = sqrt(pow(edge_normal_x[ip],2) + pow(edge_normal_y[ip],2) + pow(edge_normal_z[ip],2));
        edge_normal_x[ip] /= temp_length;
        edge_normal_y[ip] /= temp_length;
        edge_normal_z[ip] /= temp_length;
    } else {
        return TS_FAIL; // no idea how to handle this situation.
    }

    temp_length=sqrt(pow(vertex_normal_x,2)+pow(vertex_normal_y,2)+pow(vertex_normal_z,2));
    vertex_normal_x/=temp_length;
    vertex_normal_y/=temp_length;
    vertex_normal_z/=temp_length;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 1.1: update vertex normal and director based, generates the director-tangent-normal frame
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vtx->nx = vertex_normal_x;
    vtx->ny = vertex_normal_y;
    vtx->nz = vertex_normal_z;

    // d = d - (d.n)n   [or -nx(nxd)]
    dot_product = (vtx->dx*vertex_normal_x)+(vtx->dy*vertex_normal_y)+(vtx->dz*vertex_normal_z);
    director_x = vtx->dx - dot_product*vertex_normal_x;
    director_y = vtx->dy - dot_product*vertex_normal_y;
    director_z = vtx->dz - dot_product*vertex_normal_z;
    temp_length=sqrt((director_x*director_x)+(director_y*director_y)+(director_z*director_z)); // should be the same as sqrt(1-*(normal . director)^2)
    director_x = director_x / temp_length;
    director_y = director_y / temp_length;
    director_z = director_z / temp_length;
    vtx->dx = director_x; // update the vertex director
    vtx->dy = director_y;
    vtx->dz = director_z;
    // calculate the third axis tangent = normal x director
    tangent_x = vertex_normal_y*director_z - vertex_normal_z*director_y;
    tangent_y = vertex_normal_z*director_x - vertex_normal_x*director_z;
    tangent_z = vertex_normal_x*director_y - vertex_normal_y*director_x;

    
    // ###############################################
    // step 2. calculate the shape operator per edge #
    // ###############################################
    for(i=0;i<vtx->tristar_no;i++){
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.1: compute curvature contribution on the left 
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // The main equation is
        // \begin{equation} 
        //     S_{ab} = \left(\begin{matrix}
        //         \frac{1}{2}\left|e\right|\vec{N}_e\cdot\hat{\sigma} 
        //         & \frac{1}{4}\left|\sigma\right|\left(\vec{N}_e\cdot\hat{\sigma}-\vec{N}_i\cdot\hat{\sigma}\right) 
        //         \\ \frac{1}{4}\left|\sigma\right|\left(\vec{N}_e\cdot\hat{\sigma}-\vec{N}_i\cdot\hat{\sigma}\right) 
        //         & -\frac{1}{2}\left|\sigma\right|\vec{N}_i\cdot\hat{e}
        //     \end{matrix}\right)
        // \end{equation}
        // note that we didn't normalize the sigmas and es, so we need to modify our equation with "neat tricks", using the area 0.5|s||e| and the squares

        Ne_sigma = edge_normal_x[i]*sigma_L_x[i] + edge_normal_y[i]*sigma_L_y[i] + edge_normal_z[i]*sigma_L_z[i];
        Ni_sigma = vertex_normal_x*sigma_L_x[i]  + vertex_normal_y*sigma_L_y[i]  + vertex_normal_z*sigma_L_z[i];
        Ni_e = vertex_normal_x*edge_e_x[i]       + vertex_normal_y*edge_e_y[i]   + vertex_normal_z*edge_e_z[i];
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.5: get the edge shape operator Se and edge 
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Se11=0.5*edge_len[i]*Ne_sigma;
        Se21=0.25*sigma_L_len[i]*(Ne_sigma-Ni_sigma);
        Se22=-0.5*sigma_L_len[i]*Ni_e;
        Se12=Se21; //symmetric matrix: hopefully the compiler gets rid of these

        Sv11 += Se11*sigma_L_x[i]*sigma_L_x[i] + Se12*sigma_L_x[i]*edge_e_x[i] + Se21*edge_e_x[i]*sigma_L_x[i] + Se22*edge_e_x[i]*edge_e_x[i];
        Sv12 += Se11*sigma_L_x[i]*sigma_L_y[i] + Se12*sigma_L_x[i]*edge_e_y[i] + Se21*edge_e_x[i]*sigma_L_y[i] + Se22*edge_e_x[i]*edge_e_y[i];
        Sv13 += Se11*sigma_L_x[i]*sigma_L_z[i] + Se12*sigma_L_x[i]*edge_e_z[i] + Se21*edge_e_x[i]*sigma_L_z[i] + Se22*edge_e_x[i]*edge_e_z[i];
        Sv21 += Se11*sigma_L_y[i]*sigma_L_x[i] + Se12*sigma_L_y[i]*edge_e_x[i] + Se21*edge_e_y[i]*sigma_L_x[i] + Se22*edge_e_y[i]*edge_e_x[i];
        Sv22 += Se11*sigma_L_y[i]*sigma_L_y[i] + Se12*sigma_L_y[i]*edge_e_y[i] + Se21*edge_e_y[i]*sigma_L_y[i] + Se22*edge_e_y[i]*edge_e_y[i];
        Sv23 += Se11*sigma_L_y[i]*sigma_L_z[i] + Se12*sigma_L_y[i]*edge_e_z[i] + Se21*edge_e_y[i]*sigma_L_z[i] + Se22*edge_e_y[i]*edge_e_z[i];
        Sv31 += Se11*sigma_L_z[i]*sigma_L_x[i] + Se12*sigma_L_z[i]*edge_e_x[i] + Se21*edge_e_z[i]*sigma_L_x[i] + Se22*edge_e_z[i]*edge_e_x[i];
        Sv32 += Se11*sigma_L_z[i]*sigma_L_y[i] + Se12*sigma_L_z[i]*edge_e_y[i] + Se21*edge_e_z[i]*sigma_L_y[i] + Se22*edge_e_z[i]*edge_e_y[i];
        Sv33 += Se11*sigma_L_z[i]*sigma_L_z[i] + Se12*sigma_L_z[i]*edge_e_z[i] + Se21*edge_e_z[i]*sigma_L_z[i] + Se22*edge_e_z[i]*edge_e_z[i];

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.2: compute curvature contribution on the right
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ip = next_small(i,vtx->neigh_no);
        Ne_sigma = edge_normal_x[ip]*sigma_R_x[i] + edge_normal_y[ip]*sigma_R_y[i] + edge_normal_z[ip]*sigma_R_z[i];
        Ni_sigma = vertex_normal_x*sigma_R_x[i]   + vertex_normal_y*sigma_R_y[i]   + vertex_normal_z*sigma_R_z[i];
        Ni_e = vertex_normal_x*edge_e_x[ip]       + vertex_normal_y*edge_e_y[ip]   + vertex_normal_z*edge_e_z[ip];
        

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.5: get the edge shape operator Se and edge 
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Se11=0.5*edge_len[ip]*Ne_sigma;
        Se21=0.25*sigma_R_len[i]*(Ne_sigma-Ni_sigma);
        Se22=-0.5*sigma_R_len[i]*Ni_e;
        Se12=Se21; //symmetric matrix: hopefully the compiler gets rid of these

        Sv11 += Se11*sigma_R_x[i]*sigma_R_x[i] + Se12*sigma_R_x[i]*edge_e_x[ip] + Se21*edge_e_x[ip]*sigma_R_x[i] + Se22*edge_e_x[ip]*edge_e_x[ip];
        Sv12 += Se11*sigma_R_x[i]*sigma_R_y[i] + Se12*sigma_R_x[i]*edge_e_y[ip] + Se21*edge_e_x[ip]*sigma_R_y[i] + Se22*edge_e_x[ip]*edge_e_y[ip];
        Sv13 += Se11*sigma_R_x[i]*sigma_R_z[i] + Se12*sigma_R_x[i]*edge_e_z[ip] + Se21*edge_e_x[ip]*sigma_R_z[i] + Se22*edge_e_x[ip]*edge_e_z[ip];
        Sv21 += Se11*sigma_R_y[i]*sigma_R_x[i] + Se12*sigma_R_y[i]*edge_e_x[ip] + Se21*edge_e_y[ip]*sigma_R_x[i] + Se22*edge_e_y[ip]*edge_e_x[ip];
        Sv22 += Se11*sigma_R_y[i]*sigma_R_y[i] + Se12*sigma_R_y[i]*edge_e_y[ip] + Se21*edge_e_y[ip]*sigma_R_y[i] + Se22*edge_e_y[ip]*edge_e_y[ip];
        Sv23 += Se11*sigma_R_y[i]*sigma_R_z[i] + Se12*sigma_R_y[i]*edge_e_z[ip] + Se21*edge_e_y[ip]*sigma_R_z[i] + Se22*edge_e_y[ip]*edge_e_z[ip];
        Sv31 += Se11*sigma_R_z[i]*sigma_R_x[i] + Se12*sigma_R_z[i]*edge_e_x[ip] + Se21*edge_e_z[ip]*sigma_R_x[i] + Se22*edge_e_z[ip]*edge_e_x[ip];
        Sv32 += Se11*sigma_R_z[i]*sigma_R_y[i] + Se12*sigma_R_z[i]*edge_e_y[ip] + Se21*edge_e_z[ip]*sigma_R_y[i] + Se22*edge_e_z[ip]*edge_e_y[ip];
        Sv33 += Se11*sigma_R_z[i]*sigma_R_z[i] + Se12*sigma_R_z[i]*edge_e_z[ip] + Se21*edge_e_z[ip]*sigma_R_z[i] + Se22*edge_e_z[ip]*edge_e_z[ip];


    } // END FOR JJ


    // ##############################################################################
    // step 3. use vertex shape operator to obtain energy and curvature information #
    // ##############################################################################

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 3.1: project to a 2x2 matrix in the surface tangent plane
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // normalize the weighted sum
    Sv11/=Av;
    Sv12/=Av;
    Sv13/=Av;
    Sv21/=Av;
    Sv22/=Av;
    Sv23/=Av;
    Sv31/=Av;
    Sv32/=Av;
    Sv33/=Av;

    // project shape operator to the director-tangent plane
    // [a b] = [ dSd dSt ]
    // [b d]   [ tSd tSt ]
    dSd =   director_x*Sv11*director_x + director_x*Sv12*director_y + director_x*Sv13*director_z
           +director_y*Sv21*director_x + director_y*Sv22*director_y + director_y*Sv23*director_z
           +director_z*Sv31*director_x + director_z*Sv32*director_y + director_z*Sv33*director_z; 
    dSt =   director_x*Sv11*tangent_x  + director_x*Sv12*tangent_y  + director_x*Sv13*tangent_z
           +director_y*Sv21*tangent_x  + director_y*Sv22*tangent_y  + director_y*Sv23*tangent_z
           +director_z*Sv31*tangent_x  + director_z*Sv32*tangent_y  + director_z*Sv33*tangent_z ; 
    tSt =   tangent_x *Sv11*tangent_x  + tangent_x *Sv12*tangent_y  + tangent_x *Sv13*tangent_z
           +tangent_y *Sv21*tangent_x  + tangent_y *Sv22*tangent_y  + tangent_y *Sv23*tangent_z
           +tangent_z *Sv31*tangent_x  + tangent_z *Sv32*tangent_y  + tangent_z *Sv33*tangent_z ;
    error_correction_scheme(&dSd,&dSt,&tSt,vtx->mean_curvature, vtx->gaussian_curvature);
    tSd = dSt; // symmetric tensor

    //DEBUG
    vtx->S[0] = dSd;
    vtx->S[1] = dSt;
    vtx->S[2] = tSd;
    vtx->S[3] = tSt;

    return update_vertex_from_curvature_tensor(vtx,Av,dSd,dSt,tSd,tSt,
        vertex_normal_x, vertex_normal_y, vertex_normal_z,
        director_x,director_y,director_z,
        tangent_x, tangent_y, tangent_z);
}



inline ts_bool tensor_curvature_energy(ts_vesicle *vesicle, ts_vertex *vtx){
    //  step 1. calculate the area assigned to the vertex and the vertex normal
    //      step 1.1: update vertex normal and director, create normal-director-tangent frame
    //  step 2. calculate and accumulate the shape operator per edge
    //      step 2.1: calculate the normalized edge vector and edge length
    //      step 2.2: get the edge adjacent triangles lm and lp
    //      step 2.3: get the edge normal and edge binormal
    //      step 2.4: get the dihedral curvature weight he[jj]
    //      step 2.5: get the edge shape operator Se and edge weight We_Av
    //      step 2.6: accumulate contribution to the vertex shape operator
    //  step 3. use vertex shape operator to obtain energy and curvature information
    //      step 3.1: project to a 2x2 matrix in the surface tangent plane
    //      step 3.2: get curvature information from the 2x2 shape operator
    //          step 3.2.1: get eigenvectors in the director-tangent directions
    //          step 3.2.2: get eigenvectors in real space
    //      step 3.3: get energy information from the 2x2 shape difference

    // we hardcoded 10 neighbor limit!
    ts_small_idx jj, i, ip;
    ts_double edge_vector_x[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_vector_y[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_vector_z[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_normal_x[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_normal_y[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_normal_z[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_binormal_x[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_binormal_y[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double edge_binormal_z[10]={0,0,0,0,0,0,0,0,0,0};
    ts_double vertex_normal_x=0.0;
    ts_double vertex_normal_y=0.0;
    ts_double vertex_normal_z=0.0;
    ts_double director_x=0.0; // vertex director
    ts_double director_y=0.0;
    ts_double director_z=0.0;
    ts_double tangent_x=0.0; // tangent vector in director-tangent-normal axis
    ts_double tangent_y=0.0;
    ts_double tangent_z=0.0;
    ts_double edge_length=0.0;

    ts_triangle *lm=NULL, *lp=NULL; // counter-clockwise and clockwise trianlges

    ts_double temp_length; // length of a vector for normalization
    ts_double dot_product; 
    ts_double cross_x, cross_y, cross_z; // temporary variable to hold (lm->normal) x (edge_normal)

    ts_double Se11=0, Se21=0, Se22=0, Se31=0, Se32=0, Se33=0; // Edgewise shape operator
    ts_double Se12, Se13, Se23; //alias for clarity of symmetric matrices: hopefully the compiler removes them
    ts_double dSd,dSt,tSt,tSd; // vertex shape operator in the director-tangent plane

    ts_triangle *t;
    ts_double s, l_m_x,l_m_y,l_m_z,l_p_x,l_p_y,l_p_z, sigma_m_x, sigma_m_y,sigma_m_z, sigma_p_x, sigma_p_y,sigma_p_z;

    ts_double We;
    ts_double Av, We_Av;

    ts_double he[10];
    ts_double Sv[3][3]={{0,0,0},{0,0,0},{0,0,0}};

    // #########################################################################
    // step 1. calculate the area assigned to the vertex and the vertex normal #
    // #########################################################################
    Av=0;
    for(i=0; i<vtx->tristar_no; i++){
        ip = next_small(i, vtx->tristar_no);
        t = vtx->tristar[i];
        l_m_x = (vtx->neigh[i]->x - vtx->x)/2;
        l_m_y = (vtx->neigh[i]->y - vtx->y)/2;
        l_m_z = (vtx->neigh[i]->z - vtx->z)/2;
        l_p_x = (vtx->neigh[ip]->x - vtx->x)/2;
        l_p_y = (vtx->neigh[ip]->y - vtx->y)/2;
        l_p_z = (vtx->neigh[ip]->z - vtx->z)/2;
        sigma_m_x = t->xcirc - (vtx->neigh[i]->x + vtx->x)/2;
        sigma_m_y = t->ycirc - (vtx->neigh[i]->y + vtx->y)/2;
        sigma_m_z = t->zcirc - (vtx->neigh[i]->z + vtx->z)/2;
        sigma_p_x = t->xcirc - (vtx->neigh[ip]->x + vtx->x)/2;
        sigma_p_y = t->ycirc - (vtx->neigh[ip]->y + vtx->y)/2;
        sigma_p_z = t->zcirc - (vtx->neigh[ip]->z + vtx->z)/2;
        // here we do N*(lxsigma) on the left and N*(lxsigma) on the right
        // which is N*(lxsigma - lxsigma)
        cross_x = -( l_p_y*sigma_p_z - l_p_z*sigma_p_y ) + ( l_m_y*sigma_m_z - l_m_z*sigma_m_y );
        cross_y = -( l_p_z*sigma_p_x - l_p_x*sigma_p_z ) + ( l_m_z*sigma_m_x - l_m_x*sigma_m_z );
        cross_z = -( l_p_x*sigma_p_y - l_p_y*sigma_p_x ) + ( l_m_x*sigma_m_y - l_m_y*sigma_m_x );
        s = 0.5 * (t->xnorm*cross_x + t->ynorm*cross_y + t->znorm*cross_z);
        vertex_normal_x -= t->xnorm*s;
        vertex_normal_y -= t->ynorm*s;
        vertex_normal_z -= t->znorm*s;
        Av += s;
    }
    //DEBUG
    vtx->area = Av;

    temp_length=sqrt(pow(vertex_normal_x,2)+pow(vertex_normal_y,2)+pow(vertex_normal_z,2));
    vertex_normal_x=vertex_normal_x/temp_length;
    vertex_normal_y=vertex_normal_y/temp_length;
    vertex_normal_z=vertex_normal_z/temp_length;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 1.1: update vertex normal and director based, generates the director-tangent-normal frame
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vtx->nx = vertex_normal_x;
    vtx->ny = vertex_normal_y;
    vtx->nz = vertex_normal_z;

    // d = d - (d.n)n   [or -nx(nxd)]
    dot_product = (vtx->dx*vertex_normal_x)+(vtx->dy*vertex_normal_y)+(vtx->dz*vertex_normal_z);
    director_x = vtx->dx - dot_product*vertex_normal_x;
    director_y = vtx->dy - dot_product*vertex_normal_y;
    director_z = vtx->dz - dot_product*vertex_normal_z;
    // !!! this operation exclusively lowers |t|: if we do manage to avoid using the size (always taking d*A*(nxd)/d^2) we can avoid the normalization
    // and just periodically make sure the size is large enough if (d^2<0.5) d=2d
    temp_length=sqrt((director_x*director_x)+(director_y*director_y)+(director_z*director_z)); // should be the same as sqrt(1-*(normal . director)^2)
    director_x = director_x / temp_length;
    director_y = director_y / temp_length;
    director_z = director_z / temp_length;
    vtx->dx = director_x; // update the vertex director
    vtx->dy = director_y;
    vtx->dz = director_z;
    // calculate the third axis tangent = normal x director
    tangent_x = vertex_normal_y*director_z - vertex_normal_z*director_y;
    tangent_y = vertex_normal_z*director_x - vertex_normal_x*director_z;
    tangent_z = vertex_normal_x*director_y - vertex_normal_y*director_x;

    
    // ###############################################
    // step 2. calculate the shape operator per edge #
    // ###############################################

    for(jj=0;jj<vtx->neigh_no;jj++){
        // !!! We start a VERY long loop over jj !!!
        // vertex must remain ordered through initial_dist and through bondflips

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.1: calculate the normalized edge vector and edge length
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        edge_vector_x[jj]=vtx->neigh[jj]->x - vtx->x;
        edge_vector_y[jj]=vtx->neigh[jj]->y - vtx->y;
        edge_vector_z[jj]=vtx->neigh[jj]->z - vtx->z;

        edge_length=sqrt(edge_vector_x[jj]*edge_vector_x[jj]+edge_vector_y[jj]*edge_vector_y[jj]+edge_vector_z[jj]*edge_vector_z[jj]);
        edge_vector_x[jj]=edge_vector_x[jj]/edge_length;
        edge_vector_y[jj]=edge_vector_y[jj]/edge_length;
        edge_vector_z[jj]=edge_vector_z[jj]/edge_length;

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.2: get the edge adjacent triangles lm and lp
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // We can get the two triangles since everything is ordered: 
        // for edge v->i, the previous triangle lm={v,i-1,i} is in position i-1 and the next triangle lp={v,i,i+1} is in position i
        lp = vtx->tristar[jj];
        lm = vtx->tristar[prev_small(jj, vtx->tristar_no)];

        // Triangle normals are NORMALIZED!

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.3: get the edge normal and edge binormal
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // we want to get the edge normal is the average of the normal of the two adjacent triangles ne = nf1+nf2
        temp_length=sqrt( pow((lm->xnorm + lp->xnorm),2) + pow((lm->ynorm + lp->ynorm), 2) + pow((lm->znorm + lp->znorm), 2));

        edge_normal_x[jj]=-(lm->xnorm + lp->xnorm)/temp_length;
        edge_normal_y[jj]=-(lm->ynorm + lp->ynorm)/temp_length;
        edge_normal_z[jj]=-(lm->znorm + lp->znorm)/temp_length;
        // edge binormal is normal x edge vector
        edge_binormal_x[jj]= (edge_normal_y[jj]*edge_vector_z[jj])-(edge_normal_z[jj]*edge_vector_y[jj]);
        edge_binormal_y[jj]=-(edge_normal_x[jj]*edge_vector_z[jj])+(edge_normal_z[jj]*edge_vector_x[jj]);
        edge_binormal_z[jj]= (edge_normal_x[jj]*edge_vector_y[jj])-(edge_normal_y[jj]*edge_vector_x[jj]);

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.4: get the dihedral curvature weight he[jj]
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cross_x = lm->ynorm*edge_normal_z[jj] - lm->znorm*edge_normal_y[jj];
        cross_y = lm->znorm*edge_normal_x[jj] - lm->xnorm*edge_normal_z[jj];
        cross_z = lm->xnorm*edge_normal_y[jj] - lm->ynorm*edge_normal_x[jj];

        he[jj]=edge_length*(cross_x*edge_vector_x[jj] + cross_y*edge_vector_y[jj] + cross_z*edge_vector_z[jj] );
        

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.5: get the edge shape operator Se and edge weight We_Av
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Se11=edge_binormal_x[jj]*edge_binormal_x[jj]*he[jj];
        Se21=edge_binormal_x[jj]*edge_binormal_y[jj]*he[jj];
        Se22=edge_binormal_y[jj]*edge_binormal_y[jj]*he[jj];
        Se31=edge_binormal_x[jj]*edge_binormal_z[jj]*he[jj];
        Se32=edge_binormal_y[jj]*edge_binormal_z[jj]*he[jj];
        Se33=edge_binormal_z[jj]*edge_binormal_z[jj]*he[jj];
        Se12=Se21; Se13=Se31; Se23=Se32; //symmetric matrix: hopefully the compiler gets rid of these

        // weight: edge normal . veretex normal / area (check if this is per edge?!)
        We=vertex_normal_x*edge_normal_x[jj]+vertex_normal_y*edge_normal_y[jj]+vertex_normal_z*edge_normal_z[jj];
        We_Av=We/Av;


        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // step 2.6: accumulate contribution to shape operator
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Sv[0][0]+=We_Av* Se11;
        Sv[0][1]+=We_Av* Se12;
        Sv[0][2]+=We_Av* Se13;
        
        Sv[1][0]+=We_Av* Se21;
        Sv[1][1]+=We_Av* Se22;
        Sv[1][2]+=We_Av* Se23;

        Sv[2][0]+=We_Av* Se31;
        Sv[2][1]+=We_Av* Se32;
        Sv[2][2]+=We_Av* Se33;

    } // END FOR JJ


    // ##############################################################################
    // step 3. use vertex shape operator to obtain energy and curvature information #
    // ##############################################################################


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 3.1: project to a 2x2 matrix in the surface tangent plane
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // splat the shape operator back into separate variables
    Se11=Sv[0][0];
    Se12=Sv[0][1];
    Se13=Sv[0][2];
    Se21=Sv[1][0];
    Se22=Sv[1][1];
    Se23=Sv[1][2];
    Se31=Sv[2][0];
    Se32=Sv[2][1];
    Se33=Sv[2][2];

    // project shape operator to the director-tangent plane
    // [a b] = [ dSd dSt ]
    // [b d]   [ tSd tSt ]
    dSd =   director_x*Se11*director_x + director_x*Se12*director_y + director_x*Se13*director_z
           +director_y*Se21*director_x + director_y*Se22*director_y + director_y*Se23*director_z
           +director_z*Se31*director_x + director_z*Se32*director_y + director_z*Se33*director_z; 
    dSt =   director_x*Se11*tangent_x  + director_x*Se12*tangent_y  + director_x*Se13*tangent_z
           +director_y*Se21*tangent_x  + director_y*Se22*tangent_y  + director_y*Se23*tangent_z
           +director_z*Se31*tangent_x  + director_z*Se32*tangent_y  + director_z*Se33*tangent_z ; 
    tSt =   tangent_x *Se11*tangent_x  + tangent_x *Se12*tangent_y  + tangent_x *Se13*tangent_z
           +tangent_y *Se21*tangent_x  + tangent_y *Se22*tangent_y  + tangent_y *Se23*tangent_z
           +tangent_z *Se31*tangent_x  + tangent_z *Se32*tangent_y  + tangent_z *Se33*tangent_z ; 
    error_correction_scheme(&dSd,&dSt,&tSt,vtx->mean_curvature, vtx->gaussian_curvature);
    tSd = dSt; // symmetric tensor

    //DEBUG
    vtx->S[0] = dSd;
    vtx->S[1] = tSd;
    vtx->S[2] = dSt;
    vtx->S[3] = tSt;

    return update_vertex_from_curvature_tensor(vtx,Av,dSd,dSt,tSd,tSt,
        vertex_normal_x, vertex_normal_y, vertex_normal_z,
        director_x,director_y,director_z,
        tangent_x, tangent_y, tangent_z);
}


/** @brief Calculation of the bending energy of the vertex.
 *  
 * Choose the right curvature model and update the vertex energy, curvature, and normal
 * @param *vtx is a pointer to vertex at which we want to calculate the energy
 * @returns TS_SUCCESS on successful calculation.
*/
inline ts_bool vertex_curvature_energy(ts_vesicle *vesicle, ts_vertex *vtx){
    ts_double temp,nx,ny,nz,mean_curvature,gaussian_curvature,mean_energy,gaussian_energy,dx,dy,dz;
    ts_flag model=vesicle->tape->curvature_model; // control how and what model we use to calculate energy: see enum curvature_model_type in general.h
    ts_bool do_calculate_shape_op=0, do_use_shape_op_e=0;
    do_calculate_shape_op = model&to_calculate_shape_operator 
                            || (model&to_use_shape_for_anisotropy_only 
                                && (vtx->type & is_anisotropic_vtx))
                            || (model&to_disable_calculate_laplace_beltrami);  
    do_use_shape_op_e = (model&to_use_shape_operator_energy 
                            && do_calculate_shape_op);    


    // calculate the energy using the right model
    if (! (model&to_disable_calculate_laplace_beltrami)) {
        laplace_beltrami_curvature_energy(vesicle,vtx);
        
        if (model&to_update_director_shapeless){
            // t = t - (t.n)n   (or -nx(nxt)
            temp = (vtx->dx*vtx->nx)+(vtx->dy*vtx->ny)+(vtx->dz*vtx->nz);
            vtx->dx = vtx->dx - temp*vtx->nx;
            vtx->dy = vtx->dy - temp*vtx->ny;
            vtx->dz = vtx->dz - temp*vtx->nz;

            // this operation exclusively lowers |t|: if we do manage to avoid using the size (always taking t*A*(nxt)/t^2) we can avoid the normalization
            // and just periodically make sure the size is large enough if (t^2<0.5) t=2t
            temp=sqrt((vtx->dx*vtx->dx)+(vtx->dy*vtx->dy)+(vtx->dz*vtx->dz)); // should be the same as sqrt(1-*(n.t)^2)
            vtx->dx=vtx->dx/temp;
            vtx->dy=vtx->dy/temp;
            vtx->dz=vtx->dz/temp;
        }
        nx = vtx->nx;
        ny = vtx->ny;
        nz = vtx->nz;
        mean_curvature = vtx->mean_curvature;
        mean_energy = vtx->mean_energy;
        gaussian_curvature = vtx->gaussian_curvature;
        gaussian_energy = vtx->gaussian_energy;
        dx = vtx->dx;
        dy = vtx->dy;
        dz = vtx->dz;
    } 

    if (do_calculate_shape_op) {
        tensor_curvature_energy(vesicle, vtx);
        // use the old energy and old normals
        if (!do_use_shape_op_e){
            vtx->nx2 = vtx->nx;
            vtx->energy = mean_energy + gaussian_energy;
            vtx->nx = nx;
            vtx->ny2 = vtx->ny;
            vtx->ny = ny;
            vtx->nz2 = vtx->nz;
            vtx->nz = nz;
            vtx->mean_curvature2 = vtx->mean_curvature;
            vtx->mean_curvature = mean_curvature;
            vtx->mean_energy2 = vtx->mean_energy;
            vtx->mean_energy = mean_energy;
            vtx->gaussian_curvature2 = vtx->gaussian_curvature;
            vtx->gaussian_curvature = gaussian_curvature;
            vtx->gaussian_energy2 =vtx->gaussian_energy;
            vtx->gaussian_energy = gaussian_energy;
            vtx->dx=dx;
            vtx->dy=dy;
            vtx->dz=dz;
        } else {
            // use the new version, write values in debug mode
            vtx->nx2 = nx;
            vtx->ny2 = ny;
            vtx->nz2 = nz;
            vtx->mean_curvature2 = mean_curvature;
            vtx->mean_energy2 = mean_energy;
            vtx->gaussian_curvature2 = gaussian_curvature;
            vtx->gaussian_energy2 = gaussian_energy;
        }
    }

    return TS_SUCCESS;
}

ts_bool sweep_attraction_bond_energy(ts_vesicle *vesicle){
    ts_idx i;
    for(i=0;i<vesicle->blist->n;i++){
        attraction_bond_energy(vesicle, vesicle->blist->bond[i]);
    }
    return TS_SUCCESS;
}

/**
 * @brief calculate bonding energy between the two vertcies in a bond  
 * 
 * 
 * 4 things determine the bonding (in theory): bond model, bond length, bonding strength, type of the vertices, and vertex orientation  
 * bond models:  vesicle->tape->bond_model  
 *    &1 : demand vertex type be the same to bond. (0: all bond, 1: same type)
 *    &2 : anisotropic vertices have nematic order (2: )
 * bond length: 
 *   currently bonding is not affected
 * vertex properties: bond->vtx1->#  and bond->vtx2->#
 *   type: only if both vertices are of bbonding types they have bonding energy 
 *        vtx->type&is_bonding_vertex ==0 -> energy=0
 *         bond model &1 also demands type equality
 *   w: bonding strength. currently used as if universal, so for now just set to mean
 *       energy = -(w1+w2)/2
 *   orientation: anisotropic vertices have direction vtx->dx,dy,dz
 *         bond_model&2: should have nematic order E~(d1*d2)^2
 * 
 * @param vesicle pointer to primary simulation data structure, had bond model
 * @param bond pointer to bond containng two vertices
 * @return ts_bool 
 */
inline ts_bool attraction_bond_energy(ts_vesicle *vesicle, ts_bond *bond){
    ts_double energy=0;
    ts_flag bond_model = vesicle->tape->bond_model;
    // 1 bit: bond by type
    if(!(bond_model&is_bonding_type_specific) ){ 
        // all bonding type bond together
        if((bond->vtx1->type&is_bonding_vtx && bond->vtx2->type&is_bonding_vtx)){
            energy=-0.5*(bond->vtx1->w+bond->vtx2->w);
        }
        else {
            energy=0.0;
        }
    }
    else{ 
        // only bond by same type
        if((bond->vtx1->type&is_bonding_vtx && bond->vtx2->type==bond->vtx1->type)){
            energy=-0.5*(bond->vtx1->w+bond->vtx2->w);
        }
        else {
            energy=0.0;
        }
    }
    //2 bit: anisotropy
    if(bond_model&is_anisotropic_bonding_nematic){ 
        // bond by director with nematic order (arc-like proteins)
        // bonding*= (d1*d2)^2
        if((bond->vtx1->type&is_anisotropic_vtx && bond->vtx2->type&is_anisotropic_vtx)){
            energy*=pow( bond->vtx1->dx*bond->vtx2->dx
                        +bond->vtx1->dy*bond->vtx2->dy
                        +bond->vtx1->dz*bond->vtx2->dz,2);
        }
    }
    bond->energy=energy;
    return TS_SUCCESS;
}

/** @brief return the work from the vertex move under force move W = -F(vtx) * dx(vtx,vtx_old)
 *  
 *  function work differently depending on vesicle->tape->model
 *  0: force in normal direction F = vtx->f * vtx->n
 *  1: inhibation F*=No_inactive/No_neigh
 *  2: inhibation F*=No_c>=0/No_neigh
 *  3: F=0 for any neighbor with c0<0
 *  16: Vicsek model, constant weight
 *  17: Vicsek model, ~ 1/R
 *
 *  @param *vesicle: the vesicle (which includes model and vertexlist)
 *  @param *vtx pointer to moved vertex
 *  @param *vtx_old pointer to backed-up pre-move vertex
 *  @returns ts_double work -fdx
*/
ts_double direct_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
    // modified to include Vicsek Interaction and basic inhibition models



    // quit if there is no force
    if(fabs(vtx->f)<1e-15) return 0.0;

    ts_flag model=vesicle->tape->force_model;
    ts_double vicsek_strength=vesicle->tape->vicsek_strength;
    ts_double vicsek_radius= vesicle->tape->vicsek_radius;

    ts_double norml,ddp=0.0;
    //ts_double xnorm=0.0,ynorm=0.0,znorm=0.0;
    ts_double vixnorm=0.0,viynorm=0.0,viznorm=0.0;

    //stupid loop variables don't go in stupid for loop due to stupid C89 compiler
    ts_idx i, j, curr_dist;
    
    // Initial allocation size for seen_vertex.
    ts_idx size_vtx_seen=64;
    // allocate the struct for the "seen vertex" (defined in general.h, functions in vertex.c)
    ts_seen_vertex *seen_vtx; //initalize in vicsek

    ts_uint No_neigh_activating; // number of activating neighbors: actual number, not an index!
    ts_double inhibition_factor;



    // if vicsek type and vicsek model is relevant
    if ( (model&is_vicsek_model)
        && fabs(vicsek_strength)>1e-15 && fabs(vicsek_radius)>1e-15 
        && vtx->type&is_vicsek_vtx ) {
        
        //vicsek model
        //force directed by Vicsek sum-over-neighbors-normals

        //prime vertex normal
        vixnorm=vtx->nx;
        viynorm=vtx->ny;
        viznorm=vtx->nz;

        //initialize seen_vtx
        seen_vtx = init_seen_vertex(size_vtx_seen);
        // we have now seen the prime vertex
        add_vtx_to_seen(seen_vtx, vtx);


        // Breadth first search using seen_vtx, layer by layer,
        // until reaching layer that is outside the maximum radius
        for (curr_dist=1 ; curr_dist<=vicsek_radius; curr_dist++){
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
                        // Vicsek model 17: weight by 1/distance
                        if (model == model_vicsek_1overR) {
                            vixnorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nx / curr_dist;
                            viynorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->ny / curr_dist;
                            viznorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nz / curr_dist;
                        }
                        else {
                            vixnorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nx;
                            viynorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->ny;
                            viznorm += vicsek_strength * seen_vtx->vtx[i]->neigh[j]->nz;
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
        vtx->fx=vtx->f*vixnorm;
        vtx->fy=vtx->f*viynorm;
        vtx->fz=vtx->f*viznorm;

        //don't forget to free! 
        seen_vertex_free(seen_vtx);

    } //end if (!Vicsek)
    else {
        //regular "force in normal direction"
        vtx->fx = vtx->f*vtx->nx;
        vtx->fy = vtx->f*vtx->ny;
        vtx->fz = vtx->f*vtx->nz;
    }

    // now we calculate the final force and work


    //we recalculate the force based on inhibition: we want the real force to be recorded!

    // HIV_Gag inhibition model
    if (model==model_active_neigh_interfer){
        // force ~ number of active neighbors
        No_neigh_activating = 0;
        for (i=0; i<vtx->neigh_no; i++){ 
            if ( !( vtx->neigh[i]->type & is_active_vtx ) ) {
                No_neigh_activating++;
            }
        }
        inhibition_factor = 1;
        inhibition_factor *= No_neigh_activating;  //take that, integer-integer division!
        inhibition_factor /= vtx->neigh_no;
        // change force by factor
        vtx->fx *= inhibition_factor;
        vtx->fy *= inhibition_factor;
        vtx->fz *= inhibition_factor;
    }

    if (model==model_concave_neigh_interfer){
        // force ~ number of c>0 neighbors
        No_neigh_activating = 0;
        for (i=0; i<vtx->neigh_no; i++){
            if ( vtx->neigh[i]->c > -1e-15   ) {
                No_neigh_activating++;
            }
        }
        inhibition_factor = 1;
        inhibition_factor *= No_neigh_activating; 
        inhibition_factor /= vtx->neigh_no;
        // change force by factor
        vtx->fx *= inhibition_factor;
        vtx->fy *= inhibition_factor;
        vtx->fz *= inhibition_factor;
    }
    if (model==model_concave_neigh_disable){
        // force disabled for any c<0 neighbors
        for (i=0; i<vtx->neigh_no; i++){
            if ( vtx->neigh[i]->c < -1e-15   ) {
                vtx->fx = 0;
                vtx->fy = 0;
                vtx->fz = 0;
                break;
            }
        }
    }


    /*calculate ddp, normal force directed displacement*/
    ddp=vtx->fx*(vtx->x-vtx_old->x)+vtx->fy*(vtx->y-vtx_old->y)+vtx->fz*(vtx->z-vtx_old->z);

    
    /*dW=-Fdx*/
    return -ddp;
}


/** @brief get work from residual force on the vesicle due to force balance on z
 *  
 *
 *  @param *vesicle: the vesicle (which includes model and vertexlist)
 *  @param *vtx pointer to moved vertex
 *  @param *vtx_old pointer to backed-up pre-move vertex
 *  @returns ts_double work
*/
ts_double direct_force_from_Fz_balance(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
    // calculate the Fz force balance
    //     We could separate from direct_force_energy()
    //     since F_direct(old) = F_direct + Fz_balance \hat{z}
    //     W = F_direct * dX + Fz_balance * dz

    ts_double f_diff=0;

    if (vtx->type&is_active_vtx) f_diff+=vtx->fz;
    if (vtx_old->type&is_active_vtx) f_diff-=vtx_old->fz;
    vesicle->fz+=f_diff;

    if (vesicle->fz>0){
        return vesicle->fz*(vtx->z-vtx_old->z)/vesicle->vlist->n;
    }
    else{
        return 0;
    }	
    
}

ts_bool total_force_on_vesicle(ts_vesicle *vesicle){
    ts_idx i;
    vesicle->fx = 0;
    vesicle->fy = 0;
    vesicle->fz = 0;
    /*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
    for (i=0;i<vesicle->vlist->n;i++){
        if(vesicle->vlist->vtx[i]->type&is_active_vtx){
            vesicle->fx+=vesicle->vlist->vtx[i]->fx;
            vesicle->fy+=vesicle->vlist->vtx[i]->fy;
            vesicle->fz+=vesicle->vlist->vtx[i]->fz;
        }
    }
    return TS_SUCCESS;
    
}

/** @brief return adhesion energy difference E_ad(vtx)-E_ad(vtx_old)
 *  
 *  functions differently depending on vesicle->tape->adhesion_model:
 *  1: step potential
 *  2: parabolic potential
 *  4: with y anisotropy (parabolic)
 *  and vesicle->tape->geometry
 *  1: plane substrate
 *  2: spherical substrate
 *  3: cylindrical substrate
 *  4: sinosoidal (x axis)
 *  5: plane with adhesive points on a grid
 *  
 *  @param *vesicle: the vesicle (which includes model, adhesion substrate values)
 *  @param *vtx pointer to moved vertex
 *  @param *vtx_old pointer to backed-up pre-move vertex
 *  @returns ts_double energy difference in state
*/
ts_double adhesion_energy_diff(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
    ts_flag model = vesicle->tape->adhesion_model;
    if (model==0) return 0;
    ts_double delta_energy=0;
    ts_double delta=0,delta_old=0;
    ts_bool oriented=0, oriented_old=0;
    ts_double dz=vesicle->tape->adhesion_cutoff;
    ts_double adhesion_factor=1,adhesion_factor_old=1;

    delta = adhesion_geometry_distance(vesicle, vtx);
    oriented = adhesion_geometry_side(vesicle,vtx);
    adhesion_factor=adhesion_geometry_factor(vesicle, vtx);
    delta_old = adhesion_geometry_distance(vesicle, vtx_old);
    oriented_old = adhesion_geometry_side(vesicle,vtx);
    adhesion_factor_old=adhesion_geometry_factor(vesicle, vtx_old);

    if( (vtx->type&is_adhesive_vtx) && (oriented) && (delta<=dz) ){
            if (model&adhesion_step_potential) {
                delta_energy-=vtx->ad_w*adhesion_factor;
            } else if (model&adhesion_parabolic_potential){
                delta_energy-=vtx->ad_w*(1-pow(delta/dz,2))*adhesion_factor;
            }
    }
    if( (vtx_old->type&is_adhesive_vtx) && (oriented_old) &&  (delta_old<=dz) ){
            if (model&adhesion_step_potential) {
                delta_energy+=vtx_old->ad_w*adhesion_factor_old;
            } else if (model&adhesion_parabolic_potential){
                delta_energy+=vtx_old->ad_w*(1-pow(delta_old/dz,2)*adhesion_factor_old);
            }
    }


    return delta_energy;
}

// Distance of vertex from the geometry: negative "inside" wall/cylinder/sphere
ts_double adhesion_geometry_distance(ts_vesicle *vesicle, ts_vertex *vtx){
    ts_double z=vtx->z;
    ts_double z0=vesicle->tape->adhesion_z;
    ts_double r=vesicle->tape->adhesion_radius;
    ts_double c0=z0-r;
    ts_flag geometry = vesicle->tape->adhesion_geometry;

    //1 for plane potential
    if(geometry==model_plane_potential || geometry==model_plane_potential_with_spots){
        return z-z0;
    }
    //2 for spherical potential
    else if(geometry==model_spherical_potential){
        return pow(pow(c0-z,2) + pow(vtx->x,2) + pow(vtx->y,2),0.5)- r;
    }
    //3 for cylindrical adhesive substrate
    else if(geometry==model_cylindrical_potential){
        return pow(pow(c0-z,2) + pow(vtx->x,2),0.5)- r;
    }
    // 4 for sinodoidal substrate
    else if(geometry==model_x_sinosoidal_potential){
        c0 = vesicle->tape->adhesion_scale;
        return z-(z0+r*cos(vtx->x/c0));
    }
    else if (geometry==model_xy_sinosoidal_potential){
        c0 = vesicle->tape->adhesion_scale;
        return z-(z0+r*cos(vtx->x/c0)*cos(vtx->y/c0));
    }

    return 0;
}

// Check if vertex normal is oriented towards the geometry
ts_bool adhesion_geometry_side(ts_vesicle *vesicle, ts_vertex *vtx){
    ts_double z=vtx->z;
    ts_double r=vesicle->tape->adhesion_radius;
    ts_double c0=vesicle->tape->adhesion_z-r;
    ts_flag geometry = vesicle->tape->adhesion_geometry;

    //1 for plane potential, surface facing the -z direction way
    if(geometry==model_plane_potential|| geometry==model_plane_potential_with_spots){
        return vtx->nz<0;
    }
    //2 for spherical potential, surface facing into the sphere -x,-y,-(z+c0)
    else if(geometry==model_spherical_potential){
        return (vtx->nx*vtx->x + vtx->ny*vtx->y + vtx->nz*(z-c0))<0; 
    }
    //3 for cylindrical adhesive substrate, surface facing the (-x,0,-(z+c0)) way
    else if(geometry==model_cylindrical_potential){
        return (vtx->nx*vtx->x + vtx->nz*(z-c0))<0 ;
    } 
    else if(geometry==model_x_sinosoidal_potential){
        c0 = vesicle->tape->adhesion_scale;
        return (vtx->nz + r*sin(vtx->x/c0)*vtx->nx/c0) < 0;
    }
    else if(geometry==model_xy_sinosoidal_potential){
        c0 = vesicle->tape->adhesion_scale;
        return (vtx->nz 
                + r*sin(vtx->x/c0)*cos(vtx->y/c0)*vtx->nx/c0 
                + r*sin(vtx->y/c0)*cos(vtx->y/c0)*vtx->ny/c0) < 0;
    }
    return 1;
}

// make spots that are more adhesive in a lattice
ts_double adhesion_geometry_factor(ts_vesicle* vesicle, ts_vertex* vtx){
    ts_flag model = vesicle->tape->adhesion_model;
    ts_flag geometry = vesicle->tape->adhesion_geometry;
    ts_double factor = 1;
    if (model&adhesion_y_anisotropy){
        ts_double scale = vesicle->tape->adhesion_scale;
        factor *= exp(-vtx->y/scale); // Shubhadeep's version: (x-xmin)/(xmax-xmin)
    }
    if (geometry==model_plane_potential_with_spots){
        ts_double lattice = vesicle->tape->adhesion_scale;
        ts_double throwaway;
        ts_double dx = modf(vtx->x/lattice,&throwaway)-0.5;
        ts_double dy = modf(vtx->y/lattice,&throwaway)-0.5;
        if (dx*dx+dy*dy <= (1/(lattice*lattice))) {
            factor *= vesicle->tape->adhesion_factor;
        }
    }
    return factor;
}

/** @brief Calculate volume, pressure, and area constraints and energy, except triangle strech
 * return TS_SUCCESS if no constraints are violate, 
 * energy difference is added to energy *energy+=delta_energy
 * triangle stretch requires too many parts, so it is calculated outside!
*/
ts_bool volume_pressure_area_energy_constraints(ts_vesicle *vesicle, ts_double *energy, ts_double dvol, ts_double darea){
    ts_double v_sph, v_sph_old;
    ts_double delta_energy=0;

    // area
    if (vesicle->tape->area_switch==1){
        // nothing: stretch energy is per triangle so it can't be wraped neatly in the function
    } else if(vesicle->tape->area_switch==2){
        /* check whether the darea is gt epsarea */
        if((fabs(vesicle->area+darea-A0)>epsarea) && (fabs(vesicle->area+darea-A0)>fabs(vesicle->area-A0))){
            return TS_FAIL;
        }
    } else if(vesicle->tape->area_switch==3){
        delta_energy += 0.5*vesicle->tape->xkA0*(pow(vesicle->area+darea-A0,2) - pow(vesicle->area-A0,2))/A0;
    }


    // pressure energy
    if(vesicle->tape->pressure_switch==1){
         delta_energy-=vesicle->pressure*dvol;
    }

    //volume
    if(vesicle->tape->volume_switch == 1){
        //retval=constvolume(vesicle, vtx, -dvol, &delta_energy_cv, &constvol_vtx_moved,&constvol_vtx_backup);
        return TS_FAIL;
    } else if(vesicle->tape->volume_switch==2){
        /*check whether the dvol is gt than epsvol */
        if( (fabs(vesicle->volume+dvol-V0)>epsvol) && (fabs(vesicle->volume+dvol-V0)>fabs(vesicle->volume-V0))){
            return TS_FAIL;
        }
    } else if(vesicle->tape->volume_switch==3){
        /*energy difference */
        delta_energy += 0.5*vesicle->tape->xkV0*(pow(vesicle->volume+dvol-V0,2) - pow(vesicle->volume-V0,2))/V0;
    } else if(vesicle->tape->volume_switch==4){
        /*energy difference */
        v_sph=sqrt(pow(vesicle->area+darea,3)/M_PI)/6;
        v_sph_old=sqrt(pow(vesicle->area,3)/M_PI)/6;
        delta_energy += 0.5*vesicle->tape->xkV0*(pow(vesicle->tape->Vfraction-((vesicle->volume+dvol)/v_sph),2) - pow(vesicle->tape->Vfraction-(vesicle->volume/v_sph_old),2));
    }

    // update energy
    *energy+=delta_energy;
    return TS_SUCCESS;
}

void stretchenergy(ts_vesicle *vesicle, ts_triangle *triangle){
    triangle->energy=vesicle->tape->xkA0/2.0*pow((triangle->area/vesicle->tlist->a0-1.0),2);
}

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


int cmpfunc(const void *x, const void *y)
{
    double diff = fabs(*(double*)x) - fabs(*(double*)y);
    if(diff<0) return 1;
    else return -1;
}

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
        if (!(vtx[i]->type&is_ghost_vtx)) {       
            energy_vertex(vesicle, vtx[i]);
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


/** @brief anisotropic subfunction of vertex_energy.
 *
 *  This function is experimental, branch from vertex_energy for anisotropic proteins 
 *  to calculate the tensor-based curvature values c1,c2 and principle directions,
 *  calculate the energy, and save it all on the vertex
 *
 *  @returns TS_SUCCESS on successful calculation
*/
/** @brief anisotropic subfunction of vertex_energy.
 *
 *  This function is experimental, branch from vertex_energy for anisotropic proteins 
 *  to calculate the tensor-based curvature values c1,c2 and principle directions,
 *  calculate the energy, and save it all on the vertex
 *
 *  @returns TS_SUCCESS on successful calculation
*/
inline ts_bool curvature_tensor_energy_vertex(ts_vesicle *vesicle, ts_vertex *vtx){
    // we hardcoded 10 neighbor limit!
    ts_small_idx jj, i;
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
    ts_double tr, det, lambda1, lambda2, discrim_sqrt;
    ts_double eigen_vec1d, eigen_vec1t, eigen_vec2d, eigen_vec2t;


    ts_double We;
    ts_double Av, We_Av;

    ts_double he[10];
    ts_double Sv[3][3]={{0,0,0},{0,0,0},{0,0,0}};

    // #########################################################################
    // step 1. calculate the area assigned to the vertex and the vertex normal #
    // #########################################################################
    Av=0;
    for(i=0; i<vtx->tristar_no; i++){
        vertex_normal_x=(vertex_normal_x - vtx->tristar[i]->xnorm*vtx->tristar[i]->area);
        vertex_normal_y=(vertex_normal_y - vtx->tristar[i]->ynorm*vtx->tristar[i]->area);
        vertex_normal_z=(vertex_normal_z - vtx->tristar[i]->znorm*vtx->tristar[i]->area);
        Av+=vtx->tristar[i]->area/3;
    }
    temp_length=sqrt(pow(vertex_normal_x,2)+pow(vertex_normal_y,2)+pow(vertex_normal_z,2));
    vertex_normal_x=vertex_normal_x/temp_length;
    vertex_normal_y=vertex_normal_y/temp_length;
    vertex_normal_z=vertex_normal_z/temp_length;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 1.1: update vertex normal and director based on the new normal
    // and generates a director-tangent-normal frame
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vtx->nx2 = vertex_normal_x;
    vtx->ny2 = vertex_normal_y;
    vtx->nz2 = vertex_normal_z;

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

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 2.1: calculate the normalized edge vector and edge length
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    if (jj==0){
        lm = vtx->tristar[vtx->tristar_no-1];
    }
    else{
        lm = vtx->tristar[jj-1];
    } 
    //Triangle normals are NORMALIZED!

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

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 2.4: get the dihedral curvature weight he[j]
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cross_x = lm->ynorm*edge_normal_z[jj] - lm->znorm*edge_normal_y[jj];
    cross_y = lm->znorm*edge_normal_x[jj] - lm->xnorm*edge_normal_z[jj];
    cross_z = lm->xnorm*edge_normal_y[jj] - lm->ynorm*edge_normal_x[jj];

    he[jj]=edge_length*(cross_x*edge_vector_x[jj] + cross_y*edge_vector_y[jj] + cross_z*edge_vector_z[jj] );
    

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 2.5: get the edge shape operator Se and edge weight We_Av
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 2.5: add contribution to shape operator
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    dSd =   director_x*Se11*director_x + director_x*Se12*director_y + director_x*Se13*director_z
           +director_y*Se21*director_x + director_y*Se22*director_y + director_y*Se23*director_z
           +director_z*Se31*director_x + director_z*Se32*director_y + director_z*Se33*director_z; 
    dSt =   director_x*Se11*tangent_x  + director_x*Se12*tangent_y  + director_x*Se13*tangent_z
           +director_y*Se21*tangent_x  + director_y*Se22*tangent_y  + director_y*Se23*tangent_z
           +director_z*Se31*tangent_x  + director_z*Se32*tangent_y  + director_z*Se33*tangent_z ; 
    tSt =   tangent_x *Se11*tangent_x  + tangent_x *Se12*tangent_y  + tangent_x *Se13*tangent_z
           +tangent_y *Se21*tangent_x  + tangent_y *Se22*tangent_y  + tangent_y *Se23*tangent_z
           +tangent_z *Se31*tangent_x  + tangent_z *Se32*tangent_y  + tangent_z *Se33*tangent_z ; 
    tSd = dSt; // symmetric tensor

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 3.2: get curvature information from the 2x2 shape operator
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // first, the primary matrix invariants, the trace and determinant
    tr  = dSd + tSt;
    det = dSd*tSt - dSt*tSd;
    // we can immidiately get the curvatures! up to signs and factors of 2
    vtx->mean_curvature2 = tr/2;
    vtx->gaussian_curvature2 = det;
    // eigenvalues: trace determinant formula. We expect real values only, so positive discriminant
    discrim_sqrt = sqrt(tr*tr - 4*det);
    lambda1 = (tr + discrim_sqrt)/2;
    lambda2 = (tr - discrim_sqrt)/2;
    vtx->eig_v0 = lambda1;
    vtx->eig_v1 = lambda2;
    vtx->eig_v2 = 0; // we annihilate the shape operator in the normal direction

    // ---------------------------------------------------------------
    // step 3.2.1: get eigenvectors in the director-tangent directions
    // ---------------------------------------------------------------
    // construct eigenvectors in the director-tangent plane: make sure we don't have any 0 vectors!
    // based on https://math.stackexchange.com/questions/4103294/is-there-a-closed-form-expression-for-the-eigenvectors-of-a-2x2-matrix
    // a. degenerate case
    if(lambda1==lambda2){
        eigen_vec1d = 1; // we pick the director and tangent vectors as the eigenvectors
        eigen_vec1t = 0;
        eigen_vec2d = 0;
        eigen_vec2t = 1;
    }
    else{
        // b. nondegenerate case
        if(dSd>tSt){ 
            eigen_vec1d = lambda1-tSt; //tSd==0 -> lambda1=dSd, eigen_vec1d!=0
            eigen_vec1t = tSd;
            eigen_vec2d = -dSt;
            eigen_vec2t = -lambda2-dSd;
        }
        else{
            eigen_vec1d = lambda1-dSd;
            eigen_vec1t = tSd;
            eigen_vec2d = -dSt;
            eigen_vec2t = -lambda2-tSt;
        }
        // normalize the eigenvectors
        temp_length = sqrt(eigen_vec1d*eigen_vec1d + eigen_vec1t*eigen_vec1t);
        eigen_vec1d/=temp_length;
        eigen_vec1t/=temp_length;
        temp_length = sqrt(eigen_vec2d*eigen_vec2d + eigen_vec2t*eigen_vec2t);
        eigen_vec2d/=temp_length;
        eigen_vec2t/=temp_length;
    }
    // -------------------------------------------
    // step 3.2.2: get eigenvectors in real space
    // -------------------------------------------
    vtx->eig0[0] = eigen_vec1d*director_x + eigen_vec1t*tangent_x;
    vtx->eig0[1] = eigen_vec1d*director_y + eigen_vec1t*tangent_y;
    vtx->eig0[2] = eigen_vec1d*director_z + eigen_vec1t*tangent_z;
    vtx->eig1[0] = eigen_vec2d*director_x + eigen_vec2t*tangent_x;
    vtx->eig1[1] = eigen_vec2d*director_y + eigen_vec2t*tangent_y;
    vtx->eig1[2] = eigen_vec2d*director_z + eigen_vec2t*tangent_z;
    
    vtx->eig2[0] = vertex_normal_x;
    vtx->eig2[1] = vertex_normal_y;
    vtx->eig2[2] = vertex_normal_z;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // step 3.3: get energy information from the 2x2 shape difference
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // spontaneous curvature and isotropy go here!
    dSd -= vtx->c + vtx->d;
    tSt -= vtx->c - vtx->d;
    tr = dSd + tSt;
    det = dSd*tSt - tSd * dSt;
    vtx->mean_energy2 = vtx->xk*Av* pow(tr,2);
    vtx->gaussian_energy2 = vtx->xk2 * Av * det;


    return TS_SUCCESS;
}

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
inline ts_bool energy_vertex(ts_vesicle *vesicle, ts_vertex *vtx){
        
    
    ts_small_idx jj;
    ts_small_idx jjp,jjm;
    ts_vertex *j,*jp, *jm;
    ts_triangle *jt;
    ts_double s=0.0,xh=0.0,yh=0.0,zh=0.0,txn=0.0,tyn=0.0,tzn=0.0;
    ts_double x1,x2,x3,ctp,ctm,tot,xlen;
    ts_double h,ht,norml;
    ts_double angle_sum=0;
    ts_double a_dot_b, a_cross_b_x, a_cross_b_y, a_cross_b_z, mag_a_cross_b;
    ts_bool model=vesicle->tape->type_of_curvature_model;


    for(jj=0; jj<vtx->neigh_no;jj++){
        jjp=next_small(jj, vtx->neigh_no);
        jjm=prev_small(jj, vtx->neigh_no);
        j=vtx->neigh[jj];
        jp=vtx->neigh[jjp];
        jm=vtx->neigh[jjm];

        jt=vtx->tristar[jj]; // not related to the j,jp,jm, just a separate sum for txn

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



        // angle stuff for gaussian curvature
        // angle_sum += atan(ctp) + atan(ctm); // simple but slow!
        // get the angle m-vtx-j
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
    
    h=xh*xh+yh *yh+zh*zh;
    ht=txn*xh+tyn*yh + tzn*zh;
    s=s/4.0; 
#ifdef TS_DOUBLE_DOUBLE
    if(ht>=0.0) {
        vtx->mean_curvature=sqrt(h)/s;
    } else {
        vtx->mean_curvature=-sqrt(h)/s;
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
    
    // with normal: project the director on the tangent plane
    if ( model==11 || (model==10 && (vtx->type & is_anisotropic_vtx))){
        // t = t - (t.n)n   (or -nx(nxt)
        a_dot_b = (vtx->dx*vtx->nx)+(vtx->dy*vtx->ny)+(vtx->dz*vtx->nz); // temporarily used as the dot product
        vtx->dx = vtx->dx - a_dot_b*vtx->nx;
        vtx->dy = vtx->dy - a_dot_b*vtx->ny;
        vtx->dz = vtx->dz - a_dot_b*vtx->nz;
        // this operation exclusively lowers |t|: if we do manage to avoid using the size (always taking t*A*(nxt)/t^2) we can avoid the normalization
        // and just periodically make sure the size is large enough if (t^2<0.5) t=2t
        norml=sqrt((vtx->dx*vtx->dx)+(vtx->dy*vtx->dy)+(vtx->dz*vtx->dz)); // should be the same as sqrt(1-*(n.t)^2)
        vtx->dx=vtx->dx/norml;
        vtx->dy=vtx->dy/norml;
        vtx->dz=vtx->dz/norml;
    }

    // c is spontaneous curvature energy for each vertex. Should be set to zero for
    // normal circumstances.
    /* the following statement is an expression for $\frac{1}{2}\int(c_1+c_2-c_0^\prime)^2\mathrm{d}A$, where $c_0^\prime=2c_0$ (twice the spontaneous curvature)  */
    vtx->mean_energy=vtx->xk* 0.5*s*(vtx->mean_curvature-vtx->c)*(vtx->mean_curvature-vtx->c);
    
    vtx->gaussian_curvature = (2*M_PI- angle_sum)/s;

    x1 = sqrt(pow(vtx->mean_curvature,2)-vtx->gaussian_curvature); // deltaC/2 in temp variable
    vtx->new_c1 = vtx->mean_curvature + x1;
    vtx->new_c2 = vtx->mean_curvature - x1;


    vtx->gaussian_energy = vtx->xk2 * s * vtx->gaussian_curvature;

    
    if (model==10 || model==11) {
        curvature_tensor_energy_vertex(vesicle, vtx);
    }
    
    if (model==10){
        vtx->energy = vtx->mean_energy2 + vtx->gaussian_energy2;
    }
    else {
        vtx->energy = vtx->mean_energy + vtx->gaussian_energy;
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


inline ts_bool attraction_bond_energy(ts_vesicle *vesicle, ts_bond *bond){

    if (vesicle->tape->type_of_bond_model==0){ // bonding type bond together
        if((bond->vtx1->type&is_bonding_vtx && bond->vtx2->type&is_bonding_vtx)){
            // f(w1,w2)
            bond->energy=-0.5*(bond->vtx1->w+bond->vtx2->w);
        }
        else {
            bond->energy=0.0;
        }
    }
    if (vesicle->tape->type_of_bond_model==1){ // bond by same type
        if((bond->vtx1->type&is_bonding_vtx && bond->vtx2->type==bond->vtx1->type)){
            // f(w1,w2)
            bond->energy=-0.5*(bond->vtx1->w+bond->vtx2->w);
        }
        else {
            bond->energy=0.0;
        }
    }
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

    ts_bool model=vesicle->tape->type_of_force_model;
    ts_double vicsek_strength=vesicle->tape->vicsek_strength;
    ts_double vicsek_radius= vesicle->tape->vicsek_radius;

    ts_double norml,ddp=0.0;
    //ts_double xnorm=0.0,ynorm=0.0,znorm=0.0;
    ts_double vixnorm=0.0,viynorm=0.0,viznorm=0.0;

    //stupid loop variables don't go in stupid for loop due to stupid C89 compiler
    ts_idx i, j, curr_dist;
    
    // estimated mallocation size for the cluster: roughly ~pi*r^2
    ts_idx max_vtx_seen=3*(( (int) vicsek_radius)+1)*(( (int) vicsek_radius)+1);
    // allocate the struct for the "seen vertex" (defined in general.h, functions in vertex.c)
    ts_seen_vertex *seen_vtx; //initalize in vicsek

    ts_uint No_neigh_activating; // number of activating neighbors: actual number, not an index!
    ts_double inhibition_factor;



    // if vicsek type and vicsek model is relevant
    if ( vtx->type&is_vicsek_vtx && (model==16 || model==17) 
        && fabs(vicsek_strength)>1e-15 && fabs(vicsek_radius)>1e-15) {
        
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
                        if (vesicle->tape->type_of_force_model == 17) {
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
    if (model==1){
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

    if (model==2){
        // force ~ number of c>0 neighbors
        No_neigh_activating = 0;
        for (i=0; i<vtx->neigh_no; i++){
            if ( vtx->neigh[i]->c > -1e-15   ) {
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
    if (model==3){
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
    // Total force is cached once and only once! : only update at each step
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
    ts_idx i;
    ts_double fz=0;
    /*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
    for (i=0;i<vesicle->vlist->n;i++){
        if(vesicle->vlist->vtx[i]->type&is_active_vtx){
            fz+=vesicle->vlist->vtx[i]->fz;
        }
    }
    return fz;
    
}

/** @brief return adhesion energy difference E_ad(vtx)-E_ad(vtx_old)
 *  
 *  functions differently depending on vesicle->adhesion_model:
 *  1: step potential
 *  2: parabolic potential
 *  3: spherical substrate (parabolic)
 *  4: cylindrical substrate (parabolic)
 *  
 *  @param *vesicle: the vesicle (which includes model, adhesion substrate values)
 *  @param *vtx pointer to moved vertex
 *  @param *vtx_old pointer to backed-up pre-move vertex
 *  @returns ts_double energy difference in state
*/
ts_double adhesion_energy_diff(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
    ts_double delta_energy=0;
    ts_double z=vtx->z;
    ts_double z_old=vtx_old->z;
    ts_double z0=vesicle->tape->z_adhesion;
    ts_double dz=vesicle->tape->adhesion_cuttoff;
    ts_double c0=vesicle->adhesion_center;
    ts_double r=vesicle->tape->adhesion_radius;

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

/** @brief anisotropic subfunction of vertex_energy.
 *
 *  This function is experimental, branch from vertex_energy for anisotropic proteins 
 *  to calculate the tensor-based curvature values c1,c2 and principle directions,
 *  calculate the energy, and save it all on the vertex
 *
 *  @returns TS_SUCCESS on successful calculation
*/
inline ts_bool debug_curvature_tensor_energy_vertex(ts_vesicle *vesicle, ts_vertex *vtx){
    // direct copy from Samo git repository
    // ...almost
    ts_small_idx jj, i;
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

    ts_triangle *lm=NULL, *lp=NULL;
    ts_double sumnorm,a,b,c;
    ts_double temp_length;
    ts_double cross_x, cross_y, cross_z;

    ts_double Se11=0, Se21=0, Se22=0, Se31=0, Se32=0, Se33=0;
    ts_double Pv11, Pv21, Pv22, Pv31, Pv32, Pv33;
    ts_double Se12, Se13, Se23, Pv12, Pv13, Pv23; 
    //alias for clarity of symmetric matrices: hopefully the compiler removes them

    ts_double We;
    ts_double Av, We_Av;

    ts_double eigenval[3];

    gsl_matrix *gsl_Sv=gsl_matrix_alloc(3,3);
    gsl_matrix *Sv_eigenV=gsl_matrix_alloc(3,3);
    gsl_vector *Sv_eigen=gsl_vector_alloc(3);
    gsl_eigen_symmv_workspace *workspace=gsl_eigen_symmv_alloc(3);

    // ts_double mprod[7], phi[7];
    ts_double he[10];
    ts_double Sv[3][3]={{0,0,0},{0,0,0},{0,0,0}};

    ts_fprintf(stdout,"|%u|: ",vtx->idx);
    //print_tri_order(vtx);

    Av=0;
    for(i=0; i<vtx->tristar_no; i++){
        vertex_normal_x=(vertex_normal_x - vtx->tristar[i]->xnorm*vtx->tristar[i]->area);
        vertex_normal_y=(vertex_normal_y - vtx->tristar[i]->ynorm*vtx->tristar[i]->area);
        vertex_normal_z=(vertex_normal_z - vtx->tristar[i]->znorm*vtx->tristar[i]->area);
        Av+=vtx->tristar[i]->area/3;
    }
    temp_length=sqrt(pow(vertex_normal_x,2)+pow(vertex_normal_y,2)+pow(vertex_normal_z,2));
    vertex_normal_x=vertex_normal_x/temp_length;
    vertex_normal_y=vertex_normal_y/temp_length;
    vertex_normal_z=vertex_normal_z/temp_length;
    vtx->nx2 = vertex_normal_x;
    vtx->ny2 = vertex_normal_y;
    vtx->nz2 = vertex_normal_z;

    Pv11=(pow(vertex_normal_x,2)+pow(vertex_normal_y,2)+pow(vertex_normal_z,2))-vertex_normal_x*vertex_normal_x;
    Pv22=(pow(vertex_normal_x,2)+pow(vertex_normal_y,2)+pow(vertex_normal_z,2))-vertex_normal_y*vertex_normal_y;
    Pv33=(pow(vertex_normal_x,2)+pow(vertex_normal_y,2)+pow(vertex_normal_z,2))-vertex_normal_z*vertex_normal_z;
    Pv21=-vertex_normal_x*vertex_normal_y;
    Pv31=-vertex_normal_x*vertex_normal_z;
    Pv32=-vertex_normal_y*vertex_normal_z;
    Pv12=Pv21; Pv13=Pv31; Pv23=Pv32; //alias for clarity of the symmetric matrix calculation


    // vertex are ordered by initial_dist and at bondflips
    for(jj=0;jj<vtx->neigh_no;jj++){
    // !!! We start a VERY long loop over jj !!!


    edge_vector_x[jj]=vtx->neigh[jj]->x-vtx->x;
    edge_vector_y[jj]=vtx->neigh[jj]->y-vtx->y;
    edge_vector_z[jj]=vtx->neigh[jj]->z-vtx->z;

    //Here we calculate normalized edge vector

    temp_length=sqrt(edge_vector_x[jj]*edge_vector_x[jj]+edge_vector_y[jj]*edge_vector_y[jj]+edge_vector_z[jj]*edge_vector_z[jj]);
    edge_vector_x[jj]=edge_vector_x[jj]/temp_length;
    edge_vector_y[jj]=edge_vector_y[jj]/temp_length;
    edge_vector_z[jj]=edge_vector_z[jj]/temp_length;


    // vtx are ordered: for edge v->i, lm={v,i-1,i} lp={v,i,i+1}, so we can get the two triangles
    lp = vtx->tristar[jj];
    if (jj==0){
        lm = vtx->tristar[vtx->tristar_no-1];
    }
    else{
        lm = vtx->tristar[jj-1];
    } 

    //Triangle normals are NORMALIZED!

    // we want to get the edge normal ne = nf+nf

    sumnorm=sqrt( pow((lm->xnorm + lp->xnorm),2) + pow((lm->ynorm + lp->ynorm), 2) + pow((lm->znorm + lp->znorm), 2));

    edge_normal_x[jj]=-(lm->xnorm + lp->xnorm)/sumnorm;
    edge_normal_y[jj]=-(lm->ynorm + lp->ynorm)/sumnorm;
    edge_normal_z[jj]=-(lm->znorm + lp->znorm)/sumnorm;


    edge_binormal_x[jj]= (edge_normal_y[jj]*edge_vector_z[jj])-(edge_normal_z[jj]*edge_vector_y[jj]);
    edge_binormal_y[jj]=-(edge_normal_x[jj]*edge_vector_z[jj])+(edge_normal_z[jj]*edge_vector_x[jj]);
    edge_binormal_z[jj]= (edge_normal_x[jj]*edge_vector_y[jj])-(edge_normal_y[jj]*edge_vector_x[jj]);

    cross_x = lm->ynorm*edge_normal_z[jj] - lm->znorm*edge_normal_y[jj];
    cross_y = lm->znorm*edge_normal_x[jj] - lm->xnorm*edge_normal_z[jj];
    cross_z = lm->xnorm*edge_normal_y[jj] - lm->ynorm*edge_normal_x[jj];

    he[jj]=temp_length*(cross_x*edge_vector_x[jj] + cross_y*edge_vector_y[jj] + cross_z*edge_vector_z[jj] );
    
    // old style: uncomment here and their variables above
    /*
    mprod[jj]=lm->xnorm*(lp->ynorm*edge_vector_z[jj]-lp->znorm*edge_vector_y[jj]) - lm->ynorm*(lp->xnorm*edge_vector_z[jj]-lp->znorm*edge_vector_x[jj])+ lm->znorm*(lp->xnorm*edge_vector_y[jj]-lp->ynorm*edge_vector_x[jj]);

    cross_x = lm->ynorm*lp->znorm - lp->ynorm*lm->znorm;
    cross_y = lm->znorm*lp->xnorm - lp->znorm*lm->xnorm;
    cross_z = lm->xnorm*lp->ynorm - lp->xnorm*lm->ynorm;
    phi[jj]=copysign(atan2(sqrt(cross_x*cross_x+cross_y*cross_y+cross_z*cross_z),lm->xnorm*lp->xnorm+lm->ynorm*lp->ynorm+lm->znorm*lp->znorm-1e-10),-mprod[jj])+M_PI;
    We = temp_length*cos(phi[jj]/2.0); //temporarily for testing: We is set later
    
    if (fabs(We - he[jj])>1e-7){
        fprintf(stdout,"%.17e by products is not the same as %.17e by cosine\n", he[jj], We);
        fatal("not equal\n",100);
    }
    */
    
    /*
    if(vtx->idx==0){
        printf("H operator of edge vertex %d (edge %d): %f\n",vtx->idx,jj,he[jj]);
    }
    */
    Se11=edge_binormal_x[jj]*edge_binormal_x[jj]*he[jj];
    Se21=edge_binormal_x[jj]*edge_binormal_y[jj]*he[jj];
    Se22=edge_binormal_y[jj]*edge_binormal_y[jj]*he[jj];
    Se31=edge_binormal_x[jj]*edge_binormal_z[jj]*he[jj];
    Se32=edge_binormal_y[jj]*edge_binormal_z[jj]*he[jj];
    Se33=edge_binormal_z[jj]*edge_binormal_z[jj]*he[jj];
    Se12=Se21; Se13=Se31; Se23=Se32; //for clarity: hopefully compiler gets rid of these
    //ts_fprintf(stdout,"Ses[%d]={%+.17f, %+.17f, %+.17f}\n",jj,Se11,Se12,Se13);
    //ts_fprintf(stdout,"        {%+.17f, %+.17f, %+.17f}\n",   Se21,Se22,Se23);
    //ts_fprintf(stdout,"        {%+.17f, %+.17f, %+.17f}\n",   Se31,Se32,Se33);
    We=vertex_normal_x*edge_normal_x[jj]+vertex_normal_y*edge_normal_y[jj]+vertex_normal_z*edge_normal_z[jj];
    //ts_fprintf(stdout,"We=%f\n", We);
    We_Av=We/Av;


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


    //if(Sv[1][0]!=Sv[0][1])fatal("01 badness\n",3);
    //if(Sv[2][0]!=Sv[0][2])fatal("02 badness\n",3);
    //if(Sv[1][2]!=Sv[2][1])fatal("12 badness\n",3);
    Se11=Sv[0][0];
    Se12=Sv[0][1];
    Se13=Sv[0][2];
    Se21=Sv[1][0];
    Se22=Sv[1][1];
    Se23=Sv[1][2];
    Se31=Sv[2][0];
    Se32=Sv[2][1];
    Se33=Sv[2][2];
    // householder transformation: get rid of n^ components
    // matrix multiplication Sv[i,j] = Pv[i,a]*Sv[a,b]*Pv^T[b,j]
    Sv[0][0]=(   Pv11*(  Se11*Pv11 + Se12*Pv21 + Se13*Pv31 )
                +Pv12*(  Se21*Pv11 + Se22*Pv21 + Se23*Pv31 )
                +Pv13*(  Se31*Pv11 + Se32*Pv21 + Se33*Pv31 ));
    Sv[0][1]=(   Pv11*(  Se11*Pv12 + Se12*Pv22 + Se13*Pv32 )
                +Pv12*(  Se21*Pv12 + Se22*Pv22 + Se23*Pv32 )
                +Pv13*(  Se31*Pv12 + Se32*Pv22 + Se33*Pv32 ));
    Sv[0][2]=(   Pv11*(  Se11*Pv13 + Se12*Pv23 + Se13*Pv33 )
                +Pv12*(  Se21*Pv13 + Se22*Pv23 + Se23*Pv33 )
                +Pv13*(  Se31*Pv13 + Se32*Pv23 + Se33*Pv33 ));
         
    Sv[1][0]=(   Pv21*(  Se11*Pv11 + Se12*Pv21 + Se13*Pv31 )
                +Pv22*(  Se21*Pv11 + Se22*Pv21 + Se23*Pv31 )
                +Pv23*(  Se31*Pv11 + Se32*Pv21 + Se33*Pv31 ));
    Sv[1][1]=(   Pv21*(  Se11*Pv12 + Se12*Pv22 + Se13*Pv32 )
                +Pv22*(  Se21*Pv12 + Se22*Pv22 + Se23*Pv32 )
                +Pv23*(  Se31*Pv12 + Se32*Pv22 + Se33*Pv32 ));
    Sv[1][2]=(   Pv21*(  Se11*Pv13 + Se12*Pv23 + Se13*Pv33 )
                +Pv22*(  Se21*Pv13 + Se22*Pv23 + Se23*Pv33 )
                +Pv23*(  Se31*Pv13 + Se32*Pv23 + Se33*Pv33 ));
         
    Sv[2][0]=(   Pv31*(  Se11*Pv11 + Se12*Pv21 + Se13*Pv31 )
                +Pv32*(  Se21*Pv11 + Se22*Pv21 + Se23*Pv31 )
                +Pv33*(  Se31*Pv11 + Se32*Pv21 + Se33*Pv31 ));
    Sv[2][1]=(   Pv31*(  Se11*Pv12 + Se12*Pv22 + Se13*Pv32 )
                +Pv32*(  Se21*Pv12 + Se22*Pv22 + Se23*Pv32 )
                +Pv33*(  Se31*Pv12 + Se32*Pv22 + Se33*Pv32 ));
    Sv[2][2]=(   Pv31*(  Se11*Pv13 + Se12*Pv23 + Se13*Pv33 )
                +Pv32*(  Se21*Pv13 + Se22*Pv23 + Se23*Pv33 )
                +Pv33*(  Se31*Pv13 + Se32*Pv23 + Se33*Pv33 ));
    ts_fprintf(stdout, "%d: {%+f,%+f,%+f}:",vtx->idx,vtx->x,vtx->y,vtx->z);
    for (jj=0; jj<vtx->neigh_no; jj++){
        ts_fprintf(stdout,"{%+f,%+f,%+f},\n",vtx->neigh[jj]->x, vtx->neigh[jj]->y,vtx->neigh[jj]->z);
    }
    for (jj=0; jj<vtx->neigh_no; jj++){
        a=vtx->neigh[jj]->x-vtx->x;
        b=vtx->neigh[jj]->y-vtx->y;
        c=vtx->neigh[jj]->z-vtx->z;
        //ts_fprintf(stdout,"edge %d={%+f,%+f,%+f}\n",vtx->neigh[jj]->idx,a, b,c);
    }
    for (jj=0; jj<vtx->neigh_no; jj++){
        //ts_fprintf(stdout,"edge (normalized) %d={%+f,%+f,%+f}\n",jj,edge_vector_x[jj], edge_vector_y[jj],edge_vector_z[jj]);
    }
    ts_fprintf(stdout,"at %d:area=%f, position={%+f,%+f,%+f} normal={%+f,%+f,%+f}\n",vtx->idx,Av,
    vtx->x,vtx->y,vtx->z,vtx->nx2, vtx->ny2,vtx->nz2);
    
    for (jj=0; jj<vtx->neigh_no; jj++){
        //ts_fprintf(stdout,"edge normal[%d]={%+f,%+f,%+f}\n",jj,edge_normal_x[jj], edge_normal_y[jj],edge_normal_z[jj]);
    }
    for (jj=0; jj<vtx->neigh_no; jj++){
        //ts_fprintf(stdout,"edge binormal[%d]={%+f,%+f,%+f}\n",jj,edge_binormal_x[jj], edge_binormal_y[jj],edge_binormal_z[jj]);
    }
    ts_fprintf(stdout,"he[%d]={",jj);
    for (jj=0; jj<vtx->neigh_no; jj++){
        fprintf(stdout,"%f,",he[jj]);
    }
    fprintf(stdout,"}\n");
    ts_fprintf(stdout,"Sv={%+.17f, %+.17f, %+.17f}\n",Sv[0][0],Sv[0][1],Sv[0][2]);
    ts_fprintf(stdout,"   {%+.17f, %+.17f, %+.17f}\n",Sv[1][0],Sv[1][1],Sv[1][2]);
    ts_fprintf(stdout,"   {%+.17f, %+.17f, %+.17f}\n",Sv[2][0],Sv[2][1],Sv[2][2]);
    //ts_fprintf(stdout,"Pv={%+.17f, %+.17f, %+.17f}\n",Pv11, Pv12, Pv13);
    //ts_fprintf(stdout,"   {%+.17f, %+.17f, %+.17f}\n",Pv21, Pv22, Pv23);
    //ts_fprintf(stdout,"   {%+.17f, %+.17f, %+.17f}\n",Pv31, Pv32, Pv33);


    // into gsl
    gsl_matrix_set(gsl_Sv, 0,0, Sv[0][0]);
    gsl_matrix_set(gsl_Sv, 0,1, Sv[0][1]);
    gsl_matrix_set(gsl_Sv, 0,2, Sv[0][2]);
    gsl_matrix_set(gsl_Sv, 1,0, Sv[1][0]);
    gsl_matrix_set(gsl_Sv, 1,1, Sv[1][1]);
    gsl_matrix_set(gsl_Sv, 1,2, Sv[1][2]);
    gsl_matrix_set(gsl_Sv, 2,0, Sv[2][0]);
    gsl_matrix_set(gsl_Sv, 2,1, Sv[2][1]);
    gsl_matrix_set(gsl_Sv, 2,2, Sv[2][2]);

    // calculate eigenvalues and eigenvectors
    gsl_eigen_symmv(gsl_Sv, Sv_eigen, Sv_eigenV, workspace);
    gsl_eigen_symmv_sort(Sv_eigen, Sv_eigenV, GSL_EIGEN_SORT_ABS_DESC);

    // get eigenvalues and eigenvectors out
    eigenval[0]= gsl_vector_get(Sv_eigen, 0);
    eigenval[1]= gsl_vector_get(Sv_eigen, 1);
    eigenval[2]= gsl_vector_get(Sv_eigen, 2);
    vtx->eig0[0] = gsl_matrix_get(Sv_eigenV,0,0);
    vtx->eig0[1] = gsl_matrix_get(Sv_eigenV,1,0);
    vtx->eig0[2] = gsl_matrix_get(Sv_eigenV,2,0);
    vtx->eig_v0 = eigenval[0];
    vtx->eig1[0] = gsl_matrix_get(Sv_eigenV,0,1);
    vtx->eig1[1] = gsl_matrix_get(Sv_eigenV,1,1);
    vtx->eig1[2] = gsl_matrix_get(Sv_eigenV,2,1);
    vtx->eig_v1 = eigenval[1];
    vtx->eig2[0] = gsl_matrix_get(Sv_eigenV,0,2);
    vtx->eig2[1] = gsl_matrix_get(Sv_eigenV,1,2);
    vtx->eig2[2] = gsl_matrix_get(Sv_eigenV,2,2);
    vtx->eig_v2 = eigenval[2];
    ts_fprintf(stdout,"eigenvals: %f,%f,%f\n",eigenval[0],eigenval[1],eigenval[2]);
    ts_fprintf(stdout,"eigenvec1: {%f,%f,%f}\n",vtx->eig0[0],vtx->eig0[1],vtx->eig0[2]);
    ts_fprintf(stdout,"eigenvec2: {%f,%f,%f}\n",vtx->eig1[0],vtx->eig1[1],vtx->eig1[2]);
    ts_fprintf(stdout,"eigenvec3: {%f,%f,%f}\n",vtx->eig2[0],vtx->eig2[1],vtx->eig2[2]);
    //And the stuff I'm tracking
    //vtx->nx = vertex_normal_x;
    //vtx->ny = vertex_normal_y;
    //vtx->nz = vertex_normal_z;
    //vtx->curvature = (eigenval[0] + eigenval[1])/2;
    //vtx->curvature2 = eigenval[0]*eigenval[1];
    vtx->mean_curvature2 = (eigenval[0]+ eigenval[1]);
    vtx->gaussian_curvature2 = eigenval[0]*eigenval[1];
    vtx->mean_energy2 = 0.25*vtx->xk*(pow(eigenval[0]+eigenval[1]-2*vtx->c,2))*Av;
    vtx->gaussian_energy2 = vtx->xk2 * Av * eigenval[0]*eigenval[1];
    ts_fprintf(stdout, "mean curvature: %f;\tgaussian curvature %f;\tmean curvature energy %f;\tgaussian curvature energy %f\n",
                        vtx->mean_curvature2, vtx->gaussian_curvature2, vtx->mean_energy2, vtx->gaussian_energy2);
    
    gsl_matrix_free(gsl_Sv);
    gsl_vector_free(Sv_eigen);
    gsl_matrix_free(Sv_eigenV);
    gsl_eigen_symmv_free(workspace);
    return TS_SUCCESS;
}


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
inline ts_bool debug_energy_vertex(ts_vesicle *vesicle, ts_vertex *vtx){
        
    
    ts_small_idx jj;
    ts_small_idx jjp,jjm;
    ts_vertex *j,*jp, *jm;
    ts_triangle *jt;
    ts_double s=0.0,xh=0.0,yh=0.0,zh=0.0,txn=0.0,tyn=0.0,tzn=0.0;
    ts_double x1,x2,x3,ctp,ctm,tot,xlen;
    ts_double h,ht,norml;
    ts_double angle_sum=0;
    ts_double a_dot_b, a_cross_b_x, a_cross_b_y, a_cross_b_z, mag_a_cross_b;
    ts_bool model=vesicle->tape->type_of_curvature_model;

    //debugging tristar order in vertex 
    ts_small_idx li, ri;
    ts_small_idx t;
    ts_vertex *vl, *vr;
    
    for(jj=0; jj<vtx->neigh_no;jj++){
        jjp=next_small(jj, vtx->neigh_no);
        jjm=prev_small(jj, vtx->neigh_no);
        j=vtx->neigh[jj];
        jp=vtx->neigh[jjp];
        jm=vtx->neigh[jjm];

        jt=vtx->tristar[jj]; // not related to the j,jp,jm, just a separate sum for txn

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



        // angle stuff for gaussian curvature
        // angle_sum += atan(ctp) + atan(ctm); // simple but slow!
        // get the angle m-vtx-j
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
        vtx->mean_curvature=sqrt(h)/s;
    } else {
        vtx->mean_curvature=-sqrt(h)/s;
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
    // c is spontaneous curvature energy for each vertex. Should be set to zero for
    // normal circumstances.
    /* the following statement is an expression for $\frac{1}{2}\int(c_1+c_2-c_0^\prime)^2\mathrm{d}A$, where $c_0^\prime=2c_0$ (twice the spontaneous curvature)  */
    vtx->mean_energy=vtx->xk* 0.5*s*(vtx->mean_curvature-vtx->c)*(vtx->mean_curvature-vtx->c);
    
    vtx->gaussian_curvature = (2*M_PI- angle_sum)/s;

    x1 = sqrt(pow(vtx->mean_curvature,2)-vtx->gaussian_curvature); // deltaC/2 in temp variable
    vtx->new_c1 = vtx->mean_curvature + x1;
    vtx->new_c2 = vtx->mean_curvature - x1;

    vtx->gaussian_energy = vtx->xk2 * s * vtx->gaussian_curvature;

    
    if (model==10 || model==11) {
        debug_curvature_tensor_energy_vertex(vesicle, vtx);
    }
    
    if (model==10){
        vtx->energy = vtx->mean_energy2 + vtx->gaussian_energy2;
    }
    else {
        vtx->energy = vtx->mean_energy + vtx->gaussian_energy;
    }

    return TS_SUCCESS;
}

//if (0){
//    //print_vertex_ordered(vtx);
//    //fprintf(stdout,"\ndiscombobulating:");
//    for (t=0; t<vtx->tristar_no; t++){
//        jj = lrand48()%vtx->tristar_no;
//        swap_triangles(vtx, t, jj);
//        //fprintf(stdout,"swapped %d, %d, ", jj, t);
//    }
//
//    //fprintf(stdout,"\nI'm still on %u, I read triangles \n", vtx->idx);
//    //print_tri_order(vtx);
//    //fprintf(stdout,"\nreordering:\n");
//    // find first triangle
//    vl = vtx->neigh[0];
//    vr = vtx->neigh[1];
//    for (t=0; t<vtx->tristar_no; t++){
//        //fprintf(stdout,"%d,",t);
//        jt = vtx->tristar[t];
//        if (in_tri(jt,vl)){
//            if (in_tri(jt,vr)){
//                jj = t;
//            }
//            else{
//                jjp = t;
//            }
//          
//        }
//        else if (in_tri(jt,vr)){
//            jjm = t;
//        }  
//    }
//    swap_triangles(vtx, jj, 0);
//    //fprintf(stdout,"swapped primary %d, %d\n",jj,0);
//    //print_tri_order(vtx);
//    if (jjp==0) jjp=jj;
//    swap_triangles(vtx, jjp, vtx->tristar_no-1);
//    //fprintf(stdout,"swapped first left %d, %d\n",jjp,vtx->tristar_no-1);
//    //print_tri_order(vtx);
//    if (jjm==0) jjm=jj;
//    if (jjm==vtx->tristar_no-1) jjm=jjp;
//    swap_triangles(vtx, jjm, 1);
//    //fprintf(stdout,"swapped first right %d, %d\n",jjm,1);
//    //print_tri_order(vtx);
//
//    ri = 2;
//    li = vtx->neigh_no-1;
//    // now triangles can only be left of left or right of right
//    while (ri+1<li){ 
//        for (t=ri; t<li; t++){
//            //fprintf(stdout,"%d,",t);
//            vl = vtx->neigh[li];
//            vr = vtx->neigh[ri];
//            jt = vtx->tristar[t];
//            if (in_tri(jt, vl)){
//                li-=1;
//                swap_triangles(vtx, t, li);
//                //fprintf(stdout,"swapped left %d, %d\n",t,li);
//                //print_tri_order(vtx);
//                if (ri+1==li) break;
//            }
//            if (in_tri(jt, vr)){
//                swap_triangles(vtx, t, ri);
//                //fprintf(stdout,"swapped right %d, %d\n",t,ri);
//                //print_tri_order(vtx);
//                ri+=1;
//                if (ri+1==li) break;
//            }
//        }
//    }
//    //fprintf(stdout,"Done reordering: I read nodes ");
//    for (jj=0; jj<vtx->neigh_no; jj++){
//        //fprintf(stdout,"%u, ", vtx->neigh[jj]->idx);
//    }
//    //fprintf(stdout,"\n");
//    //print_vertex_ordered(vtx);
//    //fprintf(stdout,"\n");
//    assert_vtx_ordered(vtx);
//    }
    

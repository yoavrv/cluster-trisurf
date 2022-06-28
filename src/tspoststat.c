/* vim: set ts=4 sts=4 sw=4 noet : */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "general.h"
//#include "vertex.h"
//#include "bond.h"
//#include "triangle.h"
//#include "cell.h"
#include "vesicle.h"
#include "io.h"
//#include "initial_distribution.h"
//#include "frame.h"
//#include "timestep.h"
//#include "poly.h"
#include "stats.h"
#include "sh.h"
#include "shcomplex.h"
#include "dumpstate.h"
#include "restore.h"
#include "cluster.h"
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <snapshot.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdio.h>

ts_vesicle *restoreVesicle(char *filename){
    ts_vesicle *vesicle = parseDump(filename);
    return vesicle;
}


void vesicle_calculate_ulm2(ts_vesicle *vesicle){
    //complex_sph_free(vesicle->sphHarmonics);

    //vesicle->sphHarmonics=complex_sph_init(vesicle->vlist,21);
    vesicle_volume(vesicle);
    preparationSh(vesicle,getR0(vesicle));
    calculateUlmComplex(vesicle);
    ts_int i,j;
    for(i=0;i<vesicle->sphHarmonics->l;i++){
            for(j=i;j<2*i+1;j++){
            printf("%e ", gsl_complex_abs2(vesicle->sphHarmonics->ulmComplex[i][j]));
            }
    }
        printf("\n");

}


ts_uint count_bonds_with_energy(ts_bond_list *blist){
    ts_idx i;
    ts_uint cnt;
    cnt=0;
    for(i=0;i<blist->n;i++){
        if(fabs(blist->bond[i]->energy)>1e-16) cnt++;
    }
    return cnt;
}


ts_bool write_histogram_data(ts_uint timestep_no, ts_vesicle *vesicle, ts_cluster_list* cstlist){
    //printf("No clusters=%d\n",cstlist->n);
    ts_idx k,i;
    int cnt, test=0;
    int max_nvtx=0;
    char filename[255];
    sprintf(filename,"histogram_%.6u.csv",timestep_no);
    FILE *fd=fopen(filename,"w");
    fprintf(fd,"Number_of_vertices_in_cluster,Number_of_clusters,\n");
    for(k=0;k<cstlist->n;k++)
        if(cstlist->cluster[k]->nvtx>max_nvtx) max_nvtx=cstlist->cluster[k]->nvtx;
    //printf("Max. number of vertices in cluster: %d\n",max_nvtx);
    for(i=1;i<=max_nvtx;i++){
        cnt=0;
        for(k=0;k<cstlist->n;k++)
            if(cstlist->cluster[k]->nvtx==i) cnt++;
        fprintf(fd,"%d,%d,\n",i,cnt);
        test+=cnt*i;
    }
    // for(k=0;k<cstlist->n;k++){
    //     printf("*Cluster %d has %d vertices\n",k,cstlist->cluster[k]->nvtx);
    // }

    fclose(fd);
    // printf("*Sum of all vertices in clusters: %d\n", test);
    // write_vertex_xml_file(vesicle,timestep_no,cstlist);
    
    return TS_SUCCESS;
}

ts_double mean_cluster_size(ts_vesicle *vesicle, ts_cluster_list* cstlist){
    // get mean size of cluster
    ts_double s=0;
    ts_idx k=0;
    for (k=0; k<cstlist->n; k++){
        s += cstlist->cluster[k]->nvtx;
    }
    s = s/cstlist->n;
    return s;
}

ts_double var_cluster_size(ts_vesicle *vesicle, ts_cluster_list* cstlist, ts_double mean_size){
    // get mean size of cluster
    ts_double s=0;
    ts_idx k=0;
    for (k=0; k<cstlist->n; k++){
        s += (cstlist->cluster[k]->nvtx-mean_size)*(cstlist->cluster[k]->nvtx-mean_size);
    }
    if (cstlist->n>1) s = s/(cstlist->n-1); //bessel correction
    return s;
}

inline ts_bool is_active(ts_vertex *vtx){
    // simple helper function: is vertex active
    if (fabs(vtx->c)>1e-16) return 1;
    return 0;
}

inline ts_double vtx_distance_sq(ts_vertex *vtx1, ts_vertex *vtx2){
    ts_double dist;
#ifdef TS_DOUBLE_DOUBLE
    dist=pow(vtx1->x-vtx2->x,2) + pow(vtx1->y-vtx2->y,2) + pow(vtx1->z-vtx2->z,2);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    dist=powl(vtx1->x-vtx2->x,2) + powl(vtx1->y-vtx2->y,2) + powl(vtx1->z-vtx2->z,2);
#endif
#ifdef TS_DOUBLE_FLOAT
    dist=powf(vtx1->x-vtx2->x,2) + powf(vtx1->y-vtx2->y,2) + powf(vtx1->z-vtx2->z,2);
#endif
    return(dist);
}

ts_double perimeter_tvvv(ts_vertex* p, ts_vertex* a, ts_vertex* b){
    // helper function for get_perimeter_triangle()
    /*perimeter caluclation is very similar to the mean curvature: 
    similar to the mean curvature:            b
    we need the same cotan(theta),           / \
    except we take both	halfs in the        /   \
    same triangle instead the two          /\    \
    sides of the bond                     /	  0   \
                                    prime ____|_<-_first part of sigma_ij
                                              | <- second part of sigma_ij
    
    */
    ts_double pa_dot_pa,pb_dot_pb,ab_dot_ab,mix_dot,cotan_theta,perim;
    //part attached to pb, angle pab
    pa_dot_pa=vtx_distance_sq(p,a); //shouldn't be zero!
    pb_dot_pb=vtx_distance_sq(p,b); //shouldn't be zero!
    ab_dot_ab=vtx_distance_sq(a,b); // shouldn't be zero!
    mix_dot=(p->x-a->x)*(b->x-a->x)+ //carefull with the order!
            (p->y-a->y)*(b->y-a->y)+ //don't want -|pa||ab|cos(theta)
            (p->z-a->z)*(b->z-a->z);
    cotan_theta=mix_dot/sqrt(pa_dot_pa*ab_dot_ab-mix_dot*mix_dot);
    perim=0.5*cotan_theta*sqrt(pb_dot_pb);
    //part attached to pa, angle pba
    mix_dot=(p->x-b->x)*(a->x-b->x)+ //flipped
            (p->y-b->y)*(a->y-b->y)+
            (p->z-b->z)*(a->z-b->z);
    cotan_theta=mix_dot/sqrt(pb_dot_pb*ab_dot_ab-mix_dot*mix_dot);
    perim+=0.5*cotan_theta*sqrt(pa_dot_pa);
    return perim;
}

ts_double get_perimeter_from_triangle(ts_triangle* triangle){
    /* Get protein-bare vertex boundary length
    a boundary goes through the triangle if one of the vertices is different
    lots of annoying combinations to go through 000,100,110,111,010,011,001 */

    // big, disgusting if tree to get the combinations: nasty
    if (is_active(triangle->vertex[0])){
        if (is_active(triangle->vertex[1])){
            if (is_active(triangle->vertex[2])){
                // 111: inside some cluster
                return 0;
            }
            else{
                // 110: 2 is different
                return perimeter_tvvv(triangle->vertex[2],triangle->vertex[0],triangle->vertex[1]);
            }
        }
        else{ 
            if (is_active(triangle->vertex[2])){
                // 101: 1 is different
                return perimeter_tvvv(triangle->vertex[1],triangle->vertex[2],triangle->vertex[0]);
            }
            else{ 
                // 100: 0 is different
                return perimeter_tvvv(triangle->vertex[0],triangle->vertex[1],triangle->vertex[2]);
            }
        }
    }
    else{
        if (is_active(triangle->vertex[1])){
            if (is_active(triangle->vertex[2])){
                // 011: 0 is different
                return perimeter_tvvv(triangle->vertex[0],triangle->vertex[1],triangle->vertex[2]);
            }
            else{
                // 010: 1 is different
                return perimeter_tvvv(triangle->vertex[1],triangle->vertex[2],triangle->vertex[0]);
            }
        }
        else{ 
            if (is_active(triangle->vertex[2])){
                // 001: 2 is different
                return perimeter_tvvv(triangle->vertex[2],triangle->vertex[0],triangle->vertex[1]);
            }
            else{ 
                // 000: inside bare area
                return 0;
            }
        }
    }
}


ts_double get_full_perimeter(ts_vesicle* vesicle){
    // get the sum total line length 
    // of the boudaries of the clusters
    ts_idx k;
    ts_double perim=0;
    for (k=0;k<vesicle->tlist->n; k++){
        perim+=get_perimeter_from_triangle(vesicle->tlist->tria[k]);
    }
    return perim;
}



int main(){
    ts_vesicle *vesicle;
    ts_char *i,*j;
    ts_uint tstep,n;
    ts_char *number;
    struct dirent **list;
    ts_double l1,l2,l3,hbar,Nbw_Nb;
    int count;
    ts_cluster_list *cstlist;
    ts_double mean_size, var_size;
    ts_double perim;

    fprintf(stdout, "OuterLoop,Volume,Area,lamdba1,lambda2,lambda3,Nbw/Nb,hbar,mean_cluster_size,std_cluster_size,line_length,\n");


    count=scandir(".",&list,0,alphasort);
    if(count<0){
        fatal("Error, cannot open directory.",1);
    }
    tstep=0;
    for(n=0;n<count;n++){
        struct dirent *ent;
        ent=list[n];	
        i=rindex(ent->d_name,'.');
        if(i==NULL) {
            continue;
        }
        if(strcmp(i+1,"vtu")==0){
            j=rindex(ent->d_name,'_');
            if(j==NULL) continue;
            number=strndup(j+1,j-i); 
            quiet=1;
            ts_fprintf(stdout,"timestep: %u filename: %s\n",atoi(number),ent->d_name);
            
            //prepare vesicle for timestep
            vesicle=restoreVesicle(ent->d_name);
            cstlist=init_cluster_list();
            clusterize_vesicle(vesicle,cstlist);

            //calculate quantities
            //vesicle_calculate_ulm2(vesicle);
            vesicle_volume(vesicle);
            vesicle_area(vesicle);
            gyration_eigen(vesicle,&l1,&l2,&l3);
            Nbw_Nb= (ts_double)count_bonds_with_energy(vesicle->blist)/(ts_double)vesicle->blist->n;
            hbar=vesicle_meancurvature(vesicle)/vesicle->area;	
            mean_size=mean_cluster_size(vesicle,cstlist);
            var_size=var_cluster_size(vesicle,cstlist,mean_size);
            perim=get_full_perimeter(vesicle);
            fprintf(stdout,"%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,\n",
                    atoi(number),vesicle->volume, vesicle->area,l1,l2,l3, 
                    Nbw_Nb,hbar, mean_size, sqrt(var_size),perim);
            tstep++;
            write_histogram_data(atoi(number), vesicle, cstlist);
            free(number);
            cluster_list_free(cstlist);
            tape_free(vesicle->tape);
            vesicle_free(vesicle);
                }
        }
    for (n=0; n<count; n++){
        free(list[n]);
    }
    
    free(list);
    return 0;
}


/* vim: set ts=4 sts=4 sw=4 noet : */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "io.h"
#include "general.h"
#include "timestep.h"
#include "vertexmove.h"
#include "bondflip.h"
#include "frame.h"
#include "io.h"
#include "stats.h"
#include "sh.h"
#include "shcomplex.h"
#include "vesicle.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include "energy.h"
#include "vertex.h"


ts_bool run_simulation(ts_vesicle *vesicle, ts_massive_idx mcsweeps, ts_idx inititer, ts_idx iterations, ts_idx start_iteration){
    ts_idx i,k;
    ts_massive_idx j;
    ts_double kc1=0,kc2=0,kc3=0,kc4=0;
    ts_double l1,l2,l3,vmsr,bfsr, vmsrt, bfsrt; //gyration eigenvalues, statistics succsess rate
    ts_ulong epochtime;
    ts_double max_z;
    FILE *fd3=NULL;
     char filename[10000];
    //struct stat st;
    strcpy(filename,command_line_args.path);
    strcat(filename,"statistics.csv");
    //int result = stat(filename, &st);
    FILE *fd;
    if(start_iteration==0)
        fd=fopen(filename,"w");
    else
        fd=fopen(filename,"a");
    if(fd==NULL){
        fatal("Cannot open statistics.csv file for writing",1);
    }
    if(start_iteration==0)
        fprintf(fd, "Epoch,OuterLoop,VertexMoveSucessRate,BondFlipSuccessRate,Volume,Area,lamdba1,lambda2,lambda3,Kc(2-9),Kc(6-9),Kc(2-end),Kc(3-6),\n");


    /* RANDOM SEED SET BY CURRENT TIME */
    epochtime=get_epoch();	
    if (vesicle->tape->random_seed!=0){
        srand48(vesicle->tape->random_seed);
        ts_fprintf(stdout,"simulation seed = %lu\n",vesicle->tape->random_seed);
    }
    else{
        srand48(epochtime);
        ts_fprintf(stdout,"simulation seed %lu\n",epochtime);
        vesicle->tape->random_seed = epochtime;
    }

    //initialize vesicle properties
    if(vesicle->tape->allow_center_mass_movement==0){
        centermass(vesicle);
    }
    cell_occupation(vesicle);
    vesicle_volume(vesicle); //needed for constant volume at this moment
    vesicle_area(vesicle); //needed for constant area at this moment
    if (vesicle->tape->V0>1){
        // override equilibrium/constraint volume by tape
        V0=vesicle->tape->V0;
    }
    if (vesicle->tape->A0>1){
        // override equilibrium/constraint volume by tape
        A0=vesicle->tape->A0;
    }
    if(V0<0.000001) {
        V0=vesicle->volume; 
    }
    ts_fprintf(stdout,"Setting volume V0=%.17f\n",V0);
    if(A0<0.000001){
        A0=vesicle->area;
    }
    ts_fprintf(stdout,"Setting area A0=%.17f\n",A0);
    epsvol=4.0*sqrt(2.0*M_PI)/pow(3.0,3.0/4.0)*V0/pow(vesicle->tlist->n,3.0/2.0);
    // printf("epsvol=%e\n",epsvol);
    epsarea=A0/(ts_double)vesicle->tlist->n;


    if(start_iteration<inititer) ts_fprintf(stdout, "Starting simulation (first %d x %d MC sweeps will not be recorded on disk)\n", inititer, mcsweeps);
    
    //PRIMARY SIMULATION LOOP
    for(i=start_iteration;i<inititer+iterations;i++){
        vmsr=0.0;
        bfsr=0.0;


        //plane confinement
        if(vesicle->tape->plane_confinement_switch){
            max_z=-1e10;
            for(k=0;k<vesicle->vlist->n;k++){
                if(vesicle->vlist->vtx[k]->z > max_z) max_z=vesicle->vlist->vtx[k]->z;
            }
            vesicle->confinement_plane.force_switch=0;
            if(max_z>=vesicle->tape->plane_d){
                ts_fprintf(stdout, "Max vertex out of bounds (z>=%e). Plane set to max_z = %e.\n",vesicle->tape->plane_d,max_z);
                vesicle->confinement_plane.z_max = max_z;
                vesicle->confinement_plane.force_switch=1;
            } else {
                vesicle->confinement_plane.z_max=vesicle->tape->plane_d;
            }

            vesicle->confinement_plane.z_min=vesicle->tape->adhesion_z - 2*vesicle->tape->adhesion_radius;

            if(vesicle->confinement_plane.force_switch) ts_fprintf(stdout,"Squeezing with force %e.\n",vesicle->tape->plane_F);
        }
        //end plane confinement

        //total force
        total_force_on_vesicle(vesicle);

        // MAIN INNER LOOP
        // MONTE CARLO SWEEP
        for(j=0;j<mcsweeps;j++){
            single_timestep(vesicle, &vmsrt, &bfsrt);
            vmsr+=vmsrt;
            bfsr+=bfsrt;
        }
        
        //post sweep processing: statistics, save state
        vmsr/=(ts_double)mcsweeps;
        bfsr/=(ts_double)mcsweeps;
        if(vesicle->tape->allow_center_mass_movement==0){
            centermass(vesicle);
        }
        cell_occupation(vesicle);
        dump_state(vesicle,i);
        vesicle_volume(vesicle); //calculates just volume. 
        vesicle_area(vesicle); //calculates area.
        if(vesicle->tape->volume_switch==0){
            V0=vesicle->volume;
        }
        if(vesicle->tape->area_switch<=1){
            A0=vesicle->area;
        }

        //update status file
        if(i>=inititer){
            write_vertex_xml_file(vesicle,i-inititer,NULL);
            write_master_xml_file(command_line_args.output_fullfilename);
            epochtime=get_epoch();			
            gyration_eigen(vesicle, &l1, &l2, &l3);
  
            fprintf(fd, "%lu,%u,%e,%e,%1.16e,%1.16e,%1.16e,%1.16e,%1.16e,%1.16e,%1.16e,%1.16e,%1.16e,\n",epochtime,i,vmsr,bfsr,vesicle->volume, vesicle->area,l1,l2,l3,kc1, kc2, kc3,kc4);

            fflush(fd);	
            //	sprintf(filename,"timestep-%05d.pov",i-inititer);
            //	write_pov_file(vesicle,filename);
        } //end if(inititer....)
            fd3=fopen(".status","w"); //write status file when everything is written to disk.
            if(fd3==NULL){
                fatal("Cannot open .status file for writing",1);
        }
        fprintf(fd3,"%d",i);
        fclose(fd3);

        ts_fprintf(stdout,"Done %d out of %d iterations (x %d MC sweeps).\n",i+1,inititer+iterations,mcsweeps);
    }
    fclose(fd);
    //	if(fd2!=NULL) fclose(fd2);
    return TS_SUCCESS;
}

ts_bool single_timestep(ts_vesicle *vesicle,ts_double *vmsr, ts_double *bfsr){
    ts_bool retval;
    ts_double rnvec[3];
    ts_idx i,j,b;
    ts_uint vmsrcnt=0;

    for(i=0;i<vesicle->vlist->n;i++){
        //rnvec[0]=drand48();
        //rnvec[1]=drand48();
        //rnvec[2]=drand48();
        if (vesicle->vlist->vtx[i]->type == is_ghost_vtx) {
            //skip! ghost vertex should not move**
            retval=TS_FAIL;
        }
        else {
            retval=single_verticle_timestep(vesicle,vesicle->vlist->vtx[i]);
        }
        if(retval==TS_SUCCESS) vmsrcnt++;
    }
    
    ts_int bfsrcnt=0;
    for(i=0;i<3*vesicle->vlist->n;i++){
        b=rand() % vesicle->blist->n;
        //find a bond and return a pointer to a bond...
        //call single_bondflip_timestep...
        retval=single_bondflip_timestep_ordered(vesicle,vesicle->blist->bond[b],rnvec);
        // b++; retval=TS_FAIL;
        if(retval==TS_SUCCESS) bfsrcnt++;   
    }


    for(i=0;i<vesicle->poly_list->n;i++){
        for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
            if (vesicle->poly_list->poly[i]->vlist->vtx[j]->type == is_ghost_vtx) {
                //skip! ghost vertex should not move**
                retval = 0;
            }
            else {
                rnvec[0]=drand48();
                rnvec[1]=drand48();
                rnvec[2]=drand48();
                retval=single_poly_vertex_move(vesicle,vesicle->poly_list->poly[i],vesicle->poly_list->poly[i]->vlist->vtx[j],rnvec);	
            }
        }
    }


    for(i=0;i<vesicle->filament_list->n;i++){
        for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
            if (vesicle->filament_list->poly[i]->vlist->vtx[j]->type == is_ghost_vtx) {
                retval = 0;
            }
            else {
                rnvec[0]=drand48();
                rnvec[1]=drand48();
                rnvec[2]=drand48();
                retval=single_filament_vertex_move(vesicle,vesicle->filament_list->poly[i],vesicle->filament_list->poly[i]->vlist->vtx[j],rnvec);
            }	
        }
    }

 


    // printf("Bondflip success rate in one sweep: %d/%d=%e\n", cnt,3*vesicle->blist->n,(double)cnt/(double)vesicle->blist->n/3.0);
    *vmsr=(ts_double)vmsrcnt/(ts_double)vesicle->vlist->n;
    *bfsr=(ts_double)bfsrcnt/(ts_double)vesicle->vlist->n/3.0;
    // vesicle_volume(vesicle);
    // fprintf(stderr,"Volume after TS=%1.16e\n", vesicle->volume);
    return TS_SUCCESS;
}



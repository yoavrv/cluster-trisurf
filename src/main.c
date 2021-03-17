/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"
#include "cell.h"
#include "vesicle.h"
#include "io.h"
#include "initial_distribution.h"
#include "frame.h"
#include "timestep.h"
#include "poly.h"
#include "sh.h"
#include "shcomplex.h"
#include "dumpstate.h"
#include "restore.h"

#include <fcntl.h>
/** Entrance function to the program
  * @param argv is a number of parameters used in program call (including the program name
  * @param argc is a pointer to strings (character arrays) which holds the arguments
  * @returns returns 0 on success, any other number on fail.
*/

int main(int argv, char *argc[]){
	ts_vesicle *vesicle;
	ts_tape *tape;
	ts_uint start_iteration=0;
	force_from_tape=0;
	/* Area and volume for constant area and constant volume are initialized to be zero */
	A0=0;
	V0=0;
	/* create lock file */
	createPidFile("ts_trisurf",".lock",0);
	parse_args(argv,argc); // sets global variable command_line_args (defined in io.h)
	ts_fprintf(stdout,"TRISURF-NG v. %s, compiled on: %s %s.\n", TS_VERSION, __DATE__, __TIME__);
	ts_fprintf(stdout,"Programming done by: Samo Penic and Miha Fosnaric\n");
	ts_fprintf(stdout,"Released under terms of GPLv3\n");
	ts_fprintf(stdout,"Starting program...\n\n");
//	vesicle = parseDump("timestep_000000.vtu");
//		run_simulation(vesicle, vesicle->tape->mcsweeps, vesicle->tape->inititer, vesicle->tape->iterations, 1);

    if(command_line_args.dump_from_vtk[0]!=0){
		ts_fprintf(stdout,"************************************************\n");
		ts_fprintf(stdout,"**** Restoring vesicle from VTK points list ****\n");
		ts_fprintf(stdout,"************************************************\n\n");
		vesicle = parseDump(command_line_args.dump_from_vtk);
//		write_vertex_xml_file(vesicle,9999); // here you can test if restoration and rewritting results in the same dump file. Only the date od creation of dump file must differ.
		tape = vesicle->tape;
		int arguments_no;
		FILE *fd=fopen(".status","r");
		if(fd!=NULL){
			arguments_no=fscanf(fd,"%u", &start_iteration);
			if(arguments_no==0){
				ts_fprintf(stdout,"No information of start iteration in .status file\n");
				}
			fclose(fd);
			start_iteration++; 
		}
		else
			ts_fprintf(stdout,"No .status file. The iteration count will start from 0\n");
/* Here you should read new tape file, reassign some values in vertex from the tape and assign read tape to vesicle->tape */
//        tape=parsetape(command_line_args.tape_fullfilename);
  //      vesicle=vtk2vesicle(command_line_args.dump_from_vtk,tape);
    }
	else if(command_line_args.force_from_tape){
		ts_fprintf(stdout,"************************************************\n");
		ts_fprintf(stdout,"**** Generating initial geometry from tape *****\n");
		ts_fprintf(stdout,"************************************************\n\n");
		tape=parsetape(command_line_args.tape_fullfilename);
		vesicle=create_vesicle_from_tape(tape);
	} else {

		ts_fprintf(stdout,"**********************************************************************\n");
		ts_fprintf(stdout,"**** Recreating vesicle from dump file and continuing simulation *****\n");
		ts_fprintf(stdout,"**********************************************************************\n\n");
		tape=parsetape(command_line_args.tape_fullfilename);
		vesicle=restore_state(&start_iteration);
        if(vesicle==NULL){
            ts_fprintf(stderr, "Dump file does not exist or is not a regular file! Did you mean to invoke trisurf with --force-from-tape option?\n\n");
            return 1;
        }
		// nove vrednosti iz tapea...
		vesicle->bending_rigidity=tape->xk0;
		vtx_set_global_values(vesicle);
		vesicle->pswitch =tape->pswitch;
		vesicle->pressure=tape->pressure;
		vesicle->dmax=tape->dmax*tape->dmax;
		poly_assign_filament_xi(vesicle,tape);
        free(vesicle->tape);
        vesicle->tape=tape;
		vesicle->clist->dmin_interspecies = tape->dmin_interspecies*tape->dmin_interspecies;



        /* spherical harmonics */
        if(tape->shc>0){
	        vesicle->sphHarmonics=complex_sph_init(vesicle->vlist,tape->shc);
        }
        else {
            vesicle->sphHarmonics=NULL;
        }

		if(command_line_args.reset_iteration_count) start_iteration=tape->inititer;
		else start_iteration++;

		if(start_iteration>=tape->iterations){
			ts_fprintf(stdout, "Simulation already completed. if you want to rerun it try with --force-from-tape or --reset-iteration-count\n\n");
			return 0;
		}

	/* if requested in tape, we can have smaller number of polymeres attached to membrane than the number of polymeres in dump file */
		if(vesicle->tape->npoly != vesicle->poly_list->n){

		ts_fprintf(stdout,"(INFO) the number of polymeres attached to membrane in tape is different than a number of polymeres in dump file!\n");
		if(vesicle->tape->npoly > vesicle->poly_list->n){
			ts_fprintf(stdout,"(INFO) It is possible to decrease the number of polymeres on the membrane, but it is not allowed to increase its number. The maximal allowed number in tape is %d The execution of program will terminate!\n",vesicle->poly_list->n);
			fatal("Terminating due to increase of number of polymeres",1);
		} else {
			remove_random_polymeres(vesicle->poly_list, vesicle->poly_list->n - vesicle->tape->npoly);
			ts_fprintf(stdout,"(INFO)\n(INFO) The new number of polymeres from tape is %d.\n\n",vesicle->poly_list->n);

		}
		}
	}
			//printf("nucleus coords: %.17e %.17e %.17e\n",vesicle->nucleus_center[0], vesicle->nucleus_center[1], vesicle->nucleus_center[2]);
//	write_vertex_xml_file(vesicle,0);
//	exit(1);

			//write_vertex_xml_file(vesicle,1000);
	run_simulation(vesicle, tape->mcsweeps, tape->inititer, tape->iterations, start_iteration);
	write_master_xml_file(command_line_args.output_fullfilename);
	write_dout_fcompat_file(vesicle,"dout");
	vesicle_free(vesicle);
	tape_free(tape);
	return 0; //program finished perfectly ok. We return 0.
} 

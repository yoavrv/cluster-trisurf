/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _IO_H
#define _IO_H
#include <confuse.h>
#include <general.h>
/** @ Global variables for I/O operations like filenames etc. */
/*static char mastername[1024];
static char prefixname[1024];
static ts_bool restore=0;
static char tape[1024]; */
extern char tapetxt[128000]; //stores text file of the tape

typedef struct{
	ts_int force_from_tape;
	ts_int reset_iteration_count;
    char path[1024]; //path where all files should be added
    char output_fullfilename[1024]; //name of the master file
    char dump_fullfilename[1024]; //name of the dump file
    char tape_fullfilename[1024]; //name of the tape file
    char tape_templatefull[1024]; //name of the tape template file
    char tape_opts[1000]; //commandline tape options
    char dump_from_vtk[1024];
} ts_args;

extern ts_args command_line_args;

ts_bool parse_args(int argc, char **argv);



/** @brief Prints the position of vertices for the whole list
 *  
 *  The function is meant more or less as a debug tool, but can be used in production
 *  environment aswell.
 *  the output is in form of idx: x y z
 *  	@param *vlist is a structure holding information on vertex list.
 *	@returns TS_SUCCESS on successful execution, TS_FAIL otherwise.
*/

ts_bool print_vertex_list(ts_vertex_list *vlist);

/** @brief Prints the neighbours of all the vertices
 *  
 *  The function is meant more or less as a debug tool, but can be used in production
 *  environment aswell.
 *  the output is in form of idx(number of neighbours): (x1,y1,z1) (x2,y2,z2) ...
 *  	@param *vlist is a structure holding information on vertex list.
 *	@returns TS_SUCCESS on successful execution, TS_FAIL otherwise.
*/
ts_bool print_vertex_neighbours(ts_vertex_list *vlist);


/** @brief Function outputs the vetex list file to comply to old fortran format
 *  
 *  	@param *vlist is a list of vertices
 *  	@param *filename is a name of the output file to be created (note that if the file already
 *  exists it will be overwritten.
 */
ts_bool write_vertex_fcompat_file(ts_vertex_list *vlist,ts_char *filename);
ts_bool fprint_vertex_list(FILE *fh,ts_vertex_list *vlist);
ts_bool fprint_tristar(FILE *fh, ts_vesicle *vesicle);
ts_bool fprint_triangle_list(FILE *fh, ts_vesicle *vesicle);
ts_bool fprint_vertex_data(FILE *fh,ts_vertex_list *vlist);
ts_bool fprint_bonds(FILE *fh,ts_vesicle *vesicle);
ts_bool write_dout_fcompat_file(ts_vesicle *vesicle, ts_char *filename);
ts_bool read_tape_fcompat_file(ts_vesicle *vesicle, ts_char *filename);
ts_bool write_adhesion_position(FILE *fh, ts_vesicle *vesicle);


/** @brief Outputs file in vtk format, compatible with paraview.
 *
 *	@param *vlist is a list of vertices
 *	@param *filename is a name of the output file. If exists, it will be overwritten
 *	@param *text is a description line (max. 255 characters) to be included in the file
 */
ts_bool write_vertex_vtk_file(ts_vesicle *vesicle,ts_char *filename, ts_char *text);
ts_bool write_vertex_xml_file(ts_vesicle *vesicle, ts_idx timestepno, ts_cluster_list *cstlist);
ts_bool write_master_xml_file(ts_char *filename);
ts_bool write_pov_file(ts_vesicle *vesicle, char *filename);

ts_tape *parsetape(char *filename);
ts_tape *parsetapebuffer(char *buffer);
ts_bool tape_free(ts_tape *tape);
ts_bool getcmdline_tape(cfg_t *cfg, char *opts);
ts_bool cmdline_to_tape(cfg_t *cfg, char *key, char *val);
ts_bool print_help(FILE *fd);
ts_bool dump_state(ts_vesicle *vesicle, ts_idx iteration);
ts_vesicle *restore_state(ts_idx *iteration, ts_tape* tape);

ts_bool update_tapetxt(char *tape_txt, char *cmd_line_tape_args);
#endif

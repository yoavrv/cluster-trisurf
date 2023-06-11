/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _GENERAL_H
#define _GENERAL_H

#include <stdio.h>
#include <stdarg.h>
#include <gsl/gsl_complex.h>
#include <math.h> // added #ifndef M_PI to stop vscode from complaning about undefined M_PI, but this should load it anyway and avert the ifndef

/** @brief This is a header file, defining general constants and structures.
  * @file header.h
  * @author Samo Penic
  * @date 5.3.2001
  * 
  * Header file for general inclusion in all the code, defining data structures
  * and general constans. All datatypes used in the code is also defined here.
  *
  * Miha: branch trisurf-polyel
  */

/* Defines */
/** @brief Return value of type bz_bool that indiceates successful function finish 
  *
  * Function usualy return some value, which are the result of certain operation. Functions that don't
  * return any parameters can return value, that indicates if the function call finished successfully.
  * In case of successful function run, the functions should return TS_SUCCESS to the caller. This define
  * is set here to get uniformity among all the functions used in program.
  *
  * Example of usage:
  *		ts_boot somefunction(ts_int param1, ....){
  *			...
  *			return TS_SUCCESS;
  *		}
  */
#define TS_SUCCESS 0

/** @brief Return value of type bz_bool that indicates unsuccessful function finish 
  *
  * Function usualy return some value, which are the result of certain operation. Functions that don't
  * return any parameters can return value, that indicates if the function call finished successfully.
  * In case of unsuccessful function run, the functions should return TS_FAIL to the caller. This define
  * is set here to get uniformity among all the functions used in program.
  *
  * Example of usage:
  *
  *		ts_boot somefunction(ts_int param1, ....){
  *			...
  *			return TS_FAIL;
  *		}
  */
#define TS_FAIL 1

/* CONSTANTS */

#define TS_ID_FILAMENT 1

/* DATA TYPES */
/** @brief Sets the default datatype for ts_double
 *
 * Requred for some functions to work, like "pow" from math.h. If ts_double is defined as
 * float program must run with "powf". Where type dependant function is used it checks this
 * define directive to decide which version to compile in. Available options
 *
 *	TS_DOUBLE_FLOAT
 *	TS_DOUBLE_DOUBLE
 *	TS_DOUBLE_LONGDOUBLE
*/
#define TS_DOUBLE_DOUBLE

/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_int (uses int)
 */
typedef int ts_int;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_uint (uses unsigned int)
 */
typedef unsigned int ts_uint;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_long (uses long)
 */
typedef long ts_long;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_ulong (uses unsigned long)
 */
typedef unsigned long ts_ulong;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_float (uses float)
 */
typedef float ts_float;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_double (uses double)
 */
typedef double ts_double;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_char (uses char)
 */
typedef char ts_char;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_uchar (uses unsigned char)
 */
typedef unsigned char ts_uchar;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_bool (uses char)
 */
typedef char ts_bool;

/* * Yoav: Can't remember what's what with the indice
* I'm just going to typedef an index that will be used everywhere we have arr[i]
so I can stop confusing int and uint have consistent wraparound behaviors
out of bound might be simpler: it's always i greater than vec->max
*/
// currently 21/10/21: mainly excluding math-capabale things: nshell,
// ncmaxes (because heck if I know what goes in the calculations there) and cellno (by extension),
// and spherical harmonics indexes.
// also not touching the FORTRAN indexes with a ten foot pole in intial distribution
typedef signed char ts_small_idx; // for tiny indices <<64, mostly short circular neighbor vectors
typedef signed int ts_idx; // for big indices
typedef signed long ts_massive_idx; //for mcsweeps

// small index functions
extern inline ts_small_idx next_small(ts_small_idx i, ts_small_idx max);
extern inline ts_idx next_idx(ts_idx i, ts_idx max);
extern inline ts_small_idx prev_small(ts_small_idx i, ts_small_idx max);
extern inline ts_idx prev_idx(ts_idx i, ts_idx max);

// yoav: I'm going to bravely denote the cell indices with their own typedef
typedef unsigned int ts_cell_idx; // index of a cell

// for the many flags: ts_flag to distinguish from ts_bool
// with the idea that ts_bool is only 0/1 like TS_SUCCESS or TS_FAIL
typedef char ts_flag;
typedef int ts_big_flag;

/** Enums for various options
 * 
 * These are intended to work by
 * x&is_option_A  (i.e. "x and is option A")
 * y&to_B      (i.e. "y and to do B")
 * or
 * z==model_vicsek_1overR
 * 
 * */
enum vertex_type {
    is_bonding_vtx=1, // bonding vertex, form bond with other bonding vertices
    is_active_vtx=2, // active vertex under normally directed force
    is_adhesive_vtx=4, // adhesive vertex, subject to adhesion energy near surfaces
    is_anisotropic_vtx=8, // anisotropic vertex, requires computing full curvature characteristic
    is_reserved_0_vtx=16, // reserved type
    is_vicsek_vtx=32, // vertex under vicsek neighbor-influenced force
    is_edge_vtx=64, // edge vertex has unordered tristars
    is_ghost_vtx=-128, // ghost vertex can only be moved artificially
}; 
enum bond_model_type {
    is_bonding_type_specific=1,
    is_anisotropic_bonding_nematic=2,
};
enum curvature_model_type{
    to_disable_calculate_laplace_beltrami=64, // The original isotropic mean curvature calculation
    to_calculate_sum_angle=1, // "isotropic" gaussian curvature calculation
    to_calculate_shape_operator=2, // new anisotropic calculation. if to_use is not enabled, this is saved but not used!
    to_update_director_shapeless=4, // update director without shape operator
    to_use_shape_operator_energy=8, // actually use the new energy, rather than just calculate and save.
    to_use_shape_for_anisotropy_only=16, // use the shape method, but only for anisotropic vertices
    to_not_rotate_directors=32, // do not rotate directro as a monte carlo step
    to_use_sum_angle_for_kx2_only=128, // only calculate the sum angle formula for vtx with kx2!=0
    
    model_laplace_beltrami_only=0, // the original method
    model_isotropic=1, // calculate gaussian curvature energy
    model_shape_operator_only=74, // use shape operator only
    model_debug_old_energy=7, // calculate everything but use old energy
    model_debug_new_energy=15, // calculate everything but use new shape operator energy
    model_debug_parallerl_transport_directors=35, // prevent director random move
    model_debug_assumed_final=146, // anisotropic shape method only for anisotropic vertices, old method otherwise
    
};
enum force_model_type{
    model_vertex_meanW=0, // regular bonding model
    model_active_neigh_interfer=1, // force is proportional to # non-active neighbors
    model_concave_neigh_interfer=2, // force is proportional to # non-concave neighbors
    model_concave_neigh_disable=3, // force is disabled by any concave neighbor
    is_vicsek_model=16, // use a vicsek model. All vicsek model should have (model_vicsek_x&is_vicsek_model) true
    model_vicsek_1overR=17, // use vicsek model with 1/R weights

};
enum adhesion_model_type{
  adhesion_step_potential=1,
  adhesion_parabolic_potential=2,
  adhesion_y_anisotropy=4,

};
enum adhesion_geometry_type{
  model_plane_potential=1,
  model_spherical_potential=2,
  model_cylindrical_potential=3,
  model_x_sinosoidal_potential=4,
  model_plane_potential_with_spots=5,
  model_xy_sinosoidal_potential=6,
};


/* STRUCTURES */


/** @brief Data structure for keeping the coordinates in selected coordinate
 * system
 */
#define TS_COORD_CARTESIAN 0
#define TS_COORD_SPHERICAL 1
#define TS_COORD_CYLINDRICAL 2
// ^ coordinates in the vtu file?

/** @brief Prevent spikiness of triangles by imposing a minimum angle between them
 * We measure the normal between neighboring two triangles n1*n2=cos(theta)
 * And impose this is greater than 
 * */
// #define MIN_INTERTRIANGLE_ANGLE_COSINE 0

typedef struct {
    ts_double e1;
    ts_double e2;
    ts_double e3;
    ts_uint coord_type;
} ts_coord;

/** @brief Data structure of all data connected to a vertex
 *
 *  ts_vertex holds the data for one single point (bead, vertex). To understand how to use it
 *  here is a detailed description of the fields in the data structure. */
struct ts_vertex {
        // degrees of freedom
        ts_double x; // position
        ts_double y; 
        ts_double z;
        ts_double dx; // director vector
        ts_double dy;
        ts_double dz;
        // innate properties (except type, at the bottom)
        ts_double xk;  //bending modulus \kappa
        ts_double xk2; //Gaussian bending modulus \kappa_G
        ts_double c;
        ts_double d;    // spontaneous curvature deviator
        ts_double w;    //binding
        ts_double ad_w; // adhesive surface bonding
        ts_double f;    // force magnitude
        // calculated properties (dynamical)
        ts_double nx; // normal vector
        ts_double ny;
        ts_double nz;
        ts_double nx2; // normal vector
        ts_double ny2;
        ts_double nz2;
        ts_double fx; // force vector
        ts_double fy;
        ts_double fz;
        ts_double mean_curvature;
        ts_double gaussian_curvature;
        ts_double mean_energy;
        ts_double gaussian_energy;
        ts_double energy;
        ts_double eig0[3]; // principal curvature direction (largest curvature)
        ts_double eig1[3]; // principal curvature direction (smallest curvature)
        ts_double eig2[3]; // shape tensor normal (0)
        ts_double eig_v0; // c1 from shape tensor (largest curvature)
        ts_double eig_v1; // c2 from shape tensor (smallest curvature)
        ts_double eig_v2; // 0 noraml eigenvalue: 3x3 shape tensor has 0 eigenvalue normal
        ts_double mean_curvature2; // from shape tensor
        ts_double gaussian_curvature2; // from shape tensor
        ts_double mean_energy2; // from shape tensor
        ts_double gaussian_energy2; // from shape tensor
        // graph
        struct ts_vertex **neigh; /**< The pointer that holds neigh_no pointers to this structure. */
        struct ts_triangle **tristar; /**< The list of triangles this vertex belongs to. This is an array of pointers to ts_triangle structure of tristar_no length */
        struct ts_bond **bond; /**< Array of pointers of lenght bond_no that stores information on bonds. */
        struct ts_cell *cell; /**< Which cell do we belong to? */
        struct ts_poly *grafted_poly;
        struct ts_cluster *cluster;
        ts_idx idx; //vertex index
        ts_idx id;	//filament index
        // number of neighbors, of all types
        // these never really go beyond even 12, so we can squeeze to char
        // (n.n is capped at 19 at vertexmove backup[], 10 at (currently broken) constvol)
        ts_small_idx neigh_no; /**< The number of neighbours. */
        ts_small_idx tristar_no; //number of triangle-neighbors
        ts_small_idx bond_no;
        /* 1st bit: bonds, 2nd bit: active, 
        3rd bit: adhesive, 4th bit: anisotropic, 
        5th bit: reserved, 6th bit: vicsek 
        7th bit: edge, 8th bit: ghost*/
        ts_flag type; 
};

typedef struct ts_vertex ts_vertex;

typedef struct {
    ts_vertex **vtx;
    ts_idx n;

} ts_vertex_list;

struct ts_bond {
    ts_double bond_length;
    ts_double energy;
    ts_double x,y,z;
    ts_vertex *vtx1;
    ts_vertex *vtx2;
    ts_idx idx;
};
typedef struct ts_bond ts_bond;

struct ts_bond_list {
    ts_bond **bond;
    ts_idx n;
};
typedef struct ts_bond_list ts_bond_list;

struct ts_triangle {
    ts_double xnorm;
    ts_double ynorm;
    ts_double znorm;
    ts_double xcirc;
    ts_double ycirc;
    ts_double zcirc;
    ts_double area; 
    ts_double volume; // tetrahedron frrom (0,0,0)
    ts_double energy;
    ts_vertex *vertex[3];
    struct ts_triangle **neigh;
    ts_idx idx;
    ts_small_idx neigh_no;
};
typedef struct ts_triangle ts_triangle;

struct ts_triangle_list{
    ts_double a0;
    ts_triangle **tria;
    ts_idx n;
};
typedef struct ts_triangle_list ts_triangle_list;


typedef struct ts_cell {
    ts_vertex **vertex;
    ts_cell_idx idx;
    ts_small_idx nvertex;
} ts_cell; 

typedef struct ts_cell_list{
    ts_double dcell;    // density (1/size) of each cell. vtx->x*dcell+shift(+-1?) is the vtx cell coordinate
    ts_double shift[3]; // shift of x,y,z to the center
    ts_double dmin_interspecies; // ? minimum distance between non-connected vertices squared ?
    ts_cell **cell;
    ts_cell_idx ncmax[3]; // no idea what kind of indexing goes on here
    ts_cell_idx cellno;
    ts_small_idx max_occupancy;
} ts_cell_list;


typedef struct {
    ts_double *vtx_relR; //stuff taken from vertex
    ts_double *vtx_solAngle;
    ts_double **ulm;
    ts_double **co;
    ts_double ***Ylmi;
    ts_double **sumUlm2;
    gsl_complex **ulmComplex;
    ts_idx n_vtx; //vlist->n
    ts_uint l;
    ts_uint N;

} ts_spharm;



struct ts_poly {
    ts_double k;
    ts_vertex_list *vlist;
    ts_bond_list *blist;
    ts_vertex *grafted_vtx;
};
typedef struct ts_poly ts_poly;


struct ts_poly_list {
    ts_poly **poly;
    ts_idx n;
};
typedef struct ts_poly_list ts_poly_list;


typedef struct{
    ts_double z_max;
    ts_double z_min;
    ts_bool force_switch;
} ts_confinement_plane;


typedef struct {
    char* tape_text;
    ts_double R_nucleus;
    ts_double R_nucleusX;
    ts_double R_nucleusY;
    ts_double R_nucleusZ;
    ts_double xkA0; // area change modulus
    ts_double xkV0; // volume change modulus
    ts_double V0;
    ts_double A0;
    ts_double Vfraction; // equilibrium_reduced_volume
    ts_double constvolprecision;
    ts_double xk0; // bending modulus
    ts_double xk2; // second bending modulus (Gaussian/ deviatoric?)
    ts_double dmax;
    ts_double dmin_interspecies;
    ts_double stepsize;
    ts_double kspring;
    ts_double xi;
    ts_double pressure;
    ts_double c0;
    ts_double d0; // spontaneous deviator
    ts_double w;
    ts_double F;
    ts_double plane_d;
    ts_double plane_F;
    ts_double vicsek_strength;
    ts_double vicsek_radius;
    ts_double adhesion_z;
    ts_double adhesion_cutoff;
    ts_double adhesion_strength;
    ts_double adhesion_radius;
    ts_double adhesion_scale;
    ts_double adhesion_factor;
    ts_double max_dihedral_angle_cosine; // prevent spikiness of triangles by imposing a minimum dihedral angle
    ts_massive_idx mcsweeps;
    ts_ulong random_seed;
    ts_idx iterations;
    ts_idx inititer;
    ts_uint nshell;
    ts_uint ncxmax;
    ts_uint ncymax;
    ts_uint nczmax;
    ts_idx number_of_vertices_with_c0;
    ts_idx npoly;
    ts_idx nmono;
    ts_idx internal_poly;
    ts_idx nfil;
    ts_idx nfono;
    ts_uint shc; // related to max l of the spherical harmonics
    ts_bool pressure_switch;
    ts_bool volume_switch;
    ts_bool area_switch;
    ts_bool quiet;
    ts_bool plane_confinement_switch;
    ts_bool allow_center_mass_movement;
    ts_bool force_balance_along_z_axis;
    ts_flag adhesion_geometry; // geometry of adhesion (none, plane, sphere, cylinder)
    ts_flag adhesion_model; // adhesion (none, step potential, parabolic potential)
    ts_flag bond_model;
    ts_flag curvature_model;
    ts_flag force_model;
} ts_tape;




typedef struct {
    // tape
    ts_tape *tape;
    // degrees of freedom
    ts_double R_nucleusX;
    ts_double R_nucleusY;
    ts_double R_nucleusZ;
    // innate properties
    ts_double dmax;
    ts_double stepsize;
    ts_double pressure;
    ts_double R_nucleus;
    // calculated properties (dynamical)
    ts_double cm[3];
    ts_double fx;
    ts_double fy;
    ts_double fz;
    ts_double volume;
    ts_double area;
    ts_double nucleus_center[3];
    // secondary structs
    ts_spharm *sphHarmonics;
    // Polymers outside the vesicle and attached to the vesicle membrane (polymer brush):
    ts_poly_list *poly_list;
    // Filaments inside the vesicle (not attached to the vesicel membrane:
    ts_poly_list *filament_list;
    ts_vertex_list *vlist;
    ts_bond_list *blist;
    ts_triangle_list *tlist;
    ts_cell_list *clist;
    ts_confinement_plane confinement_plane;
} ts_vesicle;



struct ts_cluster{
    ts_vertex **vtx;
    ts_idx nvtx;
    ts_idx idx;
};

typedef struct ts_cluster ts_cluster;

typedef struct{
    ts_cluster **cluster;
    ts_idx n;
} ts_cluster_list;

// list of "seen" vertex for layer-based breadth-first search
// a vertex list, with 4 locations: 
// * n_top, the top of the list, where new vertices are added
// * n_next, first vertex of the layer under construction
// * n_curr, first vertex of the last completed layer, used to construct the next
// * n_prev, first vertex of the completed layer before the current one
// notice! seen_vertex->vtx[seen_vertex->n_top] is NEVER a valid vertex
typedef struct {
    ts_vertex **vtx;
    ts_idx n_prev;
    ts_idx n_curr;
    ts_idx n_next;
    ts_idx n_top;
    ts_idx size;
} ts_seen_vertex;

/* GLOBAL VARIABLES */

extern ts_bool quiet;
extern ts_double V0;
extern ts_double A0;
extern ts_double epsvol;
extern ts_double epsarea;

// global structure from io.h

typedef struct{
	ts_int force_from_tape;
	ts_int reset_iteration_count;
    char path[1024]; //path where all files should be added
    char output_fullfilename[1024]; //name of the master file
    char dump_fullfilename[1024]; //name of the dump file
    char tape_fullfilename[1024]; //name of the tape file
    char tape_templatefull[1024]; //name of the tape template file
    char tape_opts[1024]; //commandline tape options
    char dump_from_vtk[1024];
} ts_args;

extern ts_args command_line_args;

/* FUNCTIONS */

/** Non-fatal error function handler:
 *      @param text is a description of an error
 *      @returns doesn't return anything
*/
void err(char *text);

/** Fatal error function handler:
 *      @param text is a description of an error
 *      @param errcode is a (non-zero) error code
 *      @returns terminates the execution of program with errcode set
*/
void fatal(char *text, ts_int errcode);

ts_uint ts_fprintf(FILE *fd, char *fmt, ...);

#define VTX(n) &(vlist->vtx[n])
#define VTX_DATA(n) vlist->vtx[n].data


/* FOR PID GENERATION ROUTINE */
#define CPF_CLOEXEC 1

int createPidFile(const char *progName, const char *pidFile, int flags);

int lockRegion(int fd, int type, int whence, int start, int len);
char *libVersion();

// ifdefs to stop vscode from complaining
#ifndef TS_VERSION
#define TS_VERSION "whatever"
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

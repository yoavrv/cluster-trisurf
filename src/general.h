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
extern inline ts_small_idx next_small(ts_small_idx i, ts_small_idx max);
extern inline ts_idx next_idx(ts_idx i, ts_idx max);
extern inline ts_small_idx prev_small(ts_small_idx i, ts_small_idx max);
extern inline ts_idx prev_idx(ts_idx i, ts_idx max);
// for the many flags: ts_flag to distinguish from ts_bool
// with the idea that ts_bool is only 0/1 like TS_SUCCESS or TS_FAIL
typedef char ts_flag;
typedef int ts_big_flag;

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
        ts_double x; /**< The x coordinate of vertex. */
        ts_double y; /**< The y coordinate of vertex. */
        ts_double z; /**< The z coordinate of vertex. */
        ts_double mean_curvature;
        ts_double gaussian_curvature; // to be determined: mean, Gaussian curvature or c1,c2
        ts_double mean_energy;
        ts_double gaussian_energy;
        ts_double energy;
        //ts_double energy_h;
        ts_double xk; //bending rigidity
        ts_double xk2; //second bending modulus (Gaussian/ deviatoric?): should be excess (Gauss-Bonet)
        ts_double w;
        ts_double c;
        ts_double nx; // normal vector
        ts_double ny;
        ts_double nz;
        ts_double nx2; // normal vector
        ts_double ny2;
        ts_double nz2;
        ts_double f;  // force
        ts_double fx; // force vector
        ts_double fy;
        ts_double fz;
        ts_double ad_w; // adhesive surface bonding
        ts_double d;  // spontaneous curvature deviator
        ts_double dx; // director vector
        ts_double dy;
        ts_double dz;
        ts_double eig0[3];
        ts_double eig1[3];
        ts_double eig2[3];
        ts_double new_c1;
        ts_double new_c2;
        ts_double eig_v0;
        ts_double eig_v1;
        ts_double eig_v2;
        ts_double mean_curvature2;
        ts_double gaussian_curvature2; 
        ts_double mean_energy2;
        ts_double gaussian_energy2;
        struct ts_vertex **neigh; /**< The pointer that holds neigh_no pointers to this structure. */
        struct ts_triangle **tristar; /**< The list of triangles this vertex belongs to. This is an array of pointers to ts_triangle structure of tristar_no length */
        struct ts_bond **bond; /**< Array of pointers of lenght bond_no that stores information on bonds. */
        struct ts_cell *cell; /**< Which cell do we belong to? */
        struct ts_poly *grafted_poly;
        struct ts_cluster *cluster;
        ts_idx idx; //vertex index
        ts_idx id;	//filament index
        // number of neighbors, of all types
        // these never really go beyond even 12, so
        // we can squeeze to char
        // (n.n is capped at 19 at vertexmove backup[], 10 at (currently broken) constvol)
        ts_small_idx neigh_no; /**< The number of neighbours. */
        ts_small_idx tristar_no; //number of triangle-neighbors
        ts_small_idx bond_no;
        /* 1st bit: bonds, 2nd bit: active, 
        3rd bit: adhesive, 4th bit: anisotropic, 
        5th bit: reserved, 6th bit: vicsek 
        7th bit: edge, 8th bit: ghost*/
        ts_flag type; 

        // apparently using a bool for the type flag does nothing, but I don't want to reserve 
        // flags for all 32 bytes of an int

};

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
enum adhesion_type{
  adhesion_step_potential=1,
  adhesion_parabolic_potential=2,

};
enum adhesion_geometry_type{
  model_plane_potential=1,
  model_spherical_potential=2,
  model_cylindrical_potential=3,

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
    ts_double area; // firstly needed for sh.c
    ts_double volume; // firstly needed for sh.c
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
    ts_uint idx;
    ts_small_idx nvertex;
} ts_cell; 

typedef struct ts_cell_list{
    ts_double dcell;
    ts_double shift;
    ts_double dmin_interspecies; // ? minimum distance between non-connected vertices squared ?
    ts_cell **cell;
    ts_uint ncmax[3]; // no idea what kind of indexing goes on here
    ts_uint cellno;
    ts_small_idx max_occupancy;
} ts_cell_list;


typedef struct {
    ts_double *vtx_relR; //stuff taken from vertex
    ts_double* vtx_solAngle;
    ts_idx n_vtx; //vlist->n
    ts_uint l;
    ts_double **ulm;
    gsl_complex **ulmComplex;
    ts_double **sumUlm2;
    ts_uint N;
    ts_double **co;
    ts_double ***Ylmi;
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
    ts_float z_max;
    ts_float z_min;
    ts_bool force_switch;
} ts_confinement_plane;


typedef struct {
    ts_double R_nucleus;
    ts_double R_nucleusX;
    ts_double R_nucleusY;
    ts_double R_nucleusZ;
    ts_double xkA0;
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
    ts_double w;
    ts_double F;
    ts_double plane_d;
    ts_double plane_F;
    ts_double vicsek_strength;
    ts_double vicsek_radius;
    ts_double adhesion_cuttoff;
    ts_double adhesion_strength;
    ts_double z_adhesion;
    ts_double adhesion_radius;
    ts_double min_dihedral_angle_cosine; // Prevent spikiness of triangles by imposing a minimum dihedral angle
    ts_double d0; // spontaneous deviator
    //  long int brezveze0;
    //	long int brezveze1;
    //	long int brezveze2;
    ts_idx iterations;
    ts_massive_idx mcsweeps;
    ts_ulong random_seed;
    ts_idx inititer;
    ts_idx number_of_vertices_with_c0;
    ts_uint nshell;
    ts_uint ncxmax;
    ts_uint ncymax;
    ts_uint nczmax;
    ts_idx npoly;
    ts_idx nmono;
    ts_idx internal_poly;
    ts_idx nfil;
    ts_idx nfono;
    ts_uint shc; // related to max l of the spherical harmonics
    ts_bool pressure_switch;
    ts_bool constvolswitch;
    ts_bool constareaswitch;
    ts_bool stretchswitch;
    ts_bool quiet;
    ts_bool plane_confinement_switch;
    ts_flag adhesion_geometry;
    ts_bool allow_xy_plane_movement;
    ts_bool force_balance_along_z_axis;
    ts_flag adhesion_model;
    ts_flag type_of_bond_model;
    ts_flag type_of_curvature_model;
    ts_flag type_of_force_model;
    //char *multiprocessing;
} ts_tape;




typedef struct {
    //ts_double bending_rigidity;
    ts_double dmax;
    ts_double stepsize;
    ts_double cm[3];
    ts_double volume;
    ts_double spring_constant;
    ts_double pressure;
    ts_double R_nucleus;
    ts_double R_nucleusX;
    ts_double R_nucleusY;
    ts_double R_nucleusZ;
    ts_double nucleus_center[3];
    ts_double area;
    ts_double adhesion_center;
    ts_tape *tape;
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
    ts_uint nshell;
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

ts_bool quiet;
ts_double V0;
ts_double A0;
ts_double epsvol;
ts_double epsarea;
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

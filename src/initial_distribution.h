/* vim: set ts=4 sts=4 sw=4 noet : */
#include "io.h"


/** @brief initial bond length */
#define DEF_A0 1.2 

/** @brief Creates initial distribution of vertices
  *
  * @param *vlist is the pointer to vertex list
  * @param N is the number of vertices to be initialized
  * @returns TS_SUCCESS on success, TS_FAIL otherwise. If allocation fails, the execution is terminated, reporting error code to the underlying operating system.
*/
ts_vesicle *initial_distribution_dipyramid(ts_uint nshell, ts_uint ncmax1, ts_uint ncmax2, ts_uint ncmax3, ts_double stepsize);

ts_vesicle *create_vesicle_from_tape(ts_tape *tape);
ts_bool set_vesicle_values_from_tape(ts_vesicle *vesicle);

ts_bool initial_population_with_c0(ts_vesicle *vesicle, ts_tape *tape);
/** Sets the initial position of the vertexes to dipyramid
 *
 *      @param *vlist is a pointer to list of vertices
 *      @returns TS_SUCCESS on success, TS_FAIL otherwise
 */
ts_bool pentagonal_dipyramid_vertex_distribution(ts_vertex_list *vlist);

/** Finds the neighbouring vertices and add them to a list of each vertex
 *
 *      @param *vlist is a pointer to a ts_vertex_list
 *      @returns TS_SUCCESS if successful, TS_FAIL otherwise
 */
ts_bool init_vertex_neighbours(ts_vertex_list *vlist);
//ts_bool init_vertex_neighbours(ts_vertex_list *vlist);

/** interior sites and their neighbours in circ. order + the triangles they are holding together */
ts_vertex_list *init_sort_neighbours(ts_bond_list *blist,ts_vertex_list *vlist);
ts_bool init_vesicle_bonds(ts_vesicle *vesicle);
ts_bool init_triangles(ts_vesicle *vesicle);
ts_bool init_triangle_neighbours(ts_vesicle *vesicle);
ts_bool init_common_vertex_triangle_neighbours(ts_vesicle *vesicle);
ts_bool init_normal_vectors(ts_triangle_list *tlist);


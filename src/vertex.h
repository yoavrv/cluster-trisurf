/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _VERTEX_H
#define _VERTEX_H

ts_bool vertex_list_assign_id(ts_vertex_list *vlist, ts_uint id);

/** @brief Creates initial vertex list
 *  
 *  Allocates memory and initializes the vertices.
 *	@param vertex is a structure holding information about 
 *      vertices
 *	@param N is a number of vertices that are used in simulation
 *	@param zero_them is boolean value. 0 skip setting zeros to idx 
 *      and (x,y,z) coordinates for each points, 1 means to zero all 
 *      information on points > 1 requests zeroing of coordinates and 
 *      indexing the vertices 0..N-1.
 *	@returns ts_bool value 1 on success, 0 otherwise
*/
ts_vertex_list *init_vertex_list(ts_uint N);
ts_bool vtx_add_neighbour(ts_vertex *vtx, ts_vertex *nvtx);
ts_bool vtx_add_cneighbour(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2);
ts_bool vtx_add_bond(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2);
ts_bool vtx_remove_neighbour(ts_vertex *vtx, ts_vertex *nvtx);
ts_bool vtx_free(ts_vertex *vtx);
ts_bool vtx_list_free(ts_vertex_list *vlist);
inline ts_double vtx_distance_sq(ts_vertex *vtx1, ts_vertex *vtx2);
ts_bool vtx_set_global_values(ts_vesicle *vesicle);
inline ts_double vtx_direct(ts_vertex *vtx1, ts_vertex *vtx2, ts_vertex *vtx3);
inline ts_bool vertex_add_tristar(ts_vertex *vtx, ts_triangle *tristarmem);
inline ts_bool vtx_insert_neighbour(ts_vertex *vtx, ts_vertex *nvtx, ts_vertex *vtxm);
inline ts_bool vtx_remove_tristar(ts_vertex *vtx, ts_triangle *tristar);
ts_bool vtx_copy(ts_vertex *cvtx,ts_vertex *ovtx);
ts_bool vtx_duplicate(ts_vertex *cvtx, ts_vertex *ovtx);
ts_vertex **vtx_neigh_copy(ts_vertex_list *vlist,ts_vertex *ovtx);
ts_vertex_list *vertex_list_copy(ts_vertex_list *ovlist);

#endif

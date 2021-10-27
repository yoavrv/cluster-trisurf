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
ts_vertex_list *init_vertex_list(ts_idx N);
ts_seen_vertex *init_seen_vertex(ts_uint max_size);
ts_bool vtx_add_neighbour(ts_vertex *vtx, ts_vertex *nvtx);
ts_bool vtx_add_cneighbour(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2);
ts_bool vtx_add_cneighbour2(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2);
ts_bool vtx_add_bond(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2);
ts_bool vtx_remove_neighbour(ts_vertex *vtx, ts_vertex *nvtx);

ts_bool vtx_free(ts_vertex *vtx);
ts_bool vtx_list_free(ts_vertex_list *vlist);
ts_bool seen_vertex_free(ts_seen_vertex *seen_vertex);

inline ts_double vtx_distance_sq(ts_vertex *vtx1, ts_vertex *vtx2);
ts_bool vtx_set_global_values(ts_vesicle *vesicle);
inline ts_double vtx_direct(ts_vertex *vtx1, ts_vertex *vtx2, ts_vertex *vtx3);

inline ts_bool vertex_add_tristar(ts_vertex *vtx, ts_triangle *tristarmem);
inline ts_bool vtx_insert_neighbour(ts_vertex *vtx, ts_vertex *nvtx, ts_vertex *vtxm);
inline ts_bool vtx_remove_tristar(ts_vertex *vtx, ts_triangle *tristar);
ts_bool add_vtx_to_seen(ts_seen_vertex *seen_vtx, ts_vertex *vtx);
ts_bool vtx_copy(ts_vertex *cvtx,ts_vertex *ovtx);
ts_bool vtx_duplicate(ts_vertex *cvtx, ts_vertex *ovtx);
ts_vertex **vtx_neigh_copy(ts_vertex_list *vlist,ts_vertex *ovtx);
ts_vertex_list *vertex_list_copy(ts_vertex_list *ovlist);

ts_bool is_in_seen_vertex(ts_seen_vertex *seen_vertex, ts_vertex *vtx);
ts_bool advance_seen_vertex_to_next_layer(ts_seen_vertex *seen_vertex);

ts_bool swap_triangles(ts_vertex* vtx, ts_small_idx i1, ts_small_idx i2);
ts_bool order_vertex_triangles(ts_vertex* vtx);
ts_bool print_tri_order(ts_vertex* vtx);
ts_bool print_vertex_ordered(ts_vertex* vtx);
ts_bool assert_vtx_ordered(ts_vertex* vtx);

ts_bool vtx_insert_tristar_at(ts_vertex *vtx, ts_triangle *tristarmem, ts_small_idx i);
ts_bool vtx_remove_tristar_at(ts_vertex *vtx, ts_small_idx i);
ts_bool vtx_insert_neigh_at(ts_vertex *vtx, ts_vertex *vtxmem, ts_small_idx i);
ts_bool vtx_remove_neigh_at(ts_vertex *vtx, ts_small_idx i);
ts_bool vtx_insert_bond_at(ts_vertex *vtx, ts_bond *bondmem, ts_small_idx i);
ts_bool vtx_remove_bond_at(ts_vertex *vtx, ts_small_idx i);
ts_bool vtx_insert_at(ts_vertex *vtx, ts_vertex *vtx_add, ts_bond* bond_add, ts_triangle* tri_add, ts_small_idx i);
ts_bool vtx_remove_at(ts_vertex *vtx, ts_small_idx i);
#endif

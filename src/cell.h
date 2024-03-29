/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _H_CELL
#define _H_CELL
ts_cell_list  *init_cell_list(ts_cell_idx ncmax1, ts_cell_idx ncmax2, ts_cell_idx ncmax3, ts_double stepsize);
ts_bool cell_free(ts_cell* cell);
ts_bool cell_list_free(ts_cell_list *clist);
inline ts_cell_idx vertex_self_avoidance(ts_vesicle *vesicle, ts_vertex *vtx);
inline ts_bool cell_add_vertex(ts_cell *cell, ts_vertex *vtx);
inline ts_bool cell_remove_vertex(ts_cell *cell, ts_vertex *vtx);
ts_bool cell_list_cell_occupation_clear(ts_cell_list *clist);
inline ts_bool cell_occupation_number_and_internal_proximity(ts_cell_list *clist, ts_cell_idx cellidx, ts_vertex *vtx);
#endif

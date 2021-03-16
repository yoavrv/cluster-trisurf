/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _POLY_H
#define _POLY_H
#include"io.h"

#include"general.h"

ts_poly	*init_poly(ts_uint n, ts_vertex *grafted_vtx);

ts_poly_list *init_poly_list(ts_uint n_poly, ts_uint n_mono, ts_vertex_list *vlist, ts_vesicle *vesicle);

ts_bool poly_free(ts_poly *poly);

ts_bool poly_list_free(ts_poly_list *poly_list);

ts_bool poly_assign_spring_const(ts_vesicle *vesicle);

ts_bool poly_assign_filament_xi(ts_vesicle *vesicle, ts_tape *tape);

ts_poly *remove_poly_with_index(ts_poly_list *poly_list, ts_uint idx);
ts_bool remove_random_polymeres(ts_poly_list *poly_list, ts_uint number);
ts_poly_list *init_empty_poly_list(ts_uint n_poly, ts_uint n_mono);
#endif

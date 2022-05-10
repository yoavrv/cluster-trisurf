/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _TRIANGLE_H
#define _TRIANGLE_H

ts_triangle_list *init_triangle_list(void);
ts_triangle *triangle_add(ts_triangle_list *tlist, ts_vertex *vtx1, ts_vertex *vtx2, ts_vertex *vtx3);
ts_bool triangle_add_neighbour(ts_triangle *tria, ts_triangle *ntria);
ts_bool triangle_normal_vector(ts_triangle *tria);
ts_bool triangle_list_free(ts_triangle_list *tlist);
ts_bool triangle_remove_neighbour(ts_triangle *tria, ts_triangle *ntria);
ts_bool in_tri(ts_triangle* t, ts_vertex* v);
ts_double triangle_dot_normals(ts_triangle *t1, ts_triangle *t2);

#endif

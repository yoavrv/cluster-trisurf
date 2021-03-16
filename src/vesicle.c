/* vim: set ts=4 sts=4 sw=4 noet : */
#include<general.h>
#include "vesicle.h"
#include "vertex.h"
#include "triangle.h"
#include "bond.h"
#include "cell.h"
#include "stdlib.h"
#include "poly.h"
#include "sh.h"
#include "shcomplex.h"

ts_vesicle *init_vesicle(ts_uint N, ts_uint ncmax1, ts_uint ncmax2, ts_uint
ncmax3, ts_double stepsize){
    ts_vesicle *vesicle=(ts_vesicle *)calloc(1,sizeof(ts_vesicle));
    vesicle->vlist=init_vertex_list(N);
    vesicle->blist=init_bond_list();
    vesicle->tlist=init_triangle_list();
    vesicle->clist=init_cell_list(ncmax1, ncmax2, ncmax3, stepsize);
    return vesicle;
}

ts_bool vesicle_translate(ts_vesicle *vesicle,ts_double x, ts_double y, ts_double z){
	ts_uint i;
	ts_vertex **vtx=vesicle->vlist->vtx;
	ts_uint nn=vesicle->vlist->n;
	for(i=0;i<nn;i++){
		vtx[i]->x+=x;
		vtx[i]->y+=y;
		vtx[i]->z+=z;
	}
	return TS_SUCCESS;
}

ts_bool vesicle_free(ts_vesicle *vesicle){
    vtx_list_free(vesicle->vlist);
    bond_list_free(vesicle->blist);
    triangle_list_free(vesicle->tlist);
    cell_list_free(vesicle->clist);
    poly_list_free(vesicle->poly_list);
    poly_list_free(vesicle->filament_list);
    complex_sph_free(vesicle->sphHarmonics);
    free(vesicle);
    return TS_SUCCESS;
}

/* @brief Function makes a sum of partial volumes of each triangle. Volumes of
 *
 * Partial volumes are calculated when we calculate normals of triangles. It is
 * relatively easy to calculate the volume of vesicle if we take into account
 * that the volume of the whole vertex is simply sum of all partial volumes of
 * all the triangles.
 */
ts_bool vesicle_volume(ts_vesicle *vesicle){
    ts_double volume;
    ts_uint i;
    ts_triangle **tria=vesicle->tlist->tria;
    volume=0;
    for(i=0; i<vesicle->tlist->n;i++){
    volume=volume+tria[i]->volume;
    }
    vesicle->volume=volume;
    return TS_SUCCESS;
}

/* @brief Function makes a sum of partial areas of each triangle.
 *
 *
 *
 */
ts_bool vesicle_area(ts_vesicle *vesicle){
    ts_double area;
    ts_uint i;
    ts_triangle **tria=vesicle->tlist->tria;
    area=0;
    for(i=0;i<vesicle->tlist->n;i++){
        area=area+tria[i]->area;
    }
    vesicle->area=area;
    return TS_SUCCESS;
}

ts_double vesicle_meancurvature(ts_vesicle *vesicle){
// Integrates (H dA) over vesicle area A, where H=(C1+C2)/2.
// (To be devided by A outside of function)
	ts_double mc;
	ts_uint i;
	mc=0;
	for(i=0;i<vesicle->vlist->n;i++){
		mc=mc+vesicle->vlist->vtx[i]->curvature;
	}
	return mc/2.0;
}

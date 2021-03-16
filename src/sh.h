/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _H_SH
#define _H_SH
#include "general.h"
ts_bool storeUlm2(ts_vesicle *vesicle);
ts_double plgndr(ts_int l, ts_int m, ts_double x);
ts_double shY(ts_int l,ts_int m,ts_double theta,ts_double fi);

ts_bool *cart2sph(ts_coord *coord, ts_double x, ts_double y, ts_double z);
ts_bool sph2cart(ts_coord *coord);
ts_bool precomputeShCoeff(ts_spharm *sph);
ts_spharm *sph_init(ts_vertex_list *vlist, ts_uint l);
ts_bool sph_free(ts_spharm *sph);
ts_double getR0(ts_vesicle *vesicle);
ts_bool preparationSh(ts_vesicle *vesicle, ts_double r0);
ts_bool calculateYlmi(ts_vesicle *vesicle);
ts_bool calculateUlm(ts_vesicle *vesicle);
ts_bool saveAvgUlm2(ts_vesicle *vesicle);
#endif

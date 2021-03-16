/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _H_SH_COMPLEX
#define _H_SH_COMPLEX
#include "general.h"
ts_bool storeUlmComplex2(ts_vesicle *vesicle);
ts_spharm *complex_sph_init(ts_vertex_list *vlist, ts_uint l);
ts_bool complex_sph_free(ts_spharm *sph);
ts_bool calculateUlmComplex(ts_vesicle *vesicle);
char *Ulm2Complex2String(ts_vesicle *vesicle);
ts_bool freeUlm2String(char *strng);
ts_double calculateKc(ts_vesicle *vesicle, ts_int lmin, ts_int lmax);
#endif

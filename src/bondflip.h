/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _H_BONDFLIP
#define _H_BONDFLIP

ts_bool single_bondflip_timestep(ts_vesicle *vesicle, ts_bond *bond, ts_double *rn);

ts_bool ts_flip_bond(ts_vesicle *vesicle,ts_vertex *k,ts_vertex *it,ts_vertex *km, ts_vertex *kp, ts_bond *bond, ts_triangle *lm, ts_triangle *lp, ts_triangle *lm2, ts_triangle *lp1);

ts_bool single_bondflip_timestep_ordered(ts_vesicle *vesicle, ts_bond *bond, ts_double *rn);

ts_bool ts_flip_bond_ordered(ts_vesicle *vesicle, ts_bond *bond,
                             ts_vertex *it, ts_small_idx neim,
                             ts_vertex *km, ts_small_idx nei_k_at_km,
                             ts_vertex *k, ts_small_idx nei_kp_at_k,
                             ts_vertex *kp, ts_small_idx nei_it_at_kp,
                             ts_triangle *lm, ts_triangle *lp,
                             ts_triangle *lm2, ts_triangle *lp1);
#endif

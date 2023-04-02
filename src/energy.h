/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _ENERGY_H
#define _ENERGY_H

ts_bool mean_curvature_and_energy(ts_vesicle *vesicle);
inline ts_bool energy_vertex(ts_vesicle *vesicle, ts_vertex *vtx);
inline ts_bool bond_energy(ts_bond *bond,ts_poly *poly);

ts_bool sweep_attraction_bond_energy(ts_vesicle *vesicle);
inline ts_bool attraction_bond_energy(ts_vesicle *vesicle, ts_bond *bond);
ts_double direct_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
ts_double direct_force_from_Fz_balance(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
inline ts_double total_force_on_vesicle(ts_vesicle *vesicle);
ts_double adhesion_energy_diff(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
ts_double adhesion_geometry_distance(ts_vesicle *vesicle,ts_vertex *vtx);
ts_bool adhesion_geometry_side(ts_vesicle *vesicle,ts_vertex *vtx);
void stretchenergy(ts_vesicle *vesicle, ts_triangle *triangle);

#endif

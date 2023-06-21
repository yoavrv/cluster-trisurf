/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _ENERGY_H
#define _ENERGY_H

ts_bool sweep_vertex_curvature_energy(ts_vesicle *vesicle);
ts_bool sweep_vertex_forces(ts_vesicle *vesicle);
ts_bool laplace_beltrami_curvature_energy(ts_vesicle *vesicle, ts_vertex *vtx);
ts_bool update_vertex_from_curvature_tensor(ts_vertex* vtx, ts_double Av,
                             ts_double s00,ts_double s01,ts_double s10,ts_double s11,
                             ts_double nx,ts_double ny,ts_double nz,
                             ts_double dx,ts_double dy,ts_double dz,
                             ts_double tx,ts_double ty,ts_double tz);
ts_bool error_correction_scheme(ts_double *s00, ts_double *s01, ts_double *s11, ts_double H, ts_double Kg);
ts_bool tensor_curvature_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_double mean_curvature,ts_double gaussian_curvature);
ts_bool tensor_curvature_energy2(ts_vesicle *vesicle, ts_vertex *vtx);
inline ts_bool vertex_curvature_energy(ts_vesicle *vesicle, ts_vertex *vtx);
inline ts_bool bond_energy(ts_bond *bond,ts_poly *poly);

ts_bool sweep_attraction_bond_energy(ts_vesicle *vesicle);
inline ts_bool attraction_bond_energy(ts_vesicle *vesicle, ts_bond *bond);
ts_double direct_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
ts_double direct_force_from_Fz_balance(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
ts_bool total_force_on_vesicle(ts_vesicle *vesicle);
ts_double adhesion_energy_diff(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
ts_double adhesion_geometry_distance(ts_vesicle *vesicle,ts_vertex *vtx);
ts_bool adhesion_geometry_side(ts_vesicle *vesicle,ts_vertex *vtx);
ts_double adhesion_geometry_factor(ts_vesicle* vesicle, ts_vertex* vtx);

ts_bool volume_pressure_area_energy_constraints(ts_vesicle *vesicle, ts_double *delta_energy, ts_double dvol, ts_double darea);
void stretchenergy(ts_vesicle *vesicle, ts_triangle *triangle);

#endif

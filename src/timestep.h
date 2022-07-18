/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _TIMESTEP_H
#define _TIMESTEP_H
ts_bool single_timestep(ts_vesicle *vesicle, ts_double *vmsr, ts_double *bfsr);
ts_bool run_simulation(ts_vesicle *vesicle, ts_massive_idx mcsweeps, ts_idx inititer, ts_idx iterations, ts_idx start_simulation);
#endif

/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _TIMESTEP_H
#define _TIMESTEP_H
ts_bool single_timestep(ts_vesicle *vesicle, ts_double *vmsr, ts_double *bfsr);
ts_bool run_simulation(ts_vesicle *vesicle, ts_uint mcsweeps, ts_uint inititer, ts_uint iterations, ts_uint start_simulation);
#endif

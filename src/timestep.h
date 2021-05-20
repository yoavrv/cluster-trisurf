/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _TIMESTEP_H
#define _TIMESTEP_H
ts_bool single_timestep(ts_vesicle *vesicle, ts_double *vmsr, ts_double *bfsr, clock_t *time_0, clock_t *time_1, clock_t *time_2, clock_t *time_3);
ts_bool run_simulation(ts_vesicle *vesicle, ts_uint mcsweeps, ts_uint inititer, ts_uint iterations, ts_uint start_simulation);
#endif

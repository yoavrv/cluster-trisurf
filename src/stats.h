/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _H_STATS
#define _H_STATS
void gyration_eigen(ts_vesicle *vesicle, ts_double *l1, ts_double *l2, ts_double *l3);
ts_ulong get_epoch();
ts_bool get_area_volume(ts_vesicle *vesicle, ts_double *area, ts_double *volume);
#endif

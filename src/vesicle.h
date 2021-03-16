/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _VESICLE_H
#define _VESICLE_H

ts_vesicle *init_vesicle(ts_uint N, ts_uint ncmax1, ts_uint ncmax2, ts_uint ncmax3, ts_double stepsize);
ts_bool vesicle_translate(ts_vesicle *vesicle,ts_double x, ts_double y, ts_double z);
ts_bool vesicle_free(ts_vesicle *vesicle);
ts_bool vesicle_volume(ts_vesicle *vesicle);
ts_bool vesicle_area(ts_vesicle *vesicle);
ts_double vesicle_meancurvature(ts_vesicle *vesicle);
#endif

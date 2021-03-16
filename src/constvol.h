/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _H_CONSTVOL
#define _H_CONSTVOL

ts_bool constvolume(ts_vesicle *vesicle, ts_vertex *vtx_avoid, ts_double Vol, ts_double *retEnergy, ts_vertex **vtx_moved_retval, ts_vertex **vtx_backup);
ts_bool constvolConstraintCheck(ts_vesicle *vesicle, ts_vertex *vtx);
ts_bool constvolumerestore(ts_vertex *vtx_moved, ts_vertex *vtx_backup);
ts_bool constvolumeaccept(ts_vesicle *vesicle, ts_vertex *vtx_moved, ts_vertex *vtx_backup);
#endif

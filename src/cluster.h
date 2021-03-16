#ifndef _H_CLUSTER
#define _H_CLUSTER
ts_cluster_list *init_cluster_list();
ts_cluster *new_cluster(ts_cluster_list *cstlist);
ts_bool cluster_add_vertex(ts_cluster *cluster, ts_vertex *vtx);
ts_bool cluster_free(ts_cluster *cluster);
ts_bool cluster_list_free(ts_cluster_list *cstlist);
ts_bool clusterize_vesicle(ts_vesicle *vesicle, ts_cluster_list *cstlist);
#endif

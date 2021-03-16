/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include "general.h"
#include "cluster.h"
#include <math.h>




ts_cluster_list *init_cluster_list(){
	ts_cluster_list *cstlist=(ts_cluster_list *)malloc(sizeof(ts_cluster_list));
	cstlist->n=0;
	cstlist->cluster=NULL;
	return cstlist;
}

ts_cluster *new_cluster(ts_cluster_list *cstlist){
	
	cstlist->n++;
	cstlist->cluster=(ts_cluster **)realloc(cstlist->cluster,cstlist->n*sizeof(ts_cluster *));
	if(cstlist->cluster==NULL) fatal("Cannot reallocate memory for additional **ts_cluster.",100);
	cstlist->cluster[cstlist->n-1]=(ts_cluster *)calloc(1,sizeof(ts_cluster));
	if(cstlist->cluster[cstlist->n-1]==NULL) fatal("Cannot allocate memory for additional *ts_cluster.",100);
	return cstlist->cluster[cstlist->n-1];
}

ts_bool cluster_add_vertex(ts_cluster *cluster, ts_vertex *vtx){
	cluster->nvtx++;
	cluster->vtx=(ts_vertex **)realloc(cluster->vtx, cluster->nvtx*sizeof(ts_vertex *));
	cluster->vtx[cluster->nvtx-1]=vtx;
	vtx->cluster=cluster;
	return TS_SUCCESS;
}

ts_bool cluster_list_compact(ts_cluster_list *cstlist){


	ts_uint i,n=cstlist->n;
	
	for(i=0;i<cstlist->n;i++){
		if(cstlist->cluster[i]==NULL){
			if(i>=n) break;
			//printf("Have to do compacting n=%d\n ",n);
			do{
				n--;
			} while(cstlist->cluster[n]==NULL && n>i);
			cstlist->cluster[i]=cstlist->cluster[n];
			cstlist->cluster[n]=NULL;
			//printf("After compacting n=%d\n",n);
		}
	}
	cstlist->cluster=(ts_cluster **)realloc(cstlist->cluster,n*sizeof(ts_cluster *));
	cstlist->n=n;
	return TS_SUCCESS;
}


ts_bool cluster_join(ts_cluster_list *cstlist, ts_cluster *cluster1, ts_cluster *cluster2){
	ts_cluster *master_cluster,*slave_cluster;
	ts_uint i;
	if(cluster1->idx<cluster2->idx){
		master_cluster=cluster1;
		slave_cluster=cluster2;
	} else {
		master_cluster=cluster2;
		slave_cluster=cluster1;
	}
	for(i=0;i<slave_cluster->nvtx;i++){
		cluster_add_vertex(master_cluster,slave_cluster->vtx[i]);
	}
	//find cluster in the list and free the location of the cluster
	for(i=0;i<cstlist->n;i++){
		if(cstlist->cluster[i]==slave_cluster){
			cluster_free(slave_cluster);
			cstlist->cluster[i]=NULL;
		}
	}
	return TS_SUCCESS;
}

ts_bool cluster_free(ts_cluster *cluster){
	if(cluster!=NULL){
		if(cluster->vtx!=NULL)
			free(cluster->vtx);
		free(cluster);
	}
	return TS_SUCCESS;
}

ts_bool cluster_list_free(ts_cluster_list *cstlist){
	ts_uint i;
	if(cstlist!=NULL){
		for(i=0;i<cstlist->n;i++){
			cluster_free(cstlist->cluster[i]);
		}
		free(cstlist->cluster);
		free(cstlist);
	}
	return TS_SUCCESS;
}


/* This is a stub function. User should check whether the vertex is clustering or not. */
ts_bool is_clusterable(ts_vertex *vtx){
	if(fabs(vtx->c)>1e-16){
		return 1;
	}
	else {
		return 0;
	}
}


ts_cluster *cluster_vertex_neighbor(ts_vertex *vtx){
	int j;
	for(j=0;j<vtx->neigh_no;j++){
		if(vtx->neigh[j]->cluster!=NULL)
			return vtx->neigh[j]->cluster;
	}
	return NULL;
}

ts_bool cluster_vertex_neighbor_check(ts_cluster_list *cstlist, ts_vertex *vtx){

	int j;
	for(j=0;j<vtx->neigh_no;j++){
		if(vtx->neigh[j]->cluster!=NULL){
			if(vtx->neigh[j]->cluster!=vtx->cluster){
				cluster_join(cstlist, vtx->cluster, vtx->neigh[j]->cluster);
			}
		}
	}
	return TS_SUCCESS;
}


ts_bool clusterize_vesicle(ts_vesicle *vesicle, ts_cluster_list *cstlist){

	int i;
	ts_vertex *vtx;
	ts_cluster *cst;
	for(i=0;i<vesicle->vlist->n;i++){
	//for each vertex
		vtx=vesicle->vlist->vtx[i];
		if(is_clusterable(vtx)){
			if(vtx->cluster==NULL){
				//find first neigbor with cluster index
				cst=cluster_vertex_neighbor(vtx);
				if(cst==NULL){
					//no clusters are around us, vo we are probably lonely vertex or no surronding vertex has been mapped yet.
					cst=new_cluster(cstlist);
					cluster_add_vertex(cst,vtx);
				} else {
					//we are added to the first cluster found
					cluster_add_vertex(cst,vtx);
					cluster_vertex_neighbor_check(cstlist, vtx);
					cluster_list_compact(cstlist);
				}

			}
		}
	}


	return TS_SUCCESS;
}

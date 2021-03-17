/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "general.h"
//#include "vertex.h"
//#include "bond.h"
//#include "triangle.h"
//#include "cell.h"
#include "vesicle.h"
#include "io.h"
//#include "initial_distribution.h"
//#include "frame.h"
//#include "timestep.h"
//#include "poly.h"
#include "stats.h"
#include "sh.h"
#include "shcomplex.h"
#include "dumpstate.h"
#include "restore.h"
#include "cluster.h"
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <snapshot.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<stdio.h>

ts_vesicle *restoreVesicle(char *filename){
	ts_vesicle *vesicle = parseDump(filename);
	return vesicle;
}

void vesicle_calculate_ulm2(ts_vesicle *vesicle){
	//complex_sph_free(vesicle->sphHarmonics);

	//vesicle->sphHarmonics=complex_sph_init(vesicle->vlist,21);
	vesicle_volume(vesicle);
	preparationSh(vesicle,getR0(vesicle));
	calculateUlmComplex(vesicle);
	ts_int i,j;
	for(i=0;i<vesicle->sphHarmonics->l;i++){
    		for(j=i;j<2*i+1;j++){
			printf("%e ", gsl_complex_abs2(vesicle->sphHarmonics->ulmComplex[i][j]));
    		}
	}
		printf("\n");

}



int count_bonds_with_energy(ts_bond_list *blist){

	unsigned int i, cnt;
	cnt=0;
	for(i=0;i<blist->n;i++){
		if(fabs(blist->bond[i]->energy)>1e-16) cnt++;
	}
	return cnt;
}



ts_bool write_histogram_data(ts_uint timestep_no, ts_vesicle *vesicle){
	ts_cluster_list *cstlist=init_cluster_list();
	clusterize_vesicle(vesicle,cstlist);
	//printf("No clusters=%d\n",cstlist->n);
	int k,i,cnt, test=0;
	int max_nvtx=0;
	char filename[255];
	sprintf(filename,"histogram_%.6u.csv",timestep_no);
	FILE *fd=fopen(filename,"w");
	fprintf(fd,"Number_of_vertices_in_cluster,Number_of_clusters,\n");
	for(k=0;k<cstlist->n;k++)
		if(cstlist->cluster[k]->nvtx>max_nvtx) max_nvtx=cstlist->cluster[k]->nvtx;
	//printf("Max. number of vertices in cluster: %d\n",max_nvtx);
	for(i=1;i<=max_nvtx;i++){
		cnt=0;
		for(k=0;k<cstlist->n;k++)
			if(cstlist->cluster[k]->nvtx==i) cnt++;
		fprintf(fd,"%d,%d,\n",i,cnt);
		test+=cnt*i;
	}
	//for(k=0;k<cstlist->n;k++){
//		printf("*Cluster %d has %d vertices\n",k,cstlist->cluster[k]->nvtx);
//	}

	fclose(fd);
//	printf("*Sum of all vertices in clusters: %d\n", test);
//	write_vertex_xml_file(vesicle,timestep_no,cstlist);
	cluster_list_free(cstlist);
	
	return TS_SUCCESS;
}


int main(){
	ts_vesicle *vesicle;
	ts_char *i,*j;
	ts_uint tstep,n;
    	ts_char *number;
	struct dirent **list;
	ts_double l1,l2,l3,hbar;
	int count;

	fprintf(stdout, "OuterLoop Volume Area lamdba1 lambda2 lambda3 Nbw/Nb hbar\n");


	count=scandir(".",&list,0,alphasort);
	if(count<0){
		fatal("Error, cannot open directory.",1);
	}
        tstep=0;
	for(n=0;n<count;n++){
		struct dirent *ent;
		ent=list[n];	
            	i=rindex(ent->d_name,'.');
            	if(i==NULL) {
				continue;
		}
            	if(strcmp(i+1,"vtu")==0){
                    j=rindex(ent->d_name,'_');
                    if(j==NULL) continue;
                    number=strndup(j+1,j-i); 
			quiet=1;
                    ts_fprintf(stdout,"timestep: %u filename: %s\n",atoi(number),ent->d_name);
//			printf("%u ",atoi(number));
			vesicle=restoreVesicle(ent->d_name);
//			vesicle_calculate_ulm2(vesicle);
			vesicle_volume(vesicle);
			vesicle_area(vesicle);
			gyration_eigen(vesicle,&l1,&l2,&l3);
			hbar=vesicle_meancurvature(vesicle)/vesicle->area;			
			fprintf(stdout,"%d %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",atoi(number),vesicle->volume, vesicle->area,l1,l2,l3, (ts_double)count_bonds_with_energy(vesicle->blist)/(ts_double)vesicle->blist->n,hbar);
                    	tstep++;
			write_histogram_data(atoi(number), vesicle);
                    free(number);
			tape_free(vesicle->tape);
			vesicle_free(vesicle);
            	}
		}
	for (n = 0; n < count; n++)
  	{
  		free(list[n]);
  	}
	
	free(list);
	return 0;
}


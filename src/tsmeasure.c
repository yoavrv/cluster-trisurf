/* vim: set ts=4 sts=4 sw=4 noet : */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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
#include "sh.h"
#include "shcomplex.h"
#include "dumpstate.h"
#include "restore.h"
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <snapshot.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


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

int main(){
    ts_vesicle *vesicle;
    ts_char *i,*j;
    ts_uint tstep,n;
    ts_char *number;
    struct dirent **list;
    int count;
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
            printf("%u ",atoi(number));
            vesicle=restoreVesicle(ent->d_name);
            vesicle_calculate_ulm2(vesicle);
            tstep++;
            free(number);
            tape_free(vesicle->tape);
            vesicle_free(vesicle);
        }
    }
    for (n = 0; n < count; n++){
        free(list[n]);
    }
    
    free(list);
    return 0;
}


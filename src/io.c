
/* vim: set ts=4 sts=4 sw=4 noet : */
#include "general.h"
#include <stdio.h>
#include "io.h"
#include "vertex.h"
#include "bond.h"
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include "initial_distribution.h"
#include "poly.h"
#include "cell.h"
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <snapshot.h>
//#include "getter_setters.h"

#define TS_WRITE_ITERATE_VTX(string,field)\
do{\
    for(i=0;i<vesicle->vlist->n;i++){\
        fprintf(fh,string, vtx[i]->field);\
    }\
    if(poly){\
        for(i=0;i<vesicle->poly_list->n;i++){\
            for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){\
                fprintf(fh,string, vesicle->poly_list->poly[i]->vlist->vtx[j]->field);\
            }\
        }\
    }\
    if(fil){\
        for(i=0;i<vesicle->filament_list->n;i++){\
            for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){\
                fprintf(fh,string, vesicle->filament_list->poly[i]->vlist->vtx[j]->field);\
            }\
        }\
    }\
}while(0)

#define TS_WRITE_VECTOR_ITERATE_VTX(string,fieldx,fieldy,fieldz)\
do{\
    for(i=0;i<vesicle->vlist->n;i++){\
        fprintf(fh,string, vtx[i]->fieldx,vtx[i]->fieldy,vtx[i]->fieldz);\
    }\
    if(poly){\
        for(i=0;i<vesicle->poly_list->n;i++){\
            for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){\
                fprintf(fh,string, vesicle->poly_list->poly[i]->vlist->vtx[j]->fieldx, \
                                   vesicle->poly_list->poly[i]->vlist->vtx[j]->fieldy, \
                                   vesicle->poly_list->poly[i]->vlist->vtx[j]->fieldz);\
            }\
        }\
    }\
    if(fil){\
        for(i=0;i<vesicle->filament_list->n;i++){\
            for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){\
                fprintf(fh,string, vesicle->filament_list->poly[i]->vlist->vtx[j]->fieldx, \
                                   vesicle->filament_list->poly[i]->vlist->vtx[j]->fieldy, \
                                   vesicle->filament_list->poly[i]->vlist->vtx[j]->fieldz);\
            }\
        }\
    }\
}while(0)

/** DUMP STATE TO DISK DRIVE **/
ts_bool dump_state(ts_vesicle *vesicle, ts_idx iteration){

    /* save current state with wrong pointers. Will fix that later */
    ts_idx i,j,k;
    FILE *fh=fopen(command_line_args.dump_fullfilename,"wb");

    /* dump vesicle */
    fwrite(vesicle, sizeof(ts_vesicle)-sizeof(ts_double),1,fh);
    /* dump vertex list */
    fwrite(vesicle->vlist, sizeof(ts_vertex_list),1,fh);
    /* dump bond list */
    fwrite(vesicle->blist, sizeof(ts_bond_list),1,fh);
    /* dump triangle list */
    fwrite(vesicle->tlist, sizeof(ts_triangle_list),1,fh);
    /* dump cell list */
    fwrite(vesicle->clist, sizeof(ts_cell_list),1,fh);
    /* dump poly list */
    fwrite(vesicle->poly_list, sizeof(ts_poly_list),1,fh);
    /* dump filament list */
    fwrite(vesicle->filament_list, sizeof(ts_poly_list),1,fh);
    /* level 1 complete */

    /*dump vertices*/
    for(i=0;i<vesicle->vlist->n;i++){
        fwrite(vesicle->vlist->vtx[i],sizeof(ts_vertex),1,fh);
        /* dump pointer offsets for:
                    neigh
                    bond
                    tria
                    cell is ignored
        */
        for(j=0;j<vesicle->vlist->vtx[i]->neigh_no;j++){
            fwrite(&vesicle->vlist->vtx[i]->neigh[j]->idx,sizeof(ts_idx),1,fh); 
        }
        for(j=0;j<vesicle->vlist->vtx[i]->bond_no;j++){
            fwrite(&vesicle->vlist->vtx[i]->bond[j]->idx,sizeof(ts_idx),1,fh); 
        }
        for(j=0;j<vesicle->vlist->vtx[i]->tristar_no;j++){
            fwrite(&vesicle->vlist->vtx[i]->tristar[j]->idx,sizeof(ts_idx),1,fh); 
        }
    }

    /*dump bonds*/
    for(i=0;i<vesicle->blist->n;i++){
        fwrite(vesicle->blist->bond[i],sizeof(ts_bond),1,fh);
        /* dump pointer offsets for vtx1 and vtx2 */
        //off=(ts_ulong)(vesicle->blist->bond[i]->vtx1-vesicle->vlist->vtx[0]);
        fwrite(&vesicle->blist->bond[i]->vtx1->idx,sizeof(ts_idx),1,fh); 
        //off=(ts_ulong)(vesicle->blist->bond[i]->vtx2-vesicle->vlist->vtx[0]);
        fwrite(&vesicle->blist->bond[i]->vtx2->idx,sizeof(ts_idx),1,fh); 
    }

    /*dump triangles*/
    for(i=0;i<vesicle->tlist->n;i++){
        fwrite(vesicle->tlist->tria[i],sizeof(ts_triangle),1,fh);
        /* dump pointer offsets for vertex */
        fwrite(&vesicle->tlist->tria[i]->vertex[0]->idx,sizeof(ts_idx),1,fh); 
        fwrite(&vesicle->tlist->tria[i]->vertex[1]->idx,sizeof(ts_idx),1,fh); 
        fwrite(&vesicle->tlist->tria[i]->vertex[2]->idx,sizeof(ts_idx),1,fh); 
        /* dump pointer offsets for neigh */
        for(j=0;j<vesicle->tlist->tria[i]->neigh_no;j++){
            fwrite(&vesicle->tlist->tria[i]->neigh[j]->idx,sizeof(ts_idx),1,fh); 
        }
    }


    /*dump polymeres */
    for(i=0;i<vesicle->poly_list->n;i++){
        fwrite(vesicle->poly_list->poly[i],sizeof(ts_poly),1,fh);
        fwrite(vesicle->poly_list->poly[i]->vlist,sizeof(ts_vertex_list),1,fh);
        fwrite(vesicle->poly_list->poly[i]->blist,sizeof(ts_bond_list),1,fh);
    } 
     
    /* dump poly vertex(monomer) list*/
    for(i=0;i<vesicle->poly_list->n;i++){
        for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
            fwrite(vesicle->poly_list->poly[i]->vlist->vtx[j],sizeof(ts_vertex),1,fh);
            /* dump offset for neigh and bond */
            for(k=0;k<vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh_no;k++){
               // off=(ts_ulong)(vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh[k]-vesicle->poly_list->poly[i]->vlist->vtx[0]);
                fwrite(&vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh[k]->idx,sizeof(ts_idx),1,fh); 
            }
            for(k=0;k<vesicle->poly_list->poly[i]->vlist->vtx[j]->bond_no;k++){
                //off=(ts_ulong)(vesicle->poly_list->poly[i]->vlist->vtx[j]->bond[k]-vesicle->poly_list->poly[i]->blist->bond[0]);
                fwrite(&vesicle->poly_list->poly[i]->vlist->vtx[j]->bond[k]->idx,sizeof(ts_idx),1,fh); 
            }
        }
    // grafted vtx on vesicle data dump
        fwrite(&vesicle->poly_list->poly[i]->grafted_vtx->idx, sizeof(ts_idx),1,fh);
    }
    /* dump poly bonds between monomers list*/
    for(i=0;i<vesicle->poly_list->n;i++){
        for(j=0;j<vesicle->poly_list->poly[i]->blist->n;j++){
            fwrite(vesicle->poly_list->poly[i]->blist->bond[j],sizeof(ts_bond),1,fh);
            /* dump vtx1 and vtx2 offsets */
            // off=(ts_ulong)(vesicle->poly_list->poly[i]->blist->bond[j]->vtx1-vesicle->poly_list->poly[i]->vlist->vtx[0]);
            fwrite(&vesicle->poly_list->poly[i]->blist->bond[j]->vtx1->idx,sizeof(ts_idx),1,fh); 
            // off=(ts_ulong)(vesicle->poly_list->poly[i]->blist->bond[j]->vtx2-vesicle->poly_list->poly[i]->vlist->vtx[0]);
            fwrite(&vesicle->poly_list->poly[i]->blist->bond[j]->vtx2->idx,sizeof(ts_idx),1,fh); 
        }
    }


    /*dump filamentes grandes svinjas */
    for(i=0;i<vesicle->filament_list->n;i++){
        fwrite(vesicle->filament_list->poly[i],sizeof(ts_poly),1,fh);
        fwrite(vesicle->filament_list->poly[i]->vlist,sizeof(ts_vertex_list),1,fh);
        fwrite(vesicle->filament_list->poly[i]->blist,sizeof(ts_bond_list),1,fh);
    }

    /* dump filamentes vertex(monomer) list*/
    for(i=0;i<vesicle->filament_list->n;i++){
        for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
            fwrite(vesicle->filament_list->poly[i]->vlist->vtx[j],sizeof(ts_vertex),1,fh);
            /* dump offset for neigh and bond */
            for(k=0;k<vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh_no;k++){
               // off=(ts_ulong)(vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh[k]-vesicle->poly_list->poly[i]->vlist->vtx[0]);
                fwrite(&vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh[k]->idx,sizeof(ts_idx),1,fh); 
            }
            for(k=0;k<vesicle->filament_list->poly[i]->vlist->vtx[j]->bond_no;k++){
                //off=(ts_ulong)(vesicle->poly_list->poly[i]->vlist->vtx[j]->bond[k]-vesicle->poly_list->poly[i]->blist->bond[0]);
                fwrite(&vesicle->filament_list->poly[i]->vlist->vtx[j]->bond[k]->idx,sizeof(ts_idx),1,fh); 
            }
        }
    }
    /* dump poly bonds between monomers list*/
    for(i=0;i<vesicle->filament_list->n;i++){
        for(j=0;j<vesicle->filament_list->poly[i]->blist->n;j++){
            fwrite(vesicle->filament_list->poly[i]->blist->bond[j],sizeof(ts_bond),1,fh);
            /* dump vtx1 and vtx2 offsets */
            // off=(ts_ulong)(vesicle->poly_list->poly[i]->blist->bond[j]->vtx1-vesicle->poly_list->poly[i]->vlist->vtx[0]);
            fwrite(&vesicle->filament_list->poly[i]->blist->bond[j]->vtx1->idx,sizeof(ts_idx),1,fh); 
            // off=(ts_ulong)(vesicle->poly_list->poly[i]->blist->bond[j]->vtx2-vesicle->poly_list->poly[i]->vlist->vtx[0]);
            fwrite(&vesicle->filament_list->poly[i]->blist->bond[j]->vtx2->idx,sizeof(ts_idx),1,fh); 
        }
    }



    /* pointer offsets for fixing the restored pointers */
    /* need pointers for 
        vlist->vtx
        blist->bond
        tlist->tria
        clist->cell
        poly_list->poly
        and for each poly:
            poly_list->poly->vtx
            poly_list->poly->bond
    */

    //	fwrite(vesicle->clist, sizeof(ts_cell_list),1,  fh);
    /* write tape information on vesicle */
    //    fwrite(vesicle->tape,sizeof(ts_tape),1,fh);
    fwrite(&iteration, sizeof(ts_idx),1,fh);
    fclose(fh);
    return TS_SUCCESS;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/** RESTORE DUMP FROM DISK **/
ts_vesicle *restore_state(ts_idx *iteration, ts_tape* tape){
    ts_idx i,j,k;
    FILE *fh=fopen(command_line_args.dump_fullfilename,"rb");

    struct stat sb;
    if (stat(command_line_args.dump_fullfilename, &sb) == -1) {
        //dump file does not exist.
        return NULL;
    }

    //check if it is regular file
    if((sb.st_mode & S_IFMT) != S_IFREG) {
        //dump file is not a regular file.
        ts_fprintf(stderr,"Dump file is not a regular file!\n");
        return NULL;
    }

    ts_uint retval; // read from function (size_t fread())
    ts_idx idx;

    /* we restore all the data from the dump */
    /* restore vesicle */
    ts_vesicle *vesicle=(ts_vesicle *)calloc(1,sizeof(ts_vesicle));
    retval=fread(vesicle, sizeof(ts_vesicle)-sizeof(ts_double),1,fh);
    //	fprintf(stderr,"was here! %e\n",vesicle->dmax);

    /* restore vertex list */
    vesicle->vlist=(ts_vertex_list *)malloc(sizeof(ts_vertex_list));
    retval=fread(vesicle->vlist, sizeof(ts_vertex_list),1,fh);
    /* restore bond list */
    vesicle->blist=(ts_bond_list *)malloc(sizeof(ts_bond_list));
    retval=fread(vesicle->blist, sizeof(ts_bond_list),1,fh);
    /* restore triangle list */
    vesicle->tlist=(ts_triangle_list *)malloc(sizeof(ts_triangle_list));
    retval=fread(vesicle->tlist, sizeof(ts_triangle_list),1,fh);
    /* restore cell list */
    vesicle->clist=(ts_cell_list *)malloc(sizeof(ts_cell_list));
    retval=fread(vesicle->clist, sizeof(ts_cell_list),1,fh);
    /* restore poly list */
    vesicle->poly_list=(ts_poly_list *)calloc(1,sizeof(ts_poly_list));
    retval=fread(vesicle->poly_list, sizeof(ts_poly_list),1,fh);
    /* restore filament list */
    vesicle->filament_list=(ts_poly_list *)calloc(1,sizeof(ts_poly_list));
    retval=fread(vesicle->filament_list, sizeof(ts_poly_list),1,fh);
    /* level 1 complete */

    /* prerequisity. Bonds must be malloced before vertexes are recreated */
  vesicle->blist->bond=(ts_bond **)calloc(vesicle->blist->n,sizeof(ts_bond *));
    for(i=0;i<vesicle->blist->n;i++){
        vesicle->blist->bond[i]=(ts_bond *)malloc(sizeof(ts_bond));
    }
    /* prerequisity. Triangles must be malloced before vertexes are recreated */
    vesicle->tlist->tria=(ts_triangle **)calloc(vesicle->tlist->n,sizeof(ts_triangle *));
    for(i=0;i<vesicle->tlist->n;i++){
        vesicle->tlist->tria[i]=(ts_triangle *)malloc(sizeof(ts_triangle));
    }
    /* prerequisity. Vertices must be malloced before vertexes are recreated */
    vesicle->vlist->vtx=(ts_vertex **)calloc(vesicle->vlist->n,sizeof(ts_vertex *));
    for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]=(ts_vertex *)malloc(sizeof(ts_vertex));
    }
    /*restore vertices*/
    for(i=0;i<vesicle->vlist->n;i++){
        retval=fread(vesicle->vlist->vtx[i],sizeof(ts_vertex),1,fh);
        /*restore neigh, bond, tristar. Ignoring cell */
        vesicle->vlist->vtx[i]->neigh=(ts_vertex **)calloc(vesicle->vlist->vtx[i]->neigh_no, sizeof(ts_vertex *));
        for(j=0;j<vesicle->vlist->vtx[i]->neigh_no;j++){
            retval=fread(&idx,sizeof(ts_idx),1,fh);
            vesicle->vlist->vtx[i]->neigh[j]=vesicle->vlist->vtx[idx];
        }
        vesicle->vlist->vtx[i]->bond=(ts_bond **)calloc(vesicle->vlist->vtx[i]->bond_no, sizeof(ts_bond *));
        for(j=0;j<vesicle->vlist->vtx[i]->bond_no;j++){
            retval=fread(&idx,sizeof(ts_idx),1,fh);
            /* pointer can be assigned only when list of bonds is fully initialized in memory. Thus bondlist popularization must be done before vertex can reference to it */
            vesicle->vlist->vtx[i]->bond[j]=vesicle->blist->bond[idx];    
        }

        vesicle->vlist->vtx[i]->tristar=(ts_triangle **)calloc(vesicle->vlist->vtx[i]->tristar_no, sizeof(ts_triangle *));
        for(j=0;j<vesicle->vlist->vtx[i]->tristar_no;j++){
            retval=fread(&idx,sizeof(ts_idx),1,fh);
            /* same comment as above */
            vesicle->vlist->vtx[i]->tristar[j]=vesicle->tlist->tria[idx];
        }

    }

    /*restore bonds*/
    // vesicle->blist->bond=(ts_bond **)calloc(vesicle->blist->n,sizeof(ts_bond *)); // done before.
    for(i=0;i<vesicle->blist->n;i++){
     //   vesicle->blist->bond[i]=(ts_bond *)malloc(sizeof(ts_bond)); //done before.
        retval=fread(vesicle->blist->bond[i],sizeof(ts_bond),1,fh);
        /* restore vtx1 and vtx2 */
        retval=fread(&idx,sizeof(ts_idx),1,fh);
        vesicle->blist->bond[i]->vtx1=vesicle->vlist->vtx[idx];
        retval=fread(&idx,sizeof(ts_idx),1,fh);
        vesicle->blist->bond[i]->vtx2=vesicle->vlist->vtx[idx];
    }

    /*restore triangles*/
//    vesicle->tlist->tria=(ts_triangle **)calloc(vesicle->tlist->n,sizeof(ts_triangle *)); // done before
    for(i=0;i<vesicle->tlist->n;i++){
 //       vesicle->tlist->tria[i]=(ts_triangle *)malloc(sizeof(ts_triangle)); // done before
        retval=fread(vesicle->tlist->tria[i],sizeof(ts_triangle),1,fh);
        /* restore pointers for vertices */
        retval=fread(&idx,sizeof(ts_idx),1,fh);
        vesicle->tlist->tria[i]->vertex[0]=vesicle->vlist->vtx[idx];
        retval=fread(&idx,sizeof(ts_idx),1,fh);
        vesicle->tlist->tria[i]->vertex[1]=vesicle->vlist->vtx[idx];
        retval=fread(&idx,sizeof(ts_idx),1,fh);
        vesicle->tlist->tria[i]->vertex[2]=vesicle->vlist->vtx[idx];
        /* restore pointers for neigh */
     vesicle->tlist->tria[i]->neigh=(ts_triangle **)malloc(vesicle->tlist->tria[i]->neigh_no*sizeof(ts_triangle *));
        for(j=0;j<vesicle->tlist->tria[i]->neigh_no;j++){
            retval=fread(&idx,sizeof(ts_idx),1,fh);
            vesicle->tlist->tria[i]->neigh[j]=vesicle->tlist->tria[idx];
        }

    }
   
    /*restore cells */
    /*TODO: do we need to recalculate cells here? */
    /* vesicle->clist->cell=(ts_cell **)malloc(vesicle->clist->cellno*sizeof(ts_cell *));
    for(i=0;i<vesicle->clist->cellno;i++){
        vesicle->clist->cell[i]=(ts_cell *)malloc(sizeof(ts_cell));
        retval=fread(vesicle->clist->cell[i],sizeof(ts_cell),1,fh);
    }
    */
    /*restore polymeres */
    vesicle->poly_list->poly = (ts_poly **)calloc(vesicle->poly_list->n,sizeof(ts_poly *));
    for(i=0;i<vesicle->poly_list->n;i++){
        vesicle->poly_list->poly[i]=(ts_poly *)calloc(1,sizeof(ts_poly));
        retval=fread(vesicle->poly_list->poly[i],sizeof(ts_poly),1,fh);
        vesicle->poly_list->poly[i]->vlist=(ts_vertex_list *)calloc(1,sizeof(ts_vertex_list));
        retval=fread(vesicle->poly_list->poly[i]->vlist,sizeof(ts_vertex_list),1,fh);
        vesicle->poly_list->poly[i]->blist=(ts_bond_list *)calloc(1,sizeof(ts_bond_list));
        retval=fread(vesicle->poly_list->poly[i]->blist,sizeof(ts_bond_list),1,fh);
    /* initialize adress space for pointers that will hold specific vertices (monomers) and bonds */
        vesicle->poly_list->poly[i]->vlist->vtx=(ts_vertex **)calloc(vesicle->poly_list->poly[i]->vlist->n,sizeof(ts_vertex *));
        vesicle->poly_list->poly[i]->blist->bond=(ts_bond **)calloc(vesicle->poly_list->poly[i]->blist->n,sizeof(ts_bond *));
     for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
            vesicle->poly_list->poly[i]->vlist->vtx[j]=(ts_vertex *)malloc(sizeof(ts_vertex));
    }
    for(j=0;j<vesicle->poly_list->poly[i]->blist->n;j++){
            vesicle->poly_list->poly[i]->blist->bond[j]=(ts_bond *)malloc(sizeof(ts_bond));
    }

    } 

     
    /* restore poly vertex(monomer) list*/
    for(i=0;i<vesicle->poly_list->n;i++){
        for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
            retval=fread(vesicle->poly_list->poly[i]->vlist->vtx[j],sizeof(ts_vertex),1,fh);
                
            /* restore neigh and bonds */
            vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh=(ts_vertex **)calloc(vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh_no, sizeof(ts_vertex *));
            for(k=0;k<vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh_no;k++){
                retval=fread(&idx,sizeof(ts_idx),1,fh);
                vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh[k]=vesicle->poly_list->poly[i]->vlist->vtx[idx];
            }
            vesicle->poly_list->poly[i]->vlist->vtx[j]->bond=(ts_bond **)calloc(vesicle->poly_list->poly[i]->vlist->vtx[j]->bond_no, sizeof(ts_bond *));
            for(k=0;k<vesicle->poly_list->poly[i]->vlist->vtx[j]->bond_no;k++){
                retval=fread(&idx,sizeof(ts_idx),1,fh);
                vesicle->poly_list->poly[i]->vlist->vtx[j]->bond[k]=vesicle->poly_list->poly[i]->blist->bond[idx];
            }

        }
        /* restore grafted vtx on vesicle and grafted_poly */
        retval=fread(&idx,sizeof(ts_idx),1,fh);
        vesicle->vlist->vtx[idx]->grafted_poly=vesicle->poly_list->poly[i];
        vesicle->poly_list->poly[i]->grafted_vtx=vesicle->vlist->vtx[idx];	
    }

    /* restore poly bonds between monomers list*/
    for(i=0;i<vesicle->poly_list->n;i++){
        for(j=0;j<vesicle->poly_list->poly[i]->blist->n;j++){
       //     vesicle->poly_list->poly[i]->blist->bond[j]=(ts_bond *)malloc(sizeof(ts_bond));
            retval=fread(vesicle->poly_list->poly[i]->blist->bond[j],sizeof(ts_bond),1,fh);
            /* restore vtx1 and vtx2 */
                retval=fread(&idx,sizeof(ts_idx),1,fh);
                vesicle->poly_list->poly[i]->blist->bond[j]->vtx1=vesicle->poly_list->poly[i]->vlist->vtx[idx];
                retval=fread(&idx,sizeof(ts_idx),1,fh);
                vesicle->poly_list->poly[i]->blist->bond[j]->vtx2=vesicle->poly_list->poly[i]->vlist->vtx[idx];
        }
    }

    /*restore filaments */
    vesicle->filament_list->poly = (ts_poly **)calloc(vesicle->filament_list->n,sizeof(ts_poly *));
    for(i=0;i<vesicle->filament_list->n;i++){
        vesicle->filament_list->poly[i]=(ts_poly *)calloc(1,sizeof(ts_poly));
        retval=fread(vesicle->filament_list->poly[i],sizeof(ts_poly),1,fh);
        vesicle->filament_list->poly[i]->vlist=(ts_vertex_list *)calloc(1,sizeof(ts_vertex_list));
        retval=fread(vesicle->filament_list->poly[i]->vlist,sizeof(ts_vertex_list),1,fh);
        vesicle->filament_list->poly[i]->blist=(ts_bond_list *)calloc(1,sizeof(ts_bond_list));
        retval=fread(vesicle->filament_list->poly[i]->blist,sizeof(ts_bond_list),1,fh);
        /* initialize adress space for pointers that will hold specific vertices (monomers) and bonds */
        vesicle->filament_list->poly[i]->vlist->vtx=(ts_vertex **)calloc(vesicle->filament_list->poly[i]->vlist->n,sizeof(ts_vertex *));
        vesicle->filament_list->poly[i]->blist->bond=(ts_bond **)calloc(vesicle->filament_list->poly[i]->blist->n,sizeof(ts_bond *));
        for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
            vesicle->filament_list->poly[i]->vlist->vtx[j]=(ts_vertex *)malloc(sizeof(ts_vertex));
        }
        for(j=0;j<vesicle->filament_list->poly[i]->blist->n;j++){
            vesicle->filament_list->poly[i]->blist->bond[j]=(ts_bond *)malloc(sizeof(ts_bond));
        }

    } 

     
    /* restore poly vertex(monomer) list*/
    for(i=0;i<vesicle->filament_list->n;i++){
        for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
            retval=fread(vesicle->filament_list->poly[i]->vlist->vtx[j],sizeof(ts_vertex),1,fh);
                
            /* restore neigh and bonds */
            vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh=(ts_vertex **)calloc(vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh_no, sizeof(ts_vertex *));
            for(k=0;k<vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh_no;k++){
                retval=fread(&idx,sizeof(ts_idx),1,fh);
                vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh[k]=vesicle->filament_list->poly[i]->vlist->vtx[idx];
            }
            vesicle->filament_list->poly[i]->vlist->vtx[j]->bond=(ts_bond **)calloc(vesicle->filament_list->poly[i]->vlist->vtx[j]->bond_no, sizeof(ts_bond *));
            for(k=0;k<vesicle->filament_list->poly[i]->vlist->vtx[j]->bond_no;k++){
                retval=fread(&idx,sizeof(ts_idx),1,fh);
                vesicle->filament_list->poly[i]->vlist->vtx[j]->bond[k]=vesicle->filament_list->poly[i]->blist->bond[idx];
            }

        }
    }

    /* restore poly bonds between monomers list*/
    for(i=0;i<vesicle->filament_list->n;i++){
        for(j=0;j<vesicle->filament_list->poly[i]->blist->n;j++){
       //     vesicle->poly_list->poly[i]->blist->bond[j]=(ts_bond *)malloc(sizeof(ts_bond));
            retval=fread(vesicle->filament_list->poly[i]->blist->bond[j],sizeof(ts_bond),1,fh);
            /* restore vtx1 and vtx2 */
                retval=fread(&idx,sizeof(ts_idx),1,fh);
                vesicle->filament_list->poly[i]->blist->bond[j]->vtx1=vesicle->filament_list->poly[i]->vlist->vtx[idx];
                retval=fread(&idx,sizeof(ts_idx),1,fh);
                vesicle->filament_list->poly[i]->blist->bond[j]->vtx2=vesicle->filament_list->poly[i]->vlist->vtx[idx];
        }
    }
    // vesicle->tape=parsetape(command_line_args.tape_fullfilename); //done outside
    // recreating space for cells // 
    vesicle->clist=init_cell_list(tape->ncxmax, tape->ncymax, tape->nczmax, tape->stepsize);
    vesicle->clist->max_occupancy=16;
    // vesicle->tape=(ts_tape *)malloc(sizeof(ts_tape));
    // retval=fread(vesicle->tape, sizeof(ts_tape),1,fh);
    retval=fread(iteration,sizeof(ts_idx),1,fh);
    if(retval); 
    fclose(fh);
    return vesicle;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// parse command line arguments and save in command_line_args
ts_bool parse_args(int argc, char **argv){
    int c;
    int retval; // for mkdir
    struct stat sb;
    sprintf(command_line_args.path, "./"); //clear string;
    sprintf(command_line_args.output_fullfilename,"output.pvd");
    sprintf(command_line_args.dump_fullfilename,"dump.bin");
    sprintf(command_line_args.tape_fullfilename,"tape");
    sprintf(command_line_args.tape_templatefull,"./tape");
    FILE *file;

    while (1){
        static struct option long_options[] =
            {
                {"force-from-tape", no_argument, &(command_line_args.force_from_tape), 1},
                {"reset-iteration-count", no_argument, &(command_line_args.reset_iteration_count), 1},
                {"tape", required_argument, 0, 't'},
                {"version", no_argument, 0, 'v'},
                {"output-file", required_argument, 0, 'o'},
                {"directory", required_argument, 0, 'd'},
                {"dump-filename", required_argument, 0, 'f'},
                {"tape-options", required_argument, 0, 'c'},
                {"tape-template", required_argument, 0, 0},
                {"restore-from-vtk", required_argument, 0, 'r'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "d:f:o:c:t:r:v",
                        long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c){
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            /* printf ("option %s", long_options[option_index].name);
            if (optarg)
            printf (" with arg %s", optarg); 
            printf ("\n"); */
            //TODO: find a better way.
            if (strcmp(long_options[option_index].name, "tape-template") == 0){
                strcpy(command_line_args.tape_templatefull, optarg);
            }
            break;
        case 'v':
            fprintf(stdout, "TRISURF-NG v. %s, compiled on: %s %s.\n", TS_VERSION, __DATE__, __TIME__);
            fprintf(stdout, "Programming done by: Samo Penic and Miha Fosnaric\n");
            fprintf(stdout, "Released under terms of GPLv3\n");
            exit(0);

        case 'c':
            strcpy(command_line_args.tape_opts, optarg);
            break;
        case 't': //tape
            strcpy(command_line_args.tape_fullfilename, optarg);
            break;

        case 'o': //set filename of master pvd output file
            strcpy(command_line_args.output_fullfilename, optarg);
            break;

        case 'd':
            //check if directory exists. If not create one. If creation is
            //successful, set directory for output files.
            //printf ("option -d with value `%s'\n", optarg);
            if (stat(optarg, &sb) == -1)
            {
                //directory does not exist
                retval = mkdir(optarg, 0700);
                if (retval)
                {
                    fatal("Could not create requested directory. Check if you have permissions", 1);
                }
            }
            //check if is a proper directory
            else if ((sb.st_mode & S_IFMT) != S_IFDIR)
            {
                //it is not a directory. fatal error.
                ts_fprintf(stderr, "%s is not a directory!\n", optarg);
                fatal("Cannot continue", 1);
            }
            strcpy(command_line_args.path, optarg);
            break;

        case 'f':
            strcpy(command_line_args.dump_fullfilename, optarg);
            break;

        case 'r':
            strcpy(command_line_args.dump_from_vtk, optarg);
            break;

        case '?':
            /* getopt_long already printed an error message. */
            print_help(stdout);
            //ts_fprintf(stderr,"\n\nhere comes the help.\n\n");
            fatal("Ooops, read help first", 1);
            break;

        default:
            exit(1);
        }
    }

    //Here we set correct values for full filenames!
    char *buffer=(char *)malloc(10000*sizeof(char));
    //correct the path and add trailing /
    if(command_line_args.path[strlen(command_line_args.path)-1]!='/') strcat(command_line_args.path,"/");
   
    /* master pvd output file */ 
    strcpy(buffer,command_line_args.path);
    strcat(buffer,command_line_args.output_fullfilename);
    if ((file = fopen(buffer, "w")) == NULL) {
        fprintf(stderr,"Could not create output file %s!\n", buffer);
        fatal("Please specify correct output file or check permissions of the file",1);
        //there is a tape template. make a copy into desired directory
    } else {
        fclose(file);
        strcpy(command_line_args.output_fullfilename,buffer);
    }

    /* tape file */
    strcpy(buffer,command_line_args.path);
    strcat(buffer,command_line_args.tape_fullfilename);
    if (stat(buffer, &sb) == -1) {
        //tape does not exist. does tape template exist?
        if(stat(command_line_args.tape_templatefull, &sb)==-1){ 
            ts_fprintf(stderr,"Tape '%s' does not exist and no tape template was specified (or does not exist)!\n",buffer);
            fatal("Please select correct tape or check permissions of the file",1);
        } else {
            //tape template found
            fatal("Samo did not program template copy yet",1); 
        }
    } else {
        strcpy(command_line_args.tape_fullfilename,buffer);
    }


    /* dump file */
    strcpy(buffer,command_line_args.path);
    strcat(buffer,command_line_args.dump_fullfilename);
    //check if dump file exist first.
    if (stat(buffer, &sb) == -1) {
        //no dump file. check if we can create one.
        if ((file = fopen(buffer, "w")) == NULL) {
            fprintf(stderr,"Could not create dump file '%s'!\n",buffer);
            fatal("Please specify correct dump file or check permissions of the file",1);
        } else {
            fclose(file);
            //good, file is writeable. delete it for now.
            remove(buffer);
        }
    }

    strcpy(command_line_args.dump_fullfilename, buffer);
    free(buffer);
    return TS_SUCCESS;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts_bool print_help(FILE *fd){
    fprintf(fd,"TRISURF-NG v. %s, compiled on: %s %s.\n", TS_VERSION, __DATE__, __TIME__);
    fprintf(fd,"Programming done by: Samo Penic and Miha Fosnaric\n");
    fprintf(fd,"Additional Programming done by: Raj-Kumar Sadhu, Yoav Ravid\n");
    fprintf(fd,"Released under terms of GPLv3\n\n");

    fprintf(fd, "Invoking trisurf-ng without any flags results in default operation. Program reads 'tape' file and 'dump.bin' binary representation of the vesicle from disk and continues the simulation where it was aborted (as specified in 'dump.bin').\n\n");
    fprintf(fd,"If 'tape' has different values than binary dump, those are used (if possible -- some values in tape, such as nshell cannot be modified).\n\n\n");
    fprintf(fd,"However, if dump.bin is not present, user is notified to specify --force-from-tape flag. The vesicle will be created from specifications in tape.\n\n\n");
    fprintf(fd,"Flags:\n\n");
    fprintf(fd,"--force-from-tape\t\t makes initial shape of the vesicle from tape. Ignores already existing binary dump and possible simulation results.\n");
    fprintf(fd,"--restore-from-vtk (or -r)\t\t VTK's file ending with '.vtu' are preferred way to make state snapshots for restoration. With this flag the restoration of the vesicle from vtk is possible. The simulation will continue if hidden '.status' file with last iteration done is present. Otherwise it will start simulation from timestep 0.\n");
    fprintf(fd,"--reset-iteration-count\t\t starts simulation from the beginning (using binary dump).\n");
    fprintf(fd,"--tape (or -t)\t\t specifies tape filename. For --force-from-tape and restoring from binary dump. Defaults to 'tape'.\n");
    fprintf(fd,"--version (or -v)\t\t Prints version information.\n");
    fprintf(fd,"--output-file (or -o)\t\t Specifies filename of .PVD file. Defaults to 'output.pvd'\n");
    fprintf(fd,"--dump-filename (or -f)\t\t specifies filename for binary dump&restore. Defaults to 'dump.bin'\n");
    fprintf(fd,"--tape-options (or -c)\t\t specifies replacement options by 'opt1=val1,opt2=val2''\n\n\n");
    fprintf(fd,"Examples:\n\n");
    fprintf(fd,"trisurf --force-from-tape\n");
    fprintf(fd,"trisurf --reset-iteration-count\n");
    fprintf(fd,"trisurf --restore-from-vtk filename.vtu\n");
    fprintf(fd,"\n\n");

    return TS_SUCCESS;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts_bool print_vertex_list(ts_vertex_list *vlist){
    ts_idx i;
    printf("Number of vertices: %u\n",vlist->n);
    for(i=0;i<vlist->n;i++){
        printf("%u: %f %f %f\n",
        vlist->vtx[i]->idx,vlist->vtx[i]->x, vlist->vtx[i]->y, vlist->vtx[i]->z);
    }
    return TS_SUCCESS;
}

ts_bool print_vertex_neighbours(ts_vertex_list *vlist){
    ts_idx i;
    ts_small_idx j;
    ts_vertex **vtx=vlist->vtx;
    printf("Vertex id(neigh no): (neighvertex coord) (neighvertex coord) ...\n");
    for(i=0;i<vlist->n;i++){
        printf("%u(%u): ",vtx[i]->idx,vtx[i]->neigh_no);
        for(j=0;j<vtx[i]->neigh_no;j++){
            printf("(%f,%f,%f)",vtx[i]->neigh[j]->x,vtx[i]->neigh[j]->y,vtx[i]->neigh[j]->z);
        }
        printf("\n");
    }

return TS_SUCCESS;
}

ts_bool write_vertex_fcompat_file(ts_vertex_list *vlist,ts_char *filename){
    ts_vertex **vtx=vlist->vtx;
    ts_idx i;
    FILE *fh;
    
    fh=fopen(filename, "w");
    if(fh==NULL){
        err("Cannot open file %s for writing");
        return TS_FAIL;
    }
    for(i=0;i<vlist->n;i++)
        fprintf(fh," %E\t%E\t%E\n",vtx[i]->x,vtx[i]->y, vtx[i]->z);

    fclose(fh);
return TS_SUCCESS;
}


ts_bool fprint_vertex_list(FILE *fh,ts_vertex_list *vlist){
    ts_idx i;
    ts_small_idx j;
    for(i=0;i<vlist->n;i++){
        fprintf(fh," %.17E\t%.17E\t%.17E\t%u\n",vlist->vtx[i]->x,
            vlist->vtx[i]->y, vlist->vtx[i]->z,
            vlist->vtx[i]->neigh_no);
        for(j=0;j<vlist->vtx[i]->neigh_no;j++){
            fprintf(fh,"\t%u",(ts_idx)(vlist->vtx[i]->neigh[j]->idx));
        //-vlist->vtx+1));
        }
        fprintf(fh,"\n");
    }
    return TS_SUCCESS;
}

ts_bool fprint_tristar(FILE *fh, ts_vesicle *vesicle){
    ts_idx i;
    ts_small_idx j;
    for(i=0;i<vesicle->vlist->n;i++){
        fprintf(fh,"\t%u",vesicle->vlist->vtx[i]->tristar_no);
        for(j=0;j<vesicle->vlist->vtx[i]->tristar_no;j++){
            fprintf(fh,"\t%u",(ts_idx)(vesicle->vlist->vtx[i]->tristar[j]->idx));//-vesicle->tlist->tria+1));
        }
        fprintf(fh,"\n");
    }
    return TS_SUCCESS;
}

ts_bool fprint_triangle_list(FILE *fh, ts_vesicle *vesicle){
    ts_triangle_list *tlist=vesicle->tlist;
    ts_idx i;
    ts_small_idx j;
    for(i=0;i<tlist->n;i++){
        fprintf(fh,"\t%u",tlist->tria[i]->neigh_no);
        for(j=0;j<tlist->tria[i]->neigh_no;j++){
            fprintf(fh,"\t%u",(ts_idx)(tlist->tria[i]->neigh[j]->idx));//-tlist->tria+1)); 
        }
        fprintf(fh,"\n");
        for(j=0;j<3;j++){
            fprintf(fh,"\t%u",(ts_idx)(tlist->tria[i]->vertex[j]->idx));//-vesicle->vlist->vtx+1)); 
        }
        fprintf(fh,"\n");
        fprintf(fh,"%.17E\t%.17E\t%.17E\n",tlist->tria[i]->xnorm, tlist->tria[i]->ynorm,tlist->tria[i]->znorm);
        fprintf(fh,"0.00000000000000000\n0.00000000000000000\n");
    }
    return TS_SUCCESS;
}

ts_bool fprint_vertex_data(FILE *fh,ts_vertex_list *vlist){
    ts_idx i;
    ts_small_idx j;
    for(i=0;i<vlist->n;i++){
        fprintf(fh," %.17E\t%.17E\t%.17E\t%.17E\t%u\n",
        vlist->vtx[i]->xk,vlist->vtx[i]->c,vlist->vtx[i]->energy,
        vlist->vtx[i]->mean_curvature, 0);
            fprintf(fh,"\n");
        for(j=0;j<vlist->vtx[i]->bond_no;j++){
            fprintf(fh," %.17E", vlist->vtx[i]->bond[j]->bond_length);
        }
            fprintf(fh,"\n");
    }
    return TS_SUCCESS;
}

ts_bool fprint_bonds(FILE *fh,ts_vesicle *vesicle){
    ts_idx i;
    for(i=0;i<vesicle->blist->n;i++){
        fprintf(fh,"\t%u\t%u\n",(ts_idx)(vesicle->blist->bond[i]->vtx1->idx),
                                //-vesicle->vlist->vtx+1),
                                (ts_idx)(vesicle->blist->bond[i]->vtx2->idx));
                                //-vesicle->vlist.vtx+1));
    }
    return TS_SUCCESS;
}

//write adhesion surface position
ts_bool write_adhesion_position(FILE *fh, ts_vesicle *vesicle){
fprintf(fh,"\t%f\n",vesicle->tape->z_adhesion);
return TS_SUCCESS;
}

ts_bool write_dout_fcompat_file(ts_vesicle *vesicle, ts_char *filename){
    FILE *fh;
    fh=fopen(filename, "w");
    if(fh==NULL){
        err("Cannot open file %s for writing");
        return TS_FAIL;
    }
    fprintf(fh,"%.17E\n%.17E\n",vesicle->stepsize,vesicle->dmax);
    fprint_vertex_list(fh,vesicle->vlist);
    fprint_tristar(fh,vesicle);
    fprint_triangle_list(fh,vesicle);
    fprint_vertex_data(fh,vesicle->vlist);
    fprint_bonds(fh,vesicle);
    write_adhesion_position(fh,vesicle);
    fclose(fh);	
    return TS_SUCCESS;
}

ts_bool read_tape_fcompat_file(ts_vesicle *vesicle, ts_char *filename){
    FILE *fh;
    char line[255];
    fh=fopen(filename, "r");
        if(fh==NULL){
                err("Cannot open file for reading... Nonexistant file?");
                return TS_FAIL;
        }
    ts_uint retval=1; // not sure what's going on here: fscanf returns int, not uint
    while(retval!=EOF){
        retval=fscanf(fh,"%s",line);
        
        fprintf(stderr,"%s",line);
    }	
    fclose(fh);	
    return TS_SUCCESS;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts_bool write_master_xml_file(ts_char *filename){
     FILE *fh;
    ts_char *i,*j;
    ts_idx tstep;
    ts_char *number;
    fh=fopen(filename, "w");
    if(fh==NULL){
        err("Cannot open file %s for writing");
        return TS_FAIL;
    }

    fprintf(fh,"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n<Collection>");
    DIR *dir = opendir(command_line_args.path);
    if(dir){
        struct dirent *ent;
        tstep=0;
        while((ent = readdir(dir)) != NULL)
        {
            i=rindex(ent->d_name,'.');
            if(i==NULL) continue;
            if(strcmp(i+1,"vtu")==0){
                    j=rindex(ent->d_name,'_');
                    if(j==NULL) continue;
                    number=strndup(j+1,j-i); 
                    fprintf(fh,"<DataSet timestep=\"%u\" group=\"\" part=\"0\" file=\"%s\"/>\n",atoi(number),ent->d_name);
                    tstep++;
                    free(number);
            }  
        }
    }
    free(dir);
    fprintf(fh,"</Collection>\n</VTKFile>\n");
    fclose(fh);
    return TS_SUCCESS;
}

ts_bool write_vertex_xml_file(ts_vesicle *vesicle, ts_idx timestepno, ts_cluster_list *cstlist){
    ts_vertex_list *vlist=vesicle->vlist;
    ts_bond_list *blist=vesicle->blist;
    ts_vertex **vtx=vlist->vtx;
    ts_idx i,j;
    //ts_double senergy=0.0;
    char filename[10000];
    char just_name[255];
    FILE *fh;
    strcpy(filename,command_line_args.path);
    sprintf(just_name,"timestep_%.6u.vtu",timestepno);
    strcat(filename,just_name);

    fh=fopen(filename, "w");
    if(fh==NULL){
        err("Cannot open file %s for writing");
        return TS_FAIL;
    }
    /* Here comes header of the file */

    //find number of extra vtxs and bonds of polymeres
    ts_idx monono=0, polyno=0, poly_idx=0, filno=0, fonono=0; //terrified: theres i<(fonono-1)*filno
    ts_bool poly=0, fil=0;
    if(vesicle->poly_list!=NULL){
        if(vesicle->poly_list->poly[0]!=NULL){
        polyno=vesicle->poly_list->n;
        monono=vesicle->poly_list->poly[0]->vlist->n;
        poly=1;
        }
    }

    if(vesicle->filament_list!=NULL){
        if(vesicle->filament_list->poly[0]!=NULL){
        filno=vesicle->filament_list->n;
        fonono=vesicle->filament_list->poly[0]->vlist->n;
        fil=1;
        }
    }

    fprintf(fh, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
    xml_trisurf_data(fh,vesicle);
    fprintf(fh, "<UnstructuredGrid>\n");
    fprintf(fh, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\">\n",vlist->n+monono*polyno+fonono*filno, blist->n+monono*polyno+filno*(fonono-1)+vesicle->tlist->n);
    fprintf(fh,"<PointData Scalars=\"vertices_idx\" Vectors=\"force\">\n<DataArray type=\"Int64\" Name=\"vertices_idx\" format=\"ascii\">");
    for(i=0;i<vlist->n;i++){
        fprintf(fh,"%u ",vtx[i]->idx);
    }
    //polymeres
    if(poly){
        poly_idx=vlist->n;
        for(i=0;i<vesicle->poly_list->n;i++){
            for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++,poly_idx++){
                fprintf(fh,"%u ", poly_idx);
            }
        }
    }
    //filaments
    if(fil){
        poly_idx=vlist->n+monono*polyno;
        for(i=0;i<vesicle->filament_list->n;i++){
            for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++,poly_idx++){
                // fprintf(stderr,"was here\n");
                fprintf(fh,"%u ", poly_idx);
            }
        }
    }

    fprintf(fh,"</DataArray>\n");
    if(cstlist!=NULL){
        fprintf(fh,"<DataArray type=\"Int64\" Name=\"vertices_in_cluster\" format=\"ascii\">");
        for(i=0;i<vlist->n;i++){
            if(vtx[i]->cluster!=NULL){
                fprintf(fh,"%u ",vtx[i]->cluster->nvtx);
            } else {
                fprintf(fh,"-1 ");
            }
        }
        if(poly){
            poly_idx=vlist->n;
            for(i=0;i<vesicle->poly_list->n;i++){
                for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++,poly_idx++){
                    fprintf(fh,"-1 ");
                }
            }
        }
        if(fil){
            poly_idx=vlist->n+monono*polyno;
            for(i=0;i<vesicle->filament_list->n;i++){
                for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++,poly_idx++){
                    fprintf(fh,"-1 ");
                }
            }
        }

        fprintf(fh,"</DataArray>\n");

    }

    //here comes additional data as needed.
    //Yoav : additional scalar data: type, (spontaneous_curature), bonding_strength (w), direct_force (f), 
    // adhesion_strength (ad_w), spontaneous_deviator (d, which has nothing at the moment), bending energy
    // curvature, second_curvature, bending_modulus, second_bending_modulus

    fprintf(fh,"<DataArray type=\"Int64\" Name=\"type\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%u ",type);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"spontaneous_curvature\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",c);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"bonding_strength\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",w);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"direct_force\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",f);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"adhesion_strength\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",ad_w);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"spontaneous_deviator\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",d);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"mean_curvature\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",mean_curvature);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"gaussian_curvature\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",gaussian_curvature);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"mean_curvature2\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",mean_curvature2);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"gaussian_curvature2\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",gaussian_curvature2);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"new_c1\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",new_c1);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"new_c2\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",new_c2);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"eigenvalue_0\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",eig_v0);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"eigenvalue_1\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",eig_v1);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"eigenvalue_2\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",eig_v2);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"bending_modulus\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",xk);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"second_bending_modulus\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",xk2);
    fprintf(fh,"</DataArray>\n");


    //here comes additional data. Energy! different form for polymers, so no macro
    fprintf(fh,"<DataArray type=\"Float64\" Name=\"bending_energy\" format=\"ascii\">");
    for(i=0;i<vlist->n;i++){
        fprintf(fh,"%.17e ",vtx[i]->energy);
    }
    //polymeres
    if(poly){
        for(i=0;i<vesicle->poly_list->n;i++){
            for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
                fprintf(fh,"%.17e ", vesicle->poly_list->poly[i]->vlist->vtx[j]->energy* vesicle->poly_list->poly[i]->k);
            }
        }
    }
    //filaments
    if(fil){
        for(i=0;i<vesicle->filament_list->n;i++){
            for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
                fprintf(fh,"%.17e ",  vesicle->filament_list->poly[i]->vlist->vtx[j]->energy*  vesicle->filament_list->poly[i]->k);
            }
        }
    }
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"mean_energy\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",mean_energy);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"gaussian_energy\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",gaussian_energy);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"mean_energy2\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",mean_energy2);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"gaussian_energy2\" format=\"ascii\">");
    TS_WRITE_ITERATE_VTX("%.17e ",gaussian_energy2);
    fprintf(fh,"</DataArray>\n");


    // Normals
    fprintf(fh,"<DataArray type=\"Float64\" Name=\"normal\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    TS_WRITE_VECTOR_ITERATE_VTX("%.17e %.17e %.17e\n",nx,ny,nz);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"normal2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    TS_WRITE_VECTOR_ITERATE_VTX("%.17e %.17e %.17e\n",nx2,ny2,nz2);
    fprintf(fh,"</DataArray>\n");

    // Vectors: force, director (currently has nothing)
    fprintf(fh,"<DataArray type=\"Float64\" Name=\"force\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    TS_WRITE_VECTOR_ITERATE_VTX("%.17e %.17e %.17e\n",fx,fy,fz);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"director\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    TS_WRITE_VECTOR_ITERATE_VTX("%.17e %.17e %.17e\n",dx,dy,dz);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"eig0\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    TS_WRITE_VECTOR_ITERATE_VTX("%.17e %.17e %.17e\n",eig0[0],eig0[1],eig0[2]);
    fprintf(fh,"</DataArray>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"eig1\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    TS_WRITE_VECTOR_ITERATE_VTX("%.17e %.17e %.17e\n",eig1[0],eig1[1],eig1[2]);
    fprintf(fh,"</DataArray>\n");

     fprintf(fh,"<DataArray type=\"Float64\" Name=\"eig2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    TS_WRITE_VECTOR_ITERATE_VTX("%.17e %.17e %.17e\n",eig2[0],eig2[1],eig2[2]);
    fprintf(fh,"</DataArray>\n");

    //end point data: start of cell data

    fprintf(fh,"</PointData>\n<CellData>\n");

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"bonding_energy\" format=\"ascii\">");
        for(i=0;i<blist->n;i++){
            fprintf(fh,"%.17e ",vesicle->blist->bond[i]->energy);
        }
        for(i=0;i<monono*polyno+filno*(fonono-1);i++){
            fprintf(fh,"0.0 ");
        }
        for(i=0;i<vesicle->tlist->n;i++){
            fprintf(fh,"0.0 ");
        }
        fprintf(fh,"</DataArray>\n");

    if(vesicle->tape->area_switch==1){
        fprintf(fh,"<DataArray type=\"Float64\" Name=\"stretching_energy\" format=\"ascii\">");
        for(i=0;i<blist->n;i++){
            fprintf(fh, "0.0 ");
        }
        for(i=0;i<monono*polyno+filno*(fonono-1);i++){
            fprintf(fh,"0.0 ");
        }
        for(i=0;i<vesicle->tlist->n;i++){
            fprintf(fh,"%.17e ",vesicle->tlist->tria[i]->energy);
        }
        fprintf(fh,"</DataArray>\n");
    }

    fprintf(fh,"<DataArray type=\"Float64\" Name=\"face_normal\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(i=0;i<blist->n;i++){
        fprintf(fh, "0.0 0.0 0.0\n");
    }
    for(i=0;i<monono*polyno+filno*(fonono-1);i++){
        fprintf(fh,"0.0 0.0 0.0\n");
    }
    for(i=0;i<vesicle->tlist->n;i++){
        fprintf(fh,"%.17e %.17e %.17e\n",vesicle->tlist->tria[i]->xnorm, vesicle->tlist->tria[i]->ynorm, vesicle->tlist->tria[i]->znorm);
    }
    fprintf(fh,"</DataArray>\n");



    fprintf(fh,"</CellData>\n<Points>\n<DataArray type=\"Float64\" Name=\"Koordinate tock\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    TS_WRITE_VECTOR_ITERATE_VTX("%.17e %.17e %.17e\n",x,y,z);

    fprintf(fh,"</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">");
    // bonds
    for(i=0;i<blist->n;i++){
            fprintf(fh,"%u %u\n",blist->bond[i]->vtx1->idx,blist->bond[i]->vtx2->idx);
    }
    if(poly){
        for(i=0;i<vesicle->poly_list->n;i++){
            for(j=0;j<vesicle->poly_list->poly[i]->blist->n;j++){
                fprintf(fh,"%u %u\n", vesicle->poly_list->poly[i]->blist->bond[j]->vtx1->idx+vlist->n+i*monono,
                                      vesicle->poly_list->poly[i]->blist->bond[j]->vtx2->idx+vlist->n+i*monono);
            }
        //grafted bonds
        fprintf(fh,"%u %u\n", vesicle->poly_list->poly[i]->grafted_vtx->idx,
                              vesicle->poly_list->poly[i]->vlist->vtx[0]->idx+vlist->n+i*monono);
        }
    }
    if(fil){
        for(i=0;i<vesicle->filament_list->n;i++){
            for(j=0;j<vesicle->filament_list->poly[i]->blist->n;j++){
                fprintf(fh,"%u %u\n", vesicle->filament_list->poly[i]->blist->bond[j]->vtx1->idx+vlist->n+monono*polyno+i*fonono,
                                      vesicle->filament_list->poly[i]->blist->bond[j]->vtx2->idx+vlist->n+monono*polyno+i*fonono);
            }
        }
    }
    //triangles
    for(i=0;i<vesicle->tlist->n;i++){
        fprintf(fh,"%u %u %u\n", vesicle->tlist->tria[i]->vertex[0]->idx, vesicle->tlist->tria[i]->vertex[1]->idx, vesicle->tlist->tria[i]->vertex[2]->idx);
    }

    fprintf(fh,"</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">");
    for (i=2;i<(blist->n+monono*polyno+(fonono-1)*filno)*2+1;i+=2){
    fprintf(fh,"%u ",i);
    }
    for(j=i+1;j<i+3*(vesicle->tlist->n);j+=3){ //let's continue counting from where we left of
        fprintf(fh,"%u ", j);
    }
    fprintf(fh,"\n");

    fprintf(fh,"</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
     for (i=0;i<blist->n+monono*polyno+(fonono-1)*filno;i++){
        fprintf(fh,"3 ");
    }
    for(i=0;i<vesicle->tlist->n;i++){
        fprintf(fh,"5 ");
    }

    fprintf(fh,"</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    fclose(fh);
    return TS_SUCCESS;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts_bool write_vertex_vtk_file(ts_vesicle *vesicle,ts_char *filename, ts_char *text){
    ts_vertex_list *vlist=vesicle->vlist;
    ts_bond_list *blist=vesicle->blist;
    ts_vertex **vtx=vlist->vtx;
    ts_idx i;
    FILE *fh;
    
    fh=fopen(filename, "w");
    if(fh==NULL){
        err("Cannot open file %s for writing");
        return TS_FAIL;
    }
    /* Here comes header of the file */
//    fprintf(stderr,"NSHELL=%u\n",nshell);
    fprintf(fh, "# vtk DataFile Version 2.0\n");
    /* TODO: Do a sanity check on text. Max 255 char, must not me \n terminated */ 
    fprintf(fh, "%s\n", text);
    fprintf(fh,"ASCII\n");
    fprintf(fh,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fh,"POINTS %u double\n", vlist->n);
    for(i=0;i<vlist->n;i++){
        fprintf(fh,"%e %e %e\n",vtx[i]->x,vtx[i]->y, vtx[i]->z);
    }
    
    fprintf(fh,"CELLS %u %u\n",blist->n,3*blist->n);
    for(i=0;i<blist->n;i++){
            fprintf(fh,"2 %u %u\n",blist->bond[i]->vtx1->idx,blist->bond[i]->vtx2->idx);
    }
    fprintf(fh,"CELL_TYPES %u\n",blist->n);
    for(i=0;i<blist->n;i++)
        fprintf(fh,"3\n");

    fprintf(fh,"POINT_DATA %u\n", vlist->n);
    fprintf(fh,"SCALARS scalars long 1\n");
    fprintf(fh,"LOOKUP_TABLE default\n");

    for(i=0;i<vlist->n;i++)
        fprintf(fh,"%u\n",vtx[i]->idx);

    fclose(fh);
    return TS_SUCCESS;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts_bool write_pov_file(ts_vesicle *vesicle, char *filename){
    FILE *fh;
    ts_idx i;
    
    fh=fopen(filename, "w");
    if(fh==NULL){
        err("Cannot open file %s for writing");
        return TS_FAIL;
    }

    for(i=0;i<vesicle->tlist->n;i++){
    
    fprintf(fh,"\ttriangle {");
    fprintf(fh,"\t<%e,%e,%e> <%e,%e,%e> <%e,%e,%e> }\n", 
    vesicle->tlist->tria[i]->vertex[0]->x,
    vesicle->tlist->tria[i]->vertex[0]->y,
    vesicle->tlist->tria[i]->vertex[0]->z,

    vesicle->tlist->tria[i]->vertex[1]->x,
    vesicle->tlist->tria[i]->vertex[1]->y,
    vesicle->tlist->tria[i]->vertex[1]->z,

    vesicle->tlist->tria[i]->vertex[2]->x,
    vesicle->tlist->tria[i]->vertex[2]->y,
    vesicle->tlist->tria[i]->vertex[2]->z
    );
    }
        
    fclose(fh);
    return TS_SUCCESS;
}


ts_tape *parsetape(char *filename){
    FILE *fd = fopen (filename, "r");
    long length;
    size_t size;
    fseek (fd, 0, SEEK_END);
    length = ftell (fd);
    fseek (fd, 0, SEEK_SET);
    size=fread (tapetxt, 1, length, fd);
    fclose(fd);
    if(size);//?


    // fix the tape that is recorded to the .vtu
    // with the new command line arguments
    update_tapetxt(tapetxt, command_line_args.tape_opts);
    ts_tape *tape=parsetapebuffer(tapetxt);
    return tape;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts_tape *parsetapebuffer(char *buffer){
    ts_tape *tape=(ts_tape *)calloc(1,sizeof(ts_tape));
    //tape->multiprocessing=calloc(255,sizeof(char));
    
    cfg_opt_t opts[] = {
        CFG_INT("nshell", 0, CFGF_NONE),
        CFG_INT("npoly", 0, CFGF_NONE),
        CFG_INT("nmono", 0, CFGF_NONE),
        CFG_INT("nfil", 0, CFGF_NONE),
        CFG_INT("nfono", 0, CFGF_NONE),
        CFG_INT("internal_poly", 0, CFGF_NONE),
        CFG_FLOAT("R_nucleus", 0, CFGF_NONE),
        CFG_FLOAT("R_nucleusX", 0, CFGF_NONE),
        CFG_FLOAT("R_nucleusY", 0, CFGF_NONE),
        CFG_FLOAT("R_nucleusZ", 0, CFGF_NONE),
        CFG_FLOAT("dmax", 1.7, CFGF_NONE),
        CFG_FLOAT("dmin_interspecies", 1.2, CFGF_NONE),
        CFG_FLOAT("xk0", 20, CFGF_NONE),
        CFG_FLOAT("xk2", 0, CFGF_NONE),
        CFG_INT("pressure_switch", 0, CFGF_NONE),
        CFG_INT("volume_switch", 0, CFGF_NONE),
        CFG_INT("area_switch", 0, CFGF_NONE),
        CFG_FLOAT("constvolprecision", 0, CFGF_NONE),
        CFG_FLOAT("xkA0", 1.0, CFGF_NONE),
        CFG_FLOAT("xkV0", 1.0, CFGF_NONE),
        CFG_FLOAT("V0", 0.0, CFGF_NONE),
        CFG_FLOAT("A0", 0.0, CFGF_NONE),
        CFG_FLOAT("Vr", 1.0, CFGF_NONE),
        CFG_FLOAT("pressure", 0, CFGF_NONE),
        CFG_FLOAT("k_spring", 800, CFGF_NONE),
        CFG_FLOAT("xi", 100, CFGF_NONE),
        CFG_FLOAT("stepsize", 0.15, CFGF_NONE),
        CFG_INT("nxmax", 100, CFGF_NONE),
        CFG_INT("nymax", 100, CFGF_NONE),
        CFG_INT("nzmax", 100, CFGF_NONE),
        CFG_INT("iterations", 100, CFGF_NONE),
        CFG_INT("mcsweeps", 200000, CFGF_NONE),
        CFG_INT("inititer", 0, CFGF_NONE),
        CFG_BOOL("quiet", 0, CFGF_NONE),
        //CFG_SIMPLE_STR("multiprocessing",&tape->multiprocessing),
        //CFG_SIMPLE_INT("smp_cores",&tape->brezveze0),
        //CFG_SIMPLE_INT("cluster_nodes",&tape->brezveze1),
        //CFG_SIMPLE_INT("distributed_processes",&tape->brezveze2),
        CFG_INT("spherical_harmonics_coefficients", 0, CFGF_NONE),
        CFG_INT("number_of_vertices_with_c0", 50, CFGF_NONE),
        CFG_INT("adhesion_geometry", 1, CFGF_NONE),
        CFG_INT("allow_xy_plane_movement", 0, CFGF_NONE),
        CFG_INT("force_balance_along_z_axis", 0, CFGF_NONE),
        CFG_FLOAT("c0", 0.5, CFGF_NONE),
        CFG_FLOAT("w", 1.5, CFGF_NONE),
        CFG_FLOAT("F", 0.7, CFGF_NONE),
        /* Variables related to plane confinement */
        CFG_INT("plane_confinement_switch", 0, CFGF_NONE),
        CFG_FLOAT("plane_d", 15, CFGF_NONE),
        CFG_FLOAT("plane_F", 1000, CFGF_NONE),
        /* Variables related to adhesion */
        CFG_INT("adhesion_model", 0, CFGF_NONE),
        CFG_FLOAT("adhesion_cuttoff", 15, CFGF_NONE),
        CFG_FLOAT("adhesion_strength", 1000, CFGF_NONE),
        CFG_FLOAT("adhesion_radius", 1000, CFGF_NONE),
        CFG_FLOAT("z_adhesion", 1000, CFGF_NONE),
        /* variables for Vicsek interaction and general interaction modification*/
        CFG_INT("force_model", 0, CFGF_NONE),
        CFG_FLOAT("vicsek_strength", 0.1, CFGF_NONE),
        CFG_FLOAT("vicsek_radius", 1.0, CFGF_NONE),
        CFG_INT("bond_model", 0, CFGF_NONE),
        CFG_INT("curvature_model", 0, CFGF_NONE),
        /* Dihedral angle cosine constraint*/
        CFG_FLOAT("min_dihedral_angle_cosine",-1,CFGF_NONE),
        CFG_FLOAT("d0", 0.5, CFGF_NONE),
        /* random seed */
        CFG_INT("random_seed",0,CFGF_NONE),
        CFG_END()
    };
    cfg_t *cfg;    
    ts_uint retval; // not int?
    cfg = cfg_init(opts, 256); //consider using CFGF_IGNORE_UNKNOWN
    retval = cfg_parse_buf(cfg, buffer);
    tape->nshell = cfg_getint(cfg,"nshell");
    tape->npoly = cfg_getint(cfg,"npoly");
    tape->nmono = cfg_getint(cfg,"nmono");
    tape->nfil = cfg_getint(cfg,"nfil");
    tape->nshell = cfg_getint(cfg,"nshell");
    tape->nfono = cfg_getint(cfg,"nfono");
    tape->internal_poly = cfg_getint(cfg,"internal_poly");
    tape->R_nucleus = cfg_getfloat(cfg,"R_nucleus");
    tape->R_nucleusX = cfg_getfloat(cfg,"R_nucleusX");
    tape->R_nucleusY = cfg_getfloat(cfg,"R_nucleusY");
    tape->R_nucleusZ = cfg_getfloat(cfg,"R_nucleusZ");
    tape->dmax = cfg_getfloat(cfg,"dmax");
    tape->dmin_interspecies = cfg_getfloat(cfg,"dmin_interspecies");
    tape->xk0 = cfg_getfloat(cfg,"xk0");
    tape->xk2 = cfg_getfloat(cfg,"xk2");
    tape->pressure_switch = cfg_getint(cfg,"pressure_switch");
    tape->volume_switch = cfg_getint(cfg,"volume_switch");
    tape->area_switch = cfg_getint(cfg,"area_switch");
    tape->constvolprecision = cfg_getfloat(cfg,"constvolprecision");
    tape->xkA0 = cfg_getfloat(cfg,"xkA0");
    tape->xkV0 = cfg_getfloat(cfg,"xkV0");
    tape->V0 = cfg_getfloat(cfg,"V0");
    tape->A0 = cfg_getfloat(cfg,"A0");
    tape->Vfraction = cfg_getfloat(cfg,"Vr");
    tape->pressure = cfg_getfloat(cfg,"pressure");
    tape->kspring = cfg_getfloat(cfg,"k_spring");
    tape->xi = cfg_getfloat(cfg,"xi");
    tape->stepsize = cfg_getfloat(cfg,"stepsize");
    tape->ncxmax = cfg_getint(cfg,"nxmax");
    tape->ncymax = cfg_getint(cfg,"nymax");
    tape->nczmax = cfg_getint(cfg,"nzmax");
    tape->iterations = cfg_getint(cfg,"iterations");
    tape->mcsweeps = cfg_getint(cfg,"mcsweeps");
    tape->inititer = cfg_getint(cfg,"inititer");  
    tape->quiet = cfg_getbool(cfg,"quiet"); 
    tape->shc = cfg_getint(cfg,"spherical_harmonics_coefficients");
    tape->number_of_vertices_with_c0 = cfg_getint(cfg,"number_of_vertices_with_c0");
    tape->adhesion_geometry = cfg_getint(cfg,"adhesion_geometry");
    tape->allow_xy_plane_movement = cfg_getint(cfg,"allow_xy_plane_movement");
    tape->force_balance_along_z_axis = cfg_getint(cfg,"force_balance_along_z_axis");
    tape->c0 = cfg_getfloat(cfg,"c0");
    tape->w = cfg_getfloat(cfg,"w");
    tape->F = cfg_getfloat(cfg,"F");
    tape->plane_confinement_switch = cfg_getint(cfg, "plane_confinement_switch");
    tape->plane_d = cfg_getfloat(cfg, "plane_d");
    tape->plane_F = cfg_getfloat(cfg, "plane_F");
    tape->type_of_force_model = cfg_getint(cfg, "force_model");
    tape->vicsek_strength = cfg_getfloat(cfg, "vicsek_strength");
    tape->vicsek_radius = cfg_getfloat(cfg, "vicsek_radius");
    tape->adhesion_model = cfg_getint(cfg, "adhesion_model");
    tape->adhesion_cuttoff = cfg_getfloat(cfg, "adhesion_cuttoff");
    tape->adhesion_strength = cfg_getfloat(cfg, "adhesion_strength");
    tape->adhesion_radius = cfg_getfloat(cfg, "adhesion_radius");
    tape->z_adhesion = cfg_getfloat(cfg, "z_adhesion");
    tape->random_seed = cfg_getint(cfg, "random_seed");
    tape->type_of_bond_model = cfg_getint(cfg, "bond_model");
    tape->type_of_curvature_model = cfg_getint(cfg, "curvature_model");
    tape->min_dihedral_angle_cosine = cfg_getfloat(cfg,"min_dihedral_angle_cosine");
    tape->d0 = cfg_getfloat(cfg,"d0");

    if (retval==CFG_FILE_ERROR){
        fatal("No tape file.", 100);
    }
    else if (retval == CFG_PARSE_ERROR){
        fatal("Invalid tape!", 100);
    }

    // this bit is not needed, since we already re-wrote the tape directly
    // In order toi have the changes applied here and saved to .vtu

    /* here we used to override all values read from tape with values from commandline*/
    //getcmdline_tape(cfg,command_line_args.tape_opts);
    cfg_free(cfg);

    /* global variables are set automatically */
    quiet = tape->quiet;
    return tape;
}

ts_bool tape_free(ts_tape *tape){
    //free(tape->multiprocessing);
    free(tape);
    return TS_SUCCESS;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts_bool getcmdline_tape(cfg_t *cfg, char *opts){

    char *commands, *backup, *saveptr, *saveopptr, *command, *operator[2];
    operator[0]=0;
    operator[1]=0;
    ts_uint i,j;
    commands=(char *)malloc(10000*sizeof(char));
    backup=commands; //since the pointer to commands will be lost, we acquire a pointer that will serve as backup.
    strcpy(commands,opts);
    for(i=0; ;i++, commands=NULL){
        //breaks comma separated list of commands into specific commands.
        command=strtok_r(commands,",",&saveptr);	
        if(command==NULL) break;
//		fprintf(stdout,"Command %d: %s\n",i,command);	
        //extracts name of command and value of command into operator[2] array.
        for(j=0; j<2;j++,command=NULL){
            operator[j]=strtok_r(command,"=",&saveopptr);
            if(operator[j]==NULL) break;
//			fprintf(stdout," ---> Operator %d: %s\n",j,operator[j]);		
        }
        //1. check: must have 2 operators.
        if(j!=2) fatal("Error. Command line tape options are not formatted properly",1);

    //    cfg_setstr(cfg,operator[0],operator[1]);
        cmdline_to_tape(cfg,operator[0],operator[1]);
        //2. check: must be named properly.
        //3. check: must be of right format (integer, double, string, ...)
        
    }
    free(backup);
    return TS_SUCCESS;
}


ts_bool cmdline_to_tape(cfg_t *cfg, char *key, char *val){

    cfg_opt_t *cfg_opt=cfg_getopt(cfg,key);
    if(cfg_opt==NULL) fatal("Commandline tape option not recognised",1); //return TS_FAIL; 
    switch (cfg_opt->type){
        case CFGT_INT:
            cfg_setint(cfg,key,atol(val));
            break;
        case CFGT_FLOAT:
            cfg_setfloat(cfg,key,atof(val));
            break;
/*        case CFGT_BOOL:
            cfg_setbool(cfg,operator[0],operator[1]);
            break; */
        case CFGT_STR:
            cfg_setstr(cfg,key,val);
            break;
        default:
            break;

    }
    return TS_SUCCESS;
}


ts_bool read_geometry_file(char *fname, ts_vesicle *vesicle){
    FILE *fh;
    ts_idx i, nvtx,nedges,ntria;
    ts_idx vtxi1,vtxi2;
    float x,y,z;
    ts_vertex_list *vlist;
    fh=fopen(fname, "r");
        if(fh==NULL){
                err("Cannot open file for reading... Nonexistant file?");
                return TS_FAIL;
        }
    ts_uint retval;
    retval=fscanf(fh,"%u %u %u",&nvtx, &nedges, &ntria);
    vesicle->vlist=init_vertex_list(nvtx);
    vlist=vesicle->vlist;
    for(i=0;i<nvtx;i++){
   //     fscanf(fh,"%F %F %F",&vlist->vtx[i]->x,&vlist->vtx[i]->y,&vlist->vtx[i]->z);
       retval=fscanf(fh,"%F %F %F",&x,&y,&z);
        vlist->vtx[i]->x=x;
        vlist->vtx[i]->y=y;
        vlist->vtx[i]->z=z;
    }
    for(i=0;i<nedges;i++){
        retval=fscanf(fh,"%u %u",&vtxi1,&vtxi2);
        bond_add(vesicle->blist,vesicle->vlist->vtx[vtxi1-1],vesicle->vlist->vtx[vtxi2-1]);
    }
    //TODO: neighbours from bonds,
    //TODO: triangles from neigbours

//    Don't need to read triangles. Already have enough data
    /*
    for(i=0;i<ntria;i++){
        retval=fscanf(fh,"%u %u %u", &bi1, &bi2, &bi3);
        vtxi1=vesicle->blist->vertex1->idx;
        vtxi2=vesicle->blist->vertex1->idx;
        
    }
    */
    if(retval);
    fclose(fh);	



    return TS_SUCCESS;
}


/**
 * @brief update a tape pointed by tape_txt by the string in cmd_line_tape_args.
 * very hacky newbie c-code running over strings with pointers
 * now preserves newlines!
 * 
 * @param tape_txt a tape file, lines of comments "# comment" and options "optname=value"
 * @param cmd_line_tape_args a string of options in format "opt1=val1,opt2=val2"
 * @return ts_bool TS_SUCCESS if worked (always)
 */
ts_bool update_tapetxt(char* tape_txt, char* cmd_line_tape_args){
    /*
    step 1: check if there are tape arguments from the command line
    step 2: find the number of options
    step 3: go over the tape, copy to new, replacing any option with the cmd_line args
    step 4: add all remaining new options in cmd_line at the end
    step 5: copy the new tape buffer to the tape and cleanup
    */

    char tapetxt_2[128000]; //tmp buffer, used to rebuild text file of the tape
    char* arg_str=cmd_line_tape_args; //hold the new options
    char* tape_p=tape_txt;      //pointer tracing the reading of tape_txt
    char* new_tape_p=tapetxt_2; //pointer tracing the writing of tape_txt_2

    ts_uint i=0;
    ts_uint num_opts=0; // number of options
    ts_bool *is_opt_transfered;  // hold which options have been transfered to tape
    ts_bool keep_old=0;  // no update, keep the old option
    ts_uint option_name_len; // length of the option name in the string "opt=1"=>3
    ts_uint option_value_len; // length of option name and the value part of the string "opt=1"=>5
    ts_uint line_len; // length of current line of tape being read

    // step 1: check if there are tape arguments
    if (arg_str!=NULL){
        num_opts++; // assume there is at least one option
    }
    else { //null pointer
        return TS_SUCCESS;
    }

    //step 2: find number of options in ,,x=1,y=2,,z,
    // we start assuming there is 1 option
    // we count how many blocks of ',' are
    // we discard the last block of ,,, if it's at the end
    // x => +. => 1
    // x=1,y=2 => +...+... => 2
    // ,x=1,,y=2,, => +...+....+- => 2
    for ( i=1  ; arg_str[i] != '\0' ;  i++ ) { //string iterations
        if (arg_str[i]==',' && arg_str[i-1]!=','){
            num_opts++;
        }
        if (i>100000) {
            fatal("passed 100,000 chars in  --tape-options. giving up!\n",100);
        }
    }
    if (arg_str[i-1]==',') num_opts--; // if last character is a separator, we had a final ,,,, run, so we ignore it
    is_opt_transfered = (char*) calloc(num_opts,sizeof(char));

    // step 3: copy the tape to the temporary new tape buffer line by line, replace any option with the updated one
    for (tape_p = tape_txt; *tape_p != '\0' ; tape_p += line_len+1){
        line_len = strcspn(tape_p,"\n");
        keep_old = 1;
        //write spaces and tabs:
        while (*tape_p==' ' || *tape_p=='\t'){
            *new_tape_p=*tape_p;
            tape_p++;
            new_tape_p++;
            line_len--;
        }

        // go over the new options in args and see if this can be replaced
        arg_str = cmd_line_tape_args;
        for (i=0; i<num_opts; i++){
            while(*arg_str==',') arg_str++; // skip ','
            option_name_len=strcspn(arg_str,"=,");
            option_value_len=strcspn(arg_str,",");

            // if the line is one of the options to replace, write from cmd_options
            if (strncmp(tape_p, arg_str, option_name_len)==0                                            // start with option name
                && option_name_len!=0                                                                   // and the option is not empty 
                && (tape_p[option_name_len]==' ' || tape_p[option_name_len]=='=' || tape_p[option_name_len]=='\n') // and not the start of something else
                ){
                strncpy(new_tape_p, arg_str, option_value_len);
                new_tape_p[option_value_len]='\n';
                // move pointer of new tape to the end (and end the string)
                new_tape_p += option_value_len+1;
                is_opt_transfered[i]=1;
                new_tape_p[0]='\0';

                keep_old = 0;
                break;
            }
            if (i<num_opts-1) arg_str += option_value_len+1; // move to next option
        }
        // if we didn't replace, keep the old option
        if (keep_old){
            strncpy(new_tape_p, tape_p, line_len+1);
            new_tape_p += line_len;
            if (new_tape_p[0]!='\0') new_tape_p++;
            new_tape_p[0]='\0';
        }
    }

    // step 4: add all remaining new options in cmd_line at the end
    arg_str = cmd_line_tape_args;
    for (i=0; i<num_opts; i++){ // i=0, arg_str = cmd_line_tape_options...
        while(*arg_str==',') arg_str++; // skip ','
        option_name_len=strcspn(arg_str,"=,");
        option_value_len=strcspn(arg_str,",");
    
        if (is_opt_transfered[i] || option_value_len==0) {
            arg_str += option_value_len+1; // skip options that were already trasfered
            continue;
        }
        else {
            strcpy(new_tape_p,"\n#--tape-options added\n");
            new_tape_p += 23; // length of that ^
            strncpy(new_tape_p, arg_str, option_value_len);
            new_tape_p[option_value_len]='\n';
            new_tape_p += option_value_len+1;
            arg_str += option_value_len+1;
            is_opt_transfered[i]=1;
        }  
    }
    *new_tape_p = '\0'; // end string

    // step 5: copy the new tape buffere on the tape and cleanup
    strcpy(tape_txt,tapetxt_2);
    free(is_opt_transfered);

    return TS_SUCCESS;
}

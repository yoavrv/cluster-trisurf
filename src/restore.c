/* vim: set ts=4 sts=4 sw=4 noet : */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <general.h>
#include <restore.h>
#include <snapshot.h>
#include <zlib.h>
#include "vesicle.h"
#include "vertex.h"
#include "triangle.h"
#include "bond.h"
#include "energy.h"
#include "poly.h"
#include "initial_distribution.h"
#include "io.h"
#include <math.h>

ts_vesicle *parseDump(char *dumpfname) {
	xmlDocPtr doc;
	xmlNodePtr cur, cur1,cur2;
	ts_vesicle *vesicle=NULL;
	doc = xmlParseFile(dumpfname);
	int i;
	if (doc == NULL ) {
		fatal("Dump file could not be found or parsed. It is correct file?",1);
	}
	
	cur = xmlDocGetRootElement(doc);
	
	if (cur == NULL) {
		fatal("Dump file is empty.",1);
	}
	
	if (xmlStrcmp(cur->name, (const xmlChar *) "VTKFile")) {
		fatal("document of the wrong type, root node != story",1);
	}
	
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {

		if ((!xmlStrcmp(cur->name, (const xmlChar *)"tape"))){
			setGlobalTapeTXTfromTapeTag(doc, cur);
		}

		if ((!xmlStrcmp(cur->name, (const xmlChar *)"trisurf"))){
			vesicle=parseTrisurfTag(doc, cur);
		}
		// START Point Position data &  Bonds
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"UnstructuredGrid"))){
			cur1 = cur->xmlChildrenNode;
			while(cur1!=NULL){
				if ((!xmlStrcmp(cur1->name, (const xmlChar *)"Piece"))){
					cur2=cur1->xmlChildrenNode;
					while(cur2!=NULL){
						if ((!xmlStrcmp(cur2->name, (const xmlChar *)"PointData"))){
							if(vesicle!=NULL)
								parseXMLPointData(vesicle,doc,cur2);
						}
						if ((!xmlStrcmp(cur2->name, (const xmlChar *)"Points"))){
							//fprintf(stderr,"Found point data\n");
							if(vesicle!=NULL)
								//fprintf(stderr,"Fils: %ld, Nfono: %ld\n", vesicle->tape->nfil, vesicle->tape->nfono);
								parseXMLVertexPosition(vesicle, doc, cur2);
						}
						if ((!xmlStrcmp(cur2->name, (const xmlChar *)"Cells"))){
						//fprintf(stderr,"Found cell(Bonds) data\n");
							if(vesicle!=NULL)
								parseXMLBonds(vesicle, doc, cur2);
						}
						cur2=cur2->next;
					}	
				}
				cur1 = cur1->next;
			}
		}
		// END Point Position data & Bonds
	cur = cur->next;
	}
	xmlFreeDoc(doc);

//	vesicle->poly_list=init_poly_list(0, 0, vesicle->vlist, vesicle);
	set_vesicle_values_from_tape(vesicle);
	init_normal_vectors(vesicle->tlist);
	mean_curvature_and_energy(vesicle);
	sweep_attraction_bond_energy(vesicle);
	if(vesicle->tape->stretchswitch==1){
		vesicle->tlist->a0=sqrt(3)/4.0*pow((vesicle->tape->dmax+1.0)/2.0,2);  
		for(i=0;i<vesicle->tlist->n;i++){
			stretchenergy(vesicle, vesicle->tlist->tria[i]);
		}
	}
/* TODO: filaments */

//	ts_fprintf(stdout,"Restoration completed\n");
//	write_vertex_xml_file(vesicle,999);
//	vesicle_free(vesicle);
//	exit(0);
	return vesicle;
}

ts_bool setGlobalTapeTXTfromTapeTag(xmlDocPtr doc, xmlNodePtr cur){
	xmlChar *tape = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	strcpy(tapetxt,(char *)tape);
	xmlFree(tape);
	return TS_SUCCESS;
}


/* this is a parser of additional data in xml */
ts_vesicle *parseTrisurfTag(xmlDocPtr doc, xmlNodePtr cur){
	//fprintf(stderr,"Parsing trisurf tag\n");
	xmlNodePtr child;

#ifdef COMPRESS
	/* base64decode */
	size_t cLen;
	/*size_t tLen;
	const unsigned char test[]="Test";
	char *cTest=base64_encode(test, 4,&tLen);
	unsigned char *cuTest=base64_decode((char *)cTest,tLen,&tLen);
	cuTest[tLen]=0;
	fprintf(stderr,"%s\n",cuTest);
	*/
	xmlChar *b64=xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	unsigned char *compressed=base64_decode((char *)b64,strlen((char *)b64)-1,&cLen);
	/* uncompress */
	unsigned char *subtree=(unsigned char *)malloc(512000*sizeof(unsigned char)); /* TODO: again, the uncompressed string must not exceed this */
	z_stream infstream;
	infstream.zalloc = Z_NULL;
	infstream.zfree = Z_NULL;
	infstream.opaque = Z_NULL;
	infstream.avail_in = (ts_uint)cLen; // size of input
    	infstream.next_in = compressed; // input char array
    	infstream.avail_out = (ts_uint)512000; // size of output
    	infstream.next_out = subtree; // output char array
     
    	// the actual DE-compression work.
    	inflateInit(&infstream);
    	inflate(&infstream, Z_NO_FLUSH);
    	inflateEnd(&infstream);	
	//fprintf(stderr,"%lu\n",cLen);
	subtree[infstream.total_out]='\0'; //zero terminate string	
	//fprintf(stderr,"%s\n",subtree);
	
	free(subtree);
#endif
	/*parse xml subtree */
	xmlChar *nvtx, *npoly, *nmono;
	nvtx = xmlGetProp(cur, (xmlChar *)"nvtx");
	npoly=xmlGetProp(cur, (xmlChar *)"npoly");
	nmono=xmlGetProp(cur, (xmlChar *)"nmono");
	ts_tape *tape=parsetapebuffer(tapetxt);
	//fprintf(stderr,"nvtx=%u\n",atoi((char *)nvtx));
	//TODO: check if nvtx is in agreement with nshell from tape
	ts_vesicle *vesicle=init_vesicle(atoi((char *)nvtx),tape->ncxmax,tape->ncymax,tape->nczmax,tape->stepsize);
//	vesicle->poly_list=init_poly_list(atoi((char *)npoly),atoi((char *)nmono), vesicle->vlist, vesicle);
	vesicle->poly_list=init_empty_poly_list(atoi((char *)npoly),atoi((char *)nmono));
	xmlFree(nvtx);
	xmlFree(npoly);
	xmlFree(nmono);

	child = cur->xmlChildrenNode;
	while (child != NULL) {
		if ((!xmlStrcmp(child->name, (const xmlChar *)"vtxn"))){
			parseTrisurfVtxn(vesicle->vlist, doc, child);
		}
		if ((!xmlStrcmp(child->name, (const xmlChar *)"tria"))){
			parseTrisurfTria(vesicle, doc, child);
		}
		if ((!xmlStrcmp(child->name, (const xmlChar *)"trianeigh"))){
			parseTrisurfTriaNeigh(vesicle, doc, child);
		}
		 if ((!xmlStrcmp(child->name, (const xmlChar *)"tristar"))){
			parseTrisurfTristar(vesicle, doc, child);
		}
		 if ((!xmlStrcmp(child->name, (const xmlChar *)"nucleus"))){
			parseTrisurfNucleus(vesicle, doc, child);
		}
	 	 if ((!xmlStrcmp(child->name, (const xmlChar *)"constant_volume"))){
			parseTrisurfConstantVolume(doc, child);
		}
	 	 if ((!xmlStrcmp(child->name, (const xmlChar *)"constant_area"))){
			parseTrisurfConstantArea(doc, child);
		}


	child = child->next;
	}

	vesicle->tape=tape;
//	set_vesicle_values_from_tape(vesicle);
	return vesicle;
}



/* Low level tags parsers */
ts_bool parseTrisurfConstantVolume(xmlDocPtr doc, xmlNodePtr cur){
	xmlChar *cvol = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	char *n=(char *)cvol;
	V0=atof(n);
	xmlFree(cvol);
	return TS_SUCCESS;
}
ts_bool parseTrisurfConstantArea(xmlDocPtr doc, xmlNodePtr cur){
	xmlChar *carea = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	char *n=(char *)carea;
	A0=atof(n);
	xmlFree(carea);
	return TS_SUCCESS;
}

ts_bool parseTrisurfNucleus(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur){
	xmlChar *coords = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	char *n=(char *)coords;
	char *token=strtok(n," ");
	ts_uint i;
	for(i=0;i<3;i++){
		vesicle->nucleus_center[i]=atof(token);
		token=strtok(NULL," ");
	}
	xmlFree(coords);
	return TS_SUCCESS;
}

ts_bool parseTrisurfVtxn(ts_vertex_list *vlist, xmlDocPtr doc, xmlNodePtr cur){

	xmlChar *chari;
	xmlChar *neighs;
	char *n;
	char *token;
	ts_uint neighi;
	ts_uint i;
	chari = xmlGetProp(cur, (xmlChar *)"idx");
	i=atoi((char *)chari);
	xmlFree(chari);
	ts_vertex *vtx=vlist->vtx[i];
	neighs = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	//fprintf(stderr,"Found neigh for vtx %u that seems to have index %u with neighs=%s\n",i,vtx->idx,neighs);

	n=(char *)neighs;
	token=strtok(n," ");
	while(token!=NULL){
		neighi=atoi(token);
		//fprintf(stderr,"%u", neighi);
		vtx_add_neighbour(vtx,vlist->vtx[neighi]);
		token=strtok(NULL," ");
	}	
	xmlFree(neighs);
	return TS_SUCCESS;
}

ts_bool parseTrisurfTria(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur){
	xmlChar *triangles;
	char *tria;
	char *vtx[3];
	
	ts_uint i,j;
	triangles = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	tria=(char *)triangles;
	vtx[0]=strtok(tria," ");
	for(i=1;i<3;i++)	vtx[i]=strtok(NULL," ");
	j=0;
	while(vtx[2]!=NULL){
		triangle_add(vesicle->tlist, vesicle->vlist->vtx[atoi(vtx[0])],vesicle->vlist->vtx[atoi(vtx[1])],vesicle->vlist->vtx[atoi(vtx[2])]);
		for(i=0;i<3;i++)	vtx[i]=strtok(NULL," ");
		j++;
	}	
	//fprintf(stderr,"Parsing triangles %s j=%d\n",triangles,j);	

	xmlFree(triangles);
	return TS_SUCCESS;
}


ts_bool parseTrisurfTriaNeigh(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur){
	xmlChar *triangles;
	char *tria;
	char *ntria[3];
	ts_uint i,j;
	triangles = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	tria=(char *)triangles;
	ntria[0]=strtok(tria," ");
	for(i=1;i<3;i++)	ntria[i]=strtok(NULL," ");
	j=0;
	while(ntria[2]!=NULL){
		triangle_add_neighbour(vesicle->tlist->tria[j],vesicle->tlist->tria[atoi(ntria[0])]);
		triangle_add_neighbour(vesicle->tlist->tria[j],vesicle->tlist->tria[atoi(ntria[1])]);
		triangle_add_neighbour(vesicle->tlist->tria[j],vesicle->tlist->tria[atoi(ntria[2])]);
		j++;
		for(i=0;i<3;i++)	ntria[i]=strtok(NULL," ");
	}	
	//fprintf(stderr,"Parsing triangle neighbors j=%d\n",j);	

	xmlFree(triangles);
	return TS_SUCCESS;
}


ts_bool parseTrisurfTristar(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur){

	xmlChar *chari;
	xmlChar *tristar;
	char *t;
	char *token;
	ts_uint neighi;
	ts_uint i;
	chari = xmlGetProp(cur, (xmlChar *)"idx");
	i=atoi((char *)chari);
	xmlFree(chari);
	ts_vertex *vtx=vesicle->vlist->vtx[i];
	tristar = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
//	fprintf(stderr,"Found tristar for vtx %u that seems to have index %u with tristar=%s\n",i,vtx->idx,tristar);

	t=(char *)tristar;
	token=strtok(t," ");
	while(token!=NULL){
		neighi=atoi(token);
		//fprintf(stderr,"%u", neighi);
		vertex_add_tristar(vtx,vesicle->tlist->tria[neighi]);
		token=strtok(NULL," ");
	}	
	xmlFree(tristar);
	return TS_SUCCESS;
}

/* this parses the data for vertices (like spontaneous curvature, etc.) */
ts_bool parseXMLPointData(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur){
	xmlNodePtr child = cur->xmlChildrenNode;
	xmlChar *property_name;
	xmlChar *values;
	char *vals;
	char *token;
	int idx, polyidx, monoidx, filidx, fonoidx;
	while (child != NULL) {
		if ((!xmlStrcmp(child->name, (const xmlChar *)"DataArray"))){
			property_name=xmlGetProp(child, (xmlChar *)"Name");
		//	fprintf(stderr,"Name: %s\n", property_name);
			if(!xmlStrcmp(property_name,(const xmlChar *)"spontaneous_curvature")){
				values=xmlNodeListGetString(doc,child->xmlChildrenNode,1);
				vals=(char *)values;
				token=strtok(vals," ");
				idx=0;
				while(token!=NULL){
					if(idx<vesicle->vlist->n){
						vesicle->vlist->vtx[idx]->c=atof(token);
					} else if(vesicle->tape->nmono && vesicle->tape->npoly && idx<vesicle->vlist->n+vesicle->tape->nmono*vesicle->tape->npoly) {
						polyidx=(idx-vesicle->vlist->n)/vesicle->tape->nmono;
						monoidx=(idx-vesicle->vlist->n)%vesicle->tape->nmono;
						vesicle->poly_list->poly[polyidx]->vlist->vtx[monoidx]->c=atof(token);
					} else {
						filidx=(idx-vesicle->vlist->n-vesicle->tape->nmono*vesicle->tape->npoly)/vesicle->tape->nfono;
						fonoidx=(idx-vesicle->vlist->n-vesicle->tape->nmono*vesicle->tape->npoly)%vesicle->tape->nfono;
						//fprintf(stderr,"filidx=%d, fonoidx=%d, coord=%s,%s,%s\n",filidx,fonoidx,token[0],token[1],token[2]);
						vesicle->filament_list->poly[filidx]->vlist->vtx[fonoidx]->c=atof(token);
					}
					idx++;
					token=strtok(NULL," ");
				}
				xmlFree(values);		
			}
			xmlFree(property_name);
		}
		
		child=child->next;
	}
	return TS_SUCCESS;
}
/* this is a parser of vertex positions and bonds from main xml data */
ts_bool parseXMLVertexPosition(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur){
	xmlNodePtr child = cur->xmlChildrenNode;
	xmlChar *points;
	char *pts;
	int i, idx, polyidx, monoidx, filidx, fonoidx;
	char *token[3];
	while (child != NULL) {
		if ((!xmlStrcmp(child->name, (const xmlChar *)"DataArray"))){
			points = xmlNodeListGetString(doc, child->xmlChildrenNode, 1);
			pts=(char *)points;
			token[0]=strtok(pts," ");
			token[1]=strtok(NULL," ");
			token[2]=strtok(NULL,"\n");
			idx=0;
			while(token[0]!=NULL){
				if(idx<vesicle->vlist->n){
					vesicle->vlist->vtx[idx]->x=atof(token[0]);
					vesicle->vlist->vtx[idx]->y=atof(token[1]);
					vesicle->vlist->vtx[idx]->z=atof(token[2]);
				} else if(vesicle->tape->nmono && vesicle->tape->npoly && idx<vesicle->vlist->n+vesicle->tape->nmono*vesicle->tape->npoly) {
					polyidx=(idx-vesicle->vlist->n)/vesicle->tape->nmono;
					monoidx=(idx-vesicle->vlist->n)%vesicle->tape->nmono;
					vesicle->poly_list->poly[polyidx]->vlist->vtx[monoidx]->x=atof(token[0]);
					vesicle->poly_list->poly[polyidx]->vlist->vtx[monoidx]->y=atof(token[1]);
					vesicle->poly_list->poly[polyidx]->vlist->vtx[monoidx]->z=atof(token[2]);
				} else {
					filidx=(idx-vesicle->vlist->n-vesicle->tape->nmono*vesicle->tape->npoly)/vesicle->tape->nfono;
					fonoidx=(idx-vesicle->vlist->n-vesicle->tape->nmono*vesicle->tape->npoly)%vesicle->tape->nfono;
					//fprintf(stderr,"filidx=%d, fonoidx=%d, coord=%s,%s,%s\n",filidx,fonoidx,token[0],token[1],token[2]);
					vesicle->filament_list->poly[filidx]->vlist->vtx[fonoidx]->x=atof(token[0]);
					vesicle->filament_list->poly[filidx]->vlist->vtx[fonoidx]->y=atof(token[1]);
					vesicle->filament_list->poly[filidx]->vlist->vtx[fonoidx]->z=atof(token[2]);
				}
				for(i=0;i<2;i++)	token[i]=strtok(NULL," ");	
				token[2]=strtok(NULL,"\n");
				idx++;
			}
			xmlFree(points);
		}
		child=child->next;
	}
	//fprintf(stderr,"Came here\n");
	//fprintf(stderr,"Vertices position j=%d\n",idx);	

	return TS_SUCCESS;
}
ts_bool parseXMLBonds(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur){
	xmlNodePtr child = cur->xmlChildrenNode;
	xmlChar *bonds, *conname;
	char *b;
	int idx, polyidx;
	char *token[2];
	int temp_cnt=0;
	while (child != NULL) {
		conname=xmlGetProp(child, (xmlChar *)"Name");
		if ((!xmlStrcmp(child->name, (const xmlChar *)"DataArray")) && !xmlStrcmp(conname, (const xmlChar *)"connectivity") ){
			bonds = xmlNodeListGetString(doc, child->xmlChildrenNode, 1);
			b=(char *)bonds;
			token[0]=strtok(b," ");
			token[1]=strtok(NULL,"\n");
			idx=0;
			while(token[0]!=NULL){
				if(idx<3*(vesicle->vlist->n-2)){
					bond_add(vesicle->blist, vesicle->vlist->vtx[atoi(token[0])], vesicle->vlist->vtx[atoi(token[1])]);
					//fprintf(stderr,"Bonds in vesicle count idx=%d\n",idx);
				}
				else {
				//find grafted vtx
					if(vesicle->tape->npoly && vesicle->tape->nmono && (vesicle->tape->nmono-1)==(idx-3*(vesicle->vlist->n-2))%(vesicle->tape->nmono)
						&& idx<(3*vesicle->vlist->n-2+vesicle->tape->nmono*vesicle->tape->npoly)){
						temp_cnt++;
						//fprintf(stderr,"%d: Bonds in poly count idx=%d, t1=%s t2=%s\n",temp_cnt,idx, token[0], token[1]);
						
						polyidx=(idx-3*(vesicle->vlist->n-2))/(vesicle->tape->nmono);
						//fprintf(stderr,"poly=%d, vertex=%d\n",polyidx,atoi(token[0]));
						vesicle->poly_list->poly[polyidx]->grafted_vtx=vesicle->vlist->vtx[atoi(token[0])];
						vesicle->vlist->vtx[atoi(token[0])]->grafted_poly=vesicle->poly_list->poly[polyidx];
					}
				}
				token[0]=strtok(NULL," ");	
				token[1]=strtok(NULL,"\n");	
				idx++;
			}
			xmlFree(bonds);
		}
		xmlFree(conname);
		child=child->next;
	}
	//fprintf(stderr,"Bond data j=%d\n",idx);	
	return TS_SUCCESS;
}




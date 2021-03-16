/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _H_SNAPSHOT
#define _H_SNAPSHOT



typedef struct{
	char *string;
	ts_uint beg;
} ts_string;

//#define COMPRESSION

ts_bool xml_trisurf_data(FILE *fh, ts_vesicle *vesicle);
ts_bool xml_trisurf_header(FILE *fh, ts_vesicle *vesicle);
ts_bool xml_trisurf_footer(FILE *fh);
ts_bool xml_trisurf_tria(ts_string *data, ts_triangle_list *tlist);
ts_bool xml_trisurf_tria_neigh(ts_string *data, ts_triangle_list *tlist);
ts_bool xml_trisurf_vtx_neigh(ts_string *data, ts_vertex_list *vlist);
ts_bool xml_trisurf_vtx_tristar(ts_string *data, ts_vertex_list *vlist);
ts_bool xml_trisurf_nucleus(ts_string *data, ts_vesicle* vesicle);
ts_bool xml_trisurf_constvolarea(ts_string *data, ts_double volume, ts_double area);
/* UTILITIES */
char *base64_encode(const unsigned char *data, size_t input_length, size_t *output_length);
unsigned char *base64_decode(const char *data, size_t input_length, size_t *output_length);
void build_decoding_table();
void base64_cleanup();
ts_uint ts_compress_string64(char *data, ts_uint data_len, char **compressed);
#endif

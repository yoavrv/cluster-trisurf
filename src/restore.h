/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _H_RESTORE
#define _H_RESTORE

ts_vesicle *parseDump(char *dumpfname);
ts_vesicle *parseTrisurfTag(xmlDocPtr doc, xmlNodePtr cur);
ts_bool setGlobalTapeTXTfromTapeTag(xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfVtxn(ts_vertex_list *vlist, xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfTria(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfTriaNeigh(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfTristar(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseXMLPointData(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseXMLVertexPosition(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseXMLBonds(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfNucleus(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfConstantVolume(xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfConstantArea(xmlDocPtr doc, xmlNodePtr cur);
#endif

/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _BOND_H
#define _BOND_H


/** Initialize bond list with zero values
 *	@param *blist is a pointer to a ts_bond_list structure
 */
ts_bond_list *init_bond_list();

/** @brief Adds bond in the bond list
 *
 *  Function allocates space for *bond member of ts_bond_list. It then sets pointers to two vertices *vtx1 and *vtx2 that
 *  are members of the list.
 *	@param *blist is a pointer to initialized ts_bond_list.
 *	@param *vtx1 is a first vertex of a bond
 *	@param *vtx2 is a second vertex of a bond
 *	@returns TS_SUCCESS on success, TS_FAIL otherwise. If memory cannot be allocated
 *  this is considered as fatal error and execution stops, returning error code to the operating
 *  system.
 */
ts_bond *bond_add(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex *vtx2);

ts_bool bond_vector(ts_bond *bond);
ts_bool bond_list_free(ts_bond_list *blist);



#endif

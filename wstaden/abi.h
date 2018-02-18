#ifndef _seqIOABI_h
#define _seqIOABI_h


/* 
 * Title:       seqIOABI
 * 
 * File: 	 seqIOABI.h
 *Purpose:	 IO of ABI sequences
 * Last update: Mon May 28 1990
 */




/* ---- Imports ---- */


#include "Read.h"


/* ---- Exports ---- */


/*
 * Read the ABI format sequence with name `fn' into a Read structure.
 * All printing characters (as defined by ANSII C `isprint')
 * are accepted, but `N's are translated to `-'s. In this respect we
 * are adhering (more or less) to the CSET_DEFAULT uncertainty code set.
 * 
 * Returns:
 *   Read *	- Success, the Read structure read.
 *   NULLRead	- Failure.
 */
extern Read *read_abi(char *fn);
extern Read *fread_abi(FILE *fp);

#endif  /*_seqIOABI_h*/

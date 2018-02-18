#ifndef _seqIOALF_h
#define _seqIOALF_h


/* 
 * Title:       seqIOALF
 *
 * File: 	 seqIOALF.h
 * Purpose:	 IO of ABI sequences
 * Last update:  Fri 9th September 1994
 */




/* ---- Imports ---- */


#include "Read.h"


/* ---- Exports ---- */


/*
 * Read the ALF format sequence with name `fn' into a Read structure.
 * All printing characters (as defined by ANSII C `isprint')
 * are accepted, but `N's are translated to `-'s.
 *
 * Returns:
 *   Read *     - Success, the Read structure read.
 *   NULLRead   - Failure.
 */
extern Read *read_alf(char *fn);
extern Read *fread_alf(FILE *fp);

#endif  /*_seqIOALF_h*/

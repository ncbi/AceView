/*  File: querydisp.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public header for graphical query operations
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 25 13:48 1998 (fw)
 * Created: Wed Nov 25 13:46:20 1998 (fw)
 *-------------------------------------------------------------------
 */
/************************************************************/

#ifndef _QUERYDISP_H
#define _QUERYDISP_H



/* qbedisp.c */
void qbeCreate (void); 

/* querybuild.c */
void qbuildCreate (void) ;

/* BQL : 2017, new unified version of the acedb query language */
void bqlDisplayCreate (void) ;


#endif /* _QUERYDISP_H */

/************************ eof *******************************/

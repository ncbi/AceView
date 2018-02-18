/*  File: acedbgraph.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: This file describes the interface between a graphical
 *              ace application and the ace code responsible for
 *              initialising that interface. It hides away the grubby
 *              bits in setting up this interface.
 *              Currently there is not much in the interface but this
 *              will probably change as interactions between the ace
 *              code and the graph code change.
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 15 10:38 1998 (edgrif)
 * Created: Tue Oct 13 11:55:47 1998 (edgrif)
 *-------------------------------------------------------------------
 */
#ifndef DEF_ACEDBGRAPH_H
#define DEF_ACEDBGRAPH_H


/* Functions to initialise the acedb <-> graph interface.                    */
/*                                                                           */


/* This function combines graphInit and acedbGraphInit to provide a single   */
/* initialisation function for acedb applications using the graph package.   */
void acedbAppGraphInit(int *argcptr, char **argv) ;

/* For those requiring finer control they can use graphInit followed by a    */
/* call to this function.                                                    */
void acedbGraphInit(void) ;



#endif /* DEF_ACEDBGRAPH_H */

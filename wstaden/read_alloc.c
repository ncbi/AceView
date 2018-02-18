/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

/*****************************************************************/

/* 
 * File: 	read_alloc.c
 * Purpose:	Performs the allocation/freeing of Read structures
 * Last update: 01/09/94
 */


/*
    The Read data type is designed so that it can hold a varying degree
    of information about sequences, yet have a single set of calls
    to access the data.

    There are plenty of assumptions around that both the number of
    bases and the number of points will fit into an int_2, a short.

*/

/* ---- Includes ---- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "xalloc.h"
#include "misc.h"
#include "Read.h"

/*
 * Allocate a new sequence, with the given sizes.
 * Returns:
 *   "Read *" for success
 *   "NULLRead" for failure
 */
Read *read_allocate(int num_points, int num_bases) {
    Read *seq = NULLRead;


    /* Allocate the body of the sequence */
    if ((seq = (Read *)xmalloc(sizeof(Read))) == NULL)
	return(NULLRead);

    seq->NPoints = num_points;
    seq->NBases  = num_bases;

    /*   
     * Initialise the body, all pointers are set to NULL so we can
     * happily call `read_deallocate()`.
    */
    seq->leftCutoff  = 0;
    seq->rightCutoff = 0;
    seq->maxTraceVal = 0;

    seq->traceC    = NULL;
    seq->traceA    = NULL;
    seq->traceG    = NULL;
    seq->traceT    = NULL;

    seq->base      = NULL;
    seq->basePos   = NULL;

    seq->info = NULL;
    seq->format = TT_ANY;
    seq->trace_name = NULL;

    seq->prob_A = NULL;
    seq->prob_C = NULL;
    seq->prob_G = NULL;
    seq->prob_T = NULL;

    seq->orig_trace_format = TT_ANY;
    seq->orig_trace = NULL;

    /* Allocate space for the bases - 1 extra for the ->base field so
     * that we can treat it as a NULL terminated string.
     */
    if (((seq->base	 = (char *)xmalloc(num_bases+1))   == NULL) ||
	((seq->basePos   = (uint_2 *)xcalloc(num_bases,2)) == NULL) ||
	((seq->prob_A    = (char *)xmalloc(num_bases))	   == NULL) ||
	((seq->prob_C    = (char *)xmalloc(num_bases))	   == NULL) ||
	((seq->prob_G    = (char *)xmalloc(num_bases))	   == NULL) ||
	((seq->prob_T    = (char *)xmalloc(num_bases))	   == NULL))
    {
	read_deallocate(seq);
	return NULLRead;
    }

    if (num_points) {
	if (((seq->traceC   =(TRACE *)xcalloc(num_points, 2))  == NULL)||
	    ((seq->traceA   =(TRACE *)xcalloc(num_points, 2))  == NULL)||
	    ((seq->traceG   =(TRACE *)xcalloc(num_points, 2))  == NULL)||
	    ((seq->traceT   =(TRACE *)xcalloc(num_points, 2))  == NULL)
	    )
	{
	    read_deallocate(seq);
	    return NULLRead;
	}
    } else {
	seq->traceA = NULL;
	seq->traceC = NULL;
	seq->traceG = NULL;
	seq->traceT = NULL;
    }
    
    return seq;
}


/*
 * Free memory allocated to a sequence by read_allocate().
 */
void read_deallocate(Read *read)
{
    if (read == NULLRead)
	return;

    if (read->traceC  != NULL)  xfree(read->traceC);
    if (read->traceA  != NULL)  xfree(read->traceA);
    if (read->traceG  != NULL)  xfree(read->traceG);
    if (read->traceT  != NULL)  xfree(read->traceT);

    if (read->base    != NULL)  xfree(read->base);
    if (read->basePos != NULL)  xfree(read->basePos);

    if (read->info    != NULL)  xfree(read->info);

    if (read->prob_A  != NULL)  xfree(read->prob_A);
    if (read->prob_C  != NULL)  xfree(read->prob_C);
    if (read->prob_G  != NULL)  xfree(read->prob_G);
    if (read->prob_T  != NULL)  xfree(read->prob_T);

    if (read->trace_name != NULL) xfree(read->trace_name);

    if (read->orig_trace != NULL) xfree(read->orig_trace);

    xfree(read);
}

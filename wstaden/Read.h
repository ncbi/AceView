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

/* 
   This file was recovered from the anonymous ftp site
     ftp.mrc-lmb.cambridge.ac.uk/pub/staden/io_lib-1.7.tar.gz
   i copy it into the acedb include directory to be able
   to find it from the acedb make file
   but did not modify it in any way apart from the present comment

   acedb codes using this file must be linked against staden/io_lib/lib*
   compiled for the correct machine

   Jean Thierry-Mieg, octobre 98
*/

#ifndef _Read_h_
#define _Read_h_

/* 
 * Title:	Read
 *
 * File: 	Read.h
 * Purpose:	Read data type
 * Last update:	June  14 1994
 */

/*
 * This module encodes the `Read' sequence data structure.
 *
 * A `Read' contains information about bases and traces which are laid
 * out along a single dimension of points. The number of points in a
 * paricular sequence is given by `getNPoints', and these are numbered
 * 0..getNPoints-1. At each point there are four trace readings, one
 * for each base.
 *
 * The number of bases is `getNBases' which are numbered 0..N-1. 
 * Bases are represented by `char's. Every base is located at a 
 * particular point.
 *
 * The behaviour of these routines is undefined if given NULLRead or
 * an undefined sequence.
 */

#include "os.h"
#include "scf.h"

/*
 *-----------------------------------------------------------------------------
 * Macros
 *-----------------------------------------------------------------------------
 */

#define NULLRead     ((Read *)NULL)

/* Trace file formats */
#define TT_ERR -1
#define TT_UNK 0
#define TT_SCF 1
#define TT_ABI 2
#define TT_ALF 3
#define TT_PLN 4
#define TT_EXP 5
#define TT_CTF 6
#define TT_ANY TT_UNK

/*
 *-----------------------------------------------------------------------------
 * Structures and typedefs
 *-----------------------------------------------------------------------------
 */

typedef uint_2 TRACE;        /* for trace heights */

typedef struct
{
    int		format;	     /* Trace file format */
    char       *trace_name;  /* Trace file name */

    int         NPoints;     /* No. of points of data */
    int         NBases;      /* No. of bases */

    /* Traces */
    TRACE      *traceA;      /* Array of length `NPoints' */
    TRACE      *traceC;      /* Array of length `NPoints' */
    TRACE      *traceG;      /* Array of length `NPoints' */
    TRACE      *traceT;      /* Array of length `NPoints' */
    TRACE       maxTraceVal; /* The maximal value in any trace */

    /* Bases */
    char       *base;        /* Array of length `NBases' */
    uint_2     *basePos;     /* Array of length `NBases' */

    /* Cutoffs */
    int         leftCutoff;  /* Number of unwanted bases */
    int         rightCutoff; /* Number of unwanted bases */

    /* Miscellaneous Sequence Information */
    char       *info;        /* misc seq info, eg comments */

    /* Probability information */
    char       *prob_A;      /* Array of length 'NBases' */
    char       *prob_C;      /* Array of length 'NBases' */
    char       *prob_G;      /* Array of length 'NBases' */
    char       *prob_T;      /* Array of length 'NBases' */

    /* The original input format data, or NULL if inapplicable */
    int orig_trace_format;
    void *orig_trace;
} Read;


/*
 *-----------------------------------------------------------------------------
 * Function prototypes
 *-----------------------------------------------------------------------------
 */


/* ----- Main I/O routines ----- */

/*
 * Read a sequence from a file "fn" of format "format". If "format" is 0
 * (TT_ANY), we automatically determine the correct format.
 *
 * Returns:
 *   Read *   for success
 *   NULLRead for failure
 */
Read *read_reading(char *fn, int format);
Read *fread_reading(FILE *fp, char *fn, int format);


/*
 * Write a sequence to a file "fn" of format "format". If "format" is 0,
 * we choose our favourite - SCF.
 *
 * Returns:
 *   0 for success
 *  -1 for failure
 */
int write_reading(char *fn, Read *read, int format);
int fwrite_reading(FILE *fp, Read *read, int format);


/* ----- Utility routines ----- */

/*
 * Allocate a new sequence, with the given sizes.
 * Returns:
 *   "Read *" for success
 *   "NULLRead" for failure
 */
Read *read_allocate(int num_points, int num_bases);


/*
 * Free memory allocated to a sequence by read_allocate().
 */
void read_deallocate(Read *read);

#ifdef ACEDB4 
#define freeSeq(_s) (read_deallocate(_s), (_s) = 0)
#define seqMax(_seq)  ((_seq)->NPoints )
#define seqMaxBase(_seq)  ((_seq)->NBases)
#endif

/*****************************************************************/

/* unix specific file deletion routine */

int remove_file(char *fn);

Read *read_abi(char *fn);
Read *fread_abi(FILE *fp);
int write_abi(char *fn, Read *read);
int fwrite_abi(FILE *fp, Read *read);

int write_alf(char *fn, Read *read);
int fwrite_alf(FILE *fp, Read *read);
Read *read_alf(char *fn);
Read *fread_alf(FILE *fp);

int write_pln(char *fn, Read *read);
int fwrite_pln(FILE *fp, Read *read);
Read *read_pln(char *fn);
Read *fread_pln(FILE *fp);

Read *read_ctf(char *fn);
Read *fread_ctf(FILE *fp);
int write_ctf(char *fn, Read *read);
int fwrite_ctf(FILE *fp, Read *read);

#ifndef ACEDB4
#include "translate.h"
#include "filecompress.h"
#endif

#endif /* _Read_h_ */

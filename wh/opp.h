/*  Last edited: May 24 19:25 1994 (mieg) */
#ifndef _opp_h
#define _opp_h

/*
  Title:       opp

  File:        opp.h
  Purpose:     Required for complimenting a sequence
  Last update: Tue Jan 15 1991

  15.01.90 SD  Taken from seqIOEdit.h

*/

/* $Id: opp.h,v 1.1.1.1 2002/07/19 20:23:20 sienkiew Exp $ */

#include "regular.h"
#include "seq.h"

/* ---- Exports ---- */

extern char opp[256]; /* complement of any given base */

extern void oppInitialize();

/* initializes the array which stores the complement 
   of any of the Staden nucleotides or ambiguity
   codes */


void complement_seq(Seq seq);

/* complement a sequence */

#endif  /*_opp_h*/






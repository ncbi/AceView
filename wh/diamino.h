 /*  File: diamino.h
 *  Author: Danielle et Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 2001
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
   Align cDNA
 * Exported functions:
 * HISTORY:
 * Created: Mars 2001 (mieg)
 *-------------------------------------------------------------------
 */
/* %W% %G% */

#ifndef DIAMINO_H_DEF
#define DIAMINO_H_DEF

typedef BOOL (*DCF) (void *vp, int i1, int i2) ;
typedef struct diamWordStruct 
{
  int p ;              /* word length */
  BitSet bb ;          /* bitset of length lMax flagging letters in this word */
  BOOL discard ;       /* word contained in longer word, ignored in final collection */
} DIAMWORD ;

/* diaminoCreate returns an Array of DIAMWORD */
Array diaminoCreate (void *vp, int lMax, DCF isCompatible, DCF isConnected, int diamMax) ; /* explanations in diamino.c header */
void diaminoDestroy (Array aa) ;  /* pass back the result for cleaning the memory */
void diaminoTest (void) ;

#endif

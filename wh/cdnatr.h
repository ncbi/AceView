 /*  File: cdnatr.h
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

#ifndef CDNATR_H_DEF
#define CDNATR_H_DEF

typedef struct s2mStruct {
  KEY cosmid, cosmid1 ;
  Array plainHits, geneHits, chain, shadows, contigs, gmrnas, linkPos, dnaD, dnaR ;  
  int bMax, type ; void **magic ; AC_HANDLE h ; 
} S2M ;

typedef struct shadowStruct { KEYSET clones, ignoredClones ; BitSet yes, no, ghost, suspect ; int bMax ; BOOL isUp ; int begin3p, end3p ; int a1, a2 ; } SH ;
typedef struct contigStruct { S2M *s2m ; SH sh ; KEYSET ks ; Array s2mContigs ; BOOL isUp ; KEY gene, cGroup ; int a1, a2, d1, d2; } SC ;
typedef struct smrnaStruct { KEY gene, cGroup, tr ; int a1, a2, bestDna, nGaps ; Array clones, estHits, hits, dnas, orfs ; } SMRNA ;

KEYSET  mrnaCreate (int type, KEY cosmid, KEY cosmid1, Array dna, Array dnaR, Array plainHits, Array geneHits, 
		    Array linkPos, KEYSET clipTops, KEYSET clipEnds, KEYSET clonesWithBadIntron) ;
BOOL cdnaReExtendHits (S2M *s2m, SC* sc, SMRNA *gmrna, KEYSET clipTops, KEYSET clipEnds) ;
#endif

/*  File: aligntools.h
 *  Author: Ulrich Sauvage (ulrich@kaa.crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Mar  6 17:12 1996 (srk)
 * Created: Thu Jun  2 13:22:15 1994 (ulrich)
 *-------------------------------------------------------------------
 */

/*  $Id: aligntools.h,v 1.3 2016/09/03 01:20:49 mieg Exp $ */

#ifndef DEFINE_ALIGNTOOLS_h
#define DEFINE_ALIGNTOOLS_h

#define OVERLAP_MIN 20 /* overlap minimum entre 2 dna successifs */
#define NB_OCTET 3 /* dodecameres sur 3 octets */
#define DODECA 12 /* taille des oligos */
#define GO_AFT_ERREUR 5
#define GO_AFT_SINGLE 4
#define GO_AFT_DOUBLE 7

unsigned int dnaComplement(unsigned int lambda, int noct) ;
Array monDnaGet (DEFCPT look, KEY contig, KEY key) ;
void monDnaInsert(DEFCPT look, KEY key, Array dna) ;
void monDnaForget(DEFCPT look, KEY key) ;
/* void monDnaRemove (DEFCPT look) ; */

BOOL alignToolsMakeOligo(Array dna, int debut, int fin, int taille, 
			 unsigned int *oligo, int *pos) ;
Array alignToolsFindShortOligo(Array dna, int xd, int xf, 
				      unsigned int win, int taille,
				      unsigned int olig) ;
void alignToolsMakeSuccesOligo (Array dna, Associator oligo, int x1, 
				int x2, int taille) ;
void alignToolsFindColorOligo (Array dna, Associator oligo, Array color,
			       int taille) ;
void alignToolsInit (void) ;
void alignToolsDestroy (DEFCPT look) ;
Array alignToolsMatch(Array longDna, int xmin, int xmax, int x1, int x2,
		      Array shortDna, int ymin, int ymax, int y1, int y2,
		      int *debut, int *fin, int *sens, int all, int taux, 
		      int nbessai, BOOL control) ;
/*void doForceAssembleSeg(KEY seg, KEY target, int let, int cut) ;*/
void doForceAssembleSeg (DEFCPT look, KEY seg, KEY target, int min, int max,
			 int let, int cut) ;
Array alignToolsCptErreur(Array longDna, Array shortDna, int *debutp, int *finp,
			  int *shdebutp, int *shfinp) ;
Array alignToolsMakeErreur(Array longDna, int *debutp, int *finp, int pols, 
			   Array shortDna, int shdeb, int shfin, int pos, int sens) ;
void localCptErreur(Array longDna, int xl1, int xl2, int pol,
		    Array shortDna, int xs1, int xs2, int pos, int sens,
		    int *nbN, int *start, int *stop, int *recouv, Array errArray) ;

void alignToolsCheckArray (Array units, int width, int decal, int sens) ;
BOOL baseCallTrackBout (Array dna1, Array dna2, int *x1p, int *x2p, int *sensp) ;

#endif

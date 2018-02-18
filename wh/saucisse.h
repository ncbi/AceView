/* $Id: saucisse.h,v 1.1.1.1 2002/07/19 20:23:20 sienkiew Exp $ */

#ifndef SAUCISSE_DEF
#define SAUCISSE_DEF

#ifndef _SAUCISSE_DEF /* private to saucisse.c */
typedef void* Saucisse ;
#endif 

Saucisse saucisseCreate (int dim, int wordLength, int isDna) ; /* 0 scf, 1:dna, 2:peptide */

void uSaucisseDestroy (Saucisse saucisse) ;
#define saucisseDestroy(_saucisse)  (uSaucisseDestroy(_saucisse),(_saucisse) = 0)

void saucisseShow (Saucisse saucisse, int threshHold) ;
void  saucisseFill (Saucisse s, Array a) ;

#endif


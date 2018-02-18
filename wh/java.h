/*  Last edited: Mar 13 21:20 1995 (mieg) */

/* $Id: java.h,v 1.2 2006/12/16 05:05:37 mieg Exp $ */

/* file wh/java.h */

#ifndef DEF_JAVA
#define DEF_JAVA

char* freejavaprotect (char* text) ; /* in freesubs.c */
BOOL javaDumpLongText (KEY key)  ; /* in lomgtext.c */
BOOL javaDumpKeySet (KEY key) ; /* in keysetdump.c */
BOOL javaDumpDNA (KEY key) ; /* in dnasubs.c */
BOOL javaDumpPeptide(KEY key) ; /* in peptide.c */

#endif

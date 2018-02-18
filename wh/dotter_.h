/*  Last edited: Sep  6 18:34 1995 (esr) */

/* $Id: dotter_.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */
/* 
  dotter_.h - Private include file for dotter
*/

#include "dotter.h"

#define AAID_MIN	1	/* Smallest letter value in binary alphabet */
#define AAID_MAX	24	/* Maximum letter value in binary alphabet */
#define AAID_CNT	24	/* Number of letters in alphabet */
#define AAID_NAR	(AAID_MAX+1)
#define AAID_IGNORE	(AAID_MAX+3)

#define NA 24 		        /* Not A residue */
#define EL (NA+1)		/* End of Line */
#define ES (NA+2)		/* End of Sequence */
#define IC AAID_IGNORE	/* Ignore this Character */

#define UNKNOWN_AA_CHR	'X'
#define STOP_AA_CHR	'*'
#define GAP_AA_CHR '-'

#define NR 23 		        /* Not A residue */
#if !defined(NAMESIZE)
#define NAMESIZE 10
#endif

char *translate(char *seq, char **code);
extern char *stdcode1[];        /* 1-letter amino acid translation code */

int winsizeFromlambdak(int mtx[24][24], int *tob, int abetsize, char *qseq, char *sseq, 
		       double *exp_res_score, double *Lambda);

void fatal(char *format, ...);

void argvAdd(int *argc, char ***argv, char *s);


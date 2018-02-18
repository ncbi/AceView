/*  Last edited: Nov 13 15:21 1998 (fw) */

/* $Id: blxview.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

/* BLXVIEW.H

   include file for BLXVIEW.C

   Other files: iupac.h, blxview.c, translate.c

   Erik Sonnhammer, 92-02-20
*/

#ifndef _BLXVIEW_H
#define _BLXVIEW_H


#define NAMESIZE    12
#define FULLNAMESIZE 30
#define INITDBSEQLEN 50000   /* Initial estimate of max database sequence length */

#define max(a,b)        (((a) > (b)) ? (a) : (b))
#define min(a,b)        (((a) < (b)) ? (a) : (b))

/* Remember to ALWAYS update blxview.c:mspcpy() when changing this !!! */
typedef struct _MSP
{
struct _MSP *next;
    int      score;
    int      id;
    char     frame[8];
    int      qstart; 
    int      qend;
    char     sname[FULLNAMESIZE+1];
    char    *desc;
    int      sstart; 
    int      send;
    char    *sseq;
    int      box ;
#ifdef ACEDB
    KEY      key;
#endif
} MSP;

Graph blxview (char *seq, char *seqname,
	       int dispStart, int offset, MSP *msp, char *opts) ;

extern char *stdcode1[];        /* 1-letter amino acid translation code */
extern int   aa_atob[];
extern int PAM120[23][23];

char *translate(char *seq, char **code);
char *revcomp(char *comp, char *seq);
void *compl(char *seq);

Graph blxreadhsp(FILE *seqfile, FILE *exblxfile, char *qname, 
int dispstart, int qoffset, char *opts, int *argc, char **argv);

#endif /* _BLXVIEW_H */

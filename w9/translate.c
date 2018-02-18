/*  Last edited: Feb 22 15:44 1995 (esr) */
/* SCCS: %A% $Revision: 1.2 $ %D%
 * translate.c - functions for translating nucleic acid sequence
 * created Tue Jan 12 11:27:29 1993, SRE
 * 
 */

/* $Id: translate.c,v 1.2 2007/03/26 22:47:48 mieg Exp $ */

#define IUPACSYMNUM 17
#include "iupac.h"
#include "mystdlib.h"
#include "regular.h" /* for messalloc */

/*#include "squid.h"*/

/* Function: Translate(char *seq, char **code)
 * 
 * Given a ptr to the start of a nucleic acid sequence,
 * and a genetic code, translate the sequence into
 * amino acid sequence.
 * 
 * code is an array of 65 strings, representing
 * the translations of the 64 codons, arranged
 * in order AAA, AAC, AAG, AAU, ..., UUA, UUC, UUG, UUU.
 * '*' or '***' is used to represent termination
 * codons, usually. The final string, code[64],
 * is the code for an ambiguous amino acid.
 *
 * Because of the way space is allocated for the amino
 * acid sequence, the amino acid strings cannot be
 * longer than 3 letters each. (I don't foresee using
 * anything but the single- and triple- letter codes.)
 * 
 * Returns a ptr to the translation string on success,
 * or NULL on failure.
 */
char *
translate(seq, code)
     char  *seq;
     char **code;
{
  int   codon;			/* index for codon         */
  char *aaseq;                  /* RETURN: the translation */
  char *aaptr;                  /* ptr into aaseq */
  int   i;
  
  if (seq == NULL) return NULL;
  if ((aaseq = (char *) messalloc (strlen(seq) + 1)) == NULL)
      return NULL;

  aaptr = aaseq;
  for (; *seq != '\0' && *(seq+1) != '\0' && *(seq+2) != '\0'; seq += 3)
    {
				/* calculate the lookup value for
				   this codon */
      codon = 0;
      for (i = 0; i < 3; i++)
	{
	  codon *= 4;
	  switch (*(seq + i)) {
	  case 'A': case 'a':             break;
	  case 'C': case 'c': codon += 1; break;
	  case 'G': case 'g': codon += 2; break;
	  case 'T': case 't': codon += 3; break;
	  case 'U': case 'u': codon += 3; break;
	  default: codon = 64; break;
	  }
	  if (codon == 64) break;
	}

      strcpy(aaptr, code[codon]);
      aaptr += strlen(code[codon]);
    }
  
  return aaseq;
}


/* revcomp.c
 * 
 * Reverse complement of a IUPAC character string
 * 
 */

char *
revcomp(comp, seq)
     char *comp;
     char *seq;
{
  long  bases;
  char *bckp, *fwdp;
  int   idx;
  long  pos;
  int   c;

  if (comp == NULL) return NULL;
  if (seq == NULL)  return NULL;
  bases = strlen(seq);

  fwdp = comp;
  bckp = seq + bases -1;
  for (pos = 0; pos < bases; pos++)
    {
      c = *bckp;
      c = ace_upper(c);
      for (idx = 0; c != iupac[idx].sym && idx < IUPACSYMNUM; idx++);
      if (idx > IUPACSYMNUM)
	{
	  *fwdp = '\0';
	  return NULL;
	}
      else
	*fwdp = iupac[idx].symcomp;
      if (islower(*bckp)) *fwdp = ace_lower(*fwdp);
      fwdp++;
      bckp--;
    }
  *fwdp = '\0';
  return comp;
}
  
/* compl.c
 * 
 * Just complement of a IUPAC character string
 * 
 * Note that it overwrites calling string!!!! (revcomp doesn't)
 */

void compl(char *seq)
{
  char *fwdp;
  int   idx;
  long  pos;
  int   c;

  if (seq == NULL)  return;

  fwdp = seq;
  for (pos = 0; pos < strlen(seq); pos++)
  {
      c = ace_upper(*fwdp);
      for (idx = 0; c != iupac[idx].sym && idx < IUPACSYMNUM; idx++);
      if (idx > IUPACSYMNUM) {
	  *fwdp = '\0';
	  return;
      }
      else c = iupac[idx].symcomp;

      if (islower(*fwdp)) *fwdp = ace_lower(c);
      else *fwdp = c;
      fwdp++;
  }
  *fwdp = '\0';
}
  

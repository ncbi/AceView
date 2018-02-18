#include <stdio.h>
#include <fcntl.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define ACEDB4

/*#include "regular.h"
#include "keyset.h"
#include "dna.h"
*/
#include "array.h"
#include "Read.h"

/*
 * #defines to modify this file to compile easily as part of the Staden
 * Package io_lib.
 */
#define freeSeq(_s) (read_deallocate(_s), (_s) = 0)
#define seqMax(_seq)  ((_seq)->NPoints )
#define seqMaxBase(_seq)  ((_seq)->NBases)

#define xmalloc malloc
#define xcalloc calloc
#define xrealloc realloc
#define xfree free

#define arrayCreate(s,t) ArrayCreate(sizeof(t),(s))
#define array(a,n,t) ARR(t,a,n)
#define arrayMax ArrayMax
#define arrayDestroy ArrayDestroy
#undef arrp
#define arrp(a,n,t) \
    &((t*)((a)->base))[n]
#define arrayReCreate(a,n,t) arrayCreate(n,t)
#define arrayExists(a) ((a)->base != NULL)

#define BOOL int
#define mysize_t size_t
#define FALSE 0
#define TRUE 1


#define MAGIC 523747
#define PREDICTIONMODE 3  /* predictor degree */
#define COMPRESSIONMODE 3 /* compressor version */


/**********************************************************/
/* create a code for the 125 most frequent words */
static void ctfCompress3Init (Array *aap, int **lp, int **mp, int *maxCodep)
{
  short *sp ; int  i, j, k, n ;
  static int lng[128], mark[128], maxCode = 0 ;
  static Array aa = 0 ; 


  *aap = aa ; *lp = lng ; *mp = mark ; *maxCodep = maxCode ;
  if (aa) return ;
  *aap = aa = arrayCreate (512, short) ;
  array (aa, 511, short) = 0 ; /* make room */
  sp = arrp (aa, 0, short) ;
  j = 0 ;

  i = 0 ;  /* empty word */
  mark[i] = j ; n  = 0 ; j += lng[i] ; 
  
  /* single values up to +- 8 */
  for (k = 1 ; i < 126 && k < 12 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = k ;
      lng [i] = 1 ; j++ ;
      i++ ; 
      mark[i] = j ; *sp++ = -k ;
      lng [i] = 1 ; j++ ;
    }
  /* double values up to 50 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = k ; *sp++ = 0 ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = -k ; *sp++ = 0 ;
      lng [i] = 2 ; j += 2 ;
    }
  /* double values up to 51 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = k ; *sp++ = 1 ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = -k ; *sp++ = 1 ;
      lng [i] = 2 ; j += 2 ;
    }
  /* double values up to 5-1 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = k ; *sp++ = -1 ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = -k ; *sp++ = -1 ;
      lng [i] = 2 ; j += 2 ;
    }
  /* double values up to 15 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = 1 ; *sp++ = k ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = 1 ; *sp++ = -k ;
      lng [i] = 2 ; j += 2 ;
    }
  /* double values up to -15 */
  for (k = 1 ; i < 126 && k < 6 ; k++)
    { 
      i++ ; 
      mark[i] = j ; *sp++ = -1 ; *sp++ = k ;
      lng [i] = 2 ; j += 2 ;
      i++ ; 
      mark[i] = j ; *sp++ = -1 ; *sp++ = -k ;
      lng [i] = 2 ; j += 2 ;
    }
  /* triple values up to 111 */

  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;

  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 3 ; j += 3 ;

  /* quadruple values up to 1111 */

  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;

  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = -1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = -1 ; *sp++ =- 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = -1 ;
  lng [i] = 4 ; j += 4 ;

  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ;
  lng [i] = 4 ; j += 4 ;

  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = -1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = 1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 1 ; *sp++ = -1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = 1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = -1 ; *sp++ = -1 ; *sp++  = -1 ;
  lng [i] = 4 ; j += 4 ;

  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = 1 ;
  lng [i] = 5 ; j += 5 ;
  i++ ; mark[i] = j ; *sp++ = 1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 5 ; j += 5 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = 1 ;
  lng [i] = 5 ; j += 5 ;
  i++ ; mark[i] = j ; *sp++ = -1 ; *sp++ = 0 ; *sp++ = 0 ; *sp++ = 0 ; *sp++  = -1 ;
  lng [i] = 5 ; j += 5 ;

  if (i >= 126) 
    { fprintf (stderr, "FATAL ERROR, ctfCompress3Init i=%d > 126", i) ;
      exit (1) ;
    }
  if (j > 511)  
    { 
      fprintf (stderr, "FATAL ERROR, ctfCompress3Init j=%d > 511", j) ; 
      exit (1) ;
    }
  *maxCodep = maxCode = i ;
}

/**********************************************************/

static Array ctfCompress3 (Array a)
{
  int n0, n1, n10, n2, n3, n4 ;
  int i = arrayMax (a), j = 0, n ;
  Array b = arrayCreate (3 *i, char) ; /* worst case */
  unsigned char *cp, *cp0 ;
  short *sp, *sp1, *wp, z ;
  int w, *lng, *mark, maxCode ;
  Array aa = 0 ; 
  BOOL debug = FALSE ;

  ctfCompress3Init (&aa, &lng, &mark, &maxCode) ;

  array (b, 3*i - 1 , unsigned char) = 0 ;  /* make room */
  cp = cp0 = arrp (b , 0, unsigned char) ;
  sp = arrp (a, 0, short) ;

  n0 = n1 = n10 = n2 = n3 = n4 = 0 ;
  while (i--)
    {
      z = *sp++ ;
      if (!z)  /* string of zeroes */
	{
	  j =1 ;
	  while (i > 0 && j < 126 && !(*sp)) { j++ ; sp++; i-- ; } ;
	  *cp++ =  j & 0x7f ;  /* bit 1 = 0 */
	  n0 += j ; n1++ ;
	  continue ;
	}
       /* search for code word */
      for (w = maxCode ; w > 1 ; w--)  /* w-- favors long code words */
	{
	  n = lng[w] ; wp = arrp (aa, mark[w], short) ; sp1 = sp - 1 ;
	  while (n-- && *wp++ == *sp1++) ;
	  if (n == -1) break ;
	}
      if (w > 1 && lng[w] < i) /* code word found */
	{
	  n2++ ; n10 += lng[w] ;  *cp++ = 0x80 | (w & 0x7f) ;
	  sp += lng[w] - 1 ; i -= lng[w] - 1 ;
	  if (lng[w] <= 0)
	    {
	      fprintf (stderr, "FATAL ERROR in ctfCompress3 bad coding lng[w]") ;
	      exit (1) ;
	    }
	}
      else if ( z < 128 && z > -129) /* transmit byte */
	{
	  j = z + 128 ; /* range 0 ... 255 */
	  *cp++ = 0x80 | 126 ;
	  *cp++ = j & 0xff ;
	  n3++ ;
	}
      else     /* transmit short */
	{
	  j = z ; 
	  *cp++ = 0x80 | 127 ;
	  *cp++ = (j >> 8) & 0xff ;
	  *cp++ = j & 0xff ;	
	  n4++ ;
	}
    }
  arrayMax(b) = cp - cp0 ;
  if (debug) 
    printf (" // compress3:\n//  %d zeros in %d bytes, %d values coded in %d byte, %d bytes, %d short. \n// Total %d char for %d shrt\n",
	  n0, n1, n10, n2, n3, n4, arrayMax(b), arrayMax(a)) ;
  if (arrayMax(a) != n0 + n10 + n3 + n4)
    { 
      fprintf (stderr, "FATAL ERROR in ctfCompress3, codind error in compress 3") ;
      exit (1) ;
    }
  return b ;
}
/*
compress3 : 
//found   10829 zeros in 1865 bytes, 16524 values coded in 9114 byte, 183 bytes, 0 short. 
// Total 11162 char for 27536 shrt
*/

/**********************************************************/

static Array ctfDecompress3 (int dataMax, int shMax,
			    unsigned char *cp)
{
  int i = dataMax, mode, arg, n ;
  unsigned char cc, cc1, cc2 ;
  short *sp, *spMax, *wp ;
  Array b = arrayCreate (shMax, short) ;
  int *lng, *mark, maxCode ;
  Array aa = 0 ; 
  int n0, n1, n10, n2, n3, n4 ;
  BOOL debug = FALSE ;

  ctfCompress3Init (&aa, &lng, &mark, &maxCode) ;

  array (b, shMax - 1, short) = 0 ;  /* make room */
  sp = arrp (b, 0, short) ;
  spMax = sp + shMax ;
  n0 = n1 = n10 = n2 = n3 = n4 = 0 ;

  while (i-- && sp < spMax)
    {
      cc = *cp++ ;
      mode = cc & 0x80 ; arg = cc & 0x7f ;
      switch (mode)
	{
	case 0: /* initial zero = string of zero */
	  if (arg <= 0) /* should not happen */
	    { 
	      fprintf (stderr,"bad decompress3") ; 
	      goto abort ;
	    }  
	  n1++ ; n0 += arg ;
	  while (arg-- && sp < spMax) *sp++ = 0 ; 
	  break ;
	case 0x80:
	  switch (arg)
	    {
	    case 127:   /* next 2 bytes is a short */
	      i -= 2 ;  /* I need 3 bytes to code a short */
	      cc1 = *cp++ ; cc2 = *cp++ ;
	      *sp++ = (cc1 << 8) | cc2 ;
	      n4++ ;
	      break ;
	    case 126:   /* next byte is a byte */
	      i-- ;     /* I need 2 bytes to code a char */
	      cc1 = *cp++ ;
	      *sp++ = cc1 - 128 ;
	      n3++ ;
	      break ;
	    default:   /* 7 bytes is a code */
	      n = lng[arg] ; 
	      n2++ ; n10 += n ;
	      wp = arrp (aa, mark[arg], short) ;
	      while (n-- && sp < spMax) *sp++ = *wp++ ;
	      break ;
	    }
	}
    }
  if (debug)
    printf (" // compress3:\n//found   %d zeros in %d bytes, %d values coded in %d byte, %d bytes, %d short. \n// Total %d char for %d shrt\n",
	    n0, n1, n10, n2, n3, n4, n1 + n2 + n3 + n4, n0 + n10 + n3 + n4) ;
  
  if (i != -1 || sp != spMax)
    goto abort ;
  return b ;

abort:
  arrayDestroy (b) ;
  return 0 ;
}

/**********************************************************/
/**********************************************************/

static Array ctfDecompress (int compressionMode, 
			    int dataMax, int traceMax, 
			    unsigned char **cpp)
{ 
  Array a = 0 ;

  switch (compressionMode)
    {
    case 3:
      a = ctfDecompress3 (dataMax, 4*traceMax, *cpp) ;
      break ;
    default:  /* unknown compression mode */
      break ;
    }

  *cpp += dataMax ;
  return a ;
}

/**********************************************************/

static Array ctfCompress (int compressionMode, Array a)
{
  switch (compressionMode)
    {
    case 3:
      return ctfCompress3 (a) ;
    default:
      fprintf (stderr,"FATAL ERROR in ctfCompress, Non existing compression mode") ;
      exit (1) ;
      return 0 ; /* for compiler happiness */
    }
}


#define EBASE 256
double entropy(unsigned char *data, int len) {
    double E[EBASE];
    double P[EBASE];
    double e;
    int i;
    
    for (i = 0; i < EBASE; i++)
	P[i] = 0;

    for (i = 0; i < len; i++) {
	P[data[i]]++;
    }

    for (i = 0; i < EBASE; i++) {
	if (P[i]) {
	    P[i] /= len;
	    E[i] = -(log(P[i])/log(EBASE));
	} else {
	    E[i] = 0;
	}
    }

    for (e = i = 0; i < len; i++)
	e += E[data[i]];

    return e;
}

double entropy2(signed short *data, int len) {
    return entropy((unsigned char *)data, len*2);
}

#define MAXBUF 500000

void delta1(signed short *t1, signed short *t2, int len, int level) {
    int i;
    int u1 = 0, u2 = 0, u3 = 0, z;

    switch (level) {
    case 1:
	for (i = 0; i < len; i++) {
	    z = u1;
	    u1 = t1[i];
	    t2[i] = t1[i] - z;
	}
	break;

    case 2:
	for (i = 0; i < len; i++) {
	    z = 2*u1 - u2;
	    u2 = u1;
	    u1 = t1[i];
	    t2[i] = t1[i] - z;
	}
	break;

    case 3:
	for (i = 0; i < len; i++) {
	    z = 3*u1 - 3*u2 + u3;
	    u3 = u2;
	    u2 = u1;
	    u1 = t1[i];
	    t2[i] = t1[i] - z;
	}
    }
}

int shrink_16to8(signed short *i16, char *i8, int len) {
    int i, j;
    for (i = j = 0; i < len; i++) {
	if (i16[i] >= -127 && i16[i] <= 127) {
	    i8[j++] = i16[i];
	} else {
	    i8[j++] = -128;
	    i8[j++] = i16[i] >> 8;
	    i8[j++] = i16[i] &0xff;
	}
    }

    return j;
}

int table[256][256]; /* Too big for a local unless we increase stack size */
int follow(unsigned char *in, char *out, int len) {
    char next[256];
    int i, j;
 
    /* Count di-freqs */
    memset(table, 0, 256*256*sizeof(int));
    for (i = 0; i < len-1; i++)
	table[in[i]][in[i+1]]++;

    /* Pick the most frequent next byte from the preceeding byte */
    for (i = 0; i < 256; i++) {
	int bestval, bestind;

	bestval = bestind = 0;
	for (j = 0; j < 256; j++) {
	    if (table[i][j] > bestval) {
		bestval = table[i][j];
		bestind = j;
	    }
	}
	next[i] = bestind;
    }

    /* Output 'next' array */
    for (i = j = 0; i < 256; i++, j++)
	out[j] = next[i];

    /* Output new 'in' as next['in'] */
    out[j++] = in[0];
    for (i = 1; i < len; i++, j++)
	out[j] = next[in[i-1]] - ((char *)in)[i];

    return j; /* New length */
}

int main(void) {
    Read *r;
    double ent[4], ent_all;
    signed short data2[MAXBUF];
    char data1[MAXBUF], data1b[MAXBUF];
    int len;
    signed short *trace[4];
    int i;

    if (NULL == (r = fread_reading(stdin, "-", TT_ANY)))
	return 1;

    trace[0] = (signed short *)(r->traceA);
    trace[1] = (signed short *)(r->traceC);
    trace[2] = (signed short *)(r->traceG);
    trace[3] = (signed short *)(r->traceT);

    /* Raw trace */
    for (ent_all = i = 0; i < 4; i++) {
	ent[i] = entropy((unsigned char *)trace[i], r->NPoints * 2);
	ent_all += ent[i];
    }
    printf("Entropy of raw data          = %6.1f (%6.1f + %6.1f + %6.1f + %6.1f)\n",
	   ent_all, ent[0], ent[1], ent[2], ent[3]);


    /* Delta level 3 */
    for (ent_all = i = 0; i < 4; i++) {
	delta1(trace[i], data2, r->NPoints, 3);
	ent[i] = entropy((unsigned char *)(data2), r->NPoints * 2);
	ent_all += ent[i];
    }
    printf("Entropy of delta3 data       = %6.1f (%6.1f + %6.1f + %6.1f + %6.1f)\n",
	   ent_all, ent[0], ent[1], ent[2], ent[3]);

    /* Delta level 3 + 16-to-8 */
    for (ent_all = i = 0; i < 4; i++) {
	delta1(trace[i], data2, r->NPoints, 3);
	len = shrink_16to8(data2, data1, r->NPoints);
	write(2, data1, len);
	ent[i] = entropy((unsigned char *)data1, len);
	ent_all += ent[i];
    }
    printf("Entropy of delta3+16to8 data = %6.1f (%6.1f + %6.1f + %6.1f + %6.1f)\n",
	   ent_all, ent[0], ent[1], ent[2], ent[3]);

    /* Delta level 2 + 16-to-8 + follow */
    for (ent_all = i = 0; i < 4; i++) {
	delta1(trace[i], data2, r->NPoints, 3);
	len = shrink_16to8(data2, data1, r->NPoints);
	len = follow((unsigned char *)data1, data1b, len);
	ent[i] = entropy((unsigned char *)data1b, len);
	ent_all += ent[i];
    }
    printf("Entropy of d3+16to8+fol data = %6.1f (%6.1f + %6.1f + %6.1f + %6.1f)\n",
	   ent_all, ent[0], ent[1], ent[2], ent[3]);

    /* Delta level 3 + CTF */
    for (ent_all = i = 0; i < 4; i++) {
	Array a1 = arrayCreate(r->NPoints, short);
	Array a2;

	delta1(trace[i], data2, r->NPoints, 3);

	array(a1, r->NPoints-1, short) = 0;
	memcpy(a1->base, data2, r->NPoints * 2);
	a2 = ctfCompress(COMPRESSIONMODE, a1);

	ent[i] = entropy((unsigned char *)a2->base, arrayMax(a2));
	ent_all += ent[i];
    }
    printf("Entropy of d3+ctf       data = %6.1f (%6.1f + %6.1f + %6.1f + %6.1f)\n",
	   ent_all, ent[0], ent[1], ent[2], ent[3]);

    /* Delta level 3 + 16to8 + follow + CTF */
    for (ent_all = i = 0; i < 4; i++) {
	Array a1 = arrayCreate(r->NPoints, short);
	Array a2;
	int j;

	delta1(trace[i], data2, r->NPoints, 3);
	len = shrink_16to8(data2, data1, r->NPoints);
	len = follow((unsigned char *)data1, data1b, len);

	for (j = 0; j < len; j++) {
	    array(a1, j, short) = data1b[j];
	}
	/* array(a1, r->NPoints-1, short) = 0; */
	/* memcpy(a1->base, data2, r->NPoints * 2); */
	a2 = ctfCompress(COMPRESSIONMODE, a1);

	ent[i] = entropy((unsigned char *)a2->base, arrayMax(a2));
	ent_all += ent[i];
    }
    printf("Entropy of d3+16to8+fol+ctf  = %6.1f (%6.1f + %6.1f + %6.1f + %6.1f)\n",
	   ent_all, ent[0], ent[1], ent[2], ent[3]);

    /* Delta level 3 + CTF + follow */
    for (ent_all = i = 0; i < 4; i++) {
	Array a1 = arrayCreate(r->NPoints, short);
	Array a2;

	delta1(trace[i], data2, r->NPoints, 3);
	array(a1, r->NPoints-1, short) = 0;
	memcpy(a1->base, data2, r->NPoints * 2);
	a2 = ctfCompress(COMPRESSIONMODE, a1);
	len = follow((unsigned char *)a2->base, data1, arrayMax(a2));

	ent[i] = entropy((unsigned char *)data1, len);
	ent_all += ent[i];
    }
    printf("Entropy of d3+ctf+fol        = %6.1f (%6.1f + %6.1f + %6.1f + %6.1f)\n",
	   ent_all, ent[0], ent[1], ent[2], ent[3]);

    return 0;
}

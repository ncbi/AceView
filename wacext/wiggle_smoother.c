/* 
 * wiggle__smoother.c
 * author: mieg@ncbi.nlm.nih,gov
 * created 2008_11_18
 *
 * to compile:
 *   gcc -O -o wiggle_smoother wiggle_smoother.c -lm
 * purpose:
 *   to smooth a genomic distribution with given sigma
 * usage:
 *   smooth -s 100000 < wiggle_file.txt >! smoothed_wiggle_file.txt
 * wiggle files are tab delimited with 4 columns:
 *   chrom  x1  x2 value
 *   chrom: an arbitrary text shorter than 1024 char
 *   x1,x2: integral coordinates of a bin
 *   value: a real number i.e.: 45, -1000, 1.3e7
 */

#include <stdio.h> /* protptypes for input/output functions */
#include <stdlib.h> /* protptypes for std functions */
#include <string.h> /* prototypes for string functions */

static double s2 = 0 ;
extern double exp(double) ;
static int verbose = 0 ;

/* z: original value
   s: sum value
   w: weight of the sum
*/ 
typedef struct wstruct {int x1, x2 ; double z, s, w ; } WW ;

/***********************************************************************/
/***********************************************************************/
/* scan the input file, always a problem in C
 * Normally one uses a function from a tool box but i wanted to 
 * write a standalone code
 */
static int smoothGetLine (char *chrom, WW *ww, int nLine)
{
  int nn = 0 ;
  char *cp, *cq, line[4096] ;
  WW *w = ww + (nLine & 1023) ;
 
 ici:
  memset (line, 0, 4096) ;
  fgets (line, 4096, stdin) ; /* get one line of data from stdin */
  memset (chrom, 0, 1024) ;
  if (feof(stdin))   /* end of input file */
    return 1 ;  /* so the last datapoint gets processed */    

  cp = line ; 
  if (cp && *cp)
    {  /* get the chrom */
      cq = cp ; while (*cq) { if  (*cq == '\n' || *cq == '\r') { *cq = 0 ; break ; } else cq++ ; }
      cp = cq = line ; while (*cq && *cq != ' ' && *cq != '\t') cq++ ; *cq = 0 ;
      if (! *cp)  /* jump emty lines */
	goto ici ;
      if (*cp)
	{
	  if (*cp == '#' || !strncmp (cp, "track", 5))    
	    goto ici  /* jump the comment lines */ ; 
	  nn++ ;
	  strncpy (chrom, cp, 1023) ;
	  /* get x1 */
	  cp = cq + 1 ; while (*cp == ' ' || *cp == '\t') cp++ ;
	  cq = cp ; while (*cq && *cq != ' ' && *cq != '\t') cq++ ; *cq = 0 ;
	  w->x1 = atoi (cp) ; nn++ ;
	  /* get x2 */
	  cp = cq + 1 ; 
	  if (cp && *cp)
	    {
	      cp = cq + 1 ; while (*cp == ' ' || *cp == '\t') cp++ ;
	      cq = cp ; while (*cq && *cq != ' ' && *cq != '\t') cq++ ; *cq = 0 ;
	      w->x2 = atoi (cp) ; nn++ ;
	      /* get z */
	      cp = cq + 1 ;
	      if (cp && *cp)
		{
		  cp = cq + 1 ; while (*cp == ' ' || *cp == '\t') cp++ ;
		  cq = cp ; while (*cq && *cq != ' ' && *cq != '\t') cq++ ; *cq = 0 ;
		  w->z = atof (cp) ; nn++ ;
		}
	    }
	}
    }
  w->s = w->w = 0 ; /* sum and weight of the sum at this position */
  if (nn < 4)
    {
      fprintf (stderr
	       , "line %d::%s\n invalid, should be: chromosome x1 x2 value: text int int float\n"
	       , nLine+1, line
	       ) ;
      exit (1) ;
    }
  return 1 ;
} /* smoothGetLine */

/***********************************************************************/
/* average central point n1, remaining inside available datapoint [nLeft,nMax[ */
static void  smoothAccumulate (WW *ww, int n1, int nLeft, int nMax)
{
  WW *w1, *w ;
  int n ;
  double dx, weight[nMax - nLeft], we, wTotal = 0 ;

  w1 = ww + (n1 & 1023) ; /* actual position in the rotating buffer */
  wTotal = 0 ;
  for (n = nLeft ; n < nMax ; n++)
    {
      w = ww + (n & 1023) ; /* actual position in the rotating buffer */

      dx = (w->x1 + w->x2 - w1->x1 - w1->x2)/2.0 ;
      we = weight[n - nLeft] = exp (- dx * dx /s2) ;
      wTotal += we ;
    }
  if (wTotal > 0)
    for (n = nLeft ; n < nMax ; n++)
      {
	w = ww + (n & 1023) ; /* actual position in the rotating buffer */
	
	we = weight[n - nLeft]/wTotal ;
	w->w += we ;         /* acumulate the weight */
	w->s += we * w1->z ;  /* accumulate the data */
	if (verbose)
	  fprintf (stderr, "Acc n1=%d nLeft=%d nMax=%d w->x1=%d w->w=%g w->s=%g we=%g w1->z=%g\n"
		   , n1, nLeft, nMax,  w->x1, w->w, w->s, we, w1->z) ;
      }
}  /* smoothAccumulate */

/***********************************************************************/
/* export every remaining point in this chromosome */
static void  smoothExport (char *chrom, WW *ww, int nleft, int nright) 
{
  WW *w ;
  int n ;

  /* Export the results for this region */
  for (n = nleft ; n < nright ; n++)
    {
      w = ww + (n & 1023) ; /* actual position in the rotating buffer */
      if (w->w > 0) 
	{
	  if (verbose)
	    printf ("%s\t%d\t%d\t%g\t%g\n", chrom, w->x1, w->x2, w->z, w->s/w->w) ;
	  else
	    printf ("%s\t%d\t%d\t%g\n", chrom, w->x1, w->x2, w->s/w->w) ;
	}
      w->s = w->w = 0 ;
    }
}  /* smoothExport */

/***********************************************************************/
/* Master control, 
 * scan one input line at a time
 * store it in a rotating buffer large enough to contain a whole Gaussian shape
 * spread out the data point at the center of the Gaussian
 * export the points at the left of the Gaussian 
 */

static int smoothControl (double sigma)
{
  char chromOld[1024], chromNew[1024] ; /* the old and new chromosome names */
  int n1, nn = 0 ;   /* the current input line */
  int nLeft = 0, nCenter = 0 ; /* the current gaussian shape */
  double gWidth = 5*sigma ;
  WW *w, *wLeft, *wCenter ;
  WW ww[1024] ;  /* the rotating buffer, 
		  * by addressing this table as ww[nn & 1023] we never need to move the data around
		  */
  s2 = 2*sigma*sigma ; /* precompute 2*sigma^2 */
  memset (chromOld, 0, 1024) ;
  memset (chromNew, 0, 1024) ;
  memset (ww,0,sizeof(ww)) ;
  wLeft = wCenter = ww ;    /* initialise on the origin of the rotating buffer */
  for (nn = 0 ; smoothGetLine (chromNew, ww, nn) ; nn++)  /* parse one line of data */
   {
     w = ww + (nn & 1023) ;
     if (strcmp (chromOld, chromNew))
       { 
	 if (*chromOld)
	   {
	     for (n1 = nCenter ; n1 < nn ; n1++)
	       smoothAccumulate (ww, n1, nLeft, nn) ; /* average every remaining point in this chromosome */
	     smoothExport (chromOld, ww, nLeft, nn) ; /* export every remaining point in this chromosome */
	     
	     nLeft = nCenter = nn ;   /* interesting region now starts on our current data point */
	     wLeft = ww + (nLeft & 1023) ;
	     wCenter = ww + (nCenter & 1023) ;
	   }
	 if (!*chromNew)  /* the input file is fully processed */
	   break ;
	 strcpy (chromOld, chromNew) ;
       }
     while (w->x1 > wCenter->x1 + gWidth) /* all positions needed to spread wCenter are available */
       {
	 smoothAccumulate (ww, nCenter, nLeft, nn) ; /* average nleft */
	 nCenter++ ; wCenter = ww + (nCenter & 1023) ;
       }
     while (wLeft->x1 < wCenter->x1 - gWidth) /* all data going to nLeft has been accumulated */
       {
	 smoothExport (chromOld, ww, nLeft, nLeft + 1) ; /* export nLeft */
	 nLeft++ ; wLeft = ww + (nLeft & 1023) ;
       }
   }
  return nn ; /* number of scanned data points */
} /* smooth */

/***********************************************************************/
/***********************************************************************/
/* This function is called if the program is not called correcly */
static void usage (char *error)
{
  if (error)
    fprintf (stderr, "%s\n", error) ;
  fprintf (stderr,
	   "author: mieg@ncbi.nlm.nih.gov\n"
	   "  no guarantee, no restriction whatsoever\n"
	   "purpose:\n"
	   "  to smooth a genomic distribution with given sigma\n"
	   "usage:\n"
	   "  wiggle_smoother -s 100000 [-v] [-h] < wiggle_file.txt >! smoothed_wiggle_file.txt\n"
	   "      -s smoothing sigma\n"
	   "      -v verbose\n"
	   "      -h print this help\n"
	   "wiggle files are tab delimited with 4 columns:\n"
	   "  chrom  x1  x2 value\n"
	   "  chrom: an arbitrary text shorter than 1024 char\n"
	   "  x1,x2: integral coordinates of a bin\n"
	   "  value: a real number i.e.: 45, -1000, 1.3e7\n"
	   ) ;
  exit (1) ;  /* erroneous completion of the program */
} /* usage */

/***********************************************************************/
/* main is called automatically at the beginning of the program */
int main (int argc, char **argv)
{
  int ii, nn = 0 ;
  double s = -1 ;

  if (argc != 3 && argc != 4)  /* number of arguments on the command line */
    usage (0) ;
  for (ii = 1 ; ii < argc ;)
    {
      if (! strcmp (argv[ii], "-h"))
	usage (0) ;
      if (! strcmp (argv[ii], "-s"))
	{
	  ii++ ;
	  if (*argv[ii] == '-')
	    usage ("Missing smoothing value after parameter -s") ;
	  s = atol (argv[ii]) ; ii++ ;
	  if (s < 0)
	    usage ("Please give a positive value for sigma") ;
	}
      else if (! strcmp (argv[ii], "-v"))
	{
	  verbose = 1 ; ii++ ;
	}
      else
	usage ("Bad argument, I expect -s") ;
    }
      
  if (s < 0)
    usage (0) ;
  nn = smoothControl (s) ;      /* smooth the input */

  fprintf (stderr, "// Smoothed %d data point at resolution sigma = %g\n", nn, s) ;
  return 0 ;  /* successful completion of the program */
} /* main */

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

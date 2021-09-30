/*
 * authors: Danielle and Jean Thierry-Mieg, NCBI, 
 * 15 Sept 2016
 * table_histo
 *   Input: A tab delimited table with some numerical tables
 *          and comment or caption lines starting with #
 *   Output: The same caption
 *      But for each numeric column, the column is replaced by its stats and its histo
 */
#define MALLOC_CHECK  
#define ARRAY_CHECK  

#include "../wac/ac.h"

typedef struct cellStruct {
  char type ;
  double z ;
  int k ;
} CC ;

typedef struct thStruct {
  ACEIN ai ;
  ACEOUT ao ; 
  AC_HANDLE h ;
  Array table ;
  Array stats ;
  Array isText ;
  BOOL gzi, gzo ;
  const char *inFileName ;
  const char *outFileName ;
} TH ;

typedef struct statsStruct { 
  BigArray histo ; 
  int nVal ;
  BOOL isPercentage ;
  const char *title ;
  double divisor ;
  double median, average, sigma, min, max, cumul, rangeMax ; 
} STATS ;


/*************************************************************************************/
/***************************** Actual work *******************************************/

static void thParse (TH *th)
{
  ACEIN ai = th->ai ;
  const char *ccp ;
  char cutter ;
  int col = 0, line = -1 ;
  Array isText, rows = 0 ;

  th->stats = arrayHandleCreate (128, STATS, th->h) ;
  isText = th->isText = arrayHandleCreate (32, int, th->h) ;
  th->table = arrayHandleCreate (32, Array, th->h) ;

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      col = 0 ; line++ ;
      ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (! ccp)
	continue ;
      if (*ccp == '#')
	{
	  col = 0 ;
	  aceOutf (th->ao, "#%s\t%s", ccp[1]=='#' ? "#" : "", ccp) ;
	  while (cutter == '\t')
	    {
	      STATS *s = arrayp (th->stats, ++col, STATS) ;
	      ccp = aceInWordCut (ai, "\t", &cutter) ;
	      aceOutf (th->ao, "\t%s", ccp ? ccp : "") ;
	      if (0 && ccp && *ccp == '%')
		s->isPercentage = TRUE ;
	      s->title = ccp ? strnew (ccp, th->h) : 0 ;
	    }
	  aceOut (th->ao, "\n") ;
	}
      else
	{
	  while (1)
	    {
	      rows = array (th->table, col, Array) ;
	      if (! rows)
		rows = array (th->table, col, Array) = arrayHandleCreate (1024, double, th->h) ;
	      if (ccp && *ccp && strcmp (ccp, "-") && strcmp (ccp, "-nan"))
		{
		  STATS *s = arrayp (th->stats, col, STATS) ;
		  char *end ;
		  int k ;
		  double z ;

		  z = strtod (ccp, &end) ;
		  if (end && *end) /* non numeric */
		    array (isText,  col, int) = 1 ;
		  else
		    array (rows, line, double) = z ;
		  k = z ;
		  if (0 && k < -10000) 
		    invokeDebugger () ;
		  s->nVal++ ;
		}  
	      if (cutter != '\t')
		break ;
	      ccp = aceInWordCut (ai, "\t", &cutter) ;
	      col++ ;
	    }
	}
    }
  return ;
} /* thParse */

/*************************************************************************************/

static void thMakeHistos (TH *th)
{
  int col, cols = arrayMax (th->table) ;
  double *zp ;
  Array isText = th->isText ;

  for (col = 0 ; col < cols ; col++)
    {
      Array aa = array (th->table, col, Array) ;
      STATS *s = arrayp (th->stats, col, STATS) ;
      BigArray histo = s->histo ;
      int row, rowZero, rows = aa ? arrayMax (aa) : 0 ;
      double dz, zmin = 0, zmax = 0, zcumul = 0, zcumul2 = 0 ;
      double divisor ;

      if (! rows || array (isText,  col, int))
	continue ;

      divisor = s->divisor = 1 ;
      histo = s->histo = bigArrayHandleCreate (101, double, th->h) ;

      while (1)
	{
	  arraySort (aa, doubleOrder) ;
	  doubleVariance (aa, &(s->median), &(s->average), &(s->sigma)) ;
	  for (row = 0, zp = arrp (aa, row, double) ; row < rows ; zp++, row++)
	    {
	      if (!row) zmin = zmax = *zp ;
	      if (zmin > *zp) zmin = *zp ;
	      if (zmax < *zp) zmax = *zp ;
	      zcumul += *zp ;
	    }
	  while (zcumul > s->divisor * (double)1000000.0)
	    s->divisor *= 1000 ;
	  if (s->divisor > divisor)
	    {
	      divisor = s->divisor ;
	      for (row = 0 ; row < arrayMax (aa) ; row++)
		arr (aa, row, double) /= divisor ;
	    }
	  else
	    break ;
	}
	
      s->min = zmin ; s->max = zmax ; s->cumul = zcumul ;
      for (rowZero = row = 0, dz = 0, zp = arrp (aa, row, double), zcumul2 = 0 ; row < rows ; zp++, row++)
	{
	  if (*zp > 0 && (!s->isPercentage || *zp > 100) && 100 * (row -rowZero) > 95 * (rows - rowZero))
	    {
	      if (!dz)
		dz = *zp ; 
	      if (1.2 * dz < *zp)
		break ;
	    }
	  if (*zp == 0) rowZero = row ;
	  zmax = *zp ;
	  zcumul2 += *zp ;
	}
      {
	float u = zmax ;
	u = utUpperRoundPart (u) ;
	s->rangeMax = u ;
      }
      /* rowMax = row ;  clip the tail which contains at most 5% of the objects */
      if (s->isPercentage && zmax < 100) s->rangeMax = zmax = 100 ;
      dz = s->rangeMax/100 ; if (dz <= 0) dz = 1 ;
      for (row = 0, zp = arrp (aa, row, double) ; row < rows ; zp++, row++)
	{
	  double z = *zp ; if (zmin < 0) z -= zmin ; 
	  {
	    long int k = z/dz ;
	    if (k < 0) continue ;
	    if (k > 100)
	      break ;
	    double ddz1 = z/dz - k ;
	    double ddz2 = 1 - ddz1 ;
	    bigArray (histo, k, double)+= ddz2 ;
	    bigArray (histo, k+1, double)+= ddz1 ;
	  }
	}
      for (; row < rows ; zp++, row++)
	bigArray (histo, 100, double)++ ;
       
    }
  return ;
} /* thMakeHistos */

/*************************************************************************************/

static void thExport (TH *th)
{
  ACEOUT ao = th->ao ;
  int  row, col ;
  int cols = arrayMax (th->table) ;

  for (row = -8 ; row <= 100 ; row++)
    {
      switch (row)
	{
	case -8:
	  aceOutf (ao, "scale") ;
	  break ;
	case -7:
	  aceOutf (ao, "\ncumul") ;
	  break ;
	case -6:
	  aceOutf (ao, "\nnVal") ;
	  break ;
	case -5:
	  aceOutf (ao, "\nmin") ;
	  break ;
	case -4:
	  aceOutf (ao, "\nmedian") ;
	  break ;
	case -3:
	  aceOutf (ao, "\naverage") ;
	  break ;
	case -2:
	  aceOutf (ao, "\nmax") ;
	  break ;
	case -1:
	  aceOutf (ao, "\nsd") ;
	  break ;
	default:
	  if (row == 0)
	    aceOutf (ao, "\nHistogram") ;
	  aceOutf (ao, "\n%d", row) ;
	}
		       
      for (col = 0 ; col < cols ; col++)
	{
	  STATS *s = arrayp (th->stats, col, STATS) ;
	  if (s->nVal < 1)
	    { aceOutf (ao, "\t") ; continue ; }
	  if (s->cumul == 0 && row != -6)
	    { aceOutf (ao, "\t") ; continue ; }
	  double divisor = s->divisor ;
	  switch (row)
	    {
	    case -8:
	      {
		float k = 1 ; 
		int kk = 0 ;
		char *unit[] = {"", " kilo"," millions"," giga"," tera"," peta", " hexa", " zillions", "","",""} ;
		while (1000 * k < divisor * s->rangeMax && kk < 7 )
		  { k *= 1000 ; kk++ ; }
		aceOutf (ao, "\t%s [%.0f-%.0f%s]"
			 , s->title ? s->title : ""
			 , divisor * s->min
			 , divisor * s->rangeMax/k
			 , unit[kk]
			 ) ;
	      }
	      break ;
	    case -7:
	      aceOutf (ao, "\t%.2f", divisor * s->cumul) ;
	      break ;
	    case -6:
	      aceOutf (ao, "\t%d", s->nVal) ;
	      break ;
	    case -5:
	      aceOutf (ao, "\t%.0f", divisor * s->min) ;
	      break ;
	    case -4:
	      aceOutf (ao, "\t%g", divisor * s->median) ;
	      break ;
	    case -3:
	      aceOutf (ao, "\t%g", divisor * s->average) ;
	      break ;
	    case -2:
	      aceOutf (ao, "\t%g", divisor * s->max) ;
	      break ;
	    case -1:
	      aceOutf (ao, "\t%g", divisor * s->sigma) ;
	      break ;
	    default:
	      if (s->histo)
		aceOutf (ao, "\t%.1f", bigArray (s->histo, row, double)) ;
	      else
		aceOutf (ao, "\t-") ;
	      break ;
	    }
	}
    }
  aceOutf (ao, "\n") ;
      
  return ;
} /* thExport */

/*************************************************************************************/

static void thRun (TH *th)
{
  thParse (th) ;
  thMakeHistos (th) ;
  thExport (th) ;
  return ;
} /* thRun */


/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// table_histo:\n" 
	    "//   Examples\n"
	    "//   table_histo -i table.gz -o foo_bar  \n" 
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2016, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "//    Construct the histo of all nuneric columns of a tab delimited table\n"
	    "//    lines staring with # are considered as comments and simply echoed\n"
	    "// Options:\n"
	    "//   -i filename \n"
	    "//       default: parse stdin, if the file is called *.gz, it will be gunzipped\n"
	    "//   -gzi\n"
	    "//       apply gunzip on the input file, useful for stdin or arbitrary input file names\n"
	    "//   -o output_file_prefix\n"
	    "//      all exported files will be named output_file.histos.txt\n" 
            "//   -gzo : gzip all output files\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try: table_histo -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  TH th ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&th, 0, sizeof (TH)) ;
  th.h = h ;
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;
  th.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  th.gzo = getCmdLineBool (&argc, argv, "-gzo") ;
  getCmdLineOption (&argc, argv, "-i", &th.inFileName) ; 
  getCmdLineOption (&argc, argv, "-o", &th.outFileName) ; 

  if (argc > 1) usage (messprintf ("Unknown parameters %s", argv[1])) ;

  th.ai = aceInCreate (th.inFileName, th.gzi, th.h) ;
  if (! th.ai)
    messcrash ("FATAL ERROR: cannot open input file %s\n", th.inFileName ? th.inFileName : "stdin") ;
  th.ao = aceOutCreate (th.outFileName,  ".histos.txt", th.gzo, th.h) ;
  if (! th.ao)
    messcrash ("FATAL ERROR: cannot open input file %s\n", th.outFileName ? th.outFileName : "stdout") ;

  thRun (&th) ;

  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

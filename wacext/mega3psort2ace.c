
#include "ac.h"
#include "vtxt.h"
#include "mytime.h"
#include "bitset.h"

#define fixPsortCoord(_cord)		 (_cord >= xPos ? (_cord+xLen) : _cord)


typedef struct baStruct BA ;
struct baStruct
{
  AC_HANDLE h ;
  ACEIN ai ;
  ACEOUT ao ;
  BOOL gzi, gzo ;
  const char *inFileName, *outFileName ;
}  ;

/* Nakai himself says it is unreliable and it hits 90% of the big human genes, so forget it, jan 2004 */

/* motifs we do not like 2011_-2_15, 
   dan ne veut plus chercher les motifs suivants which are in PFAM or are useless like LL dileucine
 "Gavel: prediction of cleavage sites for mitochondrial preseq</A>\0\0","<A HREF=\0\0",0,"Mitochondrial_cleavage_site",4,3},
 
 
 {"Zinc finger, C2H2 type, domain (PS00028):  *** found ***\0\0","\n\n\0<A HREF=\0\0",0,"Zinc_finger",0,2},
 {"Dileucine motif in the tail:</A> found\0\0","<A HREF=\0\0",0,"Dileucine",0,2},
 {"Ribosomal protein S4e signature (PS00528):  *** found ***\0\0","<A HREF=\0\0",0,"Ribosomal_protein",0,2},
 {"Leucine zipper pattern (PS00029):  *** found ***\0\0","\n\n\0<A HREF=\0\0",0,"Leucine_zipper",0,2},
 {"Sigma-54 interaction domain ATP-binding region A signature (PS00675):  *** found ***\0\0","<A HREF=\0\0",0,"ATP_binding",0,2},
 {"Actinin-type actin-binding motif:</A>\0\0","type 2\0<A HREF=\0\0","type 1: found\n","actin_binding_1",0,2},
 {"Actinin-type actin-binding motif:</A>\0\0","<A HREF=\0\0","type 2: found\n","actin_binding_2",0,2},
 {"RNA-binding motif:</A> found\0\0","<A HREF=\0\0",0,"RNA_binding",0,2},
 {"ER Membrane Retention Signals:</A>\0\0","<A HREF=\0\0","motif in the N-terminus:","ER_membrane",0,-1}, // was mistakenly written C-terminus
 {"ER retention motif in the C-terminus:\0\0","\n\0<A HREF=\0\0",0,"ER_retention",0,-1},

 motifs we still like
 
 peroxisomal targeting signal in the C-terminus: CHL
 NMYR</A>: N-myristoylation pattern : MGILLGL

 2nd peroxisomal targeting signal:  found
      TLPN at 35
 possible vacuolar targeting motif: found
      TLPN at 35
transport motif from cell surface to Golgi: found
      YQRL at 881

 Prenylation motif:</A>       CC motif near the C-terminus: CCFN
 Prenylation motif:</A>       CaaX motif in the C-terminus: CIVA

*/


/*************************************************************************************/

static void psortOtherDomains (BA *ba, char *seq, char *content)
{
  int n, ii, x, y ;
  int ln = strlen (seq) ;
  char *cp, *cq, *cr, *motif ;
  const char *anchor[] = { "2nd peroxisomal targeting signal:  found" , "2nd_peroximal_domain" 
			   , "possible vacuolar targeting motif: found", "vacuolar_domain"
			   , "transport motif from cell surface to Golgi:", "Golgi_transport_domain"
			   , 0, 0  /* the other ones are treated in dedicated code */ 
  } ;
  const char *anchor2[] = {  "peroxisomal targeting signal in the C-terminus", "peroxisomal_domain", "Fin"
			     , "N-myristoylation pattern :", "N_myristoylation_domain", "Debut"
			     , "Prenylation motif:</A>       CaaX motif in the C-terminus: ", "prenylation_domain", "Fin"
			     , "Prenylation motif:</A>       CC motif near the C-terminus: ", "prenylation_domain", "VersLaFin"
			     , 0, 0
  } ;
  
  for (ii = 0 ; anchor[ii] ; ii+= 2)
    {      
      /* 
	 <A HREF="/psort/helpwww2.html#pox">SKL2</A>: 2nd peroxisomal targeting signal:  found
	 KLDARGIQL at 812
      */
      
      cp = strstr (content, anchor[ii]) ;
      if (! cp) continue ;
      cp = strstr (cp, "\n") ;
      cr = strstr (cp + 1, "\n") ;
      for (cp = cp + 1 ; *cp == ' ' ; cp++) ;
      cq = strstr (cp, " at ") ;
      if (! cq || cq > cr) continue ;
      *cq = 0 ; motif = cp ;
      cp = cq + 4 ;
      n = sscanf (cp, "%d", &x) ;
      if (n == 1)
	{
	  y = x + strlen (motif) - 1 ;
	  aceOutf (ba->ao, "Domain %s Psort 1. %d %d 1 1 \"%s\"\n", anchor[ii+1],x, y,  motif) ;
	}
      *cq = ' ' ;
    }
  
  /* peroxisomal targeting signal in the C-terminus: CHL */
  for (ii = 0 ; anchor2[ii] ; ii+= 3)
    {
      cp = strstr (content, anchor2[ii]) ;
      if (cp)
	{
	  cp = strstr (cp, ": ") ; 
	  if (cp) 
	    {
	      cp += 2 ; motif = cp ;
	      if (strncmp (motif, "none", 4) )
		{
		  cr = strstr (cp, "\n") ; *cr = 0 ; x = y = 1 ;
		  if (*anchor2[ii+2] == 'F')
		    { x = ln - strlen (motif) + 1 ; y = ln ; }
		  if (*anchor2[ii+2] == 'V')
		    { 
		      x = ln - strlen (motif) + 1 ; y = ln ; 
		      if (ln > 20)
			{
			  cq = strstr (seq + ln - 20, motif) ;
			  if (cq) { x = cq - seq + 1 ; y = x +  strlen (motif)  - 1 ; }
			}
		    }
		  if (*anchor2[ii+2] == 'D')
		    { x = 1 ; y = strlen (motif) ; }
		  aceOutf (ba->ao, "Domain %s Psort 1. %d %d 1 1 \"%s\"\n",  anchor2[ii+1], x, y,  motif) ;	  
		  *cr = 0 ;
		}
	    }
	}
    }
} /* psortOtherDomains */

/*************************************************************************************/

static void psortTransMembrane (BA *ba, char *seq, char *content)
{
  char *cend = 0, ccend ; /* block and restore */
  char *cp, *cq, *cr, cc ;
  int n, x, y, ln ;
  
  ln = strlen (seq) ;
  cp = strstr (content, "Tentative number of TMS(s)") ;
  if (cp)
    {
      cend = strstr (cp, "<A") ;
      if (cend) { ccend = *cend ; *cend = 0 ; }
      while ((cp = strstr(cp, "Transmembrane")))
	{
	  n = x = y = 0 ;
	  n = sscanf (cp, "Transmembrane  %d - %d\n", &x, &y) ;
	  if (n == 2 && x >= 1 && y >= 1 && x <= ln && y <= ln) /* success */
	    {
	      cq = seq + x - 1 ;
	      cr = seq + y  ;
	      cc = *cr ; *cr = 0 ;
	      aceOutf (ba->ao, "Domain Transmembrane_domain Psort 1. %d %d 1 1 \"%s\"\n"
		      , x, y, cq
		      ) ;
	      *cr = cc ;
	    }
	  cp += 10 ; /* arbitrarily move a little forward */
	}
    }
  if (cend) *cend = ccend ;
} /* psortTransMembrane */

/*************************************************************************************/

static void psortNuclearLocalization (BA *ba, char *seq, char *content)
{
  char *cend = 0, ccend ; /* block and restore */
  char *cp, *cq, *cr, cc ;
  int n, i, x, y, xmax = 0, ln ;
  BitSet bb = bitSetCreate (1000, 0) ;    

  ln = strlen (seq) ;
  /*
    <A HREF="/psort/helpwww2.html#nuc">NUCDISC: discrimination of nuclear localization signals</A>
    pat4: RPRR (4) at    4
    pat4: PRRR (4) at   69
    pat4: RRRK (5) at   70
    pat4: RRKR (5) at   71
    pat4: RKRR (5) at   72
    pat4: KRRR (5) at   73
    pat4: RRRK (5) at   74
    pat4: RRKR (5) at   75
    pat7: PRRRKRR (5) at   69
    pat7: PVATRRR (3) at   91
    bipartite: RRTQISGPRRRKRRRKR at   62
    bipartite: RRPTSTWRWSPVATRRR at   81
    content of basic residues:  24.8%
    NLS Score:  4.86
  */
  
  cp = strstr (content, "NUCDISC: discrimination of nuclear localization signals") ;
  if (cp)
    {
      cend = strstr (cp, "<A") ;
      if (cend) { ccend = *cend ; *cend = 0 ; }
      while ((cq = strstr(cp, "pat4:")))
	{
	  cp = cq + 5 ; cq = strstr (cp, "at") ;
	  cr = strstr (cp, "\n") ; 
	  if (cq && cq < cr)
	    {
	      cq += 2 ; while (*cq== ' ') cq++ ;
	      n = sscanf (cq, "%d", &x) ;
	      if (n == 1)
		{
		  for (i = 0 ; i < 4 ; i++)
		    bitSet (bb, x + i) ;
		  if (x + i > xmax) xmax = x + i ;
		}
	    }
	  cp = cr + 1 ;
	}
      while ((cq = strstr(cp, "pat7:")))
	{
	  cp = cq + 4 ; cq = strstr (cp, "at") ;
	  cr = strstr (cp, "\n") ; 
	  if (cq && cq < cr)
	    {
	      cq += 2 ; while (*cq== ' ') cq++ ;
	      n = sscanf (cq, "%d", &x) ;
	      if (n == 1)
		{
		  for (i = 0 ; i < 7 ; i++)
		    bitSet (bb, x + i) ;
		  if (x + i > xmax) xmax = x + i ;
		}
	    }
	  cp = cr + 1 ;
	}
      
      while ((cq = strstr(cp, "bipartite:")))
	{
	  cp = cq + 11;   cq = strstr (cp, "at") ;
	  cr = strstr (cp, "\n") ; 
	  if (cq && cq < cr)
	    {
	      y = cq - cp - 1 ;
	      cq += 2 ; while (*cq== ' ') cq++ ;
	      n = sscanf (cq, "%d", &x) ;
	      if (n == 1)
		{
		  for (i = 0 ; i < y ; i++)
		    bitSet (bb, x + i) ;
		  if (x + i > xmax) xmax = x + i ;
		}
	    } 
	}
    }
  /* export the merged segments */
  bitUnSet (bb, xmax) ;
  for (i = x = y = 0 ; i <= xmax ; i++)
    {
      if (bit (bb, i))
	{
	  if (x) y++ ;
	  else x = y = i ;
	}
      else
	{
	  if (x)
	    {  /* export */
	      cq = "" ; cr = 0 ;
	      if (x >= 1 && y >= 1 && x <= ln && y <= ln) /* success */
		{
		  cq = seq + x - 1 ;
		  cr = seq + y ;
		  cc = *cr ; *cr = 0 ;
		}
	      aceOutf (ba->ao, "Domain Nuclear_localization_domain Psort 1. %d %d 1 1 \"%s\"\n"
		       , x, y
		       , cq
		       ) ;
	      if (cr) *cr = cc ;
	      x = y = 0 ;
	    }
	}
    }
  bitSetDestroy (bb) ;
  if (cend) *cend = ccend ;
} /* psortNuclearLocalization */

/*************************************************************************************/

static void psortCellularLocalization (BA *ba, char *seq, char *content)
{
  char *cp, *cq, *cr, *cs ;
  int n ;
  float score = 0 ;

  cp = strstr (content, "<H1>Results of the <i>k</i>-NN Prediction</H1></A>") ;
  if (cp) cp = strstr (cp, "<PRE>\n") ;
  if (! cp) return ;
  cr = strstr (cp, ">>") ;
  cp = strstr (cp, "\n") ; cp++ ;
  while (1)
    {
      if (!cp || cp > cr) break ;
      for (++cp ; *cp == ' ' ; cp++) ;
      n = sscanf (cp, "%f %%:", &score) ;
      if (n != 1) break ;
      cp = strstr (cp, ":") ;
      if (! cp) break ;
      cp += 2 ;
      cq = strstr (cp, "\n") ;
      if (! cq) break ;
      *cq = 0 ;
      *cp = ace_upper (*cp) ;
      for (cs = cp ; *cs ; cs++)
	if (*cs == ' ') *cs = '_' ;
      if (strstr (cp, "xtracellular"))
	cp = "Secreted" ;
      aceOutf (ba->ao,"Psort Localization \"%s\" %4.1f \n", cp, score) ;
      *cq = '\n' ;
      cp = cq + 1 ;
    }
} /* psortCellularLocalization */
 
/*************************************************************************************/

static void psortCoiledCoil (BA *ba, char *seq, char *content)
{
  char *cp, *cq, *cr, cc = 0 ;
  int n, x, y ;
  int ln = strlen (seq) ;  

  cp = strstr (content, "algorithm to detect coiled-coil regions") ;
  if (cp) cp = strstr (cp, "</A>") ;
  if (! cp) return ;
  for (cp += 4 ; *cp == ' ' ; cp++) ;
  n = sscanf (cp, "%d", &x) ;    /* zero based offset */
  if (n != 1) return ;
  
  cp = strstr (cp, "total") ;
  if (! cp) return ;
  n = sscanf (cp, "total: %d residue", &y) ; /* length */
  if (n != 1) return ;
  
  if (x >= 0 && x+y >= 0 && x < ln && x+y <= ln) /* success */
    {
      /* in this output, the coordinates are zero based */
      cq = seq + x  ;
      cr = seq + x + y  ;
      cc = *cr ; *cr = 0 ;
      
      aceOutf (ba->ao, "Domain Coiled_coil_region Psort 1. %d %d 1 1 \"%s\"\n", x, x+y-1, cq) ;
      if (cc) *cr = cc ;
    }
} /* psortCoiledCoil */

/*************************************************************************************/

static void psortCleavableSignalPeptide (BA *ba, char *seq, char *content)
{
  char *cp, *cq, *cr, cc = 0 ;
  int n = 0, x, y ;
  int ln = strlen (seq) ;  

  cp = strstr (content, "Seems to have a cleavable signal peptide") ;
  if (cp) n = sscanf (cp , "Seems to have a cleavable signal peptide (%d to %d)", &x, &y) ;

  if (n == 2 && x >= 1 && y >= 1 && x <= ln && y <= ln) /* success */
    {
      /* in this output, the coordinates are zero based */
      cq = seq + x  - 1 ;
      cr = seq + y  ;
      cc = *cr ; *cr = 0 ;
      aceOutf (ba->ao, "Domain Cleavable_signal_peptide Psort 1. %d %d 1 1 \"%s\"\n", x, y, cq) ;
      if (cc) *cr = cc ;
    }
} /* psortCleavableSignalPeptide */

/*************************************************************************************/
/* break the file per protein */
static void psort2aceAnalyse (BA *ba, char *seq, char *content)
{
  if (1)
    {
      psortNuclearLocalization (ba, seq, content) ;
      psortOtherDomains (ba, seq, content) ;
      psortCleavableSignalPeptide (ba, seq, content) ;
      psortTransMembrane (ba, seq, content) ;
      psortCoiledCoil (ba, seq, content) ;
      psortCellularLocalization (ba, seq, content) ;
    }
} /* psort2aceAnalyse */

/*************************************************************************************/
/* break the file per protein */
static void psort2aceOne (BA *ba, char *content)
{
  ACEOUT ao = ba->ao ;
  char *cp, *nm = 0, *seq = 0 ;

  /* locate and export the identifier */
  cp = strstr (content, "<PRE>\n") ;
  if (cp) { cp += 6 ; nm = cp ; cp = strstr(cp, " ") ;}
  if (cp) { *cp++ = 0 ; }
  if (! nm || !cp)
    return ;
  aceOutf (ao, "\nKantor %s\n", nm) ;
  aceOutf (ao,"-D PSR\n") ;
  aceOutf (ao, "Psort_Date %s\n",timeShowNow()) ;
  aceOutf (ao, "Kantor_Date %s\n",timeShowNow());

  /* separate the sequence from the content */
  if (cp) cp = strstr(cp, "\n") ;
  if (cp) cp = strstr(cp+1, "\n") ;
  if (cp) 
    {
      seq = cp+1 ;
      content = strstr (seq, "\n") ;
      *content++ = 0 ;
      
      psort2aceAnalyse (ba, seq, content) ;
    }
} /* psort2aceOne */

/*************************************************************************************/
/* break the file per protein */
static int psort2ace (BA *ba)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = ba->ai ;
  int n = 0, nn = 0 ;
  char *ccp ;

  vTXT txt = vtxtHandleCreate (h) ;


  while (aceInCard (ai))
    {
      ccp = aceInPos (ai) ;
      if (! ccp || *ccp == '#') continue ;
      
      if (strcasecmp (ccp, "<H1>Input Sequence</H1>"))
	{ 
	  n++ ;
	  vtxtPrintf (txt, "%s\n", ccp) ;
	  continue ;
	}
      if (n)
	{
	  psort2aceOne (ba, vtxtPtr (txt)) ;
	  nn++ ; n = 0 ;
	}
      vtxtClear (txt) ;
    }

  if (n)
    {
      psort2aceOne (ba, vtxtPtr (txt)) ;
      nn++ ; n = 0 ;
    }
  
  ac_free (h) ;
  return nn ;
} /* psort2ace */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  fprintf  (stderr,
	    "// bestali.c:\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, Feb 2011, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "//   Analyse the html output of Nakai PSORT program\n"
	    "// Input:\n"
	    "//   -i file_name [-gzi] : default stdin\n"
	    "//      if the file is called .gz or if the option -gzi is specified, invokes gunzip\n"
	    "//   -o out_file_name [-gzo] : default stdout\n"
	    "//      redirect the output\n"
	    "//      if -gzo is specified, invokes gzip and add a .gz suffix to the file name\n"
	    "// Help\n"
	    "//   -help : export this on line help\n"
	    "// Caveat:\n"
	    "//   The program depends on the exact layout of the input file as of 2011_02_15\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  BA ba ;
  AC_HANDLE h = ac_new_handle () ;

  memset (&ba, 0, sizeof (BA)) ;
  ba.h = h ;

  if (argc == 0) /* for this program no argument is ok */
    usage (0) ;

  /* optional arguments */

  getCmdLineOption (&argc, argv, "-i", &(ba.inFileName)) ;
  getCmdLineOption (&argc, argv, "-o", &(ba.outFileName)) ;
  ba.gzi = getCmdLineOption (&argc, argv, "-gzi", 0) ;
  ba.gzo = getCmdLineOption (&argc, argv, "-gzo", 0) ;

  if (getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  if (argc > 1)
    {
      fprintf (stderr, "Unknown argument %s, try -help", argv[argc-1]) ;
      exit (1) ;
    }

  /* actions */

  fprintf (stderr, "Start %s\n", timeShowNow ()) ;

  ba.ai = aceInCreate (ba.inFileName, ba.gzi, h) ;
  ba.ao = aceOutCreate (ba.outFileName, 0, ba.gzi, h) ;

  if (1) /* default */
    {
      psort2ace (&ba) ;
    }
     
  /* clean up */
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;
  
  /* report the global counters */

  fprintf (stderr, "Done %s\n", timeShowNow ()) ;
  ac_free (h) ;

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/


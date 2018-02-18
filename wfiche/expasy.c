/*
  file expasy.c
 
  This program implements expasy.   Really, it just calls up somebody's
  web server and asks for the answer.
 
  Feed this program a FASTA file on stdin, get a .ace file on stdout

  If you look on www.expasy.org, you see there are several mirror
  sites.  I have good results with www.expasy.org, ca.expasy.org,
  and us.expasy.org.  The others are not so well connected from
  our site.
 
  expasyUrl = "www.expasy.org";
  expasyUrl = "ca.expasy.org";
  expasyUrl = "us.expasy.org";
*/

#define expasyUrl "ca.expasy.org"

/*
  the expasy web site wants a protein sequence as a parameter.
  The result looks like
     Theoretical pI/Mw: 6.94 / 23528.66
 
  In the acedb database, you want
 	Kantor name
 	Expasy f1 f2
  where f1= Mw / 1000 rounded to the nearest .1
  and f2 = pI
 
*/

/*************************************************************************/
#include "regular.h"
#include "mytime.h"
#include "vtxt.h"

/* receives a peptide name and sequence, printf the result in .ace format */ 
static void getPage (char *kantor, char *pep)
{
  char *cp, *page ;
  vTXT txt = 0 ;
  float pI, mw ;
  int mm ;
  char fullUrl [300] ;

  /* initialise the acedb object */
  printf ("Kantor : \"%s\"\n", kantor) ;
  printf ("-D Expasy\n") ;

  { /* format the date and always export to prevent looping */
    char *cq, buf [128] ;
    
    strcpy (buf, timeShowNow ()) ;
    if ((cq = strstr (buf, "_")))
      *cq = 0 ;
    printf ("Expasy_date %s\n", buf) ;
  }
  
  if ((cp = strstr (pep, "XXX")))
    {
      printf ("// XXX in sequence %s\n\n", cp ? cp : " pep > 150 ") ;
      return ;
    } 
  
  
  sprintf (fullUrl, "http://%s/cgi-bin/pi_tool?protein=%s", expasyUrl, pep) ;
  txt = vtxtGetURL (fullUrl, 0) ;
  if (!txt) return ;
  page = vtxtPtr (txt) ;
  	  
  cp = page ? strstr (page, "<H2>Theoretical pI/Mw: ") : 0 ;
  if (!cp && 
      (cp = strstr (page, "Its pI and Mw cannot be computed.")))
    {
      /* It is normal that some peptides do not have a definite answer.
       * For those, we do not want to have the pI/Mw numbers, but we
       * do want to know that we have processed this peptide.
       */
      printf ("// %s not computable\n\n", kantor) ;
    } 
  else if ((cp = strstr (page, "<H2>Theoretical pI/Mw: ")) &&
	   (cp += 23) &&
	   sscanf (cp, "%f / %f", &pI, &mw) == 2)
    { /*  sucess we expect: Theoretical pI/Mw: 6.94 / 23528.66 */
      mm = mw + 499 ;
      mw = mm/1000 ;
      printf ("Expasy %6.1f %6.2f\n\n", mw, pI) ;
    }
  else
    {
      /* This means that we got a bogus response from the web site.  All
       * we really know is that we do not know the correct answer.  By
       * blanking the data fields, we demand manual attention for this peptide.
       */
      printf("// %s has unrecognized response\n\n", kantor) ;
    }  
  vtxtDestroy (txt) ;
}

/*************************************************************************/

int main (int argc, char **argv)
{
  int level ;
  char *cp, kantor [1000] ;
  Stack s = stackCreate (10000) ;

  freeinit () ;

  level = freesetfile (stdin, 0) ;
  while (freecard (level))
    {
      cp = freeword () ;
      if (cp && *cp == '>') /* begin of a trace */
	{
	  if (stackMark (s))
	    getPage (kantor, stackText (s, 0)) ;
	  strcpy (kantor, cp + 1) ;
	  stackClear (s) ;
	}
      else if (cp)
	catText (s, cp) ;
    }
  
  if (stackMark (s))
    getPage (kantor, stackText (s, 0)) ;

  printf ("\n\n// Done \n\n") ;
  return 0 ;
}

/*************************************************************************/
/*************************************************************************/


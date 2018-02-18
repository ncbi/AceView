#define ARRAY_CHECK

#include "../wac/ac.h"
#include "vtxt.h"
#include "vtxt_.h"
#include "utils.h"

/*************************************************************/

static BOOL omimParseOne (char* omim, AC_HANDLE h)
{
  vTXT txt = xmlGetDDContent (omim, h) ;
  char *cp, *cq, *cr, *xml ;
  BOOL isDisease ;

  if (txt)
    {
      xml = vtxtPtr (txt) ;
      isDisease = FALSE ;
      cp = xmlGetTagContent (&xml, "Mim-entry_mimNumber", 0, h) ;
      if (cp)
	{
	  printf ("Extern \"OMIM_%s\"\nOMIM\n", cp) ;
	  cp = xml ;
	  if (1)
	    while ((cq = xmlGetNextTagContent (&cp, "Mim-entry_mimType", h)))
	      {
		/* == 2 are moved to, to be treated by hand */
		if (*cq == '2') { isDisease = TRUE ; printf ("OMIM_disease\n") ; }
		else if (*cq == '3') { isDisease = TRUE ; printf ("OMIM_disease\n") ; }
		else if (*cq == '4')  printf ("OMIM_molecular\n") ;
		else if (*cq == '5') { isDisease = TRUE ; printf ("OMIM_disease\n") ; }
		else if (*cq == '1') printf ("OMIM_molecular\n") ;
		if (0) printf ("OMIM_unkown %s\n", cq) ;
	      }
	  cp = xml ;
	  if (1)
	    while ((cq = xmlGetNextTagContent (&cp, "Mim-entry_title", h)))
	      {
		cr = strstr (cq, ";") ;
		if (cr) *cr = 0 ;
		printf ("%s %s\n"
			, "OMIM_title"
			, freeprotect (vtextLowerCase(cq))
			) ;
		if (cr) *cr = ';' ;
	      }
	  cp = xml ;
	  if (1) 
	    while ((cq = xmlGetNextTagContent (&cp, "Mim-entry_aliases_E", h)))
	      {
		printf ("%s %s\n"
			, isDisease ? "OMIM_alias" : "Properties"
			, freeprotect (vtextLowerCase(cq))
			) ;
	      }
	  cp = xml ;
	  if (0) 
	    while ((cq = xmlGetNextTagContent (&cp, "Mim-entry_included_E", h)))
	      {
		cr = strstr (cq, ", INCLUDED") ;
		if (cr) *cr = 0 ;
		printf ("%s %s\n" 
			, "Mim-entry_included"
			, freeprotect (vtextLowerCase(cq))
			) ;
	      }
	  printf ("\n") ;
	}
      return TRUE ;
    }
  
  return FALSE ;
} /* omimParseOne */

/*************************************************************/

static BOOL omimGetOne (const char* omim)
{
  BOOL ok = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  char *base = 0 ; /*"http://www.omim.org/entry/"  https://www.ncbi.nlm.nih.gov/sites/entrez?cmd=Retrieve&db=OMIM&dopt=XML&list_uids=" ; */
  char *xml, *xml2 ;

  messcrash ("we cannot access omim.org via a program ,this then kills all NLM access to OMIM.org") ;

  vTXT txt = vtxtHandleGetURL (messprintf ("%s%s", base, omim), 0, h) ;

  if (txt)
    {
      ok = TRUE ;
      xml = strnew (vtxtPtr (txt), h) ;
      xml2 = xmlGetTagContent (&xml, "title", 0, h) ;
      omimParseOne (xml2, h) ;
    }
  else
    printf ("cannot find %s%s", base, omim) ;

  ac_free (h) ;
  return ok ;
}  /* omimGetOne */

static void usage (void)
{
  printf ("// Usage : omim2ace ( -db $ACEDB  | -q omim_id | -list acedb_style_OMIM.list | file_-i xmlInputFileName\n") ;
}


int main (int argc, const char **argv)
{
  const char *listName = getArgV (&argc, argv, "-list") ;
  const char *omimId = getArgV (&argc, argv, "-q") ;
  const char *xmlInputFileName = getArgV (&argc, argv, "-i") ;
  AC_HANDLE h = 0 ;

  freeinit () ;
  
  if (!listName && !omimId && !xmlInputFileName )
    usage () ;

  h = ac_new_handle () ;
  if (omimId)
    omimGetOne (omimId) ;
  else if (listName)
    {
      char *cp ;
      ACEIN ai = aceInCreate (listName, FALSE, h) ;

      if (! ai)
	messcrash ("cannot open %s", listName) ;
      while (aceInCard (ai)) 
	{
	  cp = aceInWord (ai) ;
	  if (cp && !strcasecmp (cp, "OMIM"))
	    {
	      aceInStep (ai, ':') ;
	      cp = aceInWord (ai) ;
	      if (cp)
		omimGetOne (cp) ;
	    }
	  if (cp && !strcasecmp (cp, "Extern"))
	    {
	      aceInStep (ai, ':') ;
	      cp = aceInWord (ai) ;
	      if (cp && ! strncmp (cp, "OMIM_", 5))
		omimGetOne (cp+5) ;
	    }
	}
    }
  else if (xmlInputFileName)
    {
      char *buf = vtxtFileGetContent (xmlInputFileName, 0, 0) ;
      if (buf)
	omimParseOne (buf, 0) ;
      else
	messcrash ("cannot find xmlInputFile %s", xmlInputFileName) ;
    }

  printf ("// done\n") ;
  ac_free (h) ;
  return 0 ;
}


/*  File: embl.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: creates an EMBL submission file from ACEDB data
 	authors for submission obtained from From_Author field
	current assumptions (not necessarily valid):
	   - all CDS's are one level down, i.e. not recursive
	   - start and end are found for all embedded sequences
	   - start < end for general features: promoters etc.
 * Exported functions:
 * HISTORY:
 * Last edited: May 12 10:55 1999 (edgrif)
 * * May 12 10:51 1999 (edgrif): Added code to emblDumpKey to popdown
 *              the help dialog posted by emblDump.
 * * Jun 28 21:06 1992 (rd): added stuff for cDNA's, and keyset dumps
 * * Feb 23 01:27 1992 (rd): order() -> join()
 * Created: Sun Feb 23 01:26:18 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: embl.c,v 1.6 2020/05/30 16:50:34 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "bs.h"
#include "sysclass.h"
#include "classes.h"
#include "tags.h"
#include "systags.h"
#include "dna.h"
#include "lex.h"
#include "peptide.h"
#include "mytime.h"

static AC_HANDLE handle = 0 ;
static FILE *dumpFile ;		/* file being written to */
static FILE *logFile ;		/* log file for errors/warnings */
static int nLog ;		/* number of log messages */
static Stack text = 0 ;		/* utility buffer */
static Array a = 0 ;		/* work array for bsFlatten */
static BSMARK mark1, mark2, mark3, mark4, mark5, mark6, mark7, mark8 ;
static BOOL isDumpKeyFromDisplay = FALSE ;		    /* Was emblDumpKey called from a
							       mouse selection on a display ? */

/**** tags *****/

static KEY emblKey ;
static KEY _scRNA, _misc_RNA, _Database, _Feature, _AC_number ;
static KEY _EMBL_dump, _EMBL_feature, _EMBL_threshold, _EMBL_qualifier ;
static KEY _EMBL_dump_method, _EMBL_dump_info, _gene_from_name ;
static KEY _ID_template, _ID_division, _DE_format ;
static KEY _OS_line, _OC_line, _RL_submission, _EMBL_reference ;
static KEY _CC_line, _EMBL_chromosome, _source_organism ;
static KEY _EMBL_dump_YES, _EMBL_dump_NO ;

/*******************************/

static void initialise (void)
{
  static int isFirst = TRUE ;

  if (isFirst)
    { KEY databaseClassKey ;

      lexaddkey ("scRNA", &_scRNA, 0) ;
      lexaddkey ("misc_RNA", &_misc_RNA, 0) ;
      lexaddkey ("Database", &_Database, 0) ;
      lexaddkey ("Feature", &_Feature, 0) ;
      lexaddkey ("AC_number", &_AC_number, 0) ;

				/* for features */
      lexaddkey ("EMBL_dump", &_EMBL_dump, 0) ;
      lexaddkey ("EMBL_feature", &_EMBL_feature, 0) ;
      lexaddkey ("EMBL_threshold", &_EMBL_threshold, 0) ;
      lexaddkey ("EMBL_qualifier", &_EMBL_qualifier, 0) ;

				/* for the sequence dump method */
      lexaddkey ("EMBL_dump_info", &_EMBL_dump_info, 0) ;
      lexaddkey ("EMBL_dump_method", &_EMBL_dump_method, 0) ;
      lexaddkey ("ID_template", &_ID_template, 0) ;
      lexaddkey ("ID_division", &_ID_division, 0) ;
      lexaddkey ("DE_format", &_DE_format, 0) ;
      lexaddkey ("OS_line", &_OS_line, 0) ;
      lexaddkey ("OS_line", &_OS_line, 0) ;
      lexaddkey ("OC_line", &_OC_line, 0) ;
      lexaddkey ("RL_submission", &_RL_submission, 0) ;
      lexaddkey ("EMBL_reference", &_EMBL_reference, 0) ;
      lexaddkey ("CC_line", &_CC_line, 0) ;
      lexaddkey ("source_organism", &_source_organism, 0) ;
      lexaddkey ("gene_from_name", &_gene_from_name, 0) ;

				/* for Map objects */
      lexaddkey ("EMBL_chromosome", &_EMBL_chromosome, 0) ;

				/* for Feature_info */
      lexaddkey ("EMBL_dump_YES", &_EMBL_dump_YES, 0) ;
      lexaddkey ("EMBL_dump_NO", &_EMBL_dump_NO, 0) ;

				/* the key for EMBL itself */
      lexaddkey ("Database", &databaseClassKey, _VMainClasses) ;
      lexaddkey ("EMBL", &emblKey, KEYKEY(databaseClassKey)) ;

      isFirst = FALSE ;
    }
}

/*******************************/

static char* emblifyName (KEY key)
{
  static char work[128] ;
  char *cp = work ;
  char *nam = name(key) ;
  OBJ obj ;

  if ((obj = bsCreate (key)))
    { if (bsGetData (obj, str2tag("EMBL_name"), _Text, &cp))
	{ bsDestroy (obj) ;
	  return cp ;
	}
      bsDestroy (obj) ;
    }

  while (*nam && *nam != ' ')
    *cp++ = *nam++ ;

  while (*nam && *nam == ' ') ++nam ;

  if (*nam)			/* initial */
    { *cp++ = ' ' ;
      *cp++ = *nam ;
      *cp++ = '.' ;
    }

  *cp = 0 ;
  return work ;
}

static void emblifyReference (KEY ref)
{ 
  OBJ *obj ;
  int i, ix, year ;
  KEY key, journal ;
  char **cpp, *volume, *p1, *p2, *cp ;

  if (!(obj = bsCreate (ref)))
    { if (logFile) ++nLog, fprintf (logFile, " * Can't open paper object %s\n", name(ref)) ;
      return ;
    }

/* RA */
  if (bsFindTag (obj, _Author) && bsFlatten (obj, 1, a))
    { ix = 0 ;
      for (i = 0 ; i < arrayMax(a) ; i += 1)
	{ cp = emblifyName (arr(a,i,BSunit).k) ;
	  if (ix + strlen (cp) > 70)
	    { fprintf (dumpFile, ",\n") ; ix = 0 ; }
	  if (ix)
	    { fprintf (dumpFile, ", ") ;  ix += 2 ; }
	  else
	    { fprintf (dumpFile, "RA   ") ; ix += 5 ; }
	  fprintf (dumpFile, "%s", cp) ;
	  ix += strlen(cp) ;
	}
      fprintf (dumpFile, ";\n") ;
    }
  else
    if (logFile) ++nLog, fprintf (logFile, " * No authors for reference %s\n", name(ref)) ;

/* RT */
  if (bsGetKey (obj, _Title, &key))
    { stackClear (text) ;
      catText (text, "\"") ;
      catText (text, name(key)) ;
      catText (text, "\";") ;
      for (cpp = uBrokenLines (stackText(text,0), 67) ; *cpp ; ++cpp)
	fprintf (dumpFile, "RT   %s\n", *cpp) ;
    }
  else
    if (logFile) ++nLog, fprintf (logFile, " * No title for reference %s\n", name(ref)) ;

/* RL */
  if (bsGetKey (obj, _Journal, &journal) &&
      bsGetData (obj, _Volume, _Text, &volume) &&
      bsGetData (obj, _Year, _Int, &year) &&
      bsGetData (obj, _Page, _Text, &p1) &&
      bsGetData (obj, _bsRight, _Text, &p2))
    fprintf (dumpFile, "RL   %s %s:%s-%s (%d).\n", 
	     name(journal), volume, p1, p2, year) ;
  else
    if (logFile) ++nLog, fprintf (logFile, " * Bad RL info for reference %s\n", name(ref)) ;

  bsDestroy (obj) ;
}

/*************************************/

static char *emblDate (mytime_t dat)
{ char *cp ;
  static char buf[1000] ;

  timeShowFormat (dat, "%d-%b-%Y", buf, 1000) ;
  for (cp = buf ; *cp ; ++cp)
    *cp = ace_upper(*cp) ;
  return buf ;
}

/*************************************/

char *getAccession (KEY key)
{ 
  OBJ obj ;
  char *acc = 0 ;

  if (!(obj = bsCreate (key)))
    return 0 ;

  if (bsGetData (obj, _AC_number, _Text, &acc))
    ;
  else if (bsFindKey (obj, _Database, emblKey) &&
	   bsGetKeyTags (obj, _bsRight, 0))
    {
      if (bsGetKey (obj, _bsRight, &key))
	acc = name(key) ;
      else
	bsGetData (obj, _bsRight, _Text, &acc) ;
    }

  bsDestroy (obj) ;
  return acc ;
}

/*******************************/

	/* coordinate transformation variables */
static KEY dumpKey ;		/* key for sequence object being dumped */
static int length ;		/* length of sequence being dumped */
static int offset ;		/* start of dumped seq in superseq coords */
static BOOL offForward ;	/* is dumpseq codirectional with superseq */
static Array pieceLeft = 0 ;	/* used to build location info */
static Array pieceRight = 0 ;	/* used to build location info */
static Array piecePrefix = 0 ;	/* for the accession numbers */

/********** convert current super seq coords to dumped seq coords *******/

static void offAdjust (int *x1, int *x2)
{
  if (offForward)
    { *x1 -= offset - 1 ;
      *x2 -= offset - 1 ;
    }
  else
    { *x1 = offset + 1 - *x1 ;
      *x2 = offset + 1 - *x2 ;
    }
}

/************ stuff to maintain map of accnos for DNA fragments *********/

typedef struct {
  KEY key ;
  char *accColon ;		/* accession number followed by colon, if not self */
  int start ;
  int end ;
} ACCMAP ;

static Array accMap ;

/* this code assumes that all participating accessions are
   (i) colinear, i.e. not reverse-complemented wrt each other
   (ii) children of an ancestor of dumpKey
   OK for one-deep LINK structures as at Sanger and St Louis
*/

static void buildAccMap (KEY key, int start, int end)
{ 
  ACCMAP *am ;

  if (!key)			/* initialise */
    { accMap = arrayHandleCreate (16, ACCMAP, handle) ;
      am = arrayp (accMap, 0, ACCMAP) ;
      am->key = dumpKey ;
      am->accColon = "" ;
      am->start = 1 ;
      am->end = length ;
    }
  else
    { offAdjust (&start, &end) ;
      if (end > start)		/* don't handle reversed DNA */
	{ am = arrayp (accMap, arrayMax(accMap), ACCMAP) ;
	  am->key = key ;
	  am->accColon = 0 ;	/* defer finding accession number */
	  am->start = start ;
	  am->end = end ;
	}
    }
}

static ACCMAP *findInAccMap (int x)
{ 
  ACCMAP *am ;
  int i ;
  char *acc ;

  for (i = 0 ; i < arrayMax (accMap) ; )
    { am = arrp (accMap, i, ACCMAP) ;
      if (x >= am->start && x <= am->end)
	{ if (am->accColon)
	    return am ;
	  if ((acc = getAccession (am->key)))
	    { am->accColon = strnew (messprintf ("%s:", acc), handle) ;
	      return am ;
	    }
				/* if no acc number, throw out of accMap */
	  if (arrayMax(accMap))
	    { arr (accMap, i, ACCMAP) = arr (accMap, arrayMax(accMap)-1, ACCMAP) ;
	      --arrayMax(accMap) ;
	    }
	}
      else
	++i ;
    }

  return 0 ;
}

/***** stuff to use accMap to make location info, and write it out *****/

static BOOL makePieces (int x1, int x2)
{
  ACCMAP *am ;

  if (!(am = findInAccMap (x1)))
    { if (logFile) ++nLog, fprintf (logFile, " * can't find accession number for coord %d\n", x1) ;
      return FALSE ;
    }

  if (x2 <= am->end)
    { array (pieceLeft, arrayMax(pieceLeft), int) = x1 - am->start + 1 ;
      array (pieceRight, arrayMax(pieceRight), int) = x2 - am->start + 1 ;
      array (piecePrefix, arrayMax(piecePrefix), char*) = am->accColon ;
      return TRUE ;	/* done! */
    }
  else
    return (makePieces (x1, am->end) && makePieces (am->end+1, x2)) ;
}      

void writeLocation (BOOL isComp, BOOL startNotFound, BOOL endNotFound)
{
  int i ;
  char **cpp ;
				/* transfer piece info to text */
  stackClear (text) ;
  if (isComp) catText (text, "complement(") ;
  if (arrayMax(pieceLeft) > 1) catText (text, "join(") ;
  for (i = 0 ; i < arrayMax(pieceLeft) ; ++i)
    { if (i)
	catText (text, ",") ;

      if (startNotFound && !i)
	catText (text, "<") ;
      catText (text, arr(piecePrefix,i,char*)) ;
      catText (text, messprintf ("%d", arr(pieceLeft,i,int))) ;

      if (arr(pieceLeft,i,int) < arr(pieceRight,i,int))
	{ catText (text, "..") ;
	  if (endNotFound && i == arrayMax(pieceLeft)-1)
	    catText (text, ">") ;
	  catText (text, messprintf ("%d", arr(pieceRight,i,int))) ;
	}
    }
  if (arrayMax(pieceLeft) > 1) catText (text, ")") ;
  if (isComp) catText (text, ")") ;

  for (cpp = uBrokenLines (stackText(text,0), 56) ; *cpp ; )
    { fprintf (dumpFile, "%s\n", *cpp) ;
      if (*++cpp)
	fprintf (dumpFile, "FT                   ") ;
    }
}

/***********************************************************************/
/********** a little package to get dump method information ************/

static OBJ lastGetDumpObj = 0 ;

static Array getDumpObjects (KEY seq)
{ 
  KEY key ;
  OBJ obj ;
  static KEY lastSeq = 0 ;
  static Array a = 0 ;
  int i ;

  if (seq != lastSeq)
    { if (a)
	for (i = 0 ; i < arrayMax (a) ; ++i)
	  bsDestroy (arr(a, i, OBJ)) ;
      a = arrayReCreate (a, 16, OBJ) ;
      key = seq ;
      while ((obj = bsCreate (key)))
	{ if (!bsFindTag (obj, _EMBL_dump_info) ||
	      !bsPushObj (obj))
	    { bsDestroy (obj) ;
	      break ;
	    }
	  array (a, arrayMax(a), OBJ) = obj ;
	  if (!bsGetKey (obj, _EMBL_dump_method, &key))
	    break ;
	}
      lastSeq = seq ;
    }

  return a ;
}

static BOOL getDumpText (KEY seq, KEY tag, char **cpp)
{ 
  Array a = getDumpObjects (seq) ;
  int i ;

  lastGetDumpObj = 0 ;
  for (i = 0 ; i < arrayMax(a) ; ++i)
    if (bsGetData (arr(a, i, OBJ), tag, _Text, cpp))
      { lastGetDumpObj = arr(a, i, OBJ) ; 
	return TRUE ;
      }

  return FALSE ;
}

static BOOL getDumpKey (KEY seq, KEY tag, KEY *kp)
{ 
  Array a = getDumpObjects (seq) ;
  int i ;

  lastGetDumpObj = 0 ;
  for (i = 0 ; i < arrayMax(a) ; ++i)
    if (bsGetKey (arr(a, i, OBJ), tag, kp))
      { lastGetDumpObj = arr(a, i, OBJ) ; 
	return TRUE ;
      }

  return FALSE ;
}

static BOOL getDumpFlag (KEY seq, KEY tag)
{
  Array a = getDumpObjects (seq) ;
  int i ;

  lastGetDumpObj = 0 ;
  for (i = 0 ; i < arrayMax(a) ; ++i)
    if (bsFindTag (arr(a, i, OBJ), tag))
      { lastGetDumpObj = arr(a, i, OBJ) ; 
	return TRUE ;
      }

  return FALSE ;
}

/********* main subroutine to dump a sub sequence (e.g. CDS) *************/

static BOOL isGeneFromName ;

static void subseqFeatures (KEY key, int start, int end)
{
  OBJ	obj  ;
  BOOL	isComp ;
  int	i, pos ;
  KEY	method, predMethod, remark, cdna, gene ;
  char  **cpp ;
  char  *rnaText, *featText ;
  BOOL  pseudo = FALSE, isCDS = FALSE ;
  int	miss1, miss2 ;

  offAdjust (&start, &end) ;
  isComp = (start > end) ;
  if ((!isComp && (end < 1 || start > length)) ||
      (isComp && (start < 1 || end > length)))
    return ;
  
  if (!(obj = bsCreate (key)))
    { if (logFile) ++nLog, fprintf (logFile, " * can't open subobject %s\n", name(key)) ;
      return ;
    }

  if (!bsGetKey (obj, str2tag("Method"), &method))
    goto finish ;
  { OBJ obj2 = bsCreate (method) ;
    if (!obj2) goto finish ;
    if (!bsGetData (obj2, _EMBL_feature, _Text, &featText))
      { bsDestroy (obj2) ; goto finish ; }
    bsDestroy (obj2) ;
  }

  isCDS = FALSE ;
  if (!strcmp (featText, "CDS"))
    {
      if (bsFindTag (obj, _CDS))
	isCDS = TRUE ;
      else if (bsFindTag (obj, _Pseudogene))
	pseudo = TRUE ;
      else
	{
	  if (logFile) ++nLog, fprintf (logFile, " * %s does not have CDS set when method requires that\n", name(key)) ;
	  goto finish ;
	}
    }

  if (isCDS)
    { Array p  ;
      signed char stop = pepEncodeChar[(int)'*'] ;
      int i ;

      if ((p = peptideGet (key)))
	{ for (i = 0 ; i < arrayMax(p) ; ++i)
	    if (arr(p,i,signed char) == stop)
	      { if (logFile) ++nLog, fprintf (logFile, " * %s made a pseudogene because contains stop codon\n",
		name(key)) ;
		pseudo = TRUE ; 
		break ; 
	      }
	  arrayDestroy (p) ;
	}
      else
	{ if (logFile) ++nLog, fprintf (logFile, " * %s made a pseudogene because can't get peptide sequence\n",
		name(key)) ;
	  pseudo = TRUE ;
	}
    }

			/* make the location info in piece Arrays */
  pieceLeft = arrayHandleCreate (32, int, handle) ;
  pieceRight = arrayHandleCreate (32, int, handle) ;
  piecePrefix = arrayHandleCreate (16, char*, handle) ;

  if (bsFindTag (obj, _Source_Exons) && bsFlatten (obj, 2, a))
    { dnaExonsSort (a) ;
      if (isComp)
	for (i = arrayMax(a)-1 ; i > 0 ; i -= 2)
	  { if (!makePieces (start + 1 - arr(a,i,BSunit).i, 
			     start + 1 - arr(a,i-1,BSunit).i))
	      { if (logFile) ++nLog, fprintf (logFile, " * makePieces failed for %s when dumping %s\n",
			 name(key), name(dumpKey)) ;
		goto finish ; 
	      }
	  }
      else
	for (i = 0; i < arrayMax(a) ; i += 2)
	  { if (!makePieces (start - 1 + arr(a,i,BSunit).i,
			     start - 1 + arr(a,i+1,BSunit).i))
	      { if (logFile) ++nLog, fprintf (logFile, " * makePieces failed for %s when dumping %s\n",
			 name(key), name(dumpKey)) ;
		goto finish ;
	      }
	  }
    }
  else
    if (isComp)
      { if (!makePieces (end, start))
	  { if (logFile) ++nLog, fprintf (logFile, " * makePieces failed for %s when dumping %s\n",
		     name(key), name(dumpKey)) ;
	    goto finish ;
	  }
      }
    else
      { if (!makePieces (start, end))
	  { if (logFile) ++nLog, fprintf (logFile, " * makePieces failed for %s when dumping %s\n",
		     name(key), name(dumpKey)) ;
	    goto finish ;
	  }
      }

				/* check start and end found */
  miss1 = miss2 = 0 ;
  if (!isComp)
    { 
      if (!miss1 && bsFindTag (obj, _Start_not_found))
	{
	  if (bsGetData (obj, _bsRight, _Int, &i))
	    miss1 = 2 + i ;
	  else
	    miss1 = 3 ;		/* no frame change */
	}
      if (!miss2 && bsFindTag (obj, _End_not_found))
	miss2 = 3 ;
    }
  else
    { 
      if (!miss2 && bsFindTag (obj, _Start_not_found))
	{
	  if (bsGetData (obj, _bsRight, _Int, &i))
	    miss2 = 2 + i ;
	  else
	    miss2 = 3 ;		/* no frame change */
	}
      if (!miss1 && bsFindTag (obj, _End_not_found))
	miss1 = 3 ;
    }
      
				/* write feature key and location */
  fprintf (dumpFile, "FT   %-16s", featText) ;
  writeLocation (isComp, (miss1 > 0), (miss2 > 0)) ;

				/* write qualifiers */

  if (bsGetData (obj, str2tag("tRNA"), _Text, &rnaText))
    fprintf (dumpFile, "FT                   /note=\"%s-tRNA\"\n", rnaText) ;
  else if (bsGetData (obj, str2tag("rRNA"), _Text, &rnaText) ||
	   bsGetData (obj, str2tag("snRNA"), _Text, &rnaText) ||
	   bsGetData (obj, str2tag("scRNA"), _Text, &rnaText) ||
	   bsGetData (obj, str2tag("misc_RNA"), _Text, &rnaText))
    fprintf (dumpFile, "FT                   /note=\"%s-RNA\"\n", rnaText) ;

				/* pseudo */
  if (pseudo)
    fprintf (dumpFile, "FT                   /pseudo\n") ;

				/* prediction method */
  if (bsGetKey (obj, _CDS_predicted_by, &predMethod))
    fprintf (dumpFile, "FT                   /note=\"predicted using %s\"\n",
	     name(predMethod)) ;

	/* 970214 extra note for preliminary genes (.a etc.) */
  { char *cp ;
    for (cp = name(key) ; *cp ; ++cp) if (*cp == '.') break ;
    if (*cp == '.' && *++cp >= 'a' && *cp <= 'z')
      fprintf (dumpFile, "FT                   /note=\"preliminary prediction\"\n") ;
  }

				/* Locus info */
  if (bsGetKey (obj, _Locus, &gene))
    fprintf (dumpFile, "FT                   /gene=\"%s\"\n", name(gene)) ;

  if (isGeneFromName)
    fprintf (dumpFile, "FT                   /gene=\"%s\"\n", name(key)) ;

				/* DB_remark or brief_id info */
  stackClear (text) ;
  if (bsGetKey (obj, _DB_remark, &remark))
    { catText (text, messprintf ("/note=\"%s", name(remark))) ;
      while (bsGetKey (obj, _bsDown, &remark))
	catText (text, messprintf (", %s", name(remark))) ;
    }
  else if (bsGetKey (obj, _Brief_identification, &remark))
    catText (text, messprintf ("/note=\"similar to %s", name(remark))) ;
  if (stackMark(text))
    { catText (text, "\"") ;
      for (cpp = uBrokenLines(stackText(text,0), 56) ; *cpp ; ++cpp)
	fprintf (dumpFile, "FT                   %s\n", *cpp) ;
    }

				/* codon start for partial CDS genes */
  if (bsFindTag (obj, _CDS))
    {
      if (miss1 && !isComp)
	fprintf (dumpFile, "FT                   /codon_start=%d\n", 
		 3 - ((miss1 - 1) % 3)) ; /* was 1 + (miss1%3) */
      else if (miss2 && isComp)
	fprintf (dumpFile, "FT                   /codon_start=%d\n", 
		 3 - ((miss2 - 1) % 3)) ; /* RMD 950703 - as above */
    }

				/* matching cDNAs */
  if (bsGetKey (obj, _Matching_cDNA, &cdna))
    { int nEST = 0, nClump, n ;
      char *acc ;

      do ++nEST ; while (bsGetKey (obj, _bsDown, &cdna)) ;
      nClump = 1 + nEST/40 ;

      n = 0 ;
      bsGetKey (obj, _Matching_cDNA, &cdna) ; /* restart at first cDNA */
      do
	{ fprintf (dumpFile, "FT                   ") ;
	  if (!(n % nClump))
	    fprintf (dumpFile, "/note=\"") ;
	  if ((acc = getAccession(cdna)))
	    fprintf (dumpFile, "cDNA EST EMBL:%s comes from this gene", acc) ;
	  else
	    { fprintf (dumpFile, "cDNA EST %s comes from this gene", name(cdna)) ;
/* RD 970815: sjj request because of new Yuji cDNAs
	      if (logFile) ++nLog, fprintf (logFile, " * no accession number for EST %s\n", name(cdna)) ; */
	    }
	  ++n ;
	  if (n == nEST || !(n % nClump))
	    fprintf (dumpFile, "\"\n") ;
	  else
	    fprintf (dumpFile, ";\n") ;
	} while (bsGetKey (obj, _bsDown, &cdna)) ;
    }

				/* TSLs */
  if (bsGetData (obj, _TSL_site, _Int, &pos))
    { ACCMAP *am ;
      pos = isComp?start+1-pos:start-1+pos ;
      if ((am = findInAccMap (pos)))
	fprintf (dumpFile,
"FT                   /note=\"Possible trans-spliced leader site at %s%d\"\n",
		 am->accColon, pos - am->start) ;
    }

 finish:
  bsDestroy (obj) ;
}

/**********************************************/

typedef struct {
  char *feature ;
  float threshold ;
  BOOL isThreshold ;
  Array quals ;			/* of QUALIFIER */
} METHOD_DUMP ;

typedef struct {
  char *text ;
  char *arg ;
} QUALIFIER ;

static Associator methodAss ;
static AC_HANDLE massHandle = 0 ;

static METHOD_DUMP *methodDumpInfo (KEY method)
{
  void *v ;

  if (!massHandle) 
    { massHandle = handleCreate() ;
      methodAss = assHandleCreate (massHandle) ;
    }
  if (assFind (methodAss, assVoid(method), &v))
    { if (v == assVoid(1))
	return 0 ;
      else
	return (METHOD_DUMP*) v ;
    }
				/* if get here, a new method */
  { OBJ obj = bsCreate(method) ;
    METHOD_DUMP *md = assVoid(1) ;
    char *cp ;
    static BSMARK mark = 0 ;

    if (obj)
      { md = (METHOD_DUMP*) messalloc (sizeof(METHOD_DUMP)) ;
	if (!bsGetData (obj, _EMBL_feature, _Text, &md->feature))
	  { messfree (md) ;
	    md = assVoid(1) ;
	    goto done ;
	  }

	md->isThreshold = bsGetData (obj, _EMBL_threshold, 
				     _Float, &md->threshold) ;

	if (bsGetData (obj, _EMBL_qualifier, _Text, &cp))
	  { QUALIFIER *q ;
	    md->quals = arrayHandleCreate (8, QUALIFIER, massHandle) ;
	    do
	      { q = arrayp(md->quals, arrayMax(md->quals), QUALIFIER) ;
		q->text = strnew (cp, massHandle) ;
		mark = bsMark (obj, mark) ;
		if (bsGetData (obj, _bsRight, _Text, &cp))
		  q->arg = strnew (cp, massHandle) ;
		else
		  q->arg = 0 ;
		bsGoto (obj, mark) ;
	      } while (bsGetData (obj, _bsDown, _Text, &cp)) ;
	  }
	bsDestroy (obj) ;
      }

  done:
    assInsert (methodAss, assVoid(method), md) ;
    if ((void*)md == assVoid(1))
      return 0 ;
    else
      return md ;
  }
}

static void dumpFeatureHomol (OBJ obj, BOOL isFeature, 
			      KEY method, int start, int end, 
			      KEY target, float score, int y1, int y2)
{
  BOOL isComp ;
  char *note ;
  BOOL hasScore ;
  BOOL hasNote ;
  BOOL hasSubObj ;
  METHOD_DUMP *md ;

  md = methodDumpInfo (method) ; /* already checked this exists */

  if (isFeature)		/* a Feature line */
    { if (bsGetData (obj, _bsRight, _Float, &score))
	{ hasScore = TRUE ;
	  if (bsGetData (obj, _bsRight, _Text, &note))
	    hasNote = TRUE ;
	  else
	    hasNote = FALSE ;
	}
      else
	hasScore = hasNote = FALSE ;
    }
  else				/* a Homol line */
    { hasScore = TRUE ;
      hasNote = FALSE ;
    }

  hasSubObj = bsPushObj (obj) ;

  if (hasSubObj && bsFindTag (obj, _EMBL_dump_NO))
    return ;

				/* check score limit */
  if (md->isThreshold && !(hasSubObj && bsFindTag (obj, _EMBL_dump_YES)))
    if (!hasScore || score < md->threshold)
      return ;

			/* adjust offsets, return if out of bounds */
  offAdjust (&start, &end) ;
  isComp = (start > end) ;
  if (isComp)
    { int tmp ;
      tmp = start ; start = end ; end = tmp ;
      tmp = y1 ; y1 = y2 ; y2 = tmp ;
    }
  if (end < 1 || start > length)
    return ;
				/* write primary line */
  fprintf (dumpFile, "FT   %-16s", md->feature) ;
  if (isComp) fprintf (dumpFile, "complement(") ;
  if (start >= 1) fprintf (dumpFile, "%d..", start) ;
    else fprintf (dumpFile, "<1..") ;
  if (end <= length) fprintf (dumpFile, "%d", end) ;
    else fprintf (dumpFile, ">%d", length) ;
  if (isComp) fprintf (dumpFile, ")") ;
  fprintf (dumpFile, "\n") ;

				/* qualifiers */
				/* first from the Method */
  if (md->quals)
    { int i ;
      QUALIFIER *q ;
      char **cpp ;
      static Stack qualStack = 0 ;

      qualStack = stackReCreate (qualStack, 128) ;
      for (i = 0 ; i < arrayMax (md->quals) ; ++i)
	{ q = arrp(md->quals,i, QUALIFIER) ;
	  if (*q->text == '/' && *stackText(qualStack,0))
	    { for (cpp = uBrokenLines(stackText(qualStack,0), 56) ;
		   *cpp ; ++cpp)
		fprintf (dumpFile, "FT                   %s\n", *cpp) ;
	      stackClear (qualStack) ;
	    }
	  if (!q->arg)		/* the simple case */
	    catText (qualStack, q->text) ;
	  else if (!strcmp (q->arg, "score"))
	    { if (!hasScore) continue ;
	      catText (qualStack, messprintf(q->text, score)) ;
	    }
	  else if (!strcmp (q->arg, "note"))
	    { if (!hasNote) continue ;
	      catText (qualStack, messprintf(q->text, note)) ;
	    }
	  else if (!strcmp (q->arg, "target"))
	    { if (isFeature) continue ;
	      catText (qualStack, messprintf(q->text, name(target))) ;
	    }
	  else if (!strcmp (q->arg, "y1"))
	    { if (isFeature) continue ;
	      catText (qualStack, messprintf(q->text, y1)) ;
	    }
	  else if (!strcmp (q->arg, "y2"))
	    { if (isFeature) continue ;
	      catText (qualStack, messprintf(q->text, y2)) ;
	    }
	  else
	    { if (logFile) ++nLog, fprintf (logFile, " * bad EMBL_qualifier arg %s\n", q->arg) ;
	      q->arg = "note" ;	/* prevent reoccurence */
	    }
#ifdef BAD_IDEA
	  else if (hasSubObj && 
		   bsGetData (obj, str2tag(q->arg), _Text, &cp)) 
	    do
	      { for (cpp = uBrokenLines(messprintf (q->text, cp), 56) ; 
		     *cpp ; ++cpp)
		  fprintf (dumpFile, "FT                   %s\n", *cpp) ;
	      } while (bsGetData (obj, _bsDown, _Text, &cp)) ;
#endif
	}
      if (*stackText(qualStack,0))
	{ for (cpp = uBrokenLines(stackText(qualStack,0), 56) ;
	       *cpp ; ++cpp)
	    fprintf (dumpFile, "FT                   %s\n", *cpp) ;
	  stackClear (qualStack) ;
	}
    }
				/* then from the object */
  { char *cp, **cpp ;
    if (hasSubObj && bsGetData (obj, _EMBL_qualifier, _Text, &cp)) do
      { for (cpp = uBrokenLines(cp, 56) ; *cpp ; ++cpp)
	  fprintf (dumpFile, "FT                   %s\n", *cpp) ;
      } while (bsGetData (obj, _bsDown, _Text, &cp)) ;
  }
}

/**********************************************/

static void emblDoDump (KEY seq)
{
  KEY key, clone = 0 ;
  OBJ Seq = 0, obj ;
  Array dna = 0 ;
  char *translate ;
  char *message = 0, *cp, *cq, **cpp ;
  int i, ix, freq[5] ;
  BOOL isRNA ;

#define abort(x)	{ message = x ; goto finish ; }
#define xxout		fprintf (dumpFile, "XX\n")

  initialise() ;

  handle = handleCreate() ;

  a = arrayReCreate (a, 32, BSunit) ;
  text = stackReCreate (text, 256) ;

  dumpKey = seq ;

  fprintf (logFile, "%s\n", name(seq)) ;

  if (class(seq) != _VSequence || !(Seq = bsCreate(seq)))
    abort("Sequence object missing") ;
  if (!(dna = dnaGet(seq)))
    abort("No DNA attached to Sequence object") ;
  if (!bsGetKey (Seq, _Clone, &clone))
    if (logFile) ++nLog, fprintf (logFile, " * No clone for sequence %s\n", name(seq)) ;

  isRNA = bsFindTag (Seq, _RNA) ;
  length = arrayMax(dna) ;

				/* start dumping! */
/* ID */
  fprintf (dumpFile, "ID   ") ;
  if (bsFindKey (Seq, _Database, emblKey) &&
      bsGetData (Seq, _bsRight, _Text, &cq)) /* ID is text */
    fprintf (dumpFile, "%-10s", cq) ;
  else if (bsFindKey (Seq, _Database, emblKey) &&
	   bsGetKey (Seq, _bsRight, &key)) /* ID is key */
    fprintf (dumpFile, "%-10s", name(key)) ;
  else if (getDumpText (seq, _ID_template, &cp))
    fprintf (dumpFile, cp, name(seq)) ;
  else
    fprintf (dumpFile, "%s", name(seq)) ;
  fprintf (dumpFile, " standard; %s;", isRNA ? "RNA" : "DNA") ;
  if (getDumpText (seq, _ID_division, &cp))
    fprintf (dumpFile, " %s;", cp) ;
  else
    { fprintf (dumpFile, " XXX;") ;
      if (logFile) ++nLog, fprintf (logFile, " * division unknown for ID line - missing from method\n") ;
    }
  fprintf (dumpFile, " %d BP.\n", length) ;
  xxout ;

/* AC */
  { BOOL foundAC = FALSE ;
    if (bsFindKey (Seq, _Database, emblKey) &&
	bsGetKeyTags (Seq, _bsRight, 0))
      {
	if (bsGetKey (Seq, _bsRight, &key))
	  { 
	    do fprintf (dumpFile, "AC   %s;\n", name(key)) ;
            while (bsGetKey (Seq, _bsDown, &key)) ;
	    foundAC = TRUE ;
	  }
	else if (bsGetData (Seq, _bsRight, _Text, &cp))
	  { 
	    do fprintf (dumpFile, "AC   %s;\n", cp) ;
	    while (bsGetData (Seq, _bsDown, _Text, &cp)) ;
	    foundAC = TRUE ;
	  }
      }
    if (!foundAC && (cp = getAccession (seq))) /* maybe under _AC_number */
      { fprintf (dumpFile, "AC   %s;", cp) ;
        foundAC = TRUE ;
      }
    if (foundAC)
      xxout ;
  }

/* DE */
  if (getDumpText (seq, _DE_format, &cp))
    { for (cpp = uBrokenLines(messprintf (cp, clone ? name(clone) : name(seq)), 67) ; 
	   *cpp ; ++cpp)
	fprintf (dumpFile, "DE   %s\n", *cpp) ;
      xxout ;
    }

/* KW */
  stackClear (text) ;
  if (bsGetKey (Seq, _Keyword, &key)) do
    { if (strlen(name(key)) > 60)
	{ messout ("Keyword too long: %s", name(key)) ;
	  continue ;
	}
      if (strlen(stackText (text,0)) + strlen(name(key)) + 3 > 65)
	{ fprintf (dumpFile, "KW   %s;\n", stackText(text,0)) ;
	  stackClear (text) ;
	  catText (text, name(key)) ;
	}
      else
	{ if (stackMark (text))
	    catText (text, "; ") ;
	  catText (text, name(key)) ;
	}
    } while (bsGetKey (Seq, _bsDown, &key)) ;
  fprintf (dumpFile, "KW   %s.\n", stackText(text,0)) ;
  xxout ;

/* OS */
  if (getDumpText (seq, _OS_line, &cp))
    fprintf (dumpFile, "OS   %s\n", cp) ;
/* OC */
  if (getDumpText (seq, _OC_line, &cp)) do
    fprintf (dumpFile, "OC   %s\n", cp) ;
  while (bsGetData (lastGetDumpObj, _bsDown, _Text, &cp)) ;
  xxout ;

/* references */
  { int nref = 0 ;
				/* submission reference */
    if (getDumpText (seq,  _RL_submission, &cp))
      { fprintf (dumpFile, "RN   [%d]\n", ++nref) ;
	fprintf (dumpFile, "RP   1-%d\n", length) ;
	ix = 0 ;
	if (bsGetKey (Seq, _From_Author, &key)) do
	  { cp = emblifyName (key) ;
	    if (ix + strlen (cp) > 70)
	      { fprintf (dumpFile, ",\n") ; ix = 0 ; }
	    if (ix)
	      { fprintf (dumpFile, ", ") ;  ix += 2 ; }
	    else
	      { fprintf (dumpFile, "RA   ") ; ix += 5 ; }
	    fprintf (dumpFile, "%s", cp) ;
	    ix += strlen(cp) ;
	  }
	while (bsGetKey (Seq, _bsDown, &key)) ;
	if (bsGetKey (Seq, str2tag("Previous_Author"), &key)) do
	  { cp = emblifyName (key) ;
	    if (ix + strlen (cp) > 70)
	      { fprintf (dumpFile, ",\n") ; ix = 0 ; }
	    if (ix)
	      { fprintf (dumpFile, ", ") ;  ix += 2 ; }
	    else
	      { fprintf (dumpFile, "RA   ") ; ix += 5 ; }
	    fprintf (dumpFile, "%s", cp) ;
	    ix += strlen(cp) ;
	  }
	while (bsGetKey (Seq, _bsDown, &key)) ;
	if (ix)
	  fprintf (dumpFile, ";\n") ;
	else
	  if (logFile) ++nLog, fprintf (logFile, " * No authors for submission reference\n") ;
	fprintf (dumpFile, "RT   ;\n") ;
			/* substitute DD-MMM-YYYY in RL line by date */
	{ mytime_t date ;
	  char *dateString ;

	  if (!bsGetData (Seq, _Submitted, _DateType, &date))
	    date = timeParse ("today") ;
	  dateString = emblDate(date) ;

	  if (getDumpText (seq, _RL_submission, &cp)) do
	    { cp = strnew (cp, 0) ; /* edit copy! */
	      for (cq = cp ; *cq ; ++cq)
		if (!strncmp (cq, "DD-MMM-YYYY", 11))
		  strncpy (cq, dateString, 11) ;
	      for (cpp = uBrokenLines(cp,67) ; 
		   *cpp ; ++cpp)
		fprintf (dumpFile, "RL   %s\n", *cpp) ;
	      messfree (cp) ;
	    }
	  while (bsGetData (lastGetDumpObj, _bsDown, _Text, &cp)) ;
	}
	xxout ;
      }
				/* standard refs from Method */
    if (getDumpKey (seq, _EMBL_reference, &key)) do
      { fprintf (dumpFile, "RN   [%d]\n", ++nref) ;
	emblifyReference (key) ;
	xxout ;
      }
    while (bsGetKey (lastGetDumpObj, _bsDown, &key)) ;
  }

/* CC */
  if (getDumpText (seq, _CC_line, &cp)) do
    { for (cpp = uBrokenLines(messprintf (cp, clone?name(clone):name(seq)),
			      67) ; *cpp ; ++cpp)
        fprintf (dumpFile, "CC   %s\n", *cpp) ;
      xxout ;
    }
  while (bsGetData (lastGetDumpObj, _bsDown, _Text, &cp)) ;

				/* left and right ends of clone */
  if (clone && (obj = bsCreate (clone)))
    { KEY sseq, cclone ;
      OBJ Sseq ;

      if ((bsFindKey (Seq, _Clone_left_end, clone) ||
	   bsFindKey (Seq, _Clone_left_end, seq)) &&
	  bsGetData (Seq, _bsRight, _Int, &ix) && ix == 1 &&
	  (bsFindKey (Seq, _Clone_right_end, clone) ||
	   bsFindKey (Seq, _Clone_right_end, seq)) && 
	  bsGetData (Seq, _bsRight, _Int, &ix) && ix == length)
	fprintf (dumpFile, "CC   This sequence is the entire insert of clone %s.\n", name(clone)) ;
      else
	{
 fprintf (dumpFile, "CC   IMPORTANT: This sequence is not the entire insert of clone %s.\n", name(clone)) ;
 fprintf (dumpFile, "CC   It may be shorter because we only sequence overlapping\n") ;
 fprintf (dumpFile, "CC   sections once, or longer because we arrange for a small\n") ;
 fprintf (dumpFile, "CC   overlap between neighbouring submissions.\n") ;

				/* left end of this clone */
      if (bsFindKey (Seq, _Clone_left_end, clone) && 
	  bsGetData (Seq, _bsRight, _Int, &ix))
	fprintf (dumpFile, "CC   The true left end of clone %s is at %d in this sequence.\n", 
		 name(clone), ix) ;
      else if (bsGetKey (obj, _Clone_left_end, &sseq) && (Sseq = bsCreate (sseq)))
	{ if (bsFindKey (Sseq, _Clone_left_end, clone) &&
	      bsGetData (Sseq, _bsRight, _Int, &ix))
	  {	    
	    if ((cp = getAccession (sseq)))
	      fprintf (dumpFile, "CC   The true left end of clone %s is at %d in\n"
		       "CC   sequence %s.\n",  name(clone), ix, cp) ;
	    else
	      {
		if (logFile) ++nLog, fprintf (logFile, " * no accession for left-end containing sequence %s\n", 
					      name(sseq)) ;
	      }
	  }
	  bsDestroy (Sseq) ;
	}
				/* right end of this clone */
      if (bsFindKey (Seq, _Clone_right_end, clone) && 
	  bsGetData (Seq, _bsRight, _Int, &ix))
	fprintf (dumpFile, "CC   The true right end of clone %s is at %d in this sequence.\n", 
		 name(clone), ix) ;
      else if (bsGetKey (obj, _Clone_right_end, &sseq) && (Sseq = bsCreate (sseq)))
	{ 
	  if (bsFindKey (Sseq, _Clone_right_end, clone) &&
	      bsGetData (Sseq, _bsRight, _Int, &ix))
	    {
	      if ((cp = getAccession (sseq)))
		fprintf (dumpFile, "CC   The true right end of clone %s is at %d in\n"
			 "CC   sequence %s.\n", name(clone), ix, cp) ;
	      else
		{
		  if (logFile) ++nLog, fprintf (logFile, " * no accession for right-end containing sequence %s\n",
						name(sseq)) ;
		}
	    }
	  bsDestroy (Sseq) ;
	}
}
				/* now report other clone ends */
      if (bsFindTag (Seq, _Clone_left_end) && bsFlatten (Seq, 2, a))
	for (i = 0 ; i < arrayMax(a) ; i += 2)
	  if (arr(a, i, BSunit).k != clone && 
	      arr(a, i, BSunit).k != seq && 
	      arr(a, i+1, BSunit).i)
	    { if ((Sseq = bsCreate(arr(a, i, BSunit).k)) &&
		  bsGetKey (Sseq, _Clone, &cclone))
		fprintf (dumpFile, "CC   The true left end of clone %s is at %d in this sequence.\n", 
			 name(cclone), arr(a, i+1, BSunit).i) ;
	      else
		fprintf (dumpFile, "CC   The true left end of clone %s is at %d in this sequence.\n", 
			 name(arr(a, i, BSunit).k), arr(a, i+1, BSunit).i) ;
	      if (Sseq) bsDestroy (Sseq) ;
	    }
      if (bsFindTag (Seq, _Clone_right_end) && bsFlatten (Seq, 2, a))
	for (i = 0 ; i < arrayMax(a) ; i += 2)
	  if (arr(a, i, BSunit).k != clone && 
	      arr(a, i, BSunit).k != seq && 
	      arr(a, i+1, BSunit).i)
	    { if ((Sseq = bsCreate(arr(a, i, BSunit).k)) &&
		  bsGetKey (Sseq, _Clone, &cclone))
		fprintf (dumpFile, "CC   The true right end of clone %s is at %d in this sequence.\n", 
			 name(cclone), arr(a, i+1, BSunit).i) ;
	      else
		fprintf (dumpFile, "CC   The true right end of clone %s is at %d in this sequence.\n", 
			 name(arr(a, i, BSunit).k), arr(a, i+1, BSunit).i) ;
	      if (Sseq) bsDestroy (Sseq) ;
	    }

      bsDestroy (obj) ;
    }

				/* left and right sequence overlaps */
  if (bsGetKey (Seq, _Overlap_left, &key))
    {
      if ((cp = getAccession (key)))
	{
	  stackClear (text) ;
	  catText (text, "The start of this sequence ") ;
	  if ((obj = bsCreate (key)))
	    { 
	      int len, off ;
	      KEY dna ;
	      if (bsGetKey (obj, _DNA, &dna) && 
		  bsGetData (obj, _bsRight, _Int, &len) &&
		  bsFindKey (obj, _Overlap_right, seq) &&
		  bsGetData (obj, _bsRight, _Int, &off))
		catText (text, messprintf ("(1..%d) ", len - off + 1)) ;
	      else
		if (logFile) ++nLog, fprintf (logFile, " * no Overlap_right data for left neighbour %s\n", name(key)) ;
	      bsDestroy (obj) ;
	    }
	  catText (text, messprintf ("overlaps with the end of sequence %s.", cp)) ;
	  for (cpp = uBrokenLines (stackText(text,0), 67) ; *cpp ; ++cpp)
	    fprintf (dumpFile, "CC   %s\n", *cpp) ;
	}
      else
	if (logFile) ++nLog, fprintf (logFile, " * no accession for left overlap sequence %s\n", name(key)) ;
    }

  if (bsGetKey (Seq, _Overlap_right, &key))
    {
      if ((cp = getAccession (key)))
	{ 
	  stackClear (text) ;
	  catText (text, "The end of this sequence ") ;
	  if (bsGetData (Seq, _bsRight, _Int, &ix))
	    catText (text, messprintf ("(%d..%d) ", ix, length)) ;
	  else
	    if (logFile) ++nLog, fprintf (logFile, " * no Overlap_right amount to right neighbour %s\n", name(key)) ;
	  catText (text, messprintf ("overlaps with the start of sequence %s.", cp)) ;
	  for (cpp = uBrokenLines (stackText(text,0), 67) ; *cpp ; ++cpp)
	    fprintf (dumpFile, "CC   %s\n", *cpp) ;
	}
      else
	if (logFile) ++nLog, fprintf (logFile, " * no accession for right overlap sequence %s\n", name(key)) ;
      xxout ;
    }
				/* explicit DB remarks */
  if (bsFindTag (Seq, _DB_remark) && bsFlatten (Seq, 1, a))
    for (i = 0 ; i < arrayMax(a) ; ++i)
      { for (cpp = uBrokenLines(name(arr(a,i,BSunit).k),67) ; 
	     *cpp ; ++cpp)
	  fprintf (dumpFile, "CC   %s\n", *cpp) ;
	xxout ;
      }

/* feature table */
  fprintf (dumpFile, "FH   Key             Location/Qualifiers\n") ;
  fprintf (dumpFile, "FH\n") ;

				/* source */
  fprintf (dumpFile, "FT   source          1..%d\n", length) ;
  if (getDumpText (seq, _source_organism, &cp))
    fprintf (dumpFile, "FT                   /organism=\"%s\"\n", cp) ;
  if (clone)
    fprintf (dumpFile, "FT                   /clone=\"%s\"\n", name(clone)) ;
  if (getDumpText (seq, str2tag("EMBL_chromosome"), &cp))
    fprintf (dumpFile, "FT                   /chromosome=\"%s\"\n", cp) ;
  else if (bsGetKey (Seq, _Map, &key) && (obj = bsCreate (key)))
    { if (bsGetData (obj, _EMBL_chromosome, _Text, &cp))
	fprintf (dumpFile, "FT                   /chromosome=\"%s\"\n", cp) ;
      bsDestroy (obj) ;
    }
  if (getDumpText (seq, str2tag("EMBL_map"), &cp))
    fprintf (dumpFile, "FT                   /map=\"%s\"\n", cp) ;

		/* subseqs and features: recurse up through parents */
		/* first obj is Seq */
  key = seq ; obj = bsCreate (key) ;
  dumpKey = seq ; offset = 1 ; offForward = TRUE ;
  isGeneFromName = getDumpFlag (seq, _gene_from_name) ;
  buildAccMap (0,0,0) ;		/* initialise */
  while (TRUE)			/* see below for exit point */
    { KEY new ;
      int x1, x2 ;
      Array aa = arrayHandleCreate (1024, BSunit, handle) ;

				/* do the work */
      if (bsFindTag (obj, _Subsequence) && bsFlatten (obj, 3, aa))
	{ for (i = 0 ; i < arrayMax(aa) ; i += 3)
	    buildAccMap (arr(aa,i,BSunit).k,
			    arr(aa,i+1,BSunit).i,
			    arr(aa,i+2,BSunit).i) ;
	  for (i = 0 ; i < arrayMax(aa) ; i += 3)
	    subseqFeatures (arr(aa,i,BSunit).k,
			    arr(aa,i+1,BSunit).i,
			    arr(aa,i+2,BSunit).i) ;
	}

		/* Features - use marks because there may be a subobj */
      if (bsFindTag (obj, _Feature))
	{ KEY method ;
	  if (bsGetKey (obj, _bsRight, &method)) do {
				/* check if method is dumpable at all */
	    if (methodDumpInfo (method))
	      { mark1 = bsMark (obj, mark1) ;
		if (bsGetData (obj, _bsRight, _Int, &x1)) do {
		  mark2 = bsMark (obj, mark2) ;
		  if (bsGetData (obj, _bsRight, _Int, &x2)) do {
		    mark3 = bsMark (obj, mark3) ;
		  
		    dumpFeatureHomol (obj, TRUE, method, x1, x2, 0, 0, 0, 0) ;

		    bsGoto (obj, mark3) ;
		  } while (bsGetData (obj, _bsDown, _Int, &x2)) ;
		  bsGoto (obj, mark2) ;
		} while (bsGetData (obj, _bsDown, _Int, &x1)) ;
		bsGoto (obj, mark1) ;
	      }
	  } while (bsGetKey (obj, _bsDown, &method)) ;
	}

		/* Homols - like features but worse! */
      if (bsFindTag (obj, _Homol))
	{ KEY target, method ; float score ; int y1, y2 ;
	  if (bsGetKeyTags (obj, _bsRight, 0)) do { /* the tag2 */
	    mark1 = bsMark (obj, mark1) ;
	    if (bsGetKey (obj, _bsRight, &target)) do {
	      mark2 = bsMark (obj, mark2) ;
	      if (bsGetKey (obj, _bsRight, &method)) do {
				/* check if method is dumpable at all */
		if (methodDumpInfo (method))
		  { mark3 = bsMark (obj, mark3) ;
		    if (bsGetData (obj, _bsRight, _Float, &score)) do {
		      mark4 = bsMark (obj, mark4) ;
		      if (bsGetData (obj, _bsRight, _Int, &x1)) do {
			mark5 = bsMark (obj, mark5) ;
			if (bsGetData (obj, _bsRight, _Int, &x2)) do {
			  mark6 = bsMark (obj, mark6) ;
			  if (bsGetData (obj, _bsRight, _Int, &y1)) do {
			    mark7 = bsMark (obj, mark7) ;
			    if (bsGetData (obj, _bsRight, _Int, &y2)) do {
			      mark8 = bsMark (obj, mark8) ;
			      
			      dumpFeatureHomol (obj, FALSE, method, x1, x2, target, score, y1, y2) ;

			      bsGoto (obj, mark8) ;

			    } while (bsGetData (obj, _bsDown, _Int, &y2)) ;
			    bsGoto (obj, mark7) ;
			  } while (bsGetData (obj, _bsDown, _Int, &y1)) ;
			  bsGoto (obj, mark6) ;
			} while (bsGetData (obj, _bsDown, _Int, &x2)) ;
			bsGoto (obj, mark5) ;
		      } while (bsGetData (obj, _bsDown, _Int, &x1)) ;
		      bsGoto (obj, mark4) ;
		    } while (bsGetData (obj, _bsDown, _Float, &score)) ;
		    bsGoto (obj, mark3) ;
		  }
		/* end of the if(methodDumpInfo (method)) */
	      } while (bsGetKey (obj, _bsDown, &method)) ;
	      bsGoto (obj, mark2) ;
	    } while (bsGetKey (obj, _bsDown, &target)) ;
	    bsGoto (obj, mark1) ;
	  } while (bsGetKeyTags (obj, _bsDown, 0)) ;
	}

				/* look for a parent */
      if (!bsGetKey (obj, _Source, &new))
	{ bsDestroy (obj) ; 
	  break ;		/* exit point of loop */
	}
      bsDestroy (obj) ;
      if (!(obj = bsCreate (new)))
	{ if (logFile) ++nLog, fprintf (logFile, " * Expected source not found for %s\n", name(key)) ;
	  break ;
	}
      if (!bsFindKey (obj, _Subsequence, key) ||
	  !bsGetData (obj, _bsRight, _Int, &x1) ||
	  !bsGetData (obj, _bsRight, _Int, &x2))
	{ if (logFile) ++nLog, fprintf (logFile, " * Can't find coord of %s in parent\n", name (key)) ;
	  break ;
	}
      if (x2 > x1)
	offset += x1 - 1 ;
      else
	{ offForward = !offForward ;
	  offset = x1 + 1 - offset ;
	}
      key = new ;
    }
  xxout ;

/* SQ */
  for (i = 5 ; i-- ; )
    freq[i] = 0 ;
  for (i = arrayMax(dna) ; i-- ; )
    switch (arr(dna,i,char))
      {
      case A_: ++freq[0] ; break ;
      case C_: ++freq[1] ; break ;
      case G_: ++freq[2] ; break ;
      case T_: ++freq[3] ; break ;
      default:	++freq[4] ; break ;
      }

  fprintf (dumpFile, "SQ   Sequence  %d BP;   %d A; %d C; %d G; %d %c; %d other;",
	   length, freq[0], freq[1], freq[2], freq[3], 
	   isRNA ? 'U' : 'T', freq[4]) ;

/* the sequence itself! */
  translate = isRNA ? rnaDecodeChar : dnaDecodeChar ;
  for (i = 0 ; i < arrayMax(dna) ; ++i)
    if (!(i%60))
      fprintf (dumpFile, "\n     %c", translate[(int)arr(dna,i,char)]) ;
    else if (!(i%10))
      fprintf (dumpFile, " %c", translate[(int)arr(dna,i,char)]) ;
    else
      fputc (translate[(int)arr(dna,i,char)], dumpFile) ;
  fprintf (dumpFile, "\n//\n") ;

finish:
  if (message)
    messout ("emblDump error with %s: %s", name(seq), message) ;
  bsDestroy (Seq) ;
  arrayDestroy (dna) ;
  handleDestroy (handle) ;
}

/*******************************/

static char dname[DIR_BUFFER_SIZE], fname[FIL_BUFFER_SIZE] ;

void emblDumpKey (KEY key)
{

  /* If the dump is from a display we need to pop down the original dialog.  */
  if (isDumpKeyFromDisplay)
    {
      displayUnBlock() ;
      isDumpKeyFromDisplay = FALSE ;
    }

  handleDestroy (massHandle) ;
  if ((dumpFile = filqueryopen (dname, fname, "embl", "w",
			   "EMBL format file")))
    { logFile = filopen (messprintf ("%s/%s", dname, fname), "emblLog", "w") ;
      nLog = 0 ;

      emblDoDump (key) ;

      if (nLog)
	messout ("EMBL dump: %d errors/warnings written to %s/%s.emblLog\n",
		 nLog, dname, fname) ;
      filclose (dumpFile) ;
      filclose (logFile) ;
    }
}

void emblDumpKeySet (KEYSET kset)
{
  int i ;

  handleDestroy (massHandle) ;
  if ((dumpFile = filqueryopen (dname, fname, "embl", "w",
			   "EMBL format file")))
    { logFile = filopen (messprintf ("%s/%s", dname, fname), "emblLog", "w") ;
      nLog = 0 ;

      for (i = 0 ; i < keySetMax (kset) ; ++i)
	emblDoDump (keySet (kset,i)) ;

      if (nLog)
	messout ("EMBL dump: %d errors/warnings written to %d/%d.emblLog\n",
		 nLog, dname, fname) ;
      filclose (dumpFile) ;
      filclose (logFile) ;
    }
}

void emblDumpKeySetFile (KEYSET kset, char *fileName)
{ 
  int i ;

  handleDestroy (massHandle) ;
  if ((dumpFile = filopen (fileName, "embl", "w")))
    { logFile = filopen (fileName, "emblLog", "w") ;
      nLog = 0 ;

      for (i = 0 ; i < keySetMax (kset) ; ++i)
	emblDoDump (keySet (kset,i)) ;

      filclose (dumpFile) ;
      if (nLog)
	messout ("EMBL dump: %d errors/warnings written to %s\n",
		 nLog, filName (fileName, "emblLog", "a")) ;
      filclose (logFile) ;
    }
}

/*******************************/

void emblDump (void)
{
  isDumpKeyFromDisplay = TRUE ;				    /* Record that this dump is from display. */
  displayBlock (emblDumpKey,
		"I will embl dump the corresponding sequence ") ;
}

/********************** end of file **************************/

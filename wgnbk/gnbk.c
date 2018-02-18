/*  File: gnbkclient.c
 *  Author: Jean Thierry-Mieg
 *  Copyright (C) R Durbin and J Thierry-Mieg, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Modified from netclient to handle communication
 * with a gnbk server
 * 
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 25 21:36 1996 (mieg)
 * Created: Wed Nov 25 20:02:45 1992 (mieg)
 *-------------------------------------------------------------------
 */

#include "regular.h"
#include "dict.h"
#include <sys/types.h>
#include <sys/fcntl.h>
#include "rpcgnbk.h"

/* $Id: gnbk.c,v 1.1.1.1 2002/07/19 20:23:18 sienkiew Exp $ */

/***************************************************************************/

typedef  struct idxstruct { long n1 ; int n2 ; int master ; } IDX ;
static  Array idx = 0 ;
static  DICT* dict = 0 ;
static  Array filNames = 0 ;

/****************************************************************/

void gnbkClose (void)
{ 
  dictDestroy (dict) ;
  arrayDestroy (idx) ;
  arrayDestroy (filNames) ;
}

/****************************************************************/

static Stack html(char type, char* vp)
{ 
  Stack s = stackCreate (1000) ;

  switch (type)
    {
/* here you could process vp and for exach line
   you say 
   case 'h':   html
      cp = "............." ;
     catText (s, cp) ;
     */
    default:   
      catText (s, vp) ;
      break ;
    }

  return s ;
}

/***************************************************************************/

static Stack genBankExport (int master, int ppos, int size)
{ char *vp =0 ;
  Stack hp = 0 ;
  int fd = -1 ,  n, nretry = 0 ;
  long pp, pos = ppos;
  char *filName = master >= 0 && master < arrayMax(filNames) ? 
    arr(filNames, master, char*) : 0 ;

  if (filName)
    fd = open(messprintf("%s.1",filName), O_RDONLY, 0);

  if (fd < 0)
    { fprintf (stderr,"// Cannot open file %s.1 given on command line \n",
	       filName) ;
      return 0 ;
    }
  if (pos < 0)
    { fprintf (stderr,"// Cannot scan %ld as a positive long", pos) ;
    return 0 ;
    }

  if (size < 0)
    { fprintf (stderr, "// Cannot scan %d as a positive int", size) ;
      return 0 ;
    }
  if((pp = lseek(fd, pos, SEEK_SET)) != pos)
    { fprintf (stderr,"// Cannot seek position %d in file %s\n", size, filName) ;
      close (fd) ;
      return 0 ;
    }
  vp = messalloc(size + 1) ;
  while(vp && (n=read(fd,vp,size))!=size)
    if (nretry++ == 5)
      { fprintf (stderr,"// Cannot read %d bytes in file %s\n",  size, filName) ;
        close (fd) ;
        return 0 ;
      }
  vp[size] = 0 ;
  close (fd) ;

  hp = html('h', vp) ;
  messfree (vp) ;
  return hp ;
}


/****************************************************************/

Stack gnbkQuery (char *cp)   /* rpc query */
{ int key = -1 ;
 Stack hp = 0 ;
  IDX *idxp ;
	
  if (dictFind (dict, cp, &key) &&
      key >=0 &&
      key < arrayMax(idx))
    {
      idxp =  arrp (idx, key, IDX) ;
      hp = genBankExport (idxp->master, idxp->n1, idxp->n2) ;
    }
  else
    { 
      hp = stackCreate (200) ;
      pushText (hp, messprintf("// don t know %s\n", cp)) ;
    }
  return hp ;
}
 
/****************************************************************/

static BOOL cleanMasterFile (void)
{
  int ii, nn,  level ;
  char *cp, *filName ;
  FILE *f, *g ;

  ii = arrayMax(filNames) ;
  if (!ii)
    messcrash ("No file names provided") ;
  nn = 0 ;
  while (ii--)
    { nn = 0 ;
      filName = arr(filNames, ii, char*) ;
      f = filopen(filName, "", "r") ;
      g = filopen(filName, "1", "w")  ;
      if (!f)
	messcrash ("Cannot read file %s, sorry\n", filName) ;
      if (!g)
	messcrash ("Cannot write file %s.1, sorry\n", filName) ;

      level = freesetfile (f,"") ;
      freespecial("\n") ;
      while (freecard(level))
	{ cp = freepos() ;
	  nn++ ;
	  fprintf(g, "%s\n", cp) ;
	}
       freeclose(level) ; f = 0 ;
       filclose (g) ;
       
       messout ("Scanned %d lines,\n// %s\n// Result stored in %s.1\n",nn, 
		"Cleaned out non printable characters\n", 
		filName) ;
    }
  return TRUE ;
}

/****************************************************************/

static BOOL makeIndexFile (int phase)
{
  int iii, ii, n1, nl, ne, nby, nbyt, level, dn, n2 ;
  char *filName, *cp, buf[50][255] ;
  FILE *f, *g ;

  if (!arrayExists(filNames) || !arrayMax(filNames))
    messcrash ("No file names provided") ;

  iii = arrayMax(filNames) ;

  nbyt = 0 ;  ne = 0 ; nl = 0 ;
  while (iii--)
    {  
      filName = arr (filNames, iii, char*) ;

      f = filopen(filName, "1", "r") ;
      if (f)
	messout ("Indexing file %s.1, \n", filName) ;
      else
	{
	  f = filopen(filName, 0, "r") ;
	  if (!f)
	    messcrash ("Cannot read file %s nor %s.1, sorry\n",
		       filName, filName) ;
	}
      g = filopen(filName, "idx", "w")  ;
      if (!g)
	messcrash ("Cannot write file %s.idx, sorry\n", filName) ;

      ii = 0 ; n1 = 0 ; nby = 0 ; 
      level = freesetfile (f,"") ;
      freespecial("\n") ;
      while (freecard(level))
	{ cp = freepos() ;
	nl++ ;
	dn = strlen(cp) + 1 ;
	nby += dn ;
	switch (phase)
	  {
	  case 2:  /* creation reelle de l'index */
	    cp = freeword () ;
	    if (!cp) break ;
	    if (!strncmp(cp, "ACCESSION", 9))
	      { ii = 0 ;
	      while (ii < 50 && (cp = freeword()))
		{
		  strncpy (buf[ii++], cp, 255) ;
		}
	      }
	    else if (cp && !strncmp(cp, "//", 2))
	      { 
		n2 = nby - n1 ; ne += ii ;
		while (ii-- > 0)
		  fprintf(g, "%s %d %d\n", buf[ii], n1, n2) ;
		n1 = nby  ;
	      }
	    break ;

	  case 21:  /* creation reelle de l'index */
	    cp = freeword () ;
	    if (!cp) break ;
	    if (!strncmp(cp, "ENTRY", 5))
	      { ii = 0 ;
	      while (ii < 50 && (cp = freeword()))
		{
		  strncpy (buf[ii++], cp, 255) ;
		}
	      }
	    else if (cp && !strncmp(cp, "///", 3))
	      { 
		n2 = nby - n1  ; ne += ii ;
		while (ii-- > 0)
		fprintf(g, "%s %d %d\n", buf[ii], n1, n2) ;
		n1 = nby  ;
	      }
	    break ;
	  }
	}
      freeclose(level) ; f = 0 ;

      n2 = nby - n1 + 3 ; ne += ii ;
      while (ii-- > 0)
	fprintf(g, "%s %d %d\n", buf[ii], n1, n2) ;
      
      fprintf(g,"\n") ;
      messout("//file %s, found  %d entries %d lines %d bytes\n", filName, ne, nl, nby) ;
      filclose (g) ;
      nbyt += nby ;
    }
  messout("//Total, Indexed  %d entries pointing to %d lines, %d bytes\n",
	   ne, nl, nbyt) ;

  return TRUE ;
}

/****************************************************************/

int gnbkMax(void)
{ return dict ? dictMax(dict) : 0 ;
}

static BOOL makeHashTable (void)  /* hashing */
{ 
  int level, ii, n1, n2, key, nn, nnt = 0 ;
  char *cp, *filName ;
  FILE *f ; IDX *idxp ;

  idx = arrayReCreate (idx, 10000, IDX) ;
  dictDestroy (dict) ;
  dict = dictCreate (10000) ;

  ii = arrayMax(filNames) ;
  if (!ii)
    messcrash ("No file names provided") ;

  while (ii--)
    { nn = 0 ;
      filName = arr (filNames, ii, char*) ;
      f = filopen(filName, "idx", "r")  ;
      if (!f)
	messcrash ("Cannot read file %s.idx, sorry\n", filName) ;

      level = freesetfile (f,"") ;
      freespecial("\n") ;
      while (freecard(level))
	{ 
	  cp = freeword() ;

	  if (cp)
	    dictAdd (dict, cp, &key) ;
	  if (!cp || 
	      !freeint(&n1)  ||
	      !freeint(&n2))
	    break ;
	  idxp = arrayp (idx, key, IDX) ;
	  idxp->master = ii ;
	  idxp->n1 = n1 ;
	  idxp->n2 = n2 ;
	  nn++ ;
	}
      messout ("// Hashing index file %s.idx, %d entries ", filName, nn ) ;
      nnt += nn ;
    }
  messout ("// Hashing total %d entries ",  nnt ) ;
  return dictMax(dict) ? TRUE : FALSE ;
}

/****************************************************************/

static void interactiveQuestions (void)  /* fontionnement normal */
{ Stack hp = 0 ;
  int level = freesetfile (stdin,"") ;
  char *cp ; int key ;
  IDX *idxp ;

  fprintf(stderr, "// Ready, %d entries available in hash table\n", dictMax(dict)) ;
  
  freespecial("\n") ;
  while (freecard (level))
    {
      if ((cp = freeword()))
	{ 
	  if (!strcasecmp(cp, "quit") ||
	      !strcasecmp(cp, "exit"))
	    return ;

	  if (dictFind (dict, cp, &key) &&
	      key >=0 &&
	      key < arrayMax(idx)) 
	    {
	      idxp =  arrayp (idx, key, IDX) ;  
	      hp = genBankExport ( idxp->master, idxp->n1, idxp->n2) ;
	      printf(stackText (hp,0)) ; 
	      printf("\n") ; 
	      stackDestroy (hp) ;
	    } 
	  else
	    printf("// don t know %s\n", cp) ;
	}
      fflush (stdout) ;
    }
}

/****************************************************************/

BOOL gnbkInit (int phase, Array fNames)
{ /* register the file names as a static array */
  filNames = fNames ;

  switch (phase)
    {
    case 1: 
      messout("Phase %d, optional clean up of the data file", phase) ;
      return cleanMasterFile () ;
      break ;
    case 2: case 21: /* case 22: and so on */
      messout("Phase %d, creation of the index file", phase) ;
      return makeIndexFile (phase) ;
      break ;
    case 3: 
      if (makeHashTable ())
	interactiveQuestions() ;
      gnbkClose () ;
      break ;
    case 4: case 5:
      return makeHashTable () ;
      /* and wait for rpc connections */
      break ;
    default:
      return FALSE ;
    }
  return FALSE ;
}

/****************************************************************/
#ifdef JUNK
BOOL gnbkInitOld (int phase, Array filNames)
{ char *cp ;
  int i, ii, key, level = 0,  n2 , dn, nl , nidx = 0, master = 0 ;
  int nn, n1, nby ;
  FILE *f = 0, *g = 0 ;
  char buf[50][ 256] ;
  IDX *idxp ;

  /* register the file names as a static array */
  ii = arrayMax(filNames) ;
  while (ii--)
  switch (phase)
    {
    case 1: 
      f = filopen(cp, "", "r") ;
      g = filopen(cp, "1", "w")  ;
      break ;
    case 2: case 21: /* case 22: and so on */
      f = filopen(cp, "1", "r") ;
      g = filopen(cp, "idx", "w")  ;
      break ;
    case 3: case 4: case 5:
      f = filopen(cp, "idx", "r") ;
      idx = arrayCreate (10000, IDX) ;
      dict = dictCreate (10000) ;
      break ;
    default:
      return FALSE ;
    }

  if (!f)
    { messcrash ("Cannot open file %s, sorry\n", cp) ;
    filclose (f) ;
    filclose (g) ;
    return FALSE ;
    }
  if (phase < 3 && !g)
    { messout ("Cannot open file %s, sorry\n", cp) ;
    filclose (f) ;
    filclose (g) ;
    return FALSE ;
    }
  
  ii = 0 ; nn = 0 ; n1 = 0 ; nl = nby = 0 ; 
  level = freesetfile (f,"") ;
  freespecial("\n") ;
  while (freecard(level))
    { cp = freepos() ;
      nl++ ;
      dn = strlen(cp) + 1 ;
      nby += dn ;
      switch (phase)
	{
	case 1: /* nettoyage des char ^M */
	  fprintf(g, "%s\n", cp) ;
	  break ;
	case 2:  /* creation reelle de l'index */
	  cp = freeword () ;
	  if (!cp) break ;
	  if (!strncmp(cp, "ACCESSION", 9))
	    { ii = 0 ;
	    while (ii < 50 && (cp = freeword()))
	      {
		strncpy (buf[ii++], cp, 255) ;
	      }
	    }
	  else if (cp && !strncmp(cp, "//", 2))
	    { 
	      n2 = nn - n1 + 3 ;
	      while (ii-- > 0)
		fprintf(g, "%s %d %d\n", buf[ii], n1, n2) ;
	      n1 = nn + 3 ;
	    }
	  break ;
	case 21:  /* creation reelle de l'index */
	  cp = freeword () ;
	  if (!cp) break ;
	  if (!strncmp(cp, "ENTRY", 5))
	    { ii = 0 ;
	    while (ii < 50 && (cp = freeword()))
	      {
		strncpy (buf[ii++], cp, 255) ;
	      }
	    }
	  else if (cp && !strncmp(cp, "///", 3))
	    { 
	      n2 = nn - n1 + 3 ;
	      while (ii-- > 0)
		fprintf(g, "%s %d %d\n", buf[ii], n1, n2) ;
	      n1 = nn + 3 ;
	    }
	  break ;
	case 3: case 4: case 5: /* hashing */
	  cp = freeword () ;
	  if (cp)
	    dictAdd (dict, cp, &key) ;
	  if (!cp || 
	      !freeint(&n1)  ||
	      !freeint(&n2))
	    break ;
	  idxp = arrayp (idx, key, IDX) ;
	  idxp->n1 = n1 ;
	  idxp->n2 = n2 ;
	  break ;	  
	}
      nn += dn ;
    }
  
  freeclose(level) ; f = 0 ;

  /* fin d'initialisation */
  switch (phase)
    {
    case 2: case 21:
      n2 = nn - n1 + 3 ;
      while (ii-- > 0)
	fprintf(g, "%s %d %d\n", buf[ii], n1, n2) ;
      
      fprintf(g,"\n") ;
       break ;
    case 3:  /* fontionnement normal */
      fprintf(stderr, "// Ready, %d accession indexed\n", dictMax(dict)) ;
      level = freesetfile (stdin,"") ;
      freespecial("\n") ;
      while (freecard (level))
	{ Stack hp = 0 ;
	  if (cp = freeword())
	      { Stack hp = 0 ;

	        if (dictFind (dict, cp, &key) &&
		    key >=0 &&
		    key < arrayMax(idx)) ;
		else
		  continue ;

		idxp =  arrayp (idx, key, IDX) ;  
		hp = genBankExport ( idxp->master, idxp->n1, idxp->n2) ;
		printf(stackText (hp,0)) ; 
		printf("\n") ; 
		stackDestroy (hp) ;
	      }
	  else
	    { 
	      if (!strcasecmp(cp, "quit") ||
		  !strcasecmp(cp, "exit"))
		break ;
	      printf("// don t know %s\n", cp) ;
	    } 
	  fflush (stdout) ;
	}
      gnbkClose () ;
      break ;
    case 4: case 5:
      /* do not destroy dict/idx you will receive questions via rpc */
      break ;
    }

abort:
  filclose (f) ;
  filclose (g) ;
  return TRUE ;
}
#endif

/****************************************************************/
/****************************************************************/


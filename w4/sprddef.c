/*  File: sprdef.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
    Read write to disk
     and or
    Get Save in object
        a table definition
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 17:19 1998 (fw)
 * Created: Fri Dec 11 13:53:37 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: sprddef.c,v 1.20 2016/11/27 00:23:20 mieg Exp $ */

#include "acedb.h"
#include "ac.h"
#include "spread_.h"
#include "session.h"		/* for thisSession */
#include "freeout.h"
#include "java.h"
#include "dump.h"
#include "lex.h"
#include "systags.h"
#include "tags.h"
#include "sysclass.h"

/*****************************************************/
/******* Input Output of the definitions *************/ 

static FREEOPT showOpts[] = {
  {6, "Show options"},
  {1, "MIN"},
  {2, "MAX"},
  {3, "Average"},
  {4, "Variance"},
  {5, "Sum"},
  {6, "Count"}
} ;

/************************************************************/

void spreadDoExportBql (SPREAD spread)
{ 
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (0, 0, 0, h) ;
  int j, kk = 0 , maxCol= arrayMax (spread->colonnes) ;
  COL *c ;
 
  aceOutf (ao, "select") ;
  for (j = 0 ; j < maxCol; j++)
    { 
      c = arrp(spread->colonnes,j, COL) ;
      if (!c->type)
	continue ;
      if (! c->hidden)
	aceOutf (ao, "%s X%d", kk++ ? "," : "", j+1) ;
    }
  aceOutf (ao, " from") ;
  for (kk = j = 0 ; j < maxCol; j++)
    { 
      c = arrp(spread->colonnes,j, COL) ;
      if (!c->type)
	continue ;
      
      aceOutf (ao, "%sX%d"
	       , kk++ ? ", " : " "
	       , j+1
	       ) ;
      if (j == 0)
	{
	  if (c->type == 'k' && c->classp)
	    aceOutf (ao, " in class \"%s\"",c->classp) ;
	}
      else if (c->showType == SHOW_COPY)
	{
	  aceOutf (ao, " = X%d", c->from) ;
	}
      else if (c->showType == SHOW_COMPUTE)
	{
	  char *cr, *cq = hprintf(h, "%s", c->tagp ? c->tagp : "") ;
	  for (cr = cq ; *cr ; cr++)
	    if (*cr == '%' && (cr == cq || *(cr-1) == ' ') && cr[1] >= '0' && cr[1] <= '9')
	      *cr = 'X' ;
	  aceOutf (ao, " = %s", cq) ;
	}
      else if (c->tagp && *c->tagp && (! strcmp (c->tagp, "HERE ") ||  ! strcmp (c->tagp, "HERE # ")))
	{
	  if (c->typep && ! strcmp (c->typep, "Next_Tag"))
	    aceOutf (ao, " in X%d[2]", c->from) ;
	  else if (c->typep)
	    aceOutf (ao, " in X%d[1]", c->from) ;
	}
      else if (c->tagp && *c->tagp && c->typep && ! strcmp (c->typep, "Show_Tag"))
	{
	  aceOutf (ao, " in X%d->%s[0]", c->from, c->tagp) ;
	}
      else if (c->tagp && *c->tagp)
	{
	  char dummy ; 
	  int k = 0 ;
	  char *cr, *cq = hprintf(h, "%s", c->tagp ? c->tagp : "") ;
	  cr = strchr(cq, ':') ;
	  if (cr &&  sscanf (cr+1, "%d%c", &k, &dummy) == 1)
	    {
	      *cr = 0 ;
	      if ( ! strcmp (c->typep, "Next_Tag"))
		k++ ;
	      cq = hprintf(h, "%s[%d]", cq, k) ;
	    }
	  if (! strncmp (cq, "HERE",4))
	    aceOutf (ao, " in X%d[1]", c->from) ;
	  else
	    aceOutf (ao, " in X%d->%s", c->from, cq) ;
	  if (c->typep && k == 0 &&  ! strcmp (c->typep, "Next_Tag"))
	    aceOutf (ao, "[1]") ;	 
	}
      else
	{
	  aceOutf (ao, " BQL definition NOT IMPLEMENTED for This table type") ;
	}

      if (c->mandatory != 1 ||
	  (c->conditionBuffer && *c->conditionBuffer)
	  )
	{
	  aceOutf (ao, " where ") ;
	  if (c->mandatory == 2)
	    aceOutf (ao, "X%d  ", j+1) ;
	  if (c->mandatory == 0)
	    aceOutf (ao, " NOT X%d ", j+1) ;
	  if (c->mandatory != 1 &&
	      (c->conditionBuffer && *c->conditionBuffer)
	      )
	    aceOutf (ao, " AND ( ") ;
	  if  (c->conditionBuffer && *c->conditionBuffer)
	    {
	      char *cr, *cq = hprintf(h, c->conditionBuffer) ;
	      for (cr = cq ; *cr ; cr++)
		if (*cr == '%' && (cr == cq || *(cr-1) == ' ' || *(cr - 1) == '\"') && cr[1] >= '0' && cr[1] <= '9')
		  {
		    *cr = 'X' ;
		    if ( *(cr - 1) == '\"')  /* case "%3" == "a" */
		      {
			*(cr - 1) = ' ' ;
			cr  = strchr(cr,'\"') ;
			if (cr) *cr = ' ' ;
		      }
		  }
	      if (*cq)
		{
		  KEY tag = 0 ;
		  if (lexword2key (cq, &tag, 0) && tag)
		    aceOutf (ao, "  X%d#%s", j+1, name (tag)) ;
		  else
		    aceOutf (ao, " %s", cq) ;
		}
	    }
	  if (c->mandatory != 1 &&
	      (c->conditionBuffer && *c->conditionBuffer)
	      )
	    aceOutf (ao, " )") ;
	}
    }
      
     
  aceOutf (ao, "\n\n") ;
  ac_free (h) ;
} /*  spreadDoExportBql */

/************************************************************/

void spreadDoExportDefinitions (SPREAD spread)
{
  int j , maxCol;
  char cc;
  const char *cp ;
  COL *c ;
  
  maxCol = arrayMax(spread->colonnes) ;
  for(j = 0 ; j < maxCol; j++)
    { c = arrp(spread->colonnes,j, COL) ;
      if (!c->type)
	continue ;
      
      freeOutf("Colonne %d\n", j + 1) ;
      if (*c->subtitleBuffer)
	freeOutf("Subtitle %s\n", c->subtitleBuffer) ;
      if (*c->legendBuffer)
	freeOutf("Legend %s\n", c->legendBuffer) ;
      
      freeOutf("Width %d\n", c->width) ;
      freeOutf("%s\n",c->optionalp) ;
      freeOutf("%s\n",c->hiddenp) ;
      if (c->typep) freeOutf("%s\n",c->typep) ;
      if ((c->type == 'k') && c->classp)
	freeOutf("Class %s\n",c->classp) ;
      if ((c->type == 'K') && c->classp)
	freeOutf("Next_Key %s\n",c->classp) ;
      if (c->type == 'c')
	freeOutf("Count\n") ;
      if (c->showType != SHOW_ALL && c->showType != SHOW_COPY && c->showType != SHOW_MULTI && c->showType != SHOW_DNA && c->showType != SHOW_PEP && c->showType != SHOW_COMPUTE)
	freeOutf ("Extract %s\n", freekey2text (c->showType, showOpts)) ;
      if (c->showType == SHOW_MULTI)
	freeOutf ("MultiData\n") ;
      if (c->showType == SHOW_COMPUTE)
	freeOutf ("Compute\n") ;
      if (j > 0)
	freeOutf("%s %s\n",c->extendp, c->fromBuffer) ;
      if (c->type == 'D') 
	 { 
	   freeOutf("DNA ") ;
	   freeOutf(" %s", freeprotect (c->dna1Buffer)) ;
	   freeOutf(" %s", freeprotect (c->dna2Buffer)) ;
	   freeOutf ("\n") ;
	 }
      if (c->type == 'P') 
	 { 
	   freeOutf("Peptide ") ;
	   freeOutf(" %s", freeprotect (c->dna1Buffer)) ;
	   freeOutf(" %s", freeprotect (c->dna2Buffer)) ;
	   freeOutf ("\n") ;
	 }
      if (c->tagp && *c->tagp)
	{ 
	  freeOutf("Tag ") ;
	  cp = c->tagp ; cp-- ;
	  while (*++cp)
	    { if (*cp == '%') 
		{ if (*(cp + 1) == '%')
		    { /* export as single % */
		      cp++ ; /* i did 2 chars */
		    }
		  else
		    freeOut("\\") ;
		}
	      cc = *cp ;
	      freeOutf ("%c", cc) ;
	    }
	  freeOut ("\n") ;
	}
      if (c->conditionBuffer && *c->conditionBuffer)
	{ 
	  freeOutf("Condition ") ;
	  cp = c->conditionBuffer ; cp-- ;
	  while (*++cp)
	    { if (*cp == '%') 
		{ if (*(cp + 1) == '%')
		    { /* export as single % */
		      cp++ ; /* i did 2 chars */
		    }
		  else
		    freeOut("\\") ;
		}
	      cc = *cp ;
	      freeOutf ("%c", cc) ;
	    }
	  freeOut ("\n") ;
	}
      freeOutf("\n") ;
    }

  return;
} /* spreadDoExportDefinitions */


void spreadDoSaveDefinitions (SPREAD spread, FILE *f)
{
  int level = freeOutSetFile(f) ;
  char *cp ;

  freeOutf("// Spread sheet definition for the ACeDB software \n") ;
  if (*thisSession.name)
    freeOutf("// User: %s\n", thisSession.name) ;
  freeOutf("// Date: %s\n\n", timeShowNow ()) ;
  freeOutf("// %%n (%%%%n in the graphic)  are parameter to be given on the command line in tace\n") ;
  freeOutf("// or by default by the Parameters given in this file\n") ;
  freeOutf("// \\%%n (%%n in the graphic) are substituted by the value of column n at run time\n") ;
  freeOutf("// Line starting with // are ignored, starting with # are comments\n\n") ;
      
  if (stackExists(spread->comments))
    { stackCursor(spread->comments, 0) ;
      while ((cp = stackNextText(spread->comments)))
	freeOutf( "# %s", cp) ;
    }

  if (*spread->titleBuffer)
    freeOutf("Title %s\n\n", spread->titleBuffer) ;
  if (spread->paramBuffer && *spread->paramBuffer)
    freeOutf("Parameters %s\n\n", spread->paramBuffer) ;
  if (spread->sortColonne)
    freeOutf("Sortcolumn %d\n\n", spread->sortColonne) ;
  if (spread->precompute)
    freeOutf("Precompute\n\n") ;

  spreadDoExportDefinitions (spread) ;
  freeOutf(" \n\n// End of these definitions\n") ;
  freeOutClose(level) ;
  filclose(f) ;

  return;
} /* spreadDoSaveDefinitions */

/******************************************************/

static FREEOPT spreadReadMenu[] = 
{ 
  {51, "Spread Definitions"},
  {'#', "#"},
  {'a', "Table"}, /* first line if from object */
  {'T', "Title"},
  {'P', "Parameters"},
  {'s', "Subtitle"},
  {'l', "Legend"},
  {'c', "Colonne"}, 
  {'r', "From"},
  {'R', "Right_of"},
  {'Y', "Copy"},
  {'g', "Tag"},
  {'D', "DNA"},
  {'A', "Peptide"},
  {'v', "Visible"},
  {'h', "Hidden"},
  {'o', "Optional"},
  {'m', "Mandatory"},
  {'N', "Null"},
  {'C', "Condition"},
  {'k', "Class"},
  {'k', "A_Class"},
  {'i', "Integer"},
  {'i', "An_Int"},
  {'f', "Float"},
  {'f', "A_Float"},
  {'d', "Date"},
  {'d', "A_Date"},
  {'t', "Text"},
  {'t', "A_Text"},
  {'b', "Boolean"},
  {'B', "Show_Tag"},
  {'n', "Next_Tag"},
  {'K', "Next_Key"},
  {'S', "Sortcolumn"},
  {'E', "Extract"},
  {'w', "Width"}, 
  {'M', "FlipMap"},
  {'L', "Type"},
  {'L', "Presence"},
  {'L', "Visibility"},
  {'L', "Origin"},
  {'1', "All"},
  {'2', "Min"},
  {'3', "Max"},
  {'4', "Average"},
  {'5', "Variance"},
  {'6', "Sum"},
  {'7', "Count"},
  {'U', "Precompute"},
  {'u', "MultiData"},
  {'p', "Compute"}
} ;

/******************************************************/
/******************************************************/

BOOL spreadDoSaveInObj (SPREAD spread, KEY tableKey)
{
  int ii, j , maxCol ; char *cp ;
  COL *c ;
  KEY tag ;
  char buf[2] ;
  OBJ Table ;
  AC_HANDLE h = ac_new_handle () ;
  Stack s = stackHandleCreate (300, h) ;

  buf[1] = 0 ; /* used in stackCat */
  if (class(tableKey) != _VTable)
    return FALSE ;
  
  if (!(Table = bsUpdate(tableKey)))
    return FALSE ;
  
  while (bsGoto (Table, 0), bsGetKeyTags (Table, _bsRight, 0))
    bsRemove (Table) ;
  
  if (stackExists(spread->comments))
    { stackCursor(spread->comments, 0) ;
    while((cp = stackNextText(spread->comments)))
      bsAddData (Table, _Comment, _Text, cp) ;
    }
  
  if (spread->titleBuffer && *spread->titleBuffer)
    bsAddData (Table, _Title, _Text, spread->titleBuffer) ;

  lexaddkey ("Precompute", &tag, 0) ;
  if (spread->precompute)
    bsAddTag (Table, tag) ;
      
  lexaddkey ("Sortcolumn", &tag, 0) ;
  if (spread->sortColonne)
    bsAddData (Table, tag, _Int, &spread->sortColonne) ;
      
  lexaddkey ("Parameters", &tag, 0) ;
  if (spread->paramBuffer && *spread->paramBuffer)
    bsAddData (Table, tag, _Text, spread->paramBuffer) ;
      
  maxCol = arrayMax(spread->colonnes) ;
  for(j = 0 ; j < maxCol; j++)
    { c = arrp(spread->colonnes,j, COL) ;
      if (!c->type)
	continue ;
      
      ii = j + 1 ;
      lexaddkey ("Colonne", &tag, 0) ;
      bsAddData (Table, tag, _Int, &ii) ;
      bsPushObj(Table) ;
      
      if (c->subtitleBuffer && *c->subtitleBuffer)
	bsAddData (Table, str2tag("Subtitle"), _Text, c->subtitleBuffer) ;
      if (c->legendBuffer && *c->legendBuffer)
	bsAddData (Table, str2tag("Legend"), _Text, c->legendBuffer) ;

      switch (c->extend)
	{
	case 'f':
	  bsAddData (Table, _From, _Int, &c->from) ;
	  break ;
	case 'r':
	  lexaddkey ("Right_of", &tag, 0) ;
	  bsAddData (Table, tag, _Int, &c->from) ;
	  break ;
	case 'c':
	  lexaddkey ("Copy", &tag, 0) ;
	  bsAddData (Table, tag, _Int, &c->from) ;
	  break ;
	}
      lexaddkey ("Tag", &tag, 0) ;
      if (c->tagp && *c->tagp) 
	bsAddData (Table, tag, _Text, c->tagp) ;

      if (c->type == 'D')
	{
	  bsAddTag (Table, _DNA) ;
	  bsAddData (Table, _bsRight, _Text, c->dna1Buffer) ;
	  bsAddData (Table, _bsRight, _Text, c->dna2Buffer) ;
	  bsFindTag (Table, _DNA) ;
	  bsGetData (Table, _bsRight, _Text, &cp) ;
	  bsGetData (Table, _bsRight, _Text, &cp) ;
	}

      if (c->type == 'P')
	{
	  bsAddTag (Table, _Peptide) ;
	  bsAddData (Table, _bsRight, _Text, c->dna1Buffer) ; /* ac_protect (c->dna1Buffer, h)) ; */
	  bsAddData (Table, _bsRight, _Text, c->dna2Buffer) ; /* ac_protect (c->dna2Buffer, h)) ; */
	}

      if (c->hidden)
	bsAddTag (Table, _Hidden) ;
      else
	bsAddTag (Table, _Visible) ;
      lexaddkey ("Width", &tag, 0) ;
      bsAddData (Table, tag, _Int, &c->width) ;
  
      switch (c->mandatory)
	{
	case 0:
	  lexaddkey ("Null", &tag, 0) ;
	  bsAddTag (Table, tag) ;
	  break ;
	case 1:
	  lexaddkey ("Optional", &tag, 0) ;
	  bsAddTag (Table, tag) ;
	  break ;
	case 2:
	  lexaddkey ("Mandatory", &tag, 0) ;
	  bsAddTag (Table, tag) ;
	  break ;
	}

      if (c->conditionBuffer && *c->conditionBuffer)
	{ 
	  lexaddkey ("Condition", &tag, 0) ;
	  if (0)
	    {
	      s = stackReCreate (s,200) ;
	      cp = c->conditionBuffer ; cp-- ;
	      while (*++cp)
		{ 
		  if (*cp == '%') 
		    { if (*(cp + 1) == '%')
		      { /* export as single % */
			cp++ ; /* i did 2 chars */
		      }
		    else
		      catText(s, "%") ; /* double it */
		    }
		  buf[0] = *cp ;
		  catText (s, buf) ;
		}
	      cp = stackText(s, 0) ;
	      bsAddData (Table, tag, _Text, cp) ;
	    }
	  else
	    bsAddData (Table, tag, _Text, c->conditionBuffer) ;
	}

      switch (c->type)
	{
	case 'k':
	  lexaddkey ("A_Class", &tag, 0) ;
	  if (c->classe)
	    bsAddData (Table, tag, _Text, name(c->classe)) ;
	  break ;
	case 'K':
	  lexaddkey ("Next_Key", &tag, 0) ;
	  if (c->classe)
	    bsAddData (Table, tag, _Text, name(c->classe)) ;
	  break ;
	case 'i':
	  lexaddkey ("An_Int", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case 'f':
	  lexaddkey ("A_Float", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case 't':
	  lexaddkey ("A_Text", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case 'd':
	  lexaddkey ("A_Date", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case 'b':
	  lexaddkey ("Boolean", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case 'B':
	  lexaddkey ("Show_Tag", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case 'n':
	  lexaddkey ("Next_Tag", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case 'c':
	  lexaddkey ("Count", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	}
      
      if (c->type == 'c' || c->showType == SHOW_MULTI)
	{
	  switch (c->type == 'c' ? c->realType : c->type)
	    {
	    case 'k':
	      lexaddkey ("A_Class", &tag, 0) ;
	      if (c->classe)
		bsAddData (Table, tag, _Text, name(c->classe)) ;
	      else
		bsAddTag (Table, tag) ;
	      break ;
	    case 'i':
	      lexaddkey ("An_Int", &tag, 0) ;
	      bsAddTag (Table,tag) ;
	      break ;
	    case 'f':
	      lexaddkey ("A_Float", &tag, 0) ;
	      bsAddTag (Table,tag) ;
	      break ;
	    case 't':
	      lexaddkey ("A_Text", &tag, 0) ;
	      bsAddTag (Table,tag) ;
	      break ;
	    case 'd':
	      lexaddkey ("A_Date", &tag, 0) ;
	      bsAddTag (Table,tag) ;
	      break ;
	    case 'b':
	      lexaddkey ("Boolean", &tag, 0) ;
	      bsAddTag (Table,tag) ;
	      break ;
	    }
	}
      
      switch (c->showType)
	{
	case SHOW_ALL: /* this is the default */
	case SHOW_DNA:
	case SHOW_PEP:
	case SHOW_COPY:
	  break ;
	case SHOW_MIN:
	  bsAddTag (Table, _Min) ;
	  break ;
	case SHOW_MAX:
	  bsAddTag (Table, _Max) ;
	  break ;
	case SHOW_SUM:
	  lexaddkey ("Sum", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case SHOW_AVG:
	  lexaddkey ("Average", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case SHOW_VAR:
	  lexaddkey ("Variance", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case SHOW_MULTI:
	  lexaddkey ("MultiData", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	case SHOW_COMPUTE:
	  lexaddkey ("Compute", &tag, 0) ;
	  bsAddTag (Table,tag) ;
	  break ;
	}

      bsGoto (Table, 0) ;
    }

  bsSave (Table) ;
  ac_free (h) ;

  return TRUE ;
} /* spreadDoSaveInObj */

/******************************************************/

static BOOL spreadDoReadFromObj (SPREAD spread, KEY tableKey, const char *parms, BOOL noParms)
{ 
  char *cp, *cq, *buffer = 0, *zp, *zq ;
  int level, n ;
  Stack s = 0, s1 = 0 ;
  OBJ obj ;
  KEY _Parameters ;
  BOOL ok ;

  if (class(tableKey) != _VTable)
    return FALSE ;

  lexaddkey ("Parameters", &_Parameters, 0) ;
  s = stackCreate (1000) ;
  level = freeOutSetStack (s) ;
  dumpKeyBeautifully (tableKey, 'a', 0) ;  /* into s */
  freeOutClose (level) ;

  if (!parms) parms = "" ; cp = 0 ;
  if ((obj = bsCreate (tableKey)))
    { bsGetData (obj, _Parameters, _Text, &cp) ;
      if (cp) cp = strnew(cp, 0) ;
      bsDestroy (obj) ;
    }
  if (cp)
    cq = strnew(messprintf ("%s %s", parms, cp), 0) ;
  else
    cq = strnew(messprintf ("%s", parms), 0) ;

  /******* change the %%n run time parameters to \\\%n computer pleasing syntax */

  n = strlen(stackText(s,0)) ;
  buffer = messalloc(4*n) ;
  zp = buffer ; zq = stackText(s,0) ;
  while (*zq)
    if (!strncmp(zq,"\\%\\%",4))
      { zq +=4 ; /* gobble %% */ 
      *zp++ = '%' ;
      }
    else 
      *zp++ = *zq++ ;
  *zp++ = 0 ; 

  /******* now unprotect to be able later to substitute the parameters */

  s1 = stackCreate (1000) ;
  if(0)
    {
      level = freesettext(buffer, 0) ;
      freespecial ("\n\t\"/%\\") ;  /* special the \ (for \%1 parameters, forbid sub shells */
      
       while( freecard(level))
	{ while ((cp = freeword()))
	  { catText (s1, cp) ; catText (s1, " ") ; }
	catText (s1, "\n") ;
	}
    }
  else
    catText (s1, buffer) ;
  messfree(buffer) ;
      
  /******* now redouble the params ! **************************/

  if (noParms)
    { 
      n = strlen(stackText(s1,0)) ;
      buffer = messalloc(2*n) ;
      zq = stackText(s1,0) ; zp = buffer ;
      while (*zq)
	{
	  if (*zq == '%' && *(zq-1) != '\\')
	    *zp++ = '%' ; /* double it */
	  *zp++ = *zq++ ;
	}
      *zp++ = 0 ; 
      stackClear(s1) ; pushText(s1, buffer) ; messfree(buffer) ;
    }

  /*************** ouf ! *****************/
  
  ok = spreadDoReadDefinitions (spread, 1, 0, s1, cq, noParms) ;

  stackDestroy (s) ; 
  stackDestroy (s1) ; 
  messfree (cp) ;
  messfree (cq) ;
  messfree (buffer) ;

  return ok ;
} /* spreadDoReadFromObj */


BOOL spreadDoReadDefinitions (SPREAD spread, KEY tableKey, 
			      FILE *f, Stack s, const char *parms, 
			      BOOL noParms)
{
  int i , level = 0 ;
  KEY option ;
  char *cp ;
  COL *c ;
  Array cols = 0, t = 0 ;
  BOOL firstCol = TRUE, shift = TRUE ;
  BOOL ret = TRUE ; /*mhmp 26.11.98 */
  Stack s1 = 0, sparm = 0 ;

  if (tableKey > 2)
    return 
      spreadDoReadFromObj (spread, tableKey, parms, noParms) ;

  t = spread->tableau ;
  if (arrayExists(t))
    { i = arrayMax(t) ;
      while(i--) 
	arrayDestroy(arr(t,i,Array)) ;
      arrayDestroy(spread->tableau) ;
    }
  arrayDestroy(spread->flags) ;
  
  cols = spread->colonnes ;
  if (arrayExists(cols))
    { i = arrayMax(cols) ;
      while (i--)
	spreadDestroyCol(arrp(cols,i,COL)) ;
    }
 
  spread->comments = stackReCreate(spread->comments, 120) ;
  spread->colonnes = arrayReCreate(spread->colonnes, 512, COL) ; 
  spread->modified = TRUE ;
  /* mhmp 09.02.99 1 ---> 0 */
  spread->sortColonne = 1 ; /* default */
  memset (spread->titleBuffer, 0, 60) ;
  memset (spread->paramBuffer, 0, 180) ;
  /* read once to extract the parameters line
     then once with the params postpended 
     */

  if (f)
    level = freesetfile(f, 0) ;
  else if (stackExists(s))
    level = freesettext(stackText(s,0), 0) ;
  else
    messcrash("spreadDoReadDefinitions received neither a proper file nor a proper stack") ;
  
  
  s1 = stackCreate (2000) ; sparm = stackCreate(200) ; 
  if (parms && *parms) pushText(sparm, parms) ;

  freespecial ("\n\"/") ;  /* do not touch, just extract parameters */
  while( freecard(level))
    {
      cp = freepos() ;
      if (cp && *cp) 
	catText (s1, cp) ;
      catText (s1, "\n") ;
      if (freekey (&option, spreadReadMenu))
	switch (option)
	  {
	  case 'P':
	    if ((cp = freepos()) && *cp)
	      {
		if (tableKey)
		  cp = freeunprotect(cp) ;
		catText (sparm, "  ") ; catText (sparm, cp) ;
	      }
	  }
      else if (freeword()) /*mhmp 26.11.98 */
	{
	  ret = FALSE ;
	  goto fin ;
	}
    }


  level = freesettext(stackText(s1,0), stackText (sparm, 0)) ;
  if (noParms)
    freespecial ("\n\t\"/@\\") ;  /* special the \ (for \%1 parameters, forbid sub shells */
  else
    freespecial ("\n\t\"/@%\\") ;  /* special the \ (for \%1 parameters, forbid sub shells */

  c = arrayp(spread->colonnes,0, COL) ; 
  while( freecard(level))
    { 
    lao:
      if (freekey (&option, spreadReadMenu))
	switch (option)
	  {
	  case '#':
	    pushText(spread->comments, freepos()) ;
	    break ;
	  case 'U': /* precompute , do nothing except from command.c */
	    spread->precompute = TRUE ;
	    break ;
	  case 'c':
	    if (freeint(&i))
	      { if (firstCol && !i)
		  shift = FALSE ; /* because in acedb.1.x i saved i, now i+1 */
		if (shift)
		  i-- ;
		c = arrayp(spread->colonnes,i, COL) ; 
		if (!c->colonne && !c->from)
		  { c->colonne = i ;
		    c->width = 12 ;
		    c->mandatory = 1 ; /* optional */
		    c->hidden = FALSE ;
		    c->showType = SHOW_ALL ;
		    c->type = 0 ;
		    c->from = 1 ; /* in col 0, used as flag */
		    c->nonLocal = FALSE ;
		  }
		if (tableKey) /* reading from object */
		  goto lao ; /* loop on the inside of the def */
	      }
	    break ;
	  case 'E':
	    if (freekey (&option, showOpts) && option >= 1 && option <= 6)
	      { 
		if (option < 6)
		  c->showType = (SpShowType) ( option ) ;
		else
		  c->type = 'c' ;
	        c->nonLocal = TRUE ;
	      }
	    break ;
	  case 'w':
	    if (freeint(&i))
	      { c->width = i ;
		strcpy (c->widthBuffer, messprintf ("%d", i)) ; 
	      }
	    break ;
	  case 'N':  /* null */
	    c->mandatory = 0 ;
	    break ;
	  case 'o':  /* optional */
	    c->mandatory = 1 ;
	    break ;
	  case 'm':  /* mandatory */
	     c->mandatory = 2 ;
	    break ;
	  case 'M':
	    c->flip = TRUE ;
	    break ;
	  case 'v':
	    c->hidden = FALSE ;
	    break ;
	  case 'h':
	    c->hidden = TRUE ;
	    break ;
	  case 'i': case 'f': case 'd': case 't':
	    c->type = c->realType = option ; 
	    c->nonLocal = FALSE ;
	    break ;
	  case 'k': case 'K': 
	    c->type = c->realType = option ; 
	    c->nonLocal = FALSE ;
	    if ((cp = freepos()) && *cp)
	      { if (tableKey)
		cp = freeunprotect(cp) ;
	      lexword2key(cp, &c->classe,_VClass) ;
	      }		  
	    break ;
	  case 'b': case 'B':
	  case 'n':  
	    if (option == 'B') option = 'b' ;
	    c->type = option ; 
	    c->realType = 'b' ;
	    c->classe = _VSystem ;
	    c->nonLocal = FALSE ;
	    break ;
	  case 'r':
	    c->extend = 'f' ;
	    if (freeint(&i))
	      c->from = i ;
	    break ;
	  case 'R':
	    c->extend = 'r' ;
	    if (freeint(&i))
	      c->from = i ;
	    break ;
	  case 'Y':
	    c->extend = 'c' ;
	    c->showType = SHOW_COPY ;
	    if (freeint(&i))
	      c->from = i ;
	    break ;
	  case 'p':
	    c->extend = 'p' ;
	    c->showType = SHOW_COMPUTE ;
	    c->nonLocal = TRUE ; 
	    c->type = c->realType = 'f' ; 
	    break ;
	  case 'g':
	    if ((cp = freepos()) && *cp)
	      { 
		char *cq, *cq0 ;
		c->tagStack = stackReCreate(c->tagStack, 30) ;
		if (tableKey)
		  pushText(c->tagStack, freeunprotect(cp)) ;
		else
		  {
		    if (0) cp = freeunprotect(cp) ;
		    cq0 = cq = messalloc (2*strlen(cp)) ;
		    while (*cp)
		      { 
			if (*cp == '\\' &&
			    *(cp + 1) == '%')
			  { cp++ ; *cq++ = *cp++ ; } /* unprotect */
			else if (*cp == '%')
			  *cq++ = '%' ;  /* double it  */
			*cq++ = *cp++ ;
		      }
		    *cq = 0 ;
		    pushText(c->tagStack, cq0) ;
		  }
		c->tagp = stackText(c->tagStack, 0) ;
	      }	    
	    break ;
	  case 'D': /* DNA */
	    c->type = 'D' ;
	    c->showType = SHOW_DNA ;
	    if ((cp = freeword ()))
	      {
		strncpy (c->dna1Buffer, freeunprotect(cp), 255) ;
		if ((cp = freeword ()))
		  strncpy (c->dna2Buffer, freeunprotect(cp), 255) ;
	      }	    
	    break ;
	  case 'A': /* Peptide */
	    c->type = 'P' ;
	    c->showType = SHOW_PEP ;
	    if ((cp = freeword ()))
	      {
		strncpy (c->dna1Buffer, freeunprotect(cp), 255) ;
		if ((cp = freeword ()))
		  strncpy (c->conditionBuffer, freeunprotect(cp), 255) ;
	      }	    
	    break ;
	  case 'C':
	    if ((cp = freepos()) && *cp)
/*	      strncpy(c->conditionBuffer, cp, 359) ; */
	      { char *cq = c->conditionBuffer ;
		int i = 0 ;

		if (tableKey)
		  { cp = freeunprotect(cp) ;
		    i = 358 ; while (i-- && (*cq++ = *cp++)) ;
		  }
		else
		  while (i++ < 359 && *cp)
		    { if (*cp == '\\' &&
			  *(cp + 1) == '%')
		      { cp++ ; *cq++ = *cp++ ; } /* unprotect */
		    else if (*cp == '%')
		      *cq++ = '%' ;  /* double it  */
		    *cq++ = *cp++ ;
		    }
		*cq = 0 ;
		while (cq >  c->conditionBuffer &&
		       (*(cq - 1) == ' ' || *(cq - 1) == '\n'))
		  *cq-- = 0 ;
	      }
	    break ;
	  case 's': 
	    if ((cp = freepos()) && *cp)
	      {
		if (tableKey)
		      cp = freeunprotect(cp) ;
		strncpy(c->subtitleBuffer, cp, 59) ;
	      }
	    break ;
	  case 'l': 
	    if ((cp = freepos()) && *cp)
	      {
		if (tableKey)
		  cp = freeunprotect(cp) ;
		strncpy(c->legendBuffer, cp, 1023) ;
	      }
	    break ;
	  case 'T': 
	    if ((cp = freepos()) && *cp)
	      {
		if (tableKey)
		  cp = freeunprotect(cp) ;
		strncpy(spread->titleBuffer, cp, 59) ;
	      }
	    break ;
	  case 'S':
	    freeint (&spread->sortColonne) ;
	    break ;
	  case 'L': /* intermediate tags, loop */
	    goto lao ;
	    break ;
	  case 'P':
	    if ((cp = freepos()) && *cp)
	      {
		if (tableKey)
		  cp = freeunprotect(cp) ;
		strncpy(spread->paramBuffer, cp, 179) ;
	      }
	    break ;
	  case 'a': /* hopefully TableKey = 1 */
	    break ;
	  case '1': 
	    c->showType = SHOW_ALL ;
	    c->nonLocal = FALSE ;
	    break ;
	  case '2': 
	    c->showType = SHOW_MIN ;
	    c->nonLocal = TRUE ;
	    break ;
	  case '3': 
	    c->showType = SHOW_MAX ;
	    c->nonLocal = TRUE ;
	    break ;
	  case '4': 
	    c->showType = SHOW_AVG ;
	    c->nonLocal = TRUE ;
	    break ;
	  case '5': 
	    c->showType = SHOW_VAR ;
	    c->nonLocal = TRUE ;
	    break ;
	  case '6': 
	    c->showType = SHOW_SUM ;
	    c->nonLocal = TRUE ;
	    break ;
	  case '7': 
	    c->type = 'c' ;  c->showType = SHOW_ALL ;
	    c->nonLocal = TRUE ;
	    break ;
	  case 'u': 
	    c->showType = SHOW_MULTI ; 
	    c->nonLocal = TRUE ;
	    break ; 
	  default:
	    freeOut ("unknown option in spreadReadDefinitions") ;
	  }
      else if (freeword())/*mhmp 26.11.98 */ 
	{
	  ret = FALSE ;
	  freeclose (level) ;
	  goto fin ;
	}
    }
  spread->pos2col = arrayReCreate(spread->pos2col, 8, int) ;
  for (i = 0 ; i < arrayMax(spread->colonnes); i++)
       array(spread->pos2col,i, int) = i ;
fin: /* mhmp 26.11.98 */
  stackDestroy (sparm) ;
  stackDestroy (s1) ;
  return ret ;
} /* spreadDoReadDefinitions */

/******************************************************/
/******************************************************/

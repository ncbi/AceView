/*  File: aql.c
 *  Author: Richard Durbin (Durbin@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1997
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: aql.c,v 1.9 2014/01/14 23:53:53 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
<<<<<<< aql.c
<<<<<<< aql.c
 * Last edited: May 23 13:24 2001 (edgrif)
=======
 * Last edited: May 23 13:44 2001 (edgrif)
>>>>>>> 1.31
=======
 * Last edited: Sep  2 15:15 2004 (edgrif)
>>>>>>> 1.36
 * * Aug 14 12:25 1998 (fw): finished context-var stuff
 * * Aug  3 15:00 1998 (fw): use of debugLevel
 * * Jul 31 14:59 1998 (fw): accept context variable values
 * * Jul  8 17:23 1998 (fw): completed as a test platform.
 * Created: Thu Aug  7 08:41:22 1997 (Durbin)
 *-------------------------------------------------------------------
 */

/**********************************************************************/

#include "acedb.h"
#include "aceio.h"
#include "waql/aql_.h"

/******************************************************************/

/* This code makes use of the setjmp/longjmp recovery method via aqlError(), but is this
 * a good idea because this routine is the constructor and perhaps it should either succeed
 * or return NULL, not longjmp and return an aql object...??? */
/* mieg: 2007_05_14 passed with argc = -1, will loop as it used to, with argc>=0, looping stops after argc calls */

static AQL   aqlDoCreate (const char *queryText,
		   ACEOUT dump_out,
		   char beauty,
		   int parseDebugLevel,
		   int evalDebugLevel, 
		   int argc, 
		   va_list args)
{
  AQL	         aql;
  const char		*varName;
  int		 dictIndex=0;
  AqlLoc	*varLoc;

  if (!queryText)
    return NULL ;

  /* allocate storage for this object - cleaned up in aqlDestroy */
  aql = messalloc(sizeof(struct AqlStruct));

  /* initialise aql object */
  aql->IsError = FALSE;
  aql->errorNumber = 0;
  aql->parseDebugLevel = parseDebugLevel;
  aql->evalDebugLevel = evalDebugLevel;

  aql->dump_out = dump_out;
  aql->beauty = beauty;
  
  aql->result_keyset = NULL;
  aql->result_keyset_handle = 0;
  
  /* this handle is used for the ->query object */
  aql->handle = handleCreate();

  aql->debug_out = aceOutCreateToStderr (aql->handle);

  /* creates the query object with a cleaned up query string */
  aql->query = aqlQueryCreate (queryText, aql->handle);

  aql->errorJmpEnv = (jmp_buf*)halloc (sizeof(jmp_buf), aql->handle);

  /*************************************************/

  /* exception handler setup */
  if (setjmp (*(aql->errorJmpEnv)))
    { 
      /* coming back with an error from somewhere in this routine, aqlParse() or aqlCheck1(),
       * via aqlError() which actually does the longjmp() call. */
      if (aql->parseDebugLevel >= 2)
	{
	  aqlNodeOut (aql->debug_out, aql->query->root, 1, "query after preprocess-error", FALSE);
	  aceOut (aql->debug_out, "\n\n");
	}

      /* clear the memory, the query has used up so far */
      handleDestroy (aql->query->handle) ;

      return aql;
    }

  /* lex/yacc parsing of the string into basic parsetree */
  aqlParse (aql) ;

  if (aql->parseDebugLevel >= 2)
    {
      aqlNodeOut (aql->debug_out, aql->query->root, 0, "parsed query", FALSE);
      aceOut (aql->debug_out, "\n\n");
    }

  /* insert the names of the context variables into the scopes of the query */
  aql->query->tableScope = aqlScopeCreate (aql->query->handle);
  aql->query->scalarScope = aqlScopeCreate (aql->query->handle);
  /* passed with argc = -1, will loop as it used to, with argc>=0, looping stops after argc calls */
  while (argc-- && (varName = va_arg(args, const char*)) != (char *) 0 && varName != (char *)0x100000000 )
    {
      /* create a locator, the value and type will be set 
	 when the actual value is passed to aqlExecute */
      varLoc = (AqlLoc*) halloc (sizeof(AqlLoc), aql->query->handle) ;

      if (varName[0] == '@')
	{
	  varName++;		/* skip the @ */

	  if (!(dictAdd(aql->query->tableScope->dict, varName, &dictIndex)))
	    aqlError (aql, 850, 0, "Multiple declaration of table variable %s in context scope", varName);
	  
	  /* add the loc to the tableScope */
	  array (aql->query->tableScope->loc, dictIndex, AqlLoc*) = varLoc ;
	}
      else if (varName[0] == '$')
	{
	  varName++;		/* skip the $ */

	  if (!(dictAdd(aql->query->scalarScope->dict, varName, &dictIndex)))
	    aqlError (aql, 852, 0, "Multiple declaration of scalar variable %s in context scope", varName);

	  /* add the loc to the scalarScope */
	  array (aql->query->scalarScope->loc, dictIndex, AqlLoc*) = varLoc ;
	}
      else
	{
	  aqlError (aql, 854, 0, "Invalid context variable %s (has to start with @ or $)", varName) ;
	}

      /* we need to signal to aqlCheck1(), that this is a context-var 
	 and isn't defined in the query itself, the checking then becomes less fussy */
      varLoc->isContext = TRUE;
    }


 /* first pass query processing into executable form */
  aqlCheck1(aql) ;
  if (aql->parseDebugLevel >= 2)
    {
      aceOut (aql->debug_out, "\n\n");
      aqlNodeOut (aql->debug_out, aql->query->root, 0, "processed after check1", FALSE);
      aceOut (aql->debug_out, "\n\n");
    }


  return aql ;
} /* aqlDoCreate */

/* i split into 2 functions becuase on 64 bits machine we seem to need to provbide the
 * number of args in the query, otherwise we crash, i cant find why
 */

AQL   aqlCreate (const char *queryText,
                 ACEOUT dump_out,
                 char beauty,
                 int parseDebugLevel,
                 int evalDebugLevel, ...)
{
  AQL	         aql ;
  va_list	 args;
  va_start (args, evalDebugLevel);
  aql = aqlDoCreate (queryText, dump_out, beauty, parseDebugLevel, evalDebugLevel, -1, args) ;
  va_end (args);
  return aql ;
}

AQL   aqlCreate2 (const char *queryText,
		  ACEOUT dump_out,
		  char beauty,
		  int parseDebugLevel,
		  int evalDebugLevel,
		  int argc, ...)
{
  AQL	         aql ;
  va_list	 args;
  va_start (args, argc);
  aql = aqlDoCreate (queryText, dump_out, beauty, parseDebugLevel, evalDebugLevel, argc, args) ;
  va_end (args);
  return aql ;
}

/**********************************************************************/


/* if tableHandle is a valid handle and non-null,
 * the returned table is created on that handle and can be cleared by destroying that handle
 * if tableHandle is NULL, the table has to be explicitly destroyed
 * using tableDestroy() after the table generated by this call
 * is no longer needed
 */
TABLE* aqlExecute(AQL aql,
                  AC_HANDLE tableHandle,
                  KEYSET *result_keyset,
                  AC_HANDLE result_keyset_handle,
                  ...)
{
  va_list args;
  int     dictIndex;
  AqlLoc *varLoc;
  char   *param;
  char   *varType = NULL ;
  char   *varName;
  char   *cp;
  TABLE  *resultTable;
  TABLE  *tmpTable;
  BOOL preprocess_error = FALSE ;

  if (!aql)
    return (TABLE*)0;

  if (aql->IsError)
    return (TABLE*)0;		/* do nothing if pre-processing failed already */
  
  aql->result_keyset = result_keyset;
  aql->result_keyset_handle = result_keyset_handle;
    
  /* Process the argument list to this function -
     We can call this function to bind values to 
     context variables that have been registered in aqlCreate()
     Such argument come in PAIRS of the following form :
     "Table @t", TABLE *t
     "Int $i", int i
     "Float $f", float i
     "Text $s", char *s
     "DateType $t", mytime_t t

     NOTE : the last parameter has to be a null pointer
   */
  va_start (args, result_keyset_handle);

  while ((param = va_arg(args, char*)), param != (char *) 0 && param != (char *)0x100000000 )
    {
      varType = strdup(param); /* make our own copy to mess around with */

      for (cp = varType; *cp && *cp != ' '; cp++); /* find the space */
      if (!*cp)			/* ran past end of string, i.e. no space in it */
	{
	  messerror ("aql.c:aqlExecute() - passed invalid context var string: \"%s\"", varType) ;
	  preprocess_error = TRUE ;
	  break ;
	}


      *cp = '\0';		/* terminates varType at the space char */
      varName = ++cp;		/* start the name after the space */

      if (varName[0] == '@' || varName[0] == '$')
	varName++;		/* skip leading @ or $ in name */

      /* we've already created the AqlLoc's for these variables, we now
	 set their values */

      if (strcasecmp(varType, "Table")==0)
	{
	  if (!(dictFind(aql->query->tableScope->dict, varName, &dictIndex)))
	    {
	      messerror ("aql.c:aqlExecute() - context variable binding of "
			 "undeclared table var: \"%s\"", varName);
	      preprocess_error = TRUE ;
	      break ;
	    }


	  varLoc = arr (aql->query->tableScope->loc, dictIndex, AqlLoc*);
	  varLoc->type = 'T';
	  varLoc->value.T = va_arg(args, TABLE*);
	}
      else
	/* all other value types use the scalar scope */
	{
	  if (!(dictFind(aql->query->scalarScope->dict, varName, &dictIndex)))
	    {
	      messerror ("aql.c:aqlExecute() - context variable binding of "
			 "undeclared scalar var: \"%s\"", varName) ;
	      preprocess_error = TRUE ;
	      break ;
	    }
	  
	  varLoc = arr (aql->query->scalarScope->loc, dictIndex, AqlLoc*);

	  if (strncasecmp(varType, "Int", 3)==0) /* also accepts Integer */
	    {
	      varLoc->type = 'i';
	      varLoc->value.i = va_arg(args, int);
	    }
	  else if (strcasecmp(varType, "Float")==0)
	    {
	      varLoc->type = 'f';
	      varLoc->value.f = (float)va_arg(args, double);
	    }
	  else if (strncasecmp(varType, "DateType", 4)==0) /* allow "DateType" and "Date" */
	    {
	      varLoc->type = 't';
	      varLoc->value.t = va_arg(args, mytime_t);
	    }
	  else if ((strcasecmp(varType, "Text")==0) || 
		   (strcasecmp(varType, "String")==0))
	    {
	      varLoc->type = 's';
	      varLoc->value.s = strnew (va_arg(args, char*), aql->query->handle);
	    }
	  else
	    {
	      messerror ("aql.c:aqlExecute() - invalid value type of "
			 "context variable: \"%s\"", varType) ;
	      preprocess_error = TRUE ;
	      break ;
	    }
	}

      free(varType);
    }

  va_end (args);


  /* A bit messy but we can't use the longjmp stuff yet.... */
  if (preprocess_error)
    {
      if (varType)
	free(varType);

      return NULL ;
    }


  /*
   * exception handler setup, from here we can use longjmp for errors.
   */
  if (setjmp (*(aql->errorJmpEnv)))	/* coming back with an error from somewhere in aqlCheck2() */
    { 
      if (aql->errorNumber == 0) /* interrupted... */
	{
	  AqlNode *tableNode = aql->query->root;
	  AqlNode *resultNode = NULL;
	  TABLE *t = NULL;

	  while (tableNode)
	    {
	      resultNode = tableNode;
	      
	      tableNode = tableNode->nxt;
	    }
	  if (resultNode && resultNode->value.T)
	    t = tableCopy(resultNode->value.T, tableHandle);

	  /* clear the memory, the query has used up so far */
	  handleDestroy (aql->query->handle) ;

	  return t;
	}
      else
	{
	  /* clear the memory, the query has used up so far */
	  handleDestroy (aql->query->handle) ;

	  if (aql->parseDebugLevel >= 2)
	    {
	      aqlNodeOut (aql->debug_out, aql->query->root, 1, "query after execution-error", FALSE);
	      aceOut (aql->debug_out, "\n\n");
	    }

	  return (TABLE*)0;
	}
    }

  /* database specific query processing */
  aqlCheck2 (aql) ;
  if (aql->parseDebugLevel >= 2)
    {
      aceOut (aql->debug_out, "\n\n");
      aqlNodeOut (aql->debug_out, aql->query->root, 0, "processed after check2", FALSE);
      aceOut (aql->debug_out, "\n\n");
    }

  if (aql->parseDebugLevel >= 1)
    {
      /* pretty parsetree print after check2 */
      aqlQueryOut(aql, aql->debug_out, aql->query);
    }


  /* run the query against the database */
  tmpTable = aqlEval (aql) ;

  /* the tmpTable has been created on the query's internal store-handle
     and will be destroyed soon, so we must make a copy to return to the caller */

  resultTable = tableCopy(tmpTable, tableHandle); /* no harm if tableHandle==NULL */

  /* the results table was created on the query->resultHandle store
     and therefore won't vanish, but we clean up all the other space
     allocated during processing */
  handleDestroy (aql->query->handle) ;


  return resultTable;
} /* aqlExecute */


/******************************************************************/


void   aqlDestroy(AQL aql)
{
  if (!aql)
    return;

  /* we assume query->handle has been cleared already */


  /* clear all memory within this object */
  handleDestroy (aql->handle);

  /* kill the actual object-structure */
  messfree(aql);

  return;
} /* aqlDestroy */


/******************************************************************/


TABLE* aqlTable (char *queryText, ACEOUT error_out, AC_HANDLE handle)
     /* simple wrapper-function to create a table from a query text
      * errors are written to error_out, which may be NULL */
{
  AQL aql = aqlDoCreate (queryText, 0, 0, 0, 0, 0, 0);
  TABLE *result = aqlExecute(aql, handle, 0, 0, 0);

  if (aqlIsError(aql))
    {
      if (error_out)
	aceOutf (error_out, "%s", aqlGetErrorReport(aql));
      result = (TABLE*)0;
    }
  
  aqlDestroy (aql) ;

  return result;
} /* aqlTable */

void aqlQuick (char *queryText)
     /* quick all-in-one function */
{
  ACEOUT fo = aceOutCreateToStdout (0);
  TABLE *result = aqlTable(queryText, fo, 0);

  if (result)
    {
      tableOut (result, '\t', 'a');
      tableDestroy (result);
    }

  aceOutDestroy (fo);

  return;
} /* aqlQuick */


/****************************************************/
/*			                 	    */
/* Access to error information in the AQL object    */
/*  after execution, the IsError should be checked, */
/*  the program can then decide to enquire further  */
/*  information about the error and display it	    */
/*  in any way it likes.			    */
/*			                 	    */
/****************************************************/
BOOL  aqlIsError         (AQL aql)
{
  if (aql)
    return aql->IsError;
  
  return TRUE;
}

int   aqlGetErrorNumber  (AQL aql)
{
  if (aql)
    return aql->errorNumber;

  return 0;
}

char *aqlGetErrorMessage (AQL aql)
{
  if (aql)
    return aql->errorMessage;
  
  return NULL;
}

char *aqlGetErrorReport  (AQL aql)
{
  if (aql)
    return aql->errorReport;
  
  return NULL;
}

/******************* end of file *****************************/

/*  File: aqldebug.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
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
 * SCCS: $Id: aqldebug.c,v 1.6 2007/03/16 18:27:41 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Sep  2 15:25 2004 (edgrif)
 * * Nov  9 12:05 1999 (fw): converted to using ACEOUT
 * * Sep  7 17:12 1999 (fw): added IS_EXPR loc-type
 * * Mar 22 15:27 1999 (fw): changed debug-names for operators o reflect
 *              their original name in the query.
 * * Aug 14 18:13 1998 (fw): added LIKE wildcard operator
 * * Aug 14 15:07 1998 (fw): added the bit for LOCAL_FIELD_NAME in prettyOut
 * * Aug 12 10:08 1998 (fw): added abs() ExprFunc
 * * Jul 27 11:05 1998 (fw): added PrettyOut bit for object constructor
 * * Jul 24 11:41 1998 (fw): added tableFunc column selection code
 * * Jul 15 15:41 1998 (fw): changed varname, as tmp-vars start with '_' now
 * * Jul 10 10:21 1998 (fw): fw's version after cleanup, 
                             but still pretty much the same
 * Created: Fri Oct 25 16:04:41 1996 (rd)
 *-------------------------------------------------------------------
 */

#include "acedb.h"
#include "aceio.h"
#include "waql/aql_.h"

/************************************************************/

static void indent (ACEOUT debug_out, int ival);
static const char *varName (AqlNode *inNode);
static int nestingLevel (AqlNode *inNode);

/******************************************************************
 ************************* public functions ***********************
 ******************************************************************/

/********** print out query structures, for debugging ***********/

/* the following three string tables have to 
   correspond exactly to enum lists in aql_.h */

/********************* enum types for node->type ********************/

char* aqlNodeTypeName (AqlNodeType inType)
{
  char* aqlNodeTypeNameTable[] = 
  {
    "**type unset**", 
    "TABLE_ASSIGN", 
    "VAR_ASSIGN",		/* added (fw-980716) */
    "TABLE_VAR", 
    "TABLE_OP",
    "TABLE_SFW", 
    "TABLE_SFW_ALL", 
    "TABLE_ORDER", 
    "SORT_FIELD", 
    "SORT_FIELD_NAME", 
    "SORT_FIELD_EXPR", 
    "SORT_QUALIFIER",
    "SELECT_FIELD", 
    "FROM_LOC", 
    "FROM_TABLE",
    "CLASS", 
    "OBJECT",
    "LOC_VAR",
    "LOC_TABLE_FIELD", 
    "LOC_TABLE_FIELD_NAME",
    "LOC_LOCAL_TAG", 
    "LOC_LOCAL_TAG_NAME", 
    "LOC_LOCAL_POS", 
    "LOC_LOCAL_FIELD", 
    "LOC_LOCAL_FIELD_NAME",
    "LOC_FOLLOW_TAG", 
    "LOC_FOLLOW_TAG_NAME", 
    "LOC_FOLLOW_POS",
    "LOC_METHOD",		/* added (fw-980817) to replace EXPR_LOC_FUNC */
    "EXPR_TABLE_FUNC", 
    "EXPR_TABLE_FUNC_FIELD",	/* added (fw-980723) */
    "EXPR_TABLE_FUNC_FIELD_NAME", /* added (fw-980723) */
    "EXPR_OP",
    "BOOL_EXISTS", 
    "BOOL_EXISTS_TAG", 
    "BOOL_COMPARISON", 
    "BOOL_NOT", 
    "BOOL_OP",
    "VAR",			/* added (fw-980716) */
    "TEXT", 
    "INT",
    "FLOAT",
    "DATE", 
    "KEY",
    "BOOL"
  } ;
  return aqlNodeTypeNameTable[inType];
} /* aqlNodeTypeName */


/******************* operators for node->op ***************/

char* aqlOpTypeName (AqlOpType inType)
{
  char* aqlOpTypeNameTable[] = 
  {
    "**op unset**", 
    "union", 
    "intersect", 
    "diff",
    "asc", 
    "desc",
    "=", 
    "!=", 
    "like",
    "<", 
    ">", 
    "<=", 
    ">=",
    "not", 
    "or", 
    "and", 
    "xor",
    "yeardiff", 
    "monthdiff", 
    "weekdiff", 
    "daydiff", 
    "hourdiff", 
    "mindiff", 
    "secdiff",
    "-",			/* UMINUS */
    "+", 
    "-", 
    "*", 
    "/",
    "count",
    "first",
    "last",
    "sum",
    "min",
    "max",
    "avg",
    "abs",
    "modulo"
  } ;
  return aqlOpTypeNameTable[inType];
} /* aqlOpTypeName */

/************** enum source types for node->loc->source ********/
char* aqlLocSourceTypeName (AqlLocSourceType inType)
{
  char* aqlLocSourceTypeNameTable[] = 
  { 
    "**src unset**", 
    "IS_CLASS",
    "IS_ROW", 
    "IS_COLUMN", 
    "IS_OBJECT", 
    "IS_LOCAL_TAG", 
    "IS_FOLLOW_TAG", 
    "IS_LOCAL_POS", 
    "IS_FOLLOW_POS",
    "IS_METHOD",
    "IS_EXPR"
  };
  return aqlLocSourceTypeNameTable[inType];
} /* aqlLocSourceTypeName */

static int lisply = 1;

/************ Query **************/

void aqlQueryOut (AQL aql, ACEOUT debug_out, AqlQuery *inQuery)
{ 
  if (inQuery) 
    {
      aceOutf (debug_out, "\ncanonical query :\n\n ");
#if 1
      aqlNodeOut(debug_out,
                 inQuery->root,
                 0,
                 "canonical",
                 FALSE);
#endif
      aqlNodeOutPretty (aql, debug_out, inQuery->root, FALSE, 1) ;
      aceOutf (debug_out, "\n\n") ;
    }
  return;
} /* aqlQueryOut */

/************ recursive node dump **************/

void aqlNodeOut (ACEOUT debug_out,
                 AqlNode *inNode,
                 int ival,
                 char *rel,
                 BOOL suppress_next)
{ 
    int outermost = 0;
    char timeBuf[25] ;

    if (lisply) {
        indent (debug_out, ival*4) ;
        if (rel[strlen(rel)-1] != ' ')
        {
            aceOutf(debug_out, "%s:\n", rel);
            indent (debug_out, ival*4) ;
            outermost = 1;
            rel = "() ";
        }
        if (!inNode)
        {
            aceOutf(debug_out, "%c%c", rel[0], rel[1]);
            return;
        }

        aceOutf (debug_out, "%c%s", rel[0], aqlNodeTypeName(inNode->type)) ;

        if (inNode->name) aceOutf (debug_out, " \"%s\"", inNode->name) ;
        if (inNode->number) aceOutf (debug_out, " %d", inNode->number) ;
        if (inNode->op) aceOutf (debug_out, " %s", aqlOpTypeName(inNode->op)) ;

        if (inNode->type == nBOOL) aceOutf (debug_out, " %s", inNode->value.b?"TRUE":"FALSE") ;
        if (inNode->type == nTEXT) aceOutf (debug_out, " \"%s\"", inNode->value.s) ;
        if (inNode->type == nINT) aceOutf (debug_out, " %d", inNode->value.i) ;
        if (inNode->type == nFLOAT) aceOutf (debug_out, " %g", inNode->value.f) ;
        if (inNode->type == nDATE) aceOutf (debug_out, " %s", timeShow(inNode->value.t, timeBuf, 25)) ;

        if (inNode->loc) 
        {
            aceOutf (debug_out, " (loc");
            if (1)
            {
                aceOutf (debug_out, " %lx", (unsigned long)(inNode->loc));
                if (0)
                {
                    aceOutf (debug_out, " def=%lx ", (unsigned long)(inNode->loc->definition));
                }
            }
            if (1)
            {
                if (inNode->loc->definition == NULL)
                {
                    aceOutf (debug_out, "NULL_definition!");
                } else if ((inNode->loc->definition != inNode) &&
                    (inNode->loc->definition->type == nFROM_LOC))
                {
                    static int in_definition = 0;
                    if (! in_definition)
                    {
                        in_definition = 1;
                        aqlNodeOut (debug_out, inNode->loc->definition->left, ival+1, "<> ", TRUE);

                        in_definition = 0;
#if 0
                        indent (debug_out, ival*4) ;
#endif
                    }
                }
            }
            if (inNode->loc->source)
                aceOutf (debug_out, " %s", aqlLocSourceTypeName(inNode->loc->source));
            aceOutf (debug_out, ")");
        }


        if (inNode->left) 
            aqlNodeOut (debug_out, inNode->left, ival+1, "() ", suppress_next) ;
   
        if (inNode->right) 
            aqlNodeOut (debug_out, inNode->right, ival+1, "{} ", suppress_next) ;

        if (inNode->nxt && !suppress_next) 
            aqlNodeOut (debug_out, inNode->nxt, ival+1, "[] ", suppress_next) ;

        aceOutf (debug_out, "%c", rel[1]);
      
    } else {
        if (!inNode) return ;

        indent (debug_out, ival) ; 
        aceOutf (debug_out, "%s %s", rel, aqlNodeTypeName(inNode->type)) ;

        if (inNode->name) aceOutf (debug_out, " name \"%s\"", inNode->name) ;
        if (inNode->number) aceOutf (debug_out, " number %d", inNode->number) ;
        if (inNode->op) aceOutf (debug_out, " op %s", aqlOpTypeName(inNode->op)) ;

        if (inNode->type == nBOOL) aceOutf (debug_out, " %s", inNode->value.b?"TRUE":"FALSE") ;
        if (inNode->type == nTEXT) aceOutf (debug_out, " %s", inNode->value.s) ;
        if (inNode->type == nINT) aceOutf (debug_out, " %d", inNode->value.i) ;
        if (inNode->type == nFLOAT) aceOutf (debug_out, " %g", inNode->value.f) ;
        if (inNode->type == nDATE) aceOutf (debug_out, " %s", timeShow(inNode->value.t, timeBuf, 25)) ;

        aceOutf (debug_out, "\t\t%lx", (unsigned long)inNode);

        if (inNode->loc) 
        {
            aceOutf (debug_out, " (loc=%lx", (unsigned long)(inNode->loc));
            if (inNode->loc->source)
                aceOutf (debug_out, ", %s", aqlLocSourceTypeName(inNode->loc->source));
            aceOutf (debug_out, ")");
        }


        if (inNode->left) 
            aqlNodeOut (debug_out, inNode->left, ival+1, "left  ", suppress_next) ;
   
        if (inNode->right) 
            aqlNodeOut (debug_out, inNode->right, ival+1, "right ", suppress_next) ;

        if (inNode->nxt) 
            aqlNodeOut (debug_out, inNode->nxt, ival, "next  ", suppress_next) ;
    }

    if (outermost) {
        aceOutf(debug_out, "\n");
    }

    return;
} /* aqlNodeOut */

void aqlNodeValueOut (ACEOUT debug_out, AqlNode *inNode)
{
  char timeBuf[25] ;
  
  if (!inNode)
    {
      aceOutf (debug_out, "(NULL AqlNode) ");
      return ;
    }

  aceOutf (debug_out, "(`%s' [%c-value=", varName(inNode), inNode->vtype);
  
  if (inNode->isValue == NULLVAL)
    aceOutf (debug_out, "NULL") ;
  else if (inNode->isValue == FUZZYVAL)
    aceOutf (debug_out, "FUZZY") ;
  else
    /* TRUE */
    {
      switch (inNode->vtype)
	{
	case 'k':
	case 'g':
	  aceOutf (debug_out, "%x\"%s\"", inNode->value.k, name(inNode->value.k)) ;
	  break;
	  
	case 'i':
	  aceOutf (debug_out, "%d", inNode->value.i) ;
	  break;

	case 'f':
	  aceOutf (debug_out, "%g", inNode->value.f) ;
	  break;

	case 's':
	  aceOutf (debug_out, "%s", inNode->value.s);
	  break;

	case 't':
	  aceOutf (debug_out, "%s", timeShow(inNode->value.t, timeBuf, 25));
	  break;
	default:
	  aceOutf (debug_out, "%x", inNode->value.k) ;
	}
    }
  aceOutf (debug_out, "]) ") ;
  
  return;
} /* aqlNodeValueOut */

/**************** output a locators value ******************************/

void aqlLocOut (ACEOUT  debug_out,
		char   *text,
		AqlLoc *theLocator)
{
  char timeBuf[25] ;
  
  if (lisply) {
        aceOutf (debug_out, "%s ", text) ;
        aceOutf (debug_out, "name=%s", theLocator->definition->name);
    } else {
        aceOutf (debug_out, "%-10s ", text) ;
        aceOutf (debug_out, "loc=%lx", (unsigned long)theLocator) ;
        aceOutf (debug_out, ", (def=%lx)->name=%s", (unsigned long)theLocator->definition, theLocator->definition->name);
    }

  if (theLocator->isValue == NULLVAL)
    aceOutf (debug_out, " NULL\n") ;
  else if (theLocator->isValue == FUZZYVAL)
    aceOutf (debug_out, " FUZZY\n") ;
  else
    {
      aceOutf (debug_out, " type=%c", theLocator->type) ;
      
      switch (theLocator->type)
	{
	case 'k': case 'g':
    if (lisply) {
        aceOutf (debug_out, " value=%s", name(theLocator->value.k)) ;
    } else {
        aceOutf (debug_out, " value=%s(%x)", name(theLocator->value.k), theLocator->value.k) ;
    }
	  break ;
	  
	case 'i':
	  aceOutf (debug_out, " value=%d", theLocator->value.i) ; 
	  break ;
	  
	case 'f':
	  aceOutf (debug_out, " value=%g", theLocator->value.f) ; 
	  break ;
	  
	case 't':
	  aceOutf (debug_out, " value=%s", timeShow(theLocator->value.t, timeBuf, 25)) ; 
	  break ;
	  
	case 's':
	  aceOutf (debug_out, " value=%s", theLocator->value.s) ; 
	  break ;
	}
      aceOutf (debug_out, "\n") ;
    }

  return;
} /* aqlLocOut */

void aqlCheckEnvOut(ACEOUT debug_out,
                    char *text,
                    AqlCheckEnv *env,
                    int level)
{
    indent(debug_out, level); aceOutf(debug_out, "CheckEnv at %s:\n", text);

    if (env == NULL)
    {
        indent(debug_out, level); aceOutf(debug_out, "null_env");
    } else {

        if (env->pCurrDecl == NULL)
        {
            indent(debug_out, level); aceOutf(debug_out, "null curr_decl_pointer");
        } else {
            indent(debug_out, level); aceOutf(debug_out, "curr_decl: ");
            if (*(env->pCurrDecl) == NULL) {
                aceOutf(debug_out, "null");
            } else {
                aqlNodeOut(debug_out,
                           *(env->pCurrDecl),
                           level+1,
                           "curr_decl",
                           FALSE);
            }
        }

        if (env->pCurrWhere == NULL)
        {
            indent(debug_out, level); aceOutf(debug_out, "null curr_Where_pointer");
        } else {
            indent(debug_out, level); aceOutf(debug_out, "curr_Where: ");
            if (*(env->pCurrWhere) == NULL) {
                aceOutf(debug_out, "null");
            } else {
                aqlNodeOut(debug_out,
                           *(env->pCurrWhere),
                           level+1,
                           "curr_Where",
                           FALSE);
            }
        }

        if (env->currFromLoc == NULL) {
            indent(debug_out, level); aceOutf(debug_out, "null CurrFromLoc");
        } else {
            aqlNodeOut(debug_out,
                       env->currFromLoc,
                       level+1,
                       "curr_FromLoc",
                       FALSE);
        }
    }

}


/******************************************************************
 ************************ private functions ***********************
 ******************************************************************/

void aqlNodeOutPretty (AQL aql,
		       ACEOUT   debug_out, 
                       AqlNode *inNode,
                       BOOL     doBrackets, 
                       int      ival)
{
  char timeBuf[25] ;

  if (!inNode)
    aqlError (aql, 1000, 0, "NULL inNode in aqlNodeOutPretty()") ;

  if (doBrackets && nestingLevel(inNode) > 1)
    {
      aceOutf (debug_out, "(");
      ival += 1;
    }

  switch (inNode->type)
    {
    case nTABLE_ASSIGN:
      aceOutf (debug_out, "@%s := ", varName(inNode)) ;
      ival += strlen(varName(inNode)) + 5;
      aqlNodeOutPretty (aql, debug_out, inNode->left, FALSE, ival);
      ival -= strlen(varName(inNode)) + 5;
      
      if (inNode->nxt)
	{
	  aceOutf (debug_out, ";");
	  indent (debug_out, ival);
	  aqlNodeOutPretty (aql, debug_out, inNode->nxt, FALSE, ival);
	}
      break ;

    case nVAR_ASSIGN:
      aceOutf (debug_out, "$%s := ", varName(inNode)) ;
      ival += strlen(varName(inNode)) + 5;
      aqlNodeOutPretty (aql, debug_out, inNode->left, FALSE, ival);

      if (inNode->nxt)
	{
	  aceOutf (debug_out, ";");
	  indent (debug_out, ival);
	  aqlNodeOutPretty (aql, debug_out, inNode->nxt, FALSE, ival);
	}
      break ;

    case nTABLE_VAR:
      aceOutf (debug_out, "@%s", varName(inNode)) ;
      break ;

    case nTABLE_OP:
      {
	int oldival = ival;

	aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival);
	if (inNode->left->type != nTABLE_VAR)
	  {
	    indent (debug_out, ival);
	    ival += 1;
	  }
	aceOutf (debug_out, " %s ", aqlOpTypeName(inNode->op)) ;
	ival += (strlen(aqlOpTypeName(inNode->op)) + 2);
	
	if (inNode->left->type != nTABLE_VAR)
	  {
	    ival = oldival;
	    indent (debug_out, ival);
	  }
      
	aqlNodeOutPretty (aql, debug_out, inNode->right, TRUE, ival);

	if (inNode->nxt)
	  {
	    aceOutf (debug_out, ";");
	    oldival = ival;
	    indent (debug_out, ival);
	    aqlNodeOutPretty (aql, debug_out, inNode->nxt, FALSE, ival);
	  }
      }
      break ;

    case nTABLE_SFW:
      {
	int oldival = ival;
	AqlNode *n=0;
	BOOL isSelectAll = FALSE;
	BOOL isFirstFromClause = TRUE;
	
	if (inNode->right && inNode->right->nxt == 0 && /* only one selectfield .. */
	    strcmp(varName(inNode->right), "") == 0 && /* .. is the default locator */
	    (inNode->right->loc && inNode->right->loc->source == IS_CLASS))
	  isSelectAll = TRUE;
	
	/*	indent (debug_out, ival);*/
	if (isSelectAll)
	  {
	    aceOutf (debug_out, "select all ");
	    ival += 11;
	  }
	else
	  {
	    aceOutf (debug_out, "select ");
	    ival += 7;
	    
	    /* n loops through multiple select-fields */
	    n = inNode->right;
	    while(n)
	      {
		aqlNodeOutPretty (aql, debug_out, n, FALSE, ival); /* LOC_VARs */
		if (strcmp(varName(n), "") != 0)
		  /* named variable */
		  {
		    if (n->name && n->loc)
		      ival += strlen(varName(n->loc->definition)) ;
		    else
		      ival += strlen(varName(n));
		  }
		else
		  ival += 3;


		if (n->nxt)
		  {
		    aceOutf (debug_out, ", ");
		    ival += 2;
		  }
		n = n->nxt;
	      }
	  }
	
	ival = oldival + 2;
	n = inNode->left;	/* first of a linked list of FROM_LOCs */
	while (n)
	  {
	    if (!isSelectAll || !isFirstFromClause) 
	      indent (debug_out, ival);

	    if (isFirstFromClause)
	      {
		if (!isSelectAll) 
		  {
		    aceOutf (debug_out, "from ");
		    ival += 5;
		  }
		isFirstFromClause = FALSE;
	      }
	    aqlNodeOutPretty (aql, debug_out, n, FALSE, ival);

	    if (n->nxt)
	      aceOutf (debug_out, ", ");
	    n = n->nxt;
	  }
	
	if (inNode->nxt)
	  {
	    aceOutf (debug_out, ";");
	    ival = oldival;
	    indent (debug_out, ival);
	    aqlNodeOutPretty (aql, debug_out, inNode->nxt, FALSE, ival);
	  }
      }
      break ; 
      
    case nTABLE_ORDER:
      {
	int oldival = ival;

	aqlNodeOutPretty (aql, debug_out, inNode->left, FALSE, ival); /* table_expr */
	ival += 2;
	indent (debug_out, ival);
	aceOutf (debug_out, "order") ; ival += 5;
	if (inNode->right)
	  {
	    aceOutf (debug_out, " by ") ; ival += 4;
	    aqlNodeOutPretty (aql, debug_out, inNode->right, FALSE, ival); /* SORT_FIELDs */
	  }
      else
	{
	  if (inNode->op != oASC)
	    {
	      aceOutf (debug_out, " %s", aqlOpTypeName(inNode->op));
	      ival += (strlen(aqlOpTypeName(inNode->op)) + 2);
	    }
	}

	ival = oldival ;
	if (inNode->nxt)
	  {
	    aceOutf (debug_out, ";");
	    indent (debug_out, ival);
	    aqlNodeOutPretty (aql, debug_out, inNode->nxt, FALSE, ival); /* another table_expr */
	  }
      }
      break ; 

    case nSORT_FIELD:
      aceOutf (debug_out, ":%d", inNode->number) ;
      ival += strlen(messprintf(":%d", inNode->number));

      if (inNode->left)
	{
	  aceOutf (debug_out, " "); ival += 1;
	  aqlNodeOutPretty (aql, debug_out, inNode->left, FALSE, ival); /* qualifier */
	}
      if (inNode->nxt)
	{
	  aceOutf (debug_out, ", "); ival += 2;
	  aqlNodeOutPretty (aql, debug_out, inNode->nxt, FALSE, ival);
	}
      break ; 

    case nSORT_FIELD_EXPR:
      aceOutf (debug_out, "******************") ;
      break ; 

    case nSORT_QUALIFIER:
      if (inNode->op != oASC)
	aceOutf (debug_out, "%s", aqlOpTypeName(inNode->op)) ;
      break ; 

    case nFROM_LOC:
    case nFROM_TABLE:
      {
	AqlNode *n;

	if (strcmp(varName(inNode), "") != 0)
	  { aceOutf (debug_out, "%s in ", varName(inNode)); }

	/* output linked list of nodes hanging off the left
	   Most commonly it starts with a LOC_VAR followed
	   by nodes like LOC_FOLLOW_TAG */
	for (n = inNode->left; n; n = n->nxt)
	  aqlNodeOutPretty (aql, debug_out, n, FALSE, ival);

	if (inNode->right)
	  {
	    ival += strlen(varName(inNode)) + 1; /* position 'where' below the 'in' in the decl */
	    indent (debug_out, ival);
	    aceOutf (debug_out, "where ");
	    aqlNodeOutPretty (aql, debug_out, inNode->right, FALSE, ival+5);
	  }
      }
      break ; 
      /*
    case nFROM_TABLE:
      aceOutf (debug_out, "%s in ", varName(inNode)) ;
      aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival);

      if (inNode->right)
	{
	  indent (debug_out, ival);
	  aceOutf (debug_out, "where ") ;
	  aqlNodeOutPretty (aql, debug_out, inNode->right, FALSE, ival+2);
	}
      if (inNode->nxt)
	{
	  aceOutf (debug_out, ", ") ;
	  indent (debug_out, ival);
	  aqlNodeOutPretty (aql, debug_out, inNode->nxt, FALSE, ival);
	}
      break ; 
      */
    case nCLASS:
      aceOutf (debug_out, "class %s", inNode->name) ;
      break ; 

    case nOBJECT:
      aceOutf (debug_out, "object(") ;
      aqlNodeOutPretty (aql, debug_out, inNode->left, FALSE, ival);
      aceOutf (debug_out, ", ") ;
      aqlNodeOutPretty (aql, debug_out, inNode->right, FALSE, ival);
      aceOutf (debug_out, ")") ;
      break ;

    case nLOC_VAR:
      {
	if (strcmp(varName(inNode), "") != 0)
	  /* named variable */
	  {
	    if (inNode->name && inNode->loc)
	      aceOutf (debug_out, "%s", varName(inNode->loc->definition)) ;
	    else
	      aceOutf (debug_out, "%s", varName(inNode));
	  }
	else
	  {
	    /* default iterator */
	    if (inNode->nxt &&
		(inNode->nxt->type == nLOC_LOCAL_TAG ||
		 inNode->nxt->type == nLOC_LOCAL_POS ||
		 inNode->nxt->type == nLOC_LOCAL_FIELD ||
		 inNode->nxt->type == nLOC_FOLLOW_TAG ||
		 inNode->nxt->type == nLOC_FOLLOW_POS ||
		 inNode->nxt->type == nLOC_METHOD))
	      aceOutf (debug_out, "");
	    else
	      aceOutf (debug_out, "->0");	/* class/row based def-loc 
						 * quoted in select clause */
	  }
      }
      break ;

    case nLOC_FOLLOW_TAG:
      aceOutf (debug_out, "->%s[0]", name(inNode->value.k)) ;
      break ;

    case nLOC_LOCAL_POS:
      aceOutf (debug_out, "[%d]", inNode->number) ;
      break ;

    case nLOC_LOCAL_TAG:
      aceOutf (debug_out, "[%s]", name(inNode->value.k)) ;
      break ;

    case nLOC_LOCAL_FIELD:
      if (inNode->name)
	aceOutf (debug_out, ":%s", inNode->name) ;
      else
	aceOutf (debug_out, ":%d", inNode->number) ;
      break ;

    case nLOC_FOLLOW_POS:
      aceOutf (debug_out, "->%d", inNode->number) ;
      break ;

    case nLOC_METHOD:
      aceOutf (debug_out, ".%s", inNode->name);
      break;

    case nLOC_TABLE_FIELD:
      aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival+5); /* the safe_table */
      aceOutf (debug_out, ":%d", inNode->number);
      break;

    case nEXPR_TABLE_FUNC:
      aceOutf (debug_out, "%s", aqlOpTypeName(inNode->op)) ;
      ival += strlen(aqlOpTypeName(inNode->op));
      aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival); /* safe_table */
      break ;

    case nEXPR_TABLE_FUNC_FIELD:
      aceOutf (debug_out, "%s", aqlOpTypeName(inNode->op)) ;
      ival += strlen(aqlOpTypeName(inNode->op)) + 1;
      aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival); /* safe_table */
      ival -= strlen(aqlOpTypeName(inNode->op)) + 1;
      aceOutf (debug_out, ":%d", inNode->number) ;
      break ;


    case nEXPR_OP:
      if (inNode->op == oUMINUS)
	{
	  aceOutf (debug_out, "%s", aqlOpTypeName(inNode->op)) ;
	  aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival); /* expr */
	}
      else if (inNode->op == oPLUS || 
	    inNode->op == oMINUS || 
	    inNode->op == oTIMES ||
	    inNode->op == oDIVIDE)
	{
	  aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival);
	  aceOutf (debug_out, " %s ", aqlOpTypeName(inNode->op)) ;
	  aqlNodeOutPretty (aql, debug_out, inNode->right, TRUE, ival);
	}
      else
	{
	  aceOutf (debug_out, "%s(", aqlOpTypeName(inNode->op)) ;
	  aqlNodeOutPretty (aql, debug_out, inNode->left, FALSE, ival);
	  if (inNode->right)
	    {
	      aceOutf (debug_out, ", ");
	      aqlNodeOutPretty (aql, debug_out, inNode->right, FALSE, ival);
	    }
	  aceOutf (debug_out, ")");
	}
      break ;

    case nBOOL_EXISTS:
      aceOutf (debug_out, "exists ") ;
      aqlNodeOutPretty (aql, debug_out, inNode->left, FALSE, ival);
      break ;

    case nBOOL_COMPARISON:
      aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival);

      aceOutf (debug_out, " %s ", aqlOpTypeName(inNode->op));

      aqlNodeOutPretty (aql, debug_out, inNode->right, TRUE, ival);
      break ;

    case nBOOL_NOT:
      aceOutf (debug_out, "not ") ;
      aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival);

      break ;

    case nBOOL_OP:
      aqlNodeOutPretty (aql, debug_out, inNode->left, TRUE, ival);

      aceOutf (debug_out, " %s ", aqlOpTypeName(inNode->op)) ;

      aqlNodeOutPretty (aql, debug_out, inNode->right, TRUE, ival);
      break ;

    case nVAR:
      aceOutf (debug_out, "$%s", varName(inNode)) ;
      break ;

    case nTEXT:
      aceOutf (debug_out, "\"%s\"", inNode->value.s) ;
      break ;

    case nINT:
      aceOutf (debug_out, "%d", inNode->value.i) ;
      break ;

    case nFLOAT:
      aceOutf (debug_out, "%g", inNode->value.f) ;
      break ;

    case nDATE:
      aceOutf (debug_out, "%s", timeShow(inNode->value.t, timeBuf, 25)) ;
      break ;

    case nBOOL:
      aceOutf (debug_out, "%s", inNode->value.b ? "true" : "false") ;
      break;

    case nKEY:
      aceOutf (debug_out, "\"%s\"", name(inNode->value.k));
      break;

    case nTABLE_SFW_ALL:
    case nSELECT_FIELD:
    case nBOOL_EXISTS_TAG:
    case nLOC_LOCAL_TAG_NAME:
    case nLOC_FOLLOW_TAG_NAME:
    case nSORT_FIELD_NAME:
    case nLOC_LOCAL_FIELD_NAME:
    case nLOC_TABLE_FIELD_NAME:
    case nEXPR_TABLE_FUNC_FIELD_NAME:
      aqlError (aql, 1001, 0,
		"aqlNodeOutPretty() - %s should have been eliminated after check2",
		aqlNodeTypeName(inNode->type)) ;
      break ;

    }

  if (doBrackets && nestingLevel(inNode) > 1)
    aceOutf (debug_out, ")");

  return;
} /* aqlNodeOutPretty */

/************************************************************/

static void indent (ACEOUT debug_out, int ival)
     /** little utility for indentation **/
     /** NB writes end of previous line **/
{ 
  int i ; 

  aceOutf (debug_out, "\n") ; 
  for (i = 0 ; i < ival ; ++i) 
    aceOutf (debug_out, " ") ; 

  return;
} /* indent */

static const char *varName (AqlNode *inNode)
{
  const char *s ;

  if (inNode->scope && inNode->scope->dict)
    s = dictName(inNode->scope->dict, inNode->number) ;
  else
    s = inNode->name ;

  if (!s)
    return "(null)";

  if (strcmp(s, "_DEF") == 0)
    return "";

  if (s[0] == '_')		/* temporary variable */
    return messprintf("tmp%s", s);

  return s ;
} /* varName */

static int nestingLevel (AqlNode *inNode)
{
  int left, right;

  if (!inNode) return 0;

  left = nestingLevel (inNode->left) + 1;
  right = nestingLevel (inNode->right) + 1;

  return left > right ? left : right;
} /* nestingLevel */


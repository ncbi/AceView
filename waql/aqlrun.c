/* $Id: aqlrun.c,v 1.10 2015/09/25 17:28:55 mieg Exp $ */
/*  File: aqlrun.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
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
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Sep  2 14:30 2004 (edgrif)
 * * Oct  4 17:25 1999 (fw): added new range of timestamp methods
 *              create_time/create_session/update_time/update_session
 *              method 'timestamp' now called 'node_session'
 * * Sep  7 17:10 1999 (fw): new loc->source type IS_EXPR supported
 *              to handle declaration based on expression values.
 * * May  3 11:32 1999 (fw): new notation to empulate 'FIND <class> <templ>' queries
 *              .. from class <className>."<templ>" .. is now evaluated
 * * May  3 10:36 1999 (fw): case-sensitivity observed for comparing object-names
 * * Mar 22 11:31 1999 (fw): moved fuzzy re-run from evalTable to aqlEval
 *              direct to evalTable in fuzzy mode now possible
 * * Mar 10 10:53 1999 (fw): length method can now report the length of array objects
 * * Feb 28 14:08 1999 (fw): consistent evaluation of expression, boolean
 *              expressions and booleans used in a where-clause context.
 * * Feb 26 10:38 1999 (fw): merged evalBool into evalExpr
 * * Feb 17 17:57 1999 (fw): cheap eval of ->0 and [0]
 * * Jan 11 16:08 1999 (fw): defaults for tableSort
 * * Aug 26 18:17 1998 (fw): proper tableSort
 * * Aug 26 13:36 1998 (fw): wrote timestamp method
 * * Aug 25 16:57 1998 (fw): modified sfwReport to omit duplicate rows
 * * Aug 25 16:56 1998 (fw): nodes containing tables don't use 
 *                              node->name for their columntypes anymore
 * * Aug 25 14:44 1998 (fw): moved tableOps out to tablePackage
 * * Aug 17 16:48 1998 (fw): method 'length' for string types written
 * * Aug 17 16:21 1998 (fw): replaced EXPR_LOC_FUNC with LOC_METHOD, 
 *                              evaluated more like FOLLOW_TAG etc.
 * * Aug 17 16:21 1998 (fw): min () and max () for string types now
 * * Aug 14 18:11 1998 (fw): added wildcard string matching using LIKE operator
 * * Aug 10 16:53 1998 (fw): checked evaluation of object (), name (), class ()
 * * Aug  6 17:19 1998 (fw): able to compare tags (by name) to strings and tags to tags
 * * Aug  5 17:33 1998 (fw): introduced threadsafe AQL handle
 * * Aug  3 16:22 1998 (fw): the object constructor now evaluates its left/right expr
 * * Aug  3 14:41 1998 (fw): all debug-info in this module uses debug-level 3 and above
 * * Jul 27 15:44 1998 (fw): checked error-codes/messages
 * * Jul 24 15:11 1998 (fw): fixed count-tableFunc to count non-empty entries only
 * * Jul 24 11:37 1998 (fw): tableFunc now works on column-selections
 * * Jul 23 14:47 1998 (fw): in case of a single result value aqlEval
 *                            creates a table return-value with one element
 * * Jul 22 17:36 1998 (fw): removed tableExpr check for just one column in tableExpr
 *                             now done earlier on in aqlCheck1 ()-postTableProcess ()
 * * Jul 22 15:34 1998 (fw): included setting of ->vtype for tablenodes in evalTable ()
 * * Jul 22 15:21 1998 (fw): fixed propagation of values of nVAR-nodes
 * * Jul 16 13:52 1998 (fw): added union/diff/intersect functionality
 * * Jul 16 12:59 1998 (fw): pre and post evaluation table-type checking 
 *                             for the two tables of a TABLE_OP node
 * * Jul 15 11:26 1998 (fw): row-var, table-var stuff works now.
 *                           loopVariables:IS_COLUMN, now set value 
 *                              of loc according to value in table
 * * Jul 14 13:35 1998 (fw): comparison involving boolean types
 * * Jul 14 11:54 1998 (fw): added min,max for tables of date values
 * * Jul 13 17:54 1998 (fw): finished count,sum,min,max,avg functions for int/float
 * * Jul 13 12:17 1998 (fw): implemented avg function, yet still awkward syntax
 * * Jul 10 17:45 1998 (fw): attempted count-function
 * * Jul 10 14:32 1998 (fw): compare keys to keys by their string values
 * * Jul 10 10:06 1998 (fw): sfwReport doesn't complain about uninit column type
 * * Jul 10 10:06 1998 (fw): false evaluation now correctly sets ->isEmpty
 * * Jul  9 17:42 1998 (fw): date comparison complete
 * * Jul  8 16:38 1998 (fw): expression evaluation of numerics and dates
 *                            almost complete
 * * Jul  8 16:37 1998 (fw): boolean evaluation almost complete
 * Created: Mon Aug  4 23:45:09 1997 (Durbin)
 *-------------------------------------------------------------------
 */

#include <wh/acedb.h>
#include <wh/aceio.h>
#include <wh/bitset.h>
#include <wh/lex.h>
#include <whooks/sysclass.h>
#include <whooks/systags.h>
#include <wh/bindex.h>
#include <wh/query.h>
#include <wh/a.h>
#include <wh/pick.h>
#include <wh/utils.h>                /* for lexstrcmp */
#include <waql/aql_.h>


/******************************************************************/

static void postAssignTypesAndClasses (AqlNode *inNode, AQL aql, int level);
static void postTableTypeCheck (AqlNode *inNode, AQL aql, int level);

static AqlFuzzyType evalTable (AqlNode *tableNode, BOOL isFuzzyRun, AQL aql, int level);
static AqlFuzzyType processThisLoc (AqlNode *sfw,  AqlNode *fromFieldNode, BOOL isFuzzyRun, AQL aql, int level) ;
static AqlFuzzyType loopVariables (AqlNode *sfw, AqlNode *fromFieldNode, BOOL isFuzzyRun, AQL aql, int level) ;
static AqlFuzzyType loopVariablesSingle (AqlNode *tableSfwNode, AqlNode *fromFieldNode, BOOL isFuzzyRun, AQL aql, int level,
                                         BOOL isFuzzyStatus, BOOL isEmptyStatus, AqlFuzzyType status);
static BOOL AqlValueEqual (char at, AqlType av, char bt, AqlType bv, AQL aql);

static AqlFuzzyType sfwReport (AqlNode *sfw, BOOL isFuzzyRun, AQL aql, int level);

static AqlFuzzyType evalExpr (AqlNode *exprNode, BOOL isFuzzyRun, AQL aql, int level) ;

static BOOL stringComparison (AQL aql, AqlOpType boolOp, char *leftString, char *rightString, BOOL isCaseSensitive);
static BOOL integerComparison (AQL aql, AqlOpType boolOp, int leftInt, int rightInt);
static BOOL floatComparison (AQL aql, AqlOpType boolOp, float leftFloat, float rightFloat);
static BOOL boolComparison (AQL aql, AqlOpType boolOp, BOOL leftBool, BOOL rightBool);
static BOOL keyComparison (AQL aql, AqlOpType boolOp, KEY leftKey, KEY rightKey);

static BOOL valueAsBool (AQL aql, AqlNode *inNode);
static int evalIntArithmetic (AQL aql, AqlOpType exprOp, int leftInt, int rightInt);
static float evalFloatArithmetic (AQL aql, AqlOpType exprOp, float leftFloat, float rightFloat);

static void indent (AQL      aql,
                   int level);

/*********************************************************************/
/********************** public functions *****************************/
/*********************************************************************/

void aqlCheck2 (AQL aql)
/* do the database dependent query processing, 
   assign and check value types */
{
    aqlTraverse (aql->query->root, aql, 0, postAssignTypesAndClasses, 0) ;

    aqlTraverse (aql->query->root, aql, 0, postTableTypeCheck, 0) ; 

    return;
} /* aqlCheck2 */

/********************************************************************/

TABLE *aqlEval (AQL aql)
/* the main entry function to the query evaluation
   returned is a table with the results,
   returns 0, if the return value is of a different type */
{
  AqlNode *tableNode;
  AqlNode *resultNode;
  AqlFuzzyType result;
  BOOL doFuzzyRun = TRUE;

  messStatus ("Exec AQL, F4 to interrupt...");

 rerun:
  tableNode = resultNode = aql->query->root;
  result = NULLVAL;
  
  while (tableNode)
    {
      result = evalTable (tableNode, doFuzzyRun, aql, 0) ;
      resultNode = tableNode;
      
      tableNode = tableNode->nxt;
    }

  if (aql->evalDebugLevel >= 3)
    {
      if (result == NULLVAL)
	{
	  aceOut (aql->debug_out, "\n**** evalTable failed - table value is NULL\n") ;
	}
      else if (result == KNOWNVAL)
	{
	  aceOut (aql->debug_out, "\n**** evalTable succeeded - table value is KNOWN\n") ;
	}
    }
  
  if (result == FUZZYVAL)
    {
      if (!doFuzzyRun)
	aqlError (aql, 600, 0, "XXX unfuzzy re-run produced fuzzy result XXX") ;

      if (aql->evalDebugLevel >= 3)
	aceOut (aql->debug_out, "\n**** evalTable isn't sure yet - table value is FUZZY\n") ;
      
      doFuzzyRun = FALSE;
      goto rerun;
    }
  
  if (resultNode->vtype != 'T')
    aqlError (aql, 610, 0, "aqlEval () - tableNode not of value type 'T'") ;
  
  if (aql->evalDebugLevel >= 1)
    aceOutf (aql->debug_out, "Returning table of %d entries\n", tableMax (resultNode->value.T));
  
  /* return the result-table of the last query in the series */
  return resultNode->value.T;
} /* aqlEval */

/*********************************************************************/
/********************* private functions *****************************/
/*********************************************************************/

/*
   process nodes of type

   LOC_LOCAL_TAG_NAME, 
     LOC_FOLLOW_TAG_NAME - find tag-name in db
   CLASS                 - find class in db
   OBJECT                - find class from text-literals
   FROM_TABLE            - look at its locator
   FROM_LOC              - assign types and sources to the locator

   executed as a postFunc by aqlTraverse as part of aqlCheck2 ()
*/
static void postAssignTypesAndClasses (AqlNode *inNode, 
				       AQL      aql,
                                       int level)
{
  KEY key ;

  switch (inNode->type)
    {
    case nLOC_LOCAL_TAG_NAME:
    case nLOC_FOLLOW_TAG_NAME:
      {
	if (!lexword2key (inNode->name, &inNode->value.k, 0))
	  aqlError (aql, 944, inNode->pos, messprintf ("Unrecognized tag %s", inNode->name)) ;
	
	if (inNode->type == nLOC_LOCAL_TAG_NAME)
	  inNode->type = nLOC_LOCAL_TAG ;
	else /* nLOC_FOLLOW_TAG_NAME  */
	  inNode->type = nLOC_FOLLOW_TAG ;
      }
      break;



    case nCLASS:
      {
	if (lexword2key (inNode->name, &key, _VClass))
	  inNode->value.k = key ;
	else
	  aqlError (aql, 941, inNode->pos, messprintf ("Unrecognized class %s", inNode->name)) ;
      }
      break;


      /* transform text-literals for class and key args to key-values */
      /* values are evaluated through IS_OBJECT clause */
    case nOBJECT:
      {
	/* text literal for the class */
	if (inNode->left->type == nTEXT ||
	    (inNode->left->loc && inNode->left->loc->type == 's'&& inNode->left->loc->isValue == KNOWNVAL)) 
	  {
	    /* we find the key to the class-literal in the database,
	       then change node->left as if it had quoted this key-value literally (nKEY node) */
	    if (lexword2key (inNode->left->value.s, &key, _VClass))
	      { 
		/* change node type, now that we've found the key */
		inNode->left->type = nKEY;
		inNode->left->vtype = 'k';
		inNode->left->value.k = key;
		inNode->left->isValue = KNOWNVAL;
	      }
	    else
	      aqlError (aql, 942, inNode->pos, messprintf ("Unrecognized class %s in construct object (\"<class>\", ...)",
							   inNode->left->value.s)) ; /* used to be warning */
	  }
      }
      break;



      /* assign source and ancillary info from variable declarations */
    case nFROM_TABLE:		/* added by fw-980714 */
      {
	if (inNode->left->type == nTABLE_VAR)
	  {
	    if (inNode->loc->source != IS_ROW)
	      aqlError (aql, 620, 0,
			"postAssignTypesAndClasses () - inNode->loc->source != IS_ROW "
			"(should have been set by aqlcheck.c:preScope)") ;

	    inNode->loc->mustExist = TRUE;
	    inNode->loc->type      = inNode->loc->definition->vtype;
	  }
	else
	  aqlError (aql, 620, 0,
		    "postAssignTypesAndClasses () - type of FROM_TABLE->left is not a TABLE_VAR");
      }
      break;



    case nFROM_LOC:
      /* assign source and ancillary info from variable declarations */
      { 
	AqlNode *fromLocNode = inNode;
	
	if (fromLocNode->left->type == nCLASS)
	  { 
	    fromLocNode->loc->source = IS_CLASS ;
	    fromLocNode->loc->type = 'k' ;
	    fromLocNode->loc->mustExist = TRUE ;
	    fromLocNode->loc->key = fromLocNode->left->value.k ; /* computed above by >case nCLASS:< */
	  }
	else if (fromLocNode->left->type == nOBJECT)
	  { 
	    fromLocNode->loc->source = IS_OBJECT;
	    fromLocNode->loc->type = 'k';
	    fromLocNode->loc->mustExist = TRUE ;
	  }
	else if (fromLocNode->left->type == nVAR
		 || fromLocNode->left->type == nINT
		 || fromLocNode->left->type == nFLOAT
		 || fromLocNode->left->type == nTEXT
		 || fromLocNode->left->type == nBOOL
		 || fromLocNode->left->type == nDATE
		 || fromLocNode->left->type == nEXPR_TABLE_FUNC
		 || fromLocNode->left->type == nEXPR_TABLE_FUNC_FIELD
		 || fromLocNode->left->type == nEXPR_TABLE_FUNC_FIELD_NAME
		 || fromLocNode->left->type == nEXPR_OP
		 || fromLocNode->left->type == nBOOL
		 || fromLocNode->left->type == nBOOL_NOT
		 || fromLocNode->left->type == nBOOL_EXISTS
		 || fromLocNode->left->type == nBOOL_COMPARISON)
	  /* these should be all possible node-types for the
	   * 'expr' rule in the grammar */
	  {
	    fromLocNode->loc->source = IS_EXPR;
	  }
	else if (fromLocNode->left->nxt)
	  {
	    /******************************************************************/
	    /* by now fromLocNode->left must be a locator-variable LOC_xxx    *
	     * it is guaranteed to have a ->nxt, so it is base upon something *
	     * see aqlcheck.c:preLocCheck ()                                   */
	    /******************************************************************/
	    if (fromLocNode->left->nxt->type == nLOC_METHOD)
	      { 
		fromLocNode->loc->source = IS_METHOD;
		/* ->type depends on the method */
	      }
	    else if (fromLocNode->left->nxt->type == nLOC_LOCAL_FIELD)
	      { 
		if (fromLocNode->left->loc->source != IS_ROW)
		  aqlError (aql, 620, 0,
			    "postAssignTypesAndClasses () - nLOC_LOCAL_FIELD has to have a IS_ROW loc->source");
		
		fromLocNode->loc->source = IS_COLUMN ;
		
		/* The inNode->left      is LOC_VAR (it's loc->source = IS_ROW)
		 *     inNode->left->nxt is LOC_LOCAL_FIELD (->number is the column No)
		 */
		fromLocNode->loc->column = fromLocNode->left->nxt->number - 1 ;
		/* Note: Adjust column numbers - Whereas the user refers to 
		 * column 1 being the leftmost, this is column 0 in the program */
	      }
	    else if (fromLocNode->left->nxt->type == nLOC_LOCAL_TAG)
	      { 
		fromLocNode->loc->source = IS_LOCAL_TAG ;
		fromLocNode->loc->type = 'g' ;	/* settle on tag */
		fromLocNode->loc->value.k = fromLocNode->left->nxt->value.k ; /* computed above by >case LOCAL_TAG/FOLLOW_NAME< */
	      }
	    else if (fromLocNode->left->nxt->type == nLOC_FOLLOW_TAG)
	      { 
		fromLocNode->loc->source = IS_FOLLOW_TAG ;
		fromLocNode->loc->type = 'g' ;	/* settle on tag */
		fromLocNode->loc->value.k = fromLocNode->left->nxt->value.k ; /* computed above by >case LOCAL_TAG/FOLLOW_NAME< */
	      }
	    else if (fromLocNode->left->nxt->type == nLOC_LOCAL_POS)
	      { 
		fromLocNode->loc->source = IS_LOCAL_POS ;
	      }
	    else if (fromLocNode->left->nxt->type == nLOC_FOLLOW_POS)
	      { 
		fromLocNode->loc->source = IS_FOLLOW_POS ;
	      }
	    else if (fromLocNode->left->nxt->type == nLOC_LOCAL_FIELD)
	      { 
		AqlNode *defNode = fromLocNode->loc->definition->left;
		/* caused by e.g.
		   select a:1 from a in class movie
		   select a:2 from a in @table:1
		*/
		aqlError (aql, 930, fromLocNode->left->nxt->pos, 
			  messprintf ("Invalid column-selection on %s:%d\n The %s is not row-based.", 
				      (strcmp (defNode->name, "_DEF") == 0) ? "" : defNode->name,
				      fromLocNode->left->nxt->number,
				      (strcmp (defNode->name, "_DEF") == 0) ? "default iterator" : defNode->name)
			  ) ;
	      }
	    else if (fromLocNode->left->nxt->type == nLOC_LOCAL_FIELD_NAME)
	      { 
		AqlNode *defNode = fromLocNode->loc->definition->left;
		/* caused by e.g.
		   select a:fieldname from a in class movie  --or--
		   select a:fieldname from a in @table:1
		*/
		aqlError (aql, 931, fromLocNode->left->nxt->pos
			  , messprintf ("Invalid column-selection on %s:%s\nThe %s is not row-based.", 
					(strcmp (defNode->name, "_DEF") == 0) ? "" : defNode->name,
					fromLocNode->loc->definition->left->name,
					(strcmp (defNode->name, "_DEF") == 0) ? "default iterator" : defNode->name)
			  ) ;
	      }
	    else 
	      aqlError (aql, 620, 0,
			messprintf ("postAssignTypesAndClasses () - invalid type of declaration for %s",
				    inNode->loc->definition->left->name)
			);
	  }
	else
	  aqlError (aql, 620, 0,
		    messprintf ("postAssignTypesAndClasses () - invalid type (%s) of fromLocNode->left", 
				aqlNodeTypeName (fromLocNode->left->type))
		    );
      }	/* end-case FROM_LOC */
      break;



    default:;
    } /* end-switch inNode->type */

} /* postAssignTypesAndClasses */

/***********************************************************************************/

static void postTableTypeCheck (AqlNode *inNode,
				                AQL      aql,
                                int level)
/* 
   process nodes of type
   TABLE_OP - compatibility check after types and classes have been assigned to locators
*/
/* executed as a postFunc by aqlTraverse as part of aqlCheck2 () */
{
  switch (inNode->type)
    {
    case nTABLE_OP:
      /* the two tables have already been processed to a certain degree :-
	 (1) both tables have the same number of columns - see aqlcheck.c:postTableProcess ()
	 (2) we know table->ncol value for both tables
	 (3) the tabletypes are initialised but yet unknown
	 (4) we might be able to find out the types of the table columns 
	     by looking at the locators that define the select-fields
       */
      {
	int    col ; 
	char   leftValueType, 
	       rightValueType;
	AqlNode *selectNodeLeft, *selectNodeRight;
	
	/* NOTE, we don't bother yet with setting the columnNames, 
	   as they might harmlessly mismatch and are not used yet anyway */

	/* we can now check the two tables for type combatability */
	/* we have to get the type of each selectField, the selectFields locator is
	   defined by a FROM_LOC, whose ->left LOC_VAR has a locator of the selectField's type */
	for (selectNodeLeft = inNode->left->right, selectNodeRight = inNode->right->right, col = 0; 
	     selectNodeLeft && selectNodeRight; 
	     selectNodeLeft = selectNodeLeft->nxt, selectNodeRight = selectNodeRight->nxt, ++col)
	  {
	    if (selectNodeLeft->loc->type)
	      leftValueType  = selectNodeLeft->loc->type;
	    else
	      leftValueType  = selectNodeLeft->loc->definition->left->loc->type;

	    if (selectNodeRight->loc->type)
	      rightValueType = selectNodeRight->loc->type;
	    else
	      rightValueType = selectNodeRight->loc->definition->left->loc->type;


	    if (leftValueType == 0  || leftValueType == '0' ||
		rightValueType == 0 || rightValueType == '0')
	      {
		/* can't decide on correct types yet, we'll have to wait until tables are evaluated */
		continue;	/* try next field */
	      }

	    if (leftValueType != rightValueType)
	      /* types are known but mismatch */
	      aqlError (aql, 960, inNode->pos, messprintf ("table combining operator - error in column %d\n"
							   " value-type mismatch between the columns the two tables", 
							   col+1)
			);
	  }


	break;
      }

    default:;
    }
} /* postTableTypeCheck */

/********************************************************************/


/***************************************************/
/******* top-level query evaluation function, ******
 ***** called on query->root and all tablenodes ****/
/***************************************************/
static AqlFuzzyType evalTable (AqlNode *tableNode,
			       BOOL     isFuzzyRun,
			       AQL      aql,
                   int      level)
{
    /* NOTE that the actual TABLE datastructures have already been created
       during check1 in postTableProcess (). Here we just manipulate exisiting
       TABLEs and fill in the values or set the types. */

    switch (tableNode->type)
    {
    case nTABLE_SFW:
    {
    AqlResultConsumer * consumer = tableNode->resultConsumer;

    if (!consumer->sorting)
    {
        aceOutf (aql->debug_out, "%x->%x not sorting, should output table header here? or a bit later...\n", tableNode, consumer);
    }
        tableClear (tableNode->value.T, 1024) ;
	
        if (tableNode->left)
            /* TABLE_SFW->left is a FROM_LOC or a FROM_TABLE */
            /* i.e. evaluate all locators */
        {
            tableNode->isValue = loopVariables (tableNode, tableNode->left, isFuzzyRun, aql, level+1) ;
        }
        else
        {
            /* no FROM fields, probably just an EXPR_TABLE_FUNC as a select field */
            tableNode->isValue = sfwReport (tableNode, isFuzzyRun, aql, level+1);
        }

        if (tableNode->resultConsumer->dictStringStack) stackDestroy (tableNode->resultConsumer->dictStringStack);
#if 0
        if (tableNode->resultConsumer->resultTableDict) dictClear (tableNode->resultConsumer->resultTableDict);
#else
        if (tableNode->resultConsumer->resultTableDict)
	  {
	    dictDestroy (tableNode->resultConsumer->resultTableDict);
	    tableNode->resultConsumer->resultTableDict = dictHandleCreate (100, aql->query->handle) ;
	  }

#endif
    }
    break; /* TABLE_SFW */


    case nTABLE_ASSIGN:
    {
        /* we've assigned the table-result of an select-from-where query to a variable name */
        /* TABLE_ASSIGN->left is a TABLE_SFW */
        tableNode->isValue = evalTable (tableNode->left, isFuzzyRun, aql, level+1);

        /* propagate the table value from the TABLE_SFW to this node */
        tableNode->value.T = tableNode->left->value.T ;
    }
    break; /* TABLE_ASSIGN */


    case nVAR_ASSIGN:
    {
        char typeString[2];

        /* we've assigned the scalar value of an expression to a variable name */
        /* VAR_ASSIGN->left is an expression, e.g. an EXPR_TABLE_FUNC */
        tableNode->isValue = evalExpr (tableNode->left, isFuzzyRun, aql, level+1);

        /* copy the type and value to the locator for this node,
           so we can refer to this value */
        tableNode->loc->type = tableNode->left->vtype;
        tableNode->loc->value = tableNode->left->value;
	

        /* set the column type according to the type of the scalar value */
        typeString[0] = tableNode->left->vtype;
        typeString[1] = '\0';

        tableSetTypes (tableNode->value.T, typeString);

        /* the value of node->left is propagated to the table of the node */
        /* set the value of the element in row1,column1 */
        switch (tableNode->left->vtype)
	  {
	  case 'i': tableInt (tableNode->value.T,0,0) = tableNode->left->value.i; break; 
	  case 'f': tableFloat (tableNode->value.T,0,0) = tableNode->left->value.f; break; 
	  case 't': tableDate (tableNode->value.T,0,0) = tableNode->left->value.t; break; 
	  case 'b': tableDate (tableNode->value.T,0,0) = tableNode->left->value.b; break; 
	  case 's': tableSetString (tableNode->value.T,0,0,tableNode->left->value.s); break; 
	  case '0': tableSetEmpty (tableNode->value.T,0,0); break;
	  default:
	    aqlError (aql, 620, 0, "evalTable () - invalid scalar result value-type");
	  }
    }
    break; /* VAR_ASSIGN */


    case nTABLE_VAR:		/* needed here, if table-vars are used in expressions */
    {
        /* try to get its value from where it was defined */
        if (tableNode->loc->definition->type != nTABLE_ASSIGN)
	  aqlError (aql, 620, 0,
		    "evalTable ()->case TABLE_VAR: has to be defined by a TABLE_ASSIGN node");
	
        /* propagate values */
        tableNode->value.T = tableNode->loc->definition->value.T;
        tableNode->vtype = tableNode->loc->definition->vtype;
        tableNode->isValue = tableNode->loc->definition->isValue;
    }
    break; /* TABLE_VAR */


    case nTABLE_OP:
        /* combine two tables using UNION, INTERSECT or DIFF operators */
        /* both tables are compatible :-
           (1) matching number of columns - see aqlcheck.c:postTableProcess ()
           (2) matching column types as far as they were known prior to evaluation - see aqlrun.c:postTableTypeCheck ()
           (3) the table value of this node is already initialised - see aqlrun.c:postTableTypeCheck ()
        */
    {
        AqlFuzzyType isValueLeft, isValueRight;
        TABLE *leftTable, *rightTable;

        /* first evaulate the two tables */
        isValueLeft  = evalTable (tableNode->left, isFuzzyRun, aql, level+1);
        isValueRight = evalTable (tableNode->right, isFuzzyRun, aql, level+1);

        if (isValueLeft == KNOWNVAL && isValueRight == KNOWNVAL)
        {
            leftTable  = tableNode->left->value.T;
            rightTable = tableNode->right->value.T;

            /**********
             *** the following check are performed here, before enetering the tablePackage
             *** in order to give specific AQL errors instaed of just negative return codes
             *** from the tablePackage
             **********/

            /* we check again that the value-types of the columns match in both tables,
               as postTableTypeCheck might not have known all types yet prior to evaluation */
            if (strcmp (leftTable->type, rightTable->type) != 0)
                aqlError (aql, 962, tableNode->pos, 
                          "table combining operator -\n"
                          " value-type mismatch between the columns the two tables");

            /* the combined table value of the TABLE_OP-node depends on the operator */
            switch (tableNode->op)
            {
            case oUNION:		/* union - rows from both tables without duplicates */
                tableUnion (leftTable, rightTable, &tableNode->value.T);
                break;

            case oDIFF:		/* difference - rows that are in one buth not the other table */
                tableDiff (leftTable, rightTable, &tableNode->value.T);
                break;

            case oINTERSECT:	/* intersection - rows that are shared by both tables */
                tableIntersect (leftTable, rightTable, &tableNode->value.T);
                break;

            default:
	      aqlError (aql, 620, 0, "evalTable () - invalid operator for TABLE_OP node");
            }
            tableNode->isValue = KNOWNVAL;
        }
        else
        {
            if (isValueLeft == FUZZYVAL || isValueRight == FUZZYVAL)
                tableNode->isValue = FUZZYVAL;
            else
                tableNode->isValue = NULLVAL;
        }
    }
    break; /* TABLE_OP */


    case nTABLE_ORDER:
    {
        AqlNode *sortSpecNode;
        char    *sortSpecString;
        /* tableNode->right is a linked list of sort criterions, 
           any one of those might have an optional ->left node
           specifiying one of oASC or oDESC paremeter */
	
        /* first get the table value of the select-from-where expression */
        tableNode->isValue = evalTable (tableNode->left, isFuzzyRun, aql, level+1);
	
        /* propagate the table value that we're going to sort */
        tableNode->value.T = tableNode->left->value.T ;

        if (tableNode->isValue == KNOWNVAL)
        {
            /* make the sort specifier from the 
               sort criteria linked to this node */
            sortSpecString = 
                (char*)halloc ((tableNode->value.T->ncol*2)+1, aql->query->handle);

            if (!tableNode->right)
            {
                int i;

                /* table ordering, default columns */
                for (i = 0; i < tableNode->value.T->ncol; ++i)
                {
                    sortSpecString[i*2]   = tableNode->op == oASC ? '+' : '-';
                    sortSpecString[i*2+1] = '1'+i;
                }
            }
            else
                for (sortSpecNode = tableNode->right; 
                     sortSpecNode; sortSpecNode = sortSpecNode->nxt)
                {
                    if (!sortSpecNode->left) /* default is asc */
                        strcat (sortSpecString, 
                               messprintf ("+%1d", sortSpecNode->number));
                    else
                    {
                        switch (sortSpecNode->left->op)
                        {
                        case oASC:
                            strcat (sortSpecString, 
                                   messprintf ("+%1d", sortSpecNode->number));
                            break;
                        case oDESC:
                            strcat (sortSpecString, 
                                   messprintf ("-%1d", sortSpecNode->number));
                            break;
                        default:
			  aqlError (aql, 620, 0, "evalTable () - invalid sort spec operator");
                        }
                    }
                }
	
            if (!tableSort (tableNode->value.T, sortSpecString))
                /* although it shouldn't fail this is safer */
                tableNode->isValue = NULLVAL;
        }
    }
    break; /* TABLE_ORDER */
    default:
      aqlError (aql, 620, 0, "evalTable () - invalid type of table node") ;
    }

    tableNode->vtype = 'T';

    return tableNode->isValue;
} /* evalTable */




static AqlFuzzyType wakeUpLocal (AqlLoc *inLocator,
                                BOOL    isFuzzyRun,
                                AQL     aql,
                                int     level)
{
    if (inLocator->isValue != KNOWNVAL) /* do we have the tag in value.g ?? */
      aqlError (aql, 610, 0, "unknown inLocator value in wakeUpLocal");

    if (!inLocator->useObj)	/* open object, if necessary */
    { 
        AqlLoc *definitionLocator = inLocator->definition->left->loc ;
	    
        if (inLocator->source != IS_FOLLOW_TAG)
	  aqlError (aql, 610, 0, "wakeUpLocal () - Should have caught error 902 earlier");

        /* inLocator->value.g is the tag we're trying to find in the
         * object by which this loc is defined */
	  
        if (definitionLocator->isValue == NULLVAL)
            /* we don't have a KEY for the object */
            return NULLVAL;
	  
        if (isFuzzyRun)
        {
            /* Even in fuzzy mode it's not enough to just know that the tag is there
             * we need to open the object and position the obj-cursor on that tag
             * but in fuzzy mode we have to wait until next time round */
            return FUZZYVAL;
        }
        else
        {
            /* we open the object and position the cursor at the specified tag */
            if (!definitionLocator->trueObj)
            {
                if (aql->evalDebugLevel >= 4) {
                    indent (aql, level);
                    aceOutf (aql->debug_out, "wakeUpLocal : about to bsCreate (%x)\n",
                                 definitionLocator->value.k);
                }
                definitionLocator->trueObj = bsCreate (definitionLocator->value.k) ;
            }
	      
            if (!definitionLocator->trueObj)
                return NULLVAL; /* can't find object */
	      
            if (!bsFindTag (definitionLocator->trueObj, inLocator->value.g))
                return NULLVAL;	/* can't find tag */
      
            /* we can make it and we found the tag, so we mark the position where the tag was found
               in the object, which is the tag of inLocator, so that we can go back to it
               in the following else clause next time round */
            inLocator->useObj = definitionLocator->trueObj ;
            inLocator->mark = bsHandleMark (inLocator->useObj, inLocator->mark, aql->query->handle);
	      
            if (aql->evalDebugLevel >= 4) {
                indent (aql, level);
                aceOutf (aql->debug_out, "wakeUpLocal : locator %lx created and marked (%s) useObj = %lx (from trueObj of loc=%lx)\n",
                             (unsigned long)inLocator,
                             bsGetMarkAsString (inLocator->mark),
                             (unsigned long)inLocator->useObj, (unsigned long)definitionLocator);
            }
            return KNOWNVAL;
        }
    }
    else
        /* we already have the object to use (will never happen in fuzzy run) */
    {
        /* go back to the position where the tag was found
         * so we reposition the object cursor at the specified tag */
        bsGoto (inLocator->useObj, inLocator->mark) ;

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level);
            aceOutf (aql->debug_out, "wakeUpLocal : locator %lx jumped to mark (%s) in useObj = %lx\n",
                         (unsigned long)inLocator,
                         bsGetMarkAsString (inLocator->mark),
                         (unsigned long)inLocator->useObj);
        }

        return KNOWNVAL;
    }


    aqlError (aql, 600, 0, "wakeUpLocal - shouldn't get here");

    return 0;			/* compiler happiness */
} /* wakeUpLocal */

/***********************************************/


/* called on locators of type IS_FOLLOW_POS */
static AqlFuzzyType wakeUpFollow (AqlLoc *inLocator,
				  BOOL    isFuzzyRun,
				  AQL     aql,
                  int     level)
{
  if (!inLocator->trueObj)	
    {
      if (inLocator->type != 'k')
	/* e.g. caused by 'select a->full_name[1] */
	aqlError (aql, 900, inLocator->definition->pos, "Follow position : Trying to dereference (->) something not an object") ;

      if (isFuzzyRun)
	return FUZZYVAL;

      if (aql->evalDebugLevel >= 4) {
          indent (aql, level);
          aceOutf (aql->debug_out, "wakeUpFollow : about to bsCreate (%x)\n",
                       inLocator->value.k);
      }
      /* try to make the object whose KEY we know */
      inLocator->trueObj = bsCreate (inLocator->value.k) ;

      if (!inLocator->trueObj)	/* can't make it */
	return NULLVAL ;

      inLocator->mark = bsHandleMark (inLocator->trueObj, inLocator->mark, aql->query->handle);

      if (aql->evalDebugLevel >= 4) {
          indent (aql, level);
          aceOutf (aql->debug_out, "wakeUpFollow : locator %lx created and marked (%s) trueObj = %lx\n",
                       (unsigned long)inLocator, bsGetMarkAsString (inLocator->mark), (unsigned long)inLocator->trueObj);
      }
    }
  else
    /* we already have the object from last time round
     * so just move to the previously marked position
     we'll never have an object in fuzzy mode */
    {
      bsGoto (inLocator->trueObj, inLocator->mark);

      if (aql->evalDebugLevel >= 4) {
          indent (aql, level);
          aceOutf (aql->debug_out, "wakeUpFollow : locator %lx jumped to mark (%s) in trueObj = %lx\n", 
                       (unsigned long)inLocator, bsGetMarkAsString (inLocator->mark), (unsigned long)inLocator->trueObj);
      }
    }

  return KNOWNVAL;
} /* wakeUpFollow */

/***********************************************/

static BOOL classNextKey (KEY theClass, KEY *keyp, char *template)
     /* currently uses lexNext to iterate over the KEYs in theClass
      * To get the first value the *keyp value has to be primed with 0 */
     /* template is NULL-string to loop over whole class */
{
  BOOL isNext = FALSE;

  if (template && *template)
    {
      /* iterate with template */
      while ((isNext = lexNext (theClass, keyp)))
	{
	  if (pickMatchCaseSensitive (name (*keyp), template,
				      pickCaseSensitive (*keyp)) > 0)
	    break;
	}
    }
  else
    {
      /* standard iterator over whole class */
      isNext = lexNext (theClass, keyp);
    }

  return isNext;
} /* classNextKey */


static AqlFuzzyType loopVariables (AqlNode *tableSfwNode,
                                  AqlNode *fromFieldNode,
                                  BOOL     isFuzzyRun,
                                  AQL      aql,
                                  int      level)
/* "the guts" - the inner recursive loop for each FROM_LOC, 
 * fill fromFieldNode->loc->isValue or fromFieldNode->loc->value and fromFieldNode->loc->type, 
 * iterating if necessary, not forgetting the NULL value once set, 
 * evaluate condition if it is there, then if happy either recurse on next
 * variable, or if none, you are happy, and should report result-values
*/

/* called from evalTable with fromFieldNode being the node->left of a TABLE_SFW-node 
 * it then traverses through the FROM_LOCs in processThisLoc, which calls it on the ->nxt's*/
{
    AqlFuzzyType status = NULLVAL ;
    int counter = 0;

    if (aql->evalDebugLevel >= 3) {
        indent (aql, level); aceOutf (aql->debug_out, "loopVariables (%s):\n", isFuzzyRun ? "fuzzy" : "definite");
        if (aql->evalDebugLevel >= 3) {
            aqlNodeOut (aql->debug_out, tableSfwNode, level, "sfw", TRUE);
            aqlNodeOut (aql->debug_out, fromFieldNode, level, "from", TRUE);
        }
    }

    switch (fromFieldNode->loc->source)
    {
        /* the first few are single valued, i.e. no looping over multi-values */

    case IS_COLUMN:
    { 
        /* this is the locator of the table we're looking at */
        AqlNode *fromLocNode = fromFieldNode;
        AqlNode *fromTableNode = fromFieldNode->loc->definition->left;

        if (!fromTableNode->loc->table)
	  aqlError (aql, 620, 0,"loopVariables ():IS_COLUMN - \n"
		    "the ->table value of the locator is NULL");

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "eval IS_COLUMN in Loc %lx\n", (unsigned long)fromFieldNode->loc) ;
        }

        if ((fromTableNode->loc->isValue == NULLVAL) ||
            (fromTableNode->loc->isValue == FUZZYVAL)) {
            /* unevaluated table or fuzzy table value */
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       (fromTableNode->loc->isValue == FUZZYVAL), (fromTableNode->loc->isValue == NULLVAL), status);
        }
        
        /* check the current row in the table for a value */
        if (tabEmpty (fromTableNode->loc->table, fromTableNode->loc->row, fromLocNode->loc->column)) {
            /* no value in the table in current row */
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, TRUE, status);
        }
        
        /* we have a value in the table */
        fromLocNode->loc->isValue = KNOWNVAL ;
	
        /* set the type of the loc by taking the type of the table column it is derived from */
        fromLocNode->loc->type = fromTableNode->loc->table->type[fromLocNode->loc->column] ; /* can I predo that?,
                                                                                              * I only need to do this once !! */

        {	
            TABLE *table = fromTableNode->loc->table;
            int row = fromTableNode->loc->row;
            int column = fromLocNode->loc->column;
            /* now get the value of the loc from the current row of the table */
            /* RD says: */
            /* can not assign union: TABLETYPE is 32 bit and AqlType can be 64 bit 
               (can have ptrs).  Only values will be set this way, so this is right */
	
            switch (fromLocNode->loc->type)
            {
            case 'k':	   
            case 'g': fromLocNode->loc->value.k = tabKey (table, row, column); break ;
            case 'i': fromLocNode->loc->value.i = tabInt (table, row, column); break ;
            case 'f': fromLocNode->loc->value.f = tabFloat (table, row, column); break ;
            case 't': fromLocNode->loc->value.t = tabDate (table, row, column); break ;
            case 'b': fromLocNode->loc->value.b = tabBool (table, row, column); break ;
            case 's': fromLocNode->loc->value.s = tabString (table, row, column); break ;
	    
            case '0':		/* uninitialised type  */
                /* this might happen, if this column of the table is an expression, that doesn't 
                   evaluate in this case. The type of the column is therefore not known  but will 
                   probably be set later on, if the expression evaluates TRUE with other values later on */
                /* So we just set this value ZERO */
                fromLocNode->loc->value.i = 0;
                break;
	    
            default:  
	      aqlError (aql, 620, 0,
			messprintf ("loopVariables ():IS_COLUMN - \n"
				    "type %c not supported for loc->type",
				    fromLocNode->loc->type ? fromLocNode->loc->type : '_')
			) ;
            }

            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, FALSE, status);
        }
    }



    case IS_METHOD:
    {
        AqlNode *methodNode = fromFieldNode->loc->definition->left;

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "eval IS_METHOD in Loc %lx\n", (unsigned long)fromFieldNode->loc) ;
        }

        if ((methodNode->loc->isValue == NULLVAL) ||
            (methodNode->loc->isValue == FUZZYVAL))
        {
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       (methodNode->loc->isValue == FUZZYVAL),
                                       (methodNode->loc->isValue == NULLVAL),
                                       status);
        }

        if (strcasecmp (methodNode->nxt->name, "name") == 0)
        {
            if (methodNode->loc->type == 'k' || methodNode->loc->type == 'g')
            {
                fromFieldNode->loc->type = 's';
                fromFieldNode->loc->isValue = KNOWNVAL;
                fromFieldNode->loc->value.s = strnew (name (methodNode->loc->value.k), aql->query->handle);
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           FALSE, FALSE, status);
            }
            else
            {
                aqlError (aql, 860, methodNode->pos, "invalid valuetype of operand for method 'name'\n"
                          "(operands can only be object- or tag-values)");
            }
        }
        else if (strcasecmp (methodNode->nxt->name, "class") == 0)
        {
            if (methodNode->loc->type == 'k')
            {
                if (!methodNode->loc->value.k)
                    return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                               FALSE, TRUE, status);

                fromFieldNode->loc->value.s = strnew (className (methodNode->loc->value.k), aql->query->handle);
                fromFieldNode->loc->type = 's';
                fromFieldNode->loc->isValue = KNOWNVAL;
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           FALSE, FALSE, status);
            }
            else
            {
                aqlError (aql, 861, methodNode->nxt->pos, "invalid valuetype of operand for method 'class'\n"
                          "(operands can only be object-values)");
            }
        }
        else if (strcasecmp (methodNode->nxt->name, "length") == 0)
        {
            if (methodNode->loc->type == 's')
            {
                fromFieldNode->loc->value.i = strlen (methodNode->loc->value.s);
                fromFieldNode->loc->type = 'i';
                fromFieldNode->loc->isValue = KNOWNVAL;
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           FALSE, FALSE, status);
            }
            else if (methodNode->loc->type == 'k')
            {
                Array ar;
                int max;

                /* The key may be an object carrying an array-class object,
                 * it's length in the number of elements */
                if (!methodNode->loc->value.k)
                    return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                               FALSE, TRUE, status);

                /* we need to find out whether it is an array-class key */
                if (pickType (lexAliasOf (methodNode->loc->value.k)) != 'A')
                {
                    aqlError (aql, 863, methodNode->nxt->pos,
                              "invalid valuetype of operand for method 'length'\n"
                              "(The object-value does not belong to an A-class)");
                }


                ar = arrayGet (methodNode->loc->value.k, char, "c") ;
                if (!ar)
                    return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                               FALSE, TRUE, status);
                max = arrayMax (ar);
                arrayDestroy (ar);

                fromFieldNode->loc->value.i = max;
                fromFieldNode->loc->type = 'i';
                fromFieldNode->loc->isValue = KNOWNVAL;
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           FALSE, FALSE, status);
            }
            else
            {
                aqlError (aql, 862, methodNode->nxt->pos,
                          "invalid valuetype of operand for method 'length'\n"
                          "(operands can only be Text-values or A-class objects)");
            }
        }
        else if (strcasecmp (methodNode->nxt->name, "create_time") == 0 ||
                 strcasecmp (methodNode->nxt->name, "update_time") == 0 ||
                 strcasecmp (methodNode->nxt->name, "create_session") == 0 ||
                 strcasecmp (methodNode->nxt->name, "update_session") == 0)
        {
            if (methodNode->loc->type != 'k')
                aqlError (aql, 864, methodNode->nxt->pos,
                          "invalid valuetype of operand for method '%s'\n"
                          "(operands can only be object-values)", methodNode->nxt->name);

            /* switch on the first and 7th letter */
            switch (ace_lower (methodNode->nxt->name[0]) + ace_lower (methodNode->nxt->name[7]))
            {
#ifdef JUNK
 not implemented in mieg code
            case 'c'+'t':	/* create_time */
                fromFieldNode->loc->value.t = lexCreationStamp (methodNode->loc->value.k);
                fromFieldNode->loc->type = 't';
                break;

            case 'u'+'t':	/* update_time */
                fromFieldNode->loc->value.t = lexUpdateStamp (methodNode->loc->value.k);
                fromFieldNode->loc->type = 't';
                break;

            case 'c'+'s':	/* create_session */
                fromFieldNode->loc->value.k = lexCreationUserSession (methodNode->loc->value.k);
                fromFieldNode->loc->type = 'k';
                break;
		
            case 'u'+'s':	/* update_session */
                fromFieldNode->loc->value.k = lexUpdateUserSession (methodNode->loc->value.k);
                fromFieldNode->loc->type = 'k';
                break;
#endif		
            default:;
            }

            /* the return values of any of the lex-functions can be ZERO,
             * in which case the result will be NULL */
            if (!fromFieldNode->loc->value.k)
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           FALSE, TRUE, status);

            fromFieldNode->loc->isValue = KNOWNVAL;
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, FALSE, status);
        }
        else if (strcasecmp (methodNode->nxt->name, "node_session") == 0 ||
                 strcasecmp (methodNode->nxt->name, "timestamp") == 0 ||
                 strcasecmp (methodNode->nxt->name, "node_time") == 0)
        {
            if (methodNode->loc->source == IS_LOCAL_TAG || 
                methodNode->loc->source == IS_FOLLOW_TAG ||
                methodNode->loc->source == IS_LOCAL_POS ||
                methodNode->loc->source == IS_FOLLOW_POS)
            {
                if (!methodNode->loc->useObj)
                    return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                               FALSE, TRUE, status);

                /* get timestamp from node at obj->curr */

                switch (ace_lower (methodNode->nxt->name[0]) + ace_lower (methodNode->nxt->name[5]))
                {
                case 'n'+'s':
                case 't'+'t':
                    fromFieldNode->loc->value.k = bsGetTimeStamp (methodNode->loc->useObj);
                    fromFieldNode->loc->type = 'k';
                    break;
#ifdef BSGETNODETIME		/* bsGetNodeTime not written yet */
                case 'n'+'t':
                    fromFieldNode->loc->value.t = bsGetNodeTime (methodNode->loc->useObj);
                    fromFieldNode->loc->type = 't';
                    break;
#endif /* BSGETNODETIME */

                default:;
                }

                if (!fromFieldNode->loc->value.k) /* no value */
                    return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                               FALSE, TRUE, status);
		
                fromFieldNode->loc->isValue = KNOWNVAL;
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           FALSE, FALSE, status);
            }
            else
            {
                /* we can only get the node-timestamp of an open object,
                 * a locator-value created by means other than
                 * LOCAL_TAG etc. will not give us an object-cursor
                 * to report a node-timestamp from */
                aqlError (aql, 865, methodNode->nxt->pos,
                          "invalid valuetype of operand for method '%s'\n"
                          "(operands can only be references to object-nodes)",
                          methodNode->nxt->name);
            }
        }
        else
            aqlError (aql, 703, methodNode->nxt->pos+1, 
		      messprintf ("syntax error : unrecognised Method-name '%s'", methodNode->nxt->name)
		      );
    }


    case IS_OBJECT:
        /* try to find loc->value.k which is the KEY to the object specified in the constructor */
        /* we won't actually bsCreate the object until it's needed - we just find its KEY here */
    { 
        AqlNode *objNode = fromFieldNode->loc->definition->left ;

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "eval IS_OBJECT in Loc %lx\n", (unsigned long)fromFieldNode->loc) ;
        }

        /****** left expression is the class of the object ************/

        /*** text literals for the class have already been ************
         *** converted to nKEY nodes in postAssignTypesAndClasses () ***/

        if (evalExpr (objNode->left, isFuzzyRun, aql, level+1) == NULLVAL) 
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, TRUE, status);

        if (objNode->left->vtype != 'k')
            aqlError (aql, 952, objNode->left->pos,"invalid value-type of left expression in object constructor");

        if (!objNode->left->value.k) /* the class KEY */
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, TRUE, status);

        fromFieldNode->loc->key = objNode->left->value.k;

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "left     expr value=%s (key=%x)\n",
                                             name (objNode->left->value.k), objNode->left->value.k);
        }

        /******* right expression is the key for the object *********/

        if (evalExpr (objNode->right, isFuzzyRun, aql, level+1) == NULLVAL)
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, TRUE, status);

        if (objNode->right->vtype == 's')		/* text-literal for key */
        {
            /* set objNode->loc->value.k */
            if (!lexword2key (objNode->right->value.s, &fromFieldNode->loc->value.k, fromFieldNode->loc->key))
                /* no such object - should we output an error message ???? XXXXXX */
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           FALSE, TRUE, status);
        }
        else if (objNode->right->vtype != 'k')
            /* only strings and keys allowed */
            aqlError (aql, 953, objNode->right->pos, "invalid value-type of right expression in object constructor");


        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "right    expr value=%s (key=%x)\n",
                                             name (objNode->right->value.k), objNode->right->value.k);
        }

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "eval done loc_value=%s (key=%x)\n",
                                             name (fromFieldNode->loc->value.k), fromFieldNode->loc->value.k);
        }

        fromFieldNode->loc->isValue = KNOWNVAL;
        return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                   FALSE, FALSE, status);
    }

    case IS_EXPR:
    {
        AqlNode *exprNode = fromFieldNode->left;

        evalExpr (exprNode, isFuzzyRun, aql, level+1);

        if ((exprNode->isValue == NULLVAL) || (exprNode->isValue == FUZZYVAL))
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       (exprNode->isValue == FUZZYVAL), (exprNode->isValue == NULLVAL), status);

        fromFieldNode->loc->type = exprNode->vtype;
        fromFieldNode->loc->value = exprNode->value;
        fromFieldNode->loc->isValue = KNOWNVAL;
        return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                   FALSE, FALSE, status);
    }

    case IS_FOLLOW_TAG:
    { 
        AqlLoc *definitionLocator = fromFieldNode->loc->definition->left->loc ; /* the loc of the LOC_VAR node of this FROM_LOC */
        int     tagFindResult ;

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aqlLocOut (aql->debug_out, "eval IS_FOLLOW_TAG in Loc", fromFieldNode->loc) ;
            indent (aql, level); aceOutf (aql->debug_out, " - FROM_LOC: %s", fromFieldNode->name);
            aqlNodeOut (aql->debug_out, fromFieldNode, level, "follow from", TRUE) ;
#if 0
            aceOutf (aql->debug_out, "eval IS_FOLLOW_TAG in Loc %lx - FROM_LOC %s (%lx)\n",
                         (unsigned long)fromFieldNode->loc,
                         fromFieldNode->name, 
                         (unsigned long)fromFieldNode) ;
#endif
        }

        if ((definitionLocator->isValue == NULLVAL) ||
            (definitionLocator->isValue == FUZZYVAL)) /* do we have an object to follow up? */
        {
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       (definitionLocator->isValue == FUZZYVAL),
                                       (definitionLocator->isValue == NULLVAL),
                                       status);
        }

        /* Now we know the KEY of the object we're following up */

        if (definitionLocator->type != 'k')	/* A tag has to FOLLOW a class or a key */
            /* definitionLocator->definition is a FROM_LOC, that would usually be defined upon a class or an object or so */
            /* caused by 'select a->full_name->born from a in class author' (full_name[1] is a string value) */
            aqlError (aql, 901, fromFieldNode->pos, "Follow tag : Trying to dereference (->) something not an object") ;

        /* look for the tag -- this is the heart of Follow */
        /* does the tag exist in the definition's key ? */
        tagFindResult = bIndexFind (definitionLocator->value.k,    /* (The KEY of the object) */
                                    fromFieldNode->loc->value.g) ; /* (The KEY of the tag we're looking for in the object)  */

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "followtag looking for tag %x (\"%s\") in object %x (\"%s\") -> (%s)\n",
                                             fromFieldNode->loc->value.g, name (fromFieldNode->loc->value.g), 
                                             definitionLocator->value.k, name (definitionLocator->value.k), 
                                             tagFindResult == 2 ? "exists" : (tagFindResult == 1 ? "might exist" : "doesn't exist")) ;
        }

        switch (tagFindResult)
        {
        case 0:		/* tag doesn't exist */
            /* would have triggered an aqlError No 944 earlier in postAssignTypesAndClasses () */
            if (aql->evalDebugLevel >= 4) {
                indent (aql, level); aceOutf (aql->debug_out, "eval done by follow: tag absent\n");
            }

            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, TRUE, status);


        case 1:		/* tag might exist - open object now and check */
            if (isFuzzyRun)
            {
                if (aql->evalDebugLevel >= 4) {
                    indent (aql, level); aceOut (aql->debug_out, "eval done : fuzzy\n");
                }
                /* next time round we'll have to open the object */
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           TRUE, FALSE, status);
            }
            else
            {
                if (!definitionLocator->trueObj)
                {
                    definitionLocator->trueObj = bsCreate (definitionLocator->value.k) ;
		    
                    if (aql->evalDebugLevel >= 4) {
                        indent (aql, level);
                        aceOutf (aql->debug_out, "IS_FOLLOW_TAG opening object: definition_locator %lx created trueObj %lx (\"%s\")\n",
                                    (unsigned long)definitionLocator,
                                    (unsigned long)definitionLocator->trueObj,
                                    name (definitionLocator->value.k));
                        bsDump (definitionLocator->trueObj); /* goes to STDERR */
                        aceOut (aql->debug_out, "\n");
                    }
                }
              
                if (!definitionLocator->trueObj)
                    return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                               FALSE, TRUE, status);
              
                if (!bsFindTag (definitionLocator->trueObj, fromFieldNode->loc->value.k))
                    return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                               FALSE, TRUE, status);
              
                fromFieldNode->loc->useObj = definitionLocator->trueObj ;
              
                if (fromFieldNode->loc->mark)
                { /* we already have a mark and need to save it */
                    if (!stackExists (fromFieldNode->loc->markStack))
                        fromFieldNode->loc->markStack = stackHandleCreate (3, aql->query->handle);
                  
                    push (fromFieldNode->loc->markStack, fromFieldNode->loc->mark, BSMARK);
                    fromFieldNode->loc->mark = 0; /* so bsHandleMark will make a new one */
                }
                fromFieldNode->loc->mark = bsHandleMark (fromFieldNode->loc->useObj, fromFieldNode->loc->mark, aql->query->handle);
              
                if (aql->evalDebugLevel >= 4) {
                    indent (aql, level); aceOutf (aql->debug_out, "eval done opening object: from_locator %lx got useObj = %lx and marked (%s)\n",
                                                    (unsigned long)fromFieldNode->loc, (unsigned long)fromFieldNode->loc->useObj,
                                                    bsGetMarkAsString (fromFieldNode->loc->mark));
                }
            }

            break;


        case 2:		/* tag does exist */
            fromFieldNode->loc->useObj = 0 ;	/* we don't need an object to know that this tag is here */

            if (aql->evalDebugLevel >= 4) {
                indent (aql, level); aceOutf (aql->debug_out, "eval done : from_locator %lx zeroes useObj\n",
                                                (unsigned long)fromFieldNode->loc);
            }

            break;

        default:
	  aqlError (aql, 630, 0,
		    "loopVariables () - bIndexFind () returned a value other than 0,1,2");
        }

        fromFieldNode->loc->isValue = KNOWNVAL;	/* type and value already set */
        return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                   FALSE, FALSE, status);
    }


    case IS_LOCAL_TAG:
    { 
        AqlNode *definitionNode = fromFieldNode->loc->definition->left;
        AqlLoc *definitionLocator = definitionNode->loc ;
        AqlFuzzyType tagStatus;

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level);
            aceOutf (aql->debug_out, "eval IS_LOCAL_TAG in Loc %lx\n", (unsigned long)fromFieldNode->loc) ;
        }

        if ((definitionLocator->isValue == NULLVAL) ||
            (definitionLocator->isValue == FUZZYVAL))
        {
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       (definitionLocator->isValue == FUZZYVAL),
                                       (definitionLocator->isValue == NULLVAL),
                                       status);
        }

        if (!definitionNode->loc->useObj &&
            definitionLocator->source != IS_FOLLOW_TAG)
            /* select all class movie where exists_tag [cast] */
            aqlError (aql, 902, definitionNode->nxt->pos, "Positioning by name within an object has to follow a tag") ;
 
        tagStatus = wakeUpLocal (definitionLocator, isFuzzyRun, aql, level+1);

        if ((tagStatus == NULLVAL) ||
            (tagStatus == FUZZYVAL))
        {
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       (tagStatus == FUZZYVAL), (tagStatus == NULLVAL), status);
        }
        else
        {
            fromFieldNode->loc->useObj = definitionLocator->useObj ;
            bsPushObj (fromFieldNode->loc->useObj) ; /* changes model for #construct if necessary, otherwise harmless */
	    
            if (!bsFindTag (fromFieldNode->loc->useObj, fromFieldNode->loc->value.k))
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           FALSE, TRUE, status);
	    
            fromFieldNode->loc->mark = bsHandleMark (fromFieldNode->loc->useObj, fromFieldNode->loc->mark, aql->query->handle);
            fromFieldNode->loc->isValue = KNOWNVAL ;	/* type and value already set */
	    
            if (aql->evalDebugLevel >= 4) {
                indent (aql, level);
                aceOutf (aql->debug_out, "loc=%lx - moved to tag %s in useObj = %lx and marked (%s)\n",
                             (unsigned long)fromFieldNode->loc,
                             name (fromFieldNode->loc->value.k),
                             (unsigned long)fromFieldNode->loc->useObj,
                             bsGetMarkAsString (fromFieldNode->loc->mark));
            }
	    
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, FALSE, status);
        }
    }

    /**********************************************************************/
    /* the rest are potentially multiple valued, 
     * so if there is data they "break" into the "while (TRUE)" loop below */
    /**********************************************************************/



    case IS_ROW:		/* NB: in declarations only - not a real loc */
        /* Note, the actual table comes from an assign-statement, which 
         * has the value of the TABLE_SFW node that created the table */
    {
        AqlNode *fromTableNode = fromFieldNode;
        /*AqlNode *fromLocNode =  fromFieldNode->loc->definition->nxt;*/
        AqlNode *fromLocNode = fromFieldNode->nxt;
        AqlNode *tableVarNode = fromFieldNode->left;

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "eval IS_ROW in Loc %lx (column %d)\n",
                                            (unsigned long)fromTableNode->loc,
                                            fromLocNode->loc->column);
        }
	
        /* get the value of the actual table for our locator */
        if (fromTableNode->left->loc->isContext)
        {
            fromTableNode->loc->table = tableVarNode->loc->value.T; /* context variables have no defNode */

            if (!fromTableNode->loc->table)
                /* table is NULL - can only happen with context-vars */
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           FALSE, TRUE, status);
        }
        else
        {
            fromTableNode->loc->table = tableVarNode->loc->definition->value.T;

            if (tableVarNode->loc->definition->isValue == FUZZYVAL)
                /* This may be a TABLE_ASSIGN which can still be fuzzy , cannot happen with context-vars */
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           TRUE, FALSE, status);
        }	
	
	
        if (tableMax (fromTableNode->loc->table) == 0)
            /* all columns of the table are empty */
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, TRUE, status);

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "zero row-counter for column %d\n",
                                            fromLocNode->loc->column);
        }
	
        fromTableNode->loc->row = 0 ; /* set start row number to ZERO (that's where we the row iterator starts from) */
    }
    break ;			/* go to while (TRUE) loop */


    case IS_CLASS:		
    {
        BOOL isNext = FALSE;

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aqlLocOut (aql->debug_out, "eval IS_CLASS in Loc:", fromFieldNode->loc) ;
        }

        /* fromFieldNode->loc->type was set in postAssignTypesAndClasses */
        /* the type is 'k' and its values are found by a loop through lexNext () */
        fromFieldNode->loc->value.k = 0 ;		/* prime lexNext () with 0 */

        isNext = classNextKey (fromFieldNode->loc->key, &fromFieldNode->loc->value.k,
                               fromFieldNode->loc->definition->left->nxt ?
                               fromFieldNode->loc->definition->left->nxt->value.s : 0);

        if (!isNext)
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, TRUE, status);
    }
    break ;			/* go to while (TRUE) loop */


    case IS_LOCAL_POS: 
    case IS_FOLLOW_POS:
    { 
        AqlNode *definitionNode = fromFieldNode->loc->definition->left ;

        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "eval IS_ (LOCAL/FOLLOW)_POS in Loc %lx\n", (unsigned long)fromFieldNode->loc) ;
            indent (aql, level); aceOutf (aql->debug_out, "locname is %s\n", fromFieldNode->loc->definition->name);
        }
          
        if ((definitionNode->loc->isValue == NULLVAL) ||
            (definitionNode->loc->isValue == FUZZYVAL))
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       (definitionNode->loc->isValue == FUZZYVAL),
                                       (definitionNode->loc->isValue == NULLVAL),
                                       status);

        /* so we know the object we're talking about - 
         * it is the KEY in definitionNode->loc->value.k */
          
        /* if we derefence or position by one we need to open the object and look inside */
        if (fromFieldNode->loc->source == IS_LOCAL_POS)
        {
            AqlFuzzyType tagStatus;
            /* make the object (definitionNode->loc->useObj
               (from its definition's trueObj)
               and try to find the tag */
              
            if (!definitionNode->loc->useObj &&
                definitionNode->loc->source != IS_FOLLOW_TAG)
                /* select a[1] from a in class Movie */
                aqlError (aql, 902, definitionNode->nxt->pos, "Positioning within an object has to follow a tag") ;
              
            if (aql->evalDebugLevel >= 4) {
                indent (aql, level); aceOutf (aql->debug_out, "local/follow_pos waking up local\n") ;
            }
              
            tagStatus = wakeUpLocal (definitionNode->loc, isFuzzyRun, aql, level+1);
              
            if ((tagStatus == NULLVAL) ||
                (tagStatus == FUZZYVAL))
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           (tagStatus == FUZZYVAL),
                                           (tagStatus == NULLVAL),          /* goto empty_value; */
                                           status);
            else
            {
                /* locator's useObj internal cursor is now positioned at the specified tag */
                /* won't happen during fuzzy eval mode */
                fromFieldNode->loc->useObj = definitionNode->loc->useObj ;
                bsPushObj (fromFieldNode->loc->useObj) ; /* changes model for #construct if necessary, otherwise harmless */
            }
        }
        else /* IS_FOLLOW_POS */
        {
            AqlFuzzyType objStatus;
            /* try to make the object whose KEY we know (definitionNode->loc->value.k),
             * and mark the position within the object */
            objStatus = wakeUpFollow (definitionNode->loc, isFuzzyRun, aql, level+1);
              
            if ((objStatus == NULLVAL) ||
                (objStatus == FUZZYVAL))
                return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                           (objStatus == FUZZYVAL), (objStatus == NULLVAL),
                                           status);
            else
                /* no need to pushObj here, since at root */
                fromFieldNode->loc->useObj = definitionNode->loc->trueObj ; 
        }
        /* Can't get here in fuzzy mode - 
         * we're now guarantedd to have the fromFieldNode->loc->useObj 
         * internally positioned at the required tag (or obj-root) */
          
          
        /* move one position right and pick up the key to the data-object */
        /* this movement executes the ->1 or [1] */
        if (!bsGetKeyTags (fromFieldNode->loc->useObj, _bsRight, &fromFieldNode->loc->key)) 
            return loopVariablesSingle (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level,
                                       FALSE, TRUE, status);
          
        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "fromFieldNode->loc->key is %s (%x)\n",
                                            name (fromFieldNode->loc->key), fromFieldNode->loc->key) ;
        }
          
        if (fromFieldNode->loc->key == _Int) 
            fromFieldNode->loc->type = 'i' ;
        else if (fromFieldNode->loc->key == _Float) 
            fromFieldNode->loc->type = 'f' ;
        else if (fromFieldNode->loc->key == _DateType) 
            fromFieldNode->loc->type = 't' ;
        else if (fromFieldNode->loc->key < _LastC)
            fromFieldNode->loc->type = 's' ;
        else
            fromFieldNode->loc->type = (class (fromFieldNode->loc->key)) ? 'k' : 'g' ; 
          
        if (fromFieldNode->loc->mark)
        { /* we already have a mark and need to save it
           * before we start iterating over the multivalues in this loc */
            if (!stackExists (fromFieldNode->loc->markStack))
                fromFieldNode->loc->markStack = stackHandleCreate (3, aql->query->handle);
              
            push (fromFieldNode->loc->markStack, fromFieldNode->loc->mark, BSMARK);

            if (aql->evalDebugLevel >= 4) {
                indent (aql, level); aceOutf (aql->debug_out, "loc = %lx - pushed existing mark (%s) in useObj\n",
                                                (unsigned long)fromFieldNode->loc,
                                                bsGetMarkAsString (fromFieldNode->loc->mark),
                                                (unsigned long)fromFieldNode->loc->useObj);
            }
            fromFieldNode->loc->mark = 0; /* so bsHandleMark will make a new one */
        }

        /* marking allows you to bsGoto back to the same position in the object */
        fromFieldNode->loc->mark = bsHandleMark (fromFieldNode->loc->useObj, fromFieldNode->loc->mark, aql->query->handle);
	
        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "eval done : from_locator %lx gets value from useObj = %lx and marked (%s)\n",
                                            (unsigned long)fromFieldNode->loc, (unsigned long)fromFieldNode->loc->useObj,
                                            bsGetMarkAsString (fromFieldNode->loc->mark));
        }

        switch (fromFieldNode->loc->type)
        {
        case 'k': 
        case 'g': bsGetKeyTags (fromFieldNode->loc->useObj, _bsHere, &fromFieldNode->loc->value.k); break;
        case 'i': bsGetData (fromFieldNode->loc->useObj, _bsHere, _Int, &fromFieldNode->loc->value.i); break;
        case 'f': bsGetData (fromFieldNode->loc->useObj, _bsHere, _Float, &fromFieldNode->loc->value.f); break;
        case 't': bsGetData (fromFieldNode->loc->useObj, _bsHere, _DateType, &fromFieldNode->loc->value.t); break;
        case 's': bsGetData (fromFieldNode->loc->useObj, _bsHere, _Text, &fromFieldNode->loc->value.s); break;
        }
        /* we have found a value */
        break ;			/* go to while (TRUE) loop */

    }	/* end-case IS_LOCAL:FOLLOW_POS: */

    default:
      aqlError (aql, 620, 0,
		messprintf ("loopVariables () - invalid fromFieldNode->loc->source %s",
			    aqlLocSourceTypeName (fromFieldNode->loc->source))
		);
    } /* end-switch fromFieldNode->loc->source */


    /*********************************************************************/


    /* only get here if we have data */
    fromFieldNode->loc->isValue = KNOWNVAL ; 

    /* About to iterate over values, set up something to accept the results */

    if (aql->evalDebugLevel >= 3) {
        BOOL am_sorting = tableSfwNode->resultConsumer->sorting;
        indent (aql, level);
        aceOutf (aql->debug_out,
                    "About to iterate over values for %x, set up something to accept the results, am_sorting=%s\n",
                    tableSfwNode,
                    am_sorting ? "true" : "false");
    }

    /* this loop evaluates the multi-valued'ness of a locator */

    while (TRUE)			/* escape only by return in iteration step, no breaks */
    {
        AqlFuzzyType locStatus;
        
        if (messIsInterruptCalled ())
            aqlError (aql, 0, 0, "", 0);

        /* process the current loc and all its follow-on nodes,
           it returns once the end of the FROM_LOC list is reached and a 
           new row has been added to the result-table */
        if (aql->evalDebugLevel >= 3) {
            indent (aql, level); aceOutf (aql->debug_out, "multi value on iteration %d", counter);
            aqlNodeOut (aql->debug_out, fromFieldNode, level, "on loop \"from\" field", TRUE);
#if 0
            aqlNodeOut (aql->debug_out, tableSfwNode, level, "and loop body", TRUE);
#endif
            indent (aql, level); aqlLocOut (aql->debug_out, "with locator assignment", fromFieldNode->loc) ;
        }
        counter++;
        locStatus = processThisLoc (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level+1);

        /* incorporate the locator's result in the overall status
           NULL status can become KNOWN, but FUZZY makes and stays FUZZY */
        if (!fromFieldNode->nxt)
        {
            if (! (status == FUZZYVAL && locStatus == NULLVAL)) /* don't allow FUZZY to be unset to NULL */
                status = locStatus;		/* we've just reported to table - it's status is important */
        }
        else
        {
            if ((locStatus == FUZZYVAL) ||
                (status == NULLVAL && locStatus == KNOWNVAL))
                status = locStatus;
        }

        if (isFuzzyRun && (status == FUZZYVAL))
        {
            /* there is no need trying any other value, because although we had a value
             * to work with, the processThisLoc () on it went fuzzy, so we rather return early */
            return status;
        }


        /* we've reported one row, so the current value is exhausted and we grab another value for evalution, 
           a new key from the class, the next row in the table or the next item in the tag */

        switch (fromFieldNode->loc->source)
        {
            /***** get a next key from the class */
        case IS_CLASS:
        {
            BOOL isNext;

            if (aql->evalDebugLevel >= 4) {
                indent (aql, level); aqlLocOut (aql->debug_out, "next CLASS in loc", fromFieldNode->loc) ;
            }

            isNext = classNextKey (fromFieldNode->loc->key, &fromFieldNode->loc->value.k,
                                   fromFieldNode->loc->definition->left->nxt ?
                                   fromFieldNode->loc->definition->left->nxt->value.s : 0);

            if (!isNext)
            {
                if (aql->evalDebugLevel >= 3) {
                    indent (aql, level); aceOut (aql->debug_out, "end CLASS multi value\n");
                }
                return status ;
            }
        }
        break ;


        /**** advance to the next row in the table */
        case IS_ROW:
        {
            AqlNode *fromTableNode = fromFieldNode;
            /*AqlNode *fromLocNode =  fromFieldNode->loc->definition->nxt;*/
            AqlNode *fromLocNode =  fromFieldNode->nxt;

            if (aql->evalDebugLevel >= 4) {
                indent (aql, level); aceOutf (aql->debug_out, "next ROW in loc=%lx column %d - current at row %d\n",
                                                (unsigned long)fromTableNode->loc,
                                                fromLocNode->loc->column,
                                                fromTableNode->loc->row);
            }

            /* go to next row in table and stop when the table-column is exhausted */
            if (++fromTableNode->loc->row >= tabMax (fromTableNode->loc->table, fromLocNode->loc->column)) 
            {
                if (aql->evalDebugLevel >= 3) {
                    indent (aql, level); aceOutf (aql->debug_out, "end ROW multi value in column %d\n",
                                                    fromLocNode->loc->column);
                }

                return status ;
            }
        }
        break ;


        /**** get the next item in the tag of the object */
        case IS_LOCAL_POS:
        case IS_FOLLOW_POS:
            if (aql->evalDebugLevel >= 4) {
                indent (aql, level); aceOutf (aql->debug_out, "next POS in loc=%lx - jumped to mark (%s) useObj = %lx\n",
                                                (unsigned long)fromFieldNode->loc,
                                                bsGetMarkAsString (fromFieldNode->loc->mark),
                                                (unsigned long)fromFieldNode->loc->useObj);
            }

            bsGoto (fromFieldNode->loc->useObj, fromFieldNode->loc->mark) ;

            if (stackExists (fromFieldNode->loc->markStack))
            { /* we had a mark before that one, so we restore it */
                fromFieldNode->loc->mark = pop (fromFieldNode->loc->markStack, BSMARK);
                if (stackEmpty (fromFieldNode->loc->markStack))
                    stackDestroy (fromFieldNode->loc->markStack);
            }

            /* RD says: the following would be much better written accessing obj->curr->n  directly, as in bsGetData () */

            /* look down in the object for another value of this tag */
            if (!bsGetKeyTags (fromFieldNode->loc->useObj, _bsDown, 0))
            {
                /* no more values for this tag */
                if (aql->evalDebugLevel >= 3) {
                    indent (aql, level); aceOutf (aql->debug_out, "end (LOCAL/FOLLOW)_POS multi value\n");
                }
	      
                /* go back to the position where we were before we started iterating over the multivalues */
                if (stackExists (fromFieldNode->loc->markStack))
                {
                    fromFieldNode->loc->mark = pop (fromFieldNode->loc->markStack, BSMARK);
                    bsGoto (fromFieldNode->loc->useObj, fromFieldNode->loc->mark);
                    if (stackEmpty (fromFieldNode->loc->markStack))
                        stackDestroy (fromFieldNode->loc->markStack);

                    if (aql->evalDebugLevel >= 4) {
                        indent (aql, level); aceOutf (aql->debug_out, "loc=%lx - popped markStack useObj = %lx and jumped to mark (%s)\n",
                                                        (unsigned long)fromFieldNode->loc,
                                                        (unsigned long)fromFieldNode->loc->useObj,
                                                        bsGetMarkAsString (fromFieldNode->loc->mark));
                    }
                }
                return status ;
            }

            fromFieldNode->loc->mark = bsHandleMark (fromFieldNode->loc->useObj, fromFieldNode->loc->mark, aql->query->handle);

            if (aql->evalDebugLevel >= 4) {
                indent (aql, level); aceOutf (aql->debug_out, "loc=%lx - moved down useObj = %lx and marked (%s)\n",
                                                (unsigned long)fromFieldNode->loc, (unsigned long)fromFieldNode->loc->useObj,
                                                bsGetMarkAsString (fromFieldNode->loc->mark));
            }

            switch (fromFieldNode->loc->type)
            {
            case 'k': 
            case 'g': bsGetKeyTags (fromFieldNode->loc->useObj, _bsHere, &fromFieldNode->loc->value.k); break;
            case 'i': bsGetData (fromFieldNode->loc->useObj, _bsHere, _Int, &fromFieldNode->loc->value.i); break;
            case 'f': bsGetData (fromFieldNode->loc->useObj, _bsHere, _Float, &fromFieldNode->loc->value.f); break;
            case 't': bsGetData (fromFieldNode->loc->useObj, _bsHere, _DateType, &fromFieldNode->loc->value.t); break;
            case 's': bsGetData (fromFieldNode->loc->useObj, _bsHere, _Text, &fromFieldNode->loc->value.s); break;
            }
            break ;

        default:
	  aqlError (aql, 620, 0,
		    messprintf ("loopVariables () : should not loop on source type %d",
				fromFieldNode->loc->source)
		    ) ;
        } /* end-switch loc->source */

        /* if we get here (as a result of a "break" above, we will have found value */
        fromFieldNode->loc->isValue = KNOWNVAL ;

        /* can't "break" out of here into the code below, just return when loop-list exhausted */
    } /* end while (TRUE) */

    {
        BOOL am_sorting = tableSfwNode->resultConsumer->sorting;
        aceOutf (aql->debug_out,
                    "Finished iterating over values, finishing accepting the results, am_sorting=%s\n",
                    am_sorting ? "true" : "false");
    }
    return status ;
} /* loopVariables */

static AqlFuzzyType loopVariablesSingle (AqlNode *tableSfwNode,
                                        AqlNode *fromFieldNode,
                                        BOOL     isFuzzyRun,
                                        AQL      aql,
                                        int      level,
                                        BOOL     isFuzzyStatus,
                                        BOOL     isEmptyStatus,
                                        AqlFuzzyType status)
{
    AqlFuzzyType locStatus;

    if (isFuzzyStatus) {
        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOut (aql->debug_out, "eval fuzzy\n");
        }
        fromFieldNode->loc->isValue = FUZZYVAL; 
        status = FUZZYVAL;
    } else {
        if (isEmptyStatus) {
            if (aql->evalDebugLevel >= 4) {
                indent (aql, level); aceOut (aql->debug_out, "eval empty\n");
            }
            fromFieldNode->loc->isValue = NULLVAL; 
        } else {
            if (fromFieldNode->loc->isValue == KNOWNVAL &&
                aql->evalDebugLevel >= 3) {
                indent (aql, level); aceOut (aql->debug_out, "single value:\n");
            }
        }
    }

    locStatus = processThisLoc (tableSfwNode, fromFieldNode, isFuzzyRun, aql, level+1);
  
    if (isFuzzyRun == FALSE && locStatus == FUZZYVAL) /* impossible */
      aqlError (aql, 610, 0,
		"loopVariables () - processThisLoc produced FUZZYVAL in non-fuzzy re-run");

    if (!fromFieldNode->nxt)
    {
        /* we are at the bottom level of the recursion,
           therefore we can simply set the status; the status
           from the other levels will be combined with this on
           the way back up the stack
        */
        status = locStatus;		/* we've just reported to table - it's status is important */
    }
    else
    {
        /* we are on the way up the recursion stack
           the status so far reflects the status of all the deeper
           levels
        */
        /* incorporate the locator's result in the overall status
           NULL status can become KNOWN, but FUZZY never gets unset */
        /* according to Jean (2001) there is a problem with KNOWN getting
           changed to FUZZY on the way up, thus spuriously provoking a re-run
           and hence a disk access
        */
        if ((locStatus == FUZZYVAL) ||
            (/*status != FUZZYVAL && */locStatus == KNOWNVAL))
        {
            if (aql->evalDebugLevel >= 3) {
                indent (aql, level); aceOutf (aql->debug_out,
                                                "locStatus is %s; status has been %s; assigning status\n",
                                                locStatus == KNOWNVAL ? "known" :
                                                locStatus == FUZZYVAL ? "fuzzy" : "null",
                                                status == KNOWNVAL ? "known" :
                                                status == FUZZYVAL ? "fuzzy" : "null"
                    );
          
           
            }
            status = locStatus;
        }
    }

    if (fromFieldNode->loc->isValue == KNOWNVAL &&
        aql->evalDebugLevel >= 3) {
        indent (aql, level); aceOutf (aql->debug_out,
                                        "loopVariablesSingle returning status \"%s\"\n",
                                        status == KNOWNVAL ? "known" :
                                        status == FUZZYVAL ? "fuzzy" : "null"
            );
    }

    return status ;
}

static AqlFuzzyType processThisLoc (AqlNode *tableSfwNode,  
                                   AqlNode *fromFieldNode, 
                                   BOOL     isFuzzyRun,
                                   AQL      aql,
                                   int      level)
{
  BOOL condition = TRUE;
  AqlFuzzyType status = NULLVAL;
  
  if (aql->evalDebugLevel >= 3) {
      indent (aql, level);
      aqlLocOut (aql->debug_out, "processLoc:", fromFieldNode->loc) ;
      if (aql->evalDebugLevel >= 3) {
          aqlNodeOut (aql->debug_out, tableSfwNode, level, "sfw", TRUE);
          aqlNodeOut (aql->debug_out, fromFieldNode, level, "from", TRUE);
      
      }
  }
  
  if (fromFieldNode->loc->mustExist)
    condition = (fromFieldNode->loc->isValue != NULLVAL);
  
  if (condition)
  {
      if (fromFieldNode->right)
          /* check the WHERE clause */
      {
          if (aql->evalDebugLevel >= 4) {
              indent (aql, level); aqlNodeOut (aql->debug_out, fromFieldNode->right, level, "evaluating where", TRUE);
      
          }

          status = evalExpr (fromFieldNode->right, isFuzzyRun, aql, level+1);

          if (status == NULLVAL)
          {
              if (fromFieldNode->right->vtype == 'b')
                  condition = fromFieldNode->right->value.b;
              else
                  condition = FALSE;
          }
          else if (isFuzzyRun && status == FUZZYVAL)
              ;/* leave it */
          else
          {
              /* result is known in either fuzzyrun or exactrun */
              if (fromFieldNode->right->vtype == 'b')
                  condition = fromFieldNode->right->value.b;
              else
                  condition = TRUE;	/* KNOWNVAL counts as true where-clause */
          }
      }
      else
      {
          status = KNOWNVAL;
      }
      if (aql->evalDebugLevel >= 4) {
          indent (aql, level); aceOutf (aql->debug_out, "where status is %s and result is %s\n",
                                          status == FUZZYVAL ? "fuzzy" :
                                          status == NULLVAL ? "null" :
                                          "known",
                                          condition ? "true" : "false");
      }
      
  }
  else
    status = isFuzzyRun ? FUZZYVAL : NULLVAL; /* default to worst case */
  
  if (condition ||		/* continue eval if where-clause was true .. */
      status == FUZZYVAL)	/* .. or if it was still fuzzy */
  {
      if (aql->evalDebugLevel >= 4) {
          indent (aql, level); aceOutf (aql->debug_out, "either, condition is true, or status fuzzy: condition %s; status  %s\n",
                      condition ? "true" : "false",
                      status == FUZZYVAL ? "fuzzy" :
                      status == NULLVAL ? "null" :
                      "known");
      }
      if (fromFieldNode->nxt) {
          /* process next variable */
          if (aql->evalDebugLevel >= 4) {
              indent (aql, level); aceOut (aql->debug_out, "There are deeper variables; processing the next one\n");
          }
          status = loopVariables (tableSfwNode, fromFieldNode->nxt, isFuzzyRun, aql, level+1);
      }
      else
      {
          if (aql->evalDebugLevel >= 4) {
              indent (aql, level); aceOutf (aql->debug_out, "There are no deeper variables;%s\n",
                                              (status != FUZZYVAL) ? "Status is definite, adding to result" : "status is still fuzzy");
          }
          /* all locators evaluated - report to result table 
           * (unless we had a fuzzy where-clause) */
          if (status != FUZZYVAL)
          {
              status = sfwReport (tableSfwNode, isFuzzyRun, aql, level+1);
          }
      }
  }
  else
  {
      if (aql->evalDebugLevel >= 4) {
          indent (aql, level); aceOutf (aql->debug_out,
                                          "either, condition is false, or status not fuzzy (condition %s; status %s) -- therefore dropping this one\n",
                                          condition ? "true" : "false",
                                          status == FUZZYVAL ? "fuzzy" :
                                          status == NULLVAL ? "null" :
                                          "known");
      }
  }

  if (fromFieldNode->loc->trueObj) /* never true in fuzzy mode */
    {
        if (aql->evalDebugLevel >= 4) {
            indent (aql, level); aceOutf (aql->debug_out, "processLoc: from_locator %lx lost trueObj = %lx\n",
                                            (unsigned long)fromFieldNode->loc, (unsigned long)fromFieldNode->loc->trueObj);
        }

      bsDestroy (fromFieldNode->loc->trueObj) ;
    }

  return status ;
} /* processThisLoc */

static BOOL AqlValueEqual (char at, AqlType av,
                          char bt, AqlType bv, AQL aql)
{
    if (at != bt) return FALSE;
    switch (at)
    {
    case 's':
        if (av.s == bv.s) return TRUE;
        if (av.s == NULL || bv.s == NULL) return FALSE;
        return (strcmp (av.s, bv.s) == 0);
    default:
        return (av.i == bv.i);
    }
}

static AqlFuzzyType sfwReport (AqlNode *sfw,
                              BOOL     isFuzzyRun,
                              AQL      aql,
                              int      level)
/* evaluate selectField expressions and append a new row of values
   at the end of the results-table of this sfw-node.

   returns rowStatus - 
   NULLVAL  - all fields are NULL (row has not been appended)
   FUZZYVAL - some fields are still fuzzy (row has not been appended)
   KNOWNVAL - append row if we haven't got it already.
*/
/* XXXXXXXX rd says: inefficient - should only recalc when necessary! XXXXXXXX */
{
    AqlNode *selectFieldNode ;
    int newRowNumber, colNum;
    int selectFieldCounter ;	/* also the number of the table-column for the select-field */
    BOOL isRowFound;
    TABLE *theSfwTable = sfw->value.T ; /* where results are gathered */
    AqlFuzzyType rowStatus;

    BOOL am_sorting = sfw->resultConsumer->sorting;

    if (aql->evalDebugLevel >= 3) {
        indent (aql, level);
        aceOutf (aql->debug_out, "report to results-%s:\n",
                     am_sorting ? "table" : "stream") ;
        if (aql->evalDebugLevel >= 3) {
            aqlNodeOut (aql->debug_out,
                       sfw, level, "sfw", TRUE);
        }
    }

    /* TABLE_SFW->right is a linked list of select fields,
       they can be LOC_VARs, EXPR_OP (combining LOC_VARs) etc */

    /* init - assume all select-field values are just NULL,
     * the row can be upgraded by a KNOWN value, but gets downgraded by FUZZY values */
    rowStatus = NULLVAL;

    /* loop through selectFields of the table */
    for (selectFieldNode = sfw->right, selectFieldCounter = 0 ; selectFieldNode ; selectFieldNode = selectFieldNode->nxt, ++selectFieldCounter)
    { 
        evalExpr (selectFieldNode, isFuzzyRun, aql, level+1);

        if (aql->evalDebugLevel >= 3) {
            indent (aql, level);
            aqlNodeValueOut (aql->debug_out, selectFieldNode);
            aceOut (aql->debug_out, "\n");
        }

        if (isFuzzyRun)
        {
            if (rowStatus == NULLVAL && selectFieldNode->isValue == KNOWNVAL)
                rowStatus = KNOWNVAL; /* upgrade status */
            else if (/*rowStatus == KNOWNVAL && */selectFieldNode->isValue == FUZZYVAL)
                rowStatus = FUZZYVAL; /* taint row as fuzzy */
        }
        else
        {
            /* mark row as KNOWN if at least one field is KNOWN */
            if (selectFieldNode->isValue == KNOWNVAL)
            {
                rowStatus = KNOWNVAL;
            }
        }
    }

    /**********************************/

    if (rowStatus == NULLVAL)
    {
        /* couldn't evaluate any of the select-field expressions
         * the row is thereby completely blank */
        if (aql->evalDebugLevel >= 3) {
            indent (aql, level);
            aceOut (aql->debug_out, "(not adding row - all fields NULL)\n") ;
        }
    }
    else if (rowStatus == FUZZYVAL)
    {
        if (aql->evalDebugLevel >= 3) {
            indent (aql, level);
            aceOut (aql->debug_out, "(not adding row - some fields fuzzy)\n") ;
        }
    }
    else if (rowStatus == KNOWNVAL)
    {
        /* at least one select-field has a KNOWN value
         * so we can add the row to the results table */
      
        if (sfw->resultConsumer->sorting)
        {
            if (aql->evalDebugLevel >= 3) {
                indent (aql, level);
                aceOutf (aql->debug_out, "outputting via table for %x, using all rows\n", sfw);
            }
            newRowNumber = tableMax (theSfwTable);	/* add at the end of table */
        } else {
            if (aql->evalDebugLevel >= 3) {
                indent (aql, level);
                aceOutf (aql->debug_out, "outputting directly for %x, only ever using one row\n", sfw);
            }

            newRowNumber = 0;	/* if outputting directly, we only ever have one row */
        }            

        if (!sfw->resultConsumer->dictStringStack)
            sfw->resultConsumer->dictStringStack = stackHandleCreate (50, aql->query->handle);
        else
            stackClear (sfw->resultConsumer->dictStringStack);

        pushText (sfw->resultConsumer->dictStringStack, " ");
      
        for (selectFieldNode = sfw->right, selectFieldCounter = 0 ; selectFieldNode ;
             selectFieldNode = selectFieldNode->nxt, ++selectFieldCounter)
        {
            if (selectFieldNode->isValue == NULLVAL)
                /* set the table field to be NULL */
            {
                tableSetEmpty (theSfwTable, newRowNumber, selectFieldCounter) ;
                catText (sfw->resultConsumer->dictStringStack, "(NULL) ");
            }
            else if (selectFieldNode->isValue == KNOWNVAL)
                /* copy the value to the table-field */
            {
                if (selectFieldCounter == 0)
                {
                    /* use the leftmost field to decide whether we can clear the uniquifying hash */
                    if (!sfw->resultConsumer->sorting)
                    {
                        if (!AqlValueEqual (selectFieldNode->vtype, selectFieldNode->value,
                                           sfw->resultConsumer->vtype, sfw->resultConsumer->value, aql))
                        {
                          /*  dictClear (sfw->resultConsumer->resultTableDict); */
			  dictDestroy (sfw->resultConsumer->resultTableDict) ;
			  sfw->resultConsumer->resultTableDict = dictHandleCreate (100, aql->query->handle);
                          sfw->resultConsumer->vtype = selectFieldNode->vtype;
                          sfw->resultConsumer->value = selectFieldNode->value;
                        }
                    }
                }

                /* set the type of the table column */
                if (theSfwTable->type[selectFieldCounter] == '0')
                {
                    /* this happens the first time a value is set for this column;
                       we set the type of the values for this table column to be the type of the select-field value */
                    theSfwTable->type[selectFieldCounter] = selectFieldNode->vtype ;
                }
                else if (selectFieldNode->vtype != theSfwTable->type[selectFieldCounter])
                {
                    if (aql->evalDebugLevel >= 3) {
                        indent (aql, level);
                        aceOutf (aql->debug_out, "(type mismatch in column %d - value = NULL)",
                                     selectFieldCounter+1) ;
                    }
		  
                    /*
                      tableSetEmpty (theSfwTable, newRowNumber, selectFieldCounter) ;
                      strcat (dictString, "(NULL) ");
                      continue;
                    */
		  
                    aqlError (aql, 810, selectFieldNode->pos
			      , messprintf ("Inconsistent value types in table column %d\n"
					    "Value type was '%c' and now the field evaluates to type '%c' in row %d",
					    selectFieldCounter+1, /* number of column we're talking about (C-numbers ---> usernumbers) */
					    selectFieldNode->vtype, /* type of value we want to set in the table */
					    theSfwTable->type[selectFieldCounter], /* type of the column values so far */
					    newRowNumber)
			      ) ;
                }
	      

                /* set the table value depending on the type of selectField-value */
                switch (selectFieldNode->vtype)
                    /* type of column values */
                {
                case 'k':		/* key */
                case 'g':		/* tag */
                    tableKey (theSfwTable, newRowNumber, selectFieldCounter) = selectFieldNode->value.k ; 
                    tableParent (theSfwTable, newRowNumber, selectFieldCounter) = selectFieldNode->value.k ; 
                    catText (sfw->resultConsumer->dictStringStack, messprintf ("%d ", selectFieldNode->value.k));
                    break ;
		  
                case 'i':		/* integer */
                    tableInt (theSfwTable, newRowNumber, selectFieldCounter) = selectFieldNode->value.i ; 
                    catText (sfw->resultConsumer->dictStringStack, messprintf ("%d ", selectFieldNode->value.i));
                    break ;
		  
                case 'f':		/* float */
                    tableFloat (theSfwTable, newRowNumber, selectFieldCounter) = selectFieldNode->value.f ; 
                    catText (sfw->resultConsumer->dictStringStack, messprintf ("%g ", selectFieldNode->value.f));
                    break ;
		  
                case 't':		/* date-time */
                    tableDate (theSfwTable, newRowNumber, selectFieldCounter) = selectFieldNode->value.t ; 
                    catText (sfw->resultConsumer->dictStringStack, messprintf ("%d ", selectFieldNode->value.t));
                    break ;
		  
                case 'b':		/* boolean */
                    tableBool (theSfwTable, newRowNumber, selectFieldCounter) = selectFieldNode->value.b ; 
                    catText (sfw->resultConsumer->dictStringStack, messprintf ("%d ", selectFieldNode->value.b));
                    break ;
		  
                case 's':		/* string */
                    tableSetString (theSfwTable, newRowNumber, selectFieldCounter, selectFieldNode->value.s) ; 
                    catText (sfw->resultConsumer->dictStringStack, selectFieldNode->value.s);
                    break ;
		  
                case '0':		/* uninitialised type  */
                    /* this might happen, if this column of the table is an expression, that doesn't 
                       evaluate in this case. The type of the column is therefore not known  but will 
                       probably be set later on, if the expression evaluates TRUE with other values later on */
                    /* So we just set this table entry blank */
                    tableSetEmpty (theSfwTable, newRowNumber, selectFieldCounter);
                    catText (sfw->resultConsumer->dictStringStack, "(NULL) ");
                    break;
		  
                default:
		  aqlError (aql, 620, 0,
			    messprintf ("sfwReport () - type %c not supported for table expr",
					selectFieldNode->vtype ? selectFieldNode->vtype : '_')
			    ) ;
                }
	      
            }
	  
        } /* end for */

        /* we have now got a dictString for the complete row
         * try to find the newly added row amongst the
         * rows before that, and then behave as if it hadn't been added */
        isRowFound = FALSE;
#if 0
        if (sfw->resultConsumer->nRows == 0)
        {
            dictAdd (sfw->resultConsumer->resultTableDict, stackText (sfw->resultConsumer->dictStringStack,0), 0); /* add first item to dict-hash */
        }
        else
        {
#endif
            if (aql->evalDebugLevel >= 3) {
                indent (aql, level);
                aceOutf (aql->debug_out,
                             "Do we already have a row like \"%s\"?\n", stackText (sfw->resultConsumer->dictStringStack,0)) ;
            }

            if (!dictAdd (sfw->resultConsumer->resultTableDict,
                         stackText (sfw->resultConsumer->dictStringStack,0),
                         0)) /* check of dictString is already in hash */
                isRowFound = TRUE;
#if 0
        }
#endif
        if (isRowFound)
        {
            if (aql->evalDebugLevel >= 3) {
                indent (aql, level);
                aceOutf (aql->debug_out, "(not adding row - exists already)\n") ;
            }
            if (sfw->resultConsumer->sorting)
            {
                /* remove the row we just added, because we had it already */
                for (colNum = 0; colNum < theSfwTable->ncol; colNum++)
                {
                    if ( arrayMax (theSfwTable->col[colNum])-1 == newRowNumber )
                        /* we've just added to this column, so get rid of it */
                        arrayMax (theSfwTable->col[colNum])--;
                }
            }
        }
        else
        {
            sfw->resultConsumer->nRows++;
            if (aql->evalDebugLevel >= 3) {
                indent (aql, level);
                aceOutf (aql->debug_out, "(adding row at line %d (%d in consumer))\n", newRowNumber, sfw->resultConsumer->nRows) ;
            }

            /* now output the row */

#if 0
            tableRowOut (ACEOUT dump_out, 
                         int row, 
                         TABLE *t, 
                         char separator, 
                         char style,
                         BOOL showEmpty,
                         Array oldClass,
                         DICT **colFlagDict,
                         char **colFlagSet,
                         int *javatypeWatermark);
#endif






        }
    }

    return rowStatus ;
} /* sfwReport */

/***********************************************/

static AqlFuzzyType evalExpr (AqlNode *exprNode,
			      BOOL     isFuzzyRun,
			      AQL      aql,
                   int      level) 
     /* if evaluation was successful, it will set
      *   ->value.x to the result
      *   ->vtype to 'x'
      *   ->isValue to KNOWNVAL
      * if the result is still unknown, it will set
      *   ->vtype to the predicted value type
      *   ->isValue to FUZZYVAL
      * if evaluation doesn't succeed, it will set
      *   ->isValue to NULLVAL
      *   the value and vtype is irrelevant for NULL results
      *   NOTE: boolean types will carry a value and vtype, 
      *         even when empty - used for where-clause evaluation.
      */
{
  int tableFuncColumn = 0;	/* set it to the first column 
				 * in case we get to the code via EXPR_TABLE_FUNC, 
				 * which has no column selection */
  if (aql->evalDebugLevel >= 4) {
      aqlNodeOut (aql->debug_out,
                 exprNode, level, "evalExpr", TRUE);
      
  }

  switch (exprNode->type)
    {
      /*************************************************************/
      /* these first few are literals (i.e. constants) and have their values 
	 already set, all we need to do here is set the type of the expression.
	 This type will then be used in type checking by the expression operators */

    case nVAR:			/* a scalar variable that has hopefuly been set earlier */
      {
	/* the value and type has been copied to the locator,
	 *  see evalTable ()->case VAR_ASSIGN: */
	exprNode->vtype = exprNode->loc->type;
	exprNode->value = exprNode->loc->value;
	exprNode->isValue = exprNode->loc->isValue;
      }
      break; /* VAR */


    case nTEXT:			/* a string literal */
      {
	exprNode->vtype = 's';
	exprNode->isValue = KNOWNVAL;
	/* value.s was set by the yacc-parser */
      }
      break; /* TEXT */


    case nINT:			/* a number literal */
      {
	/* NOTE: this is a positive integer as negative integers 
	 * are nINT nodes which are attached to a nUMINUS node */
	exprNode->vtype = 'i';
	exprNode->isValue = KNOWNVAL;
	/* value.i was set by the yacc-parser */
      }
      break; /* INT */


    case nFLOAT:		/* a float number literal */
      {
	exprNode->vtype = 'f';
	exprNode->isValue = KNOWNVAL;
	/* value.f was set by the yacc-parser */
      }
      break; /* FLOAT */


    case nDATE:
      {
	exprNode->vtype = 't';	/* a date/time literal */
	exprNode->isValue = KNOWNVAL;
	/* value.t was set by the yacc-parser */
      }
      break; /* DATE */

    case nBOOL:			/* boolean literal */
      {
	exprNode->vtype = 'b';
	exprNode->isValue = KNOWNVAL;
	/* value.tb was set by the yacc-parser */
      }
      break; /* BOOL */


      
    case nKEY:			/* derived from a CLASS text-literal in object (classExpr,keyExpr) */
      {
	/* vtype/value.k/isValue already set in postAssignTypesAndClasses ():case (nOBJECT): */
      }
      break; /* KEY */



      /************************************************************/
      /* the variable gets the type and value of its locator */
    case nLOC_VAR:
      {
	exprNode->isValue = exprNode->loc->isValue ;
	
	if (exprNode->isValue == KNOWNVAL)
	  {
	    /* if the locator wasn't empty */
	    /* the expression node gets the type and value of its locator */
	    exprNode->vtype = exprNode->loc->type ; 
	    exprNode->value = exprNode->loc->value ;
	  }
	else if (exprNode->isValue == FUZZYVAL)
	  {
	    /* we'd already know the type, but don't have a value yet */
	    exprNode->vtype = exprNode->loc->type ;
	  }
	else			/* NULLVAL */
	  {
	    /* value and type are irrelevant */
	    exprNode->vtype = '0';
	    exprNode->value.i = 0;
	  }
      }
      break; /* LOC_VAR */



      /****************************************************************/
      /* the expression operators have to evaluate both sides first
	 before combining their results */
    case nEXPR_OP:
      {
	AqlFuzzyType isValueLeft;
	AqlFuzzyType isValueRight;

	/*** these are the nodes that create an expression of type EXPR_OP ***/
	/*
	  | '-' expr %prec UMINUS	    op = oUMINUS ; z->left = $<Ptr>2
	  | ExprFunc ' (' expr ')'	    op = $<Int>1 ; z->left = $<Ptr>3 ; 
	  | ExprExprFunc ' (' expr ',' expr ')'	    op = $<Int>1 ; z->left = $<Ptr>3 ; z->right = $<Ptr>5
	  | DateFunc ' (' expr ',' expr ')'  op = $<Int>1 ; z->left = $<Ptr>3 ; z->right = $<Ptr>5
	  | expr '+' expr		    op = oPLUS ; z->left = $<Ptr>1 ; z->right = $<Ptr>3
	  | expr '-' expr		    op = oMINUS ; z->left = $<Ptr>1 ; z->right = $<Ptr>3
	  | expr '*' expr		    op = oTIMES ; z->left = $<Ptr>1 ; z->right = $<Ptr>3
	  | expr '/' expr		    op = oDIVIDE ; z->left = $<Ptr>1 ; z->right = $<Ptr>3
	  | expr '%' expr		    op = oMOD ; z->left = $<Ptr>1 ; z->right = $<Ptr>3
	*/

	/* a UMINUS node is dealt with specially, because it doesn't have two expressions (hence UNARY minus :-)*/
	if (exprNode->op == oUMINUS)
	  {
	    evalExpr (exprNode->left, isFuzzyRun, aql, level+1);

	    if (exprNode->left->isValue == NULLVAL)
	      exprNode->isValue = NULLVAL;
	    else if (exprNode->left->isValue == FUZZYVAL)
	      exprNode->isValue = FUZZYVAL;
	    else
	      /* KNOWNVAL */
	      {
		/* we can only "minus" numeric expressions, i.e. a negative date or string is nonsense ! */
		if (! (exprNode->left->vtype == 'i' || exprNode->left->vtype == 'f'))
		  aqlError (aql, 898, exprNode->left->pos
			    , messprintf ("incorrect expression type (%c) for unary MINUS operator", exprNode->left->vtype)
			    );
		
		exprNode->vtype = exprNode->left->vtype;
		
		/* value is TRUE - do the minus op */
		if (exprNode->left->vtype == 'i')
		  exprNode->value.i = - (exprNode->left->value.i);
		else if (exprNode->left->vtype == 'f')
		  exprNode->value.f = - (exprNode->left->value.f);

		exprNode->isValue = KNOWNVAL;
	      }
	  }
	else if (exprNode->op == oABS)
	  {
	    evalExpr (exprNode->left, isFuzzyRun, aql, level+1);

	    if (exprNode->left->isValue == NULLVAL)
	      exprNode->isValue = NULLVAL;
	    else if (exprNode->left->isValue == FUZZYVAL)
	      exprNode->isValue = FUZZYVAL;
	    else
	      /* KNOWNVAL */
	      {
		/* we can only "abs ()" numeric expressions, absolute value of a date or string is nonsense ! */
		if (! (exprNode->left->vtype == 'i' || exprNode->left->vtype == 'f'))
		  aqlError (aql, 899, exprNode->left->pos
			    , messprintf ("incorrect expression type (%c) for abs () ExprFunc", exprNode->left->vtype)
			    );
		
		exprNode->vtype = exprNode->left->vtype;

		if (exprNode->left->vtype == 'i')
		  exprNode->value.i = (exprNode->left->value.i > 0) ? exprNode->left->value.i : - (exprNode->left->value.i);
		else if (exprNode->left->vtype == 'f')
		  exprNode->value.f = (exprNode->left->value.f > 0) ? exprNode->left->value.f : - (exprNode->left->value.f);
	    
		exprNode->isValue = KNOWNVAL;
	      }
	  }
	else  /* func (expr,expr) */
	  {
	    /****** all the remaining expressiontype have two operands - left/right ******/
	    /***** we first need to evaluate both sides of the expression, which thereby get a type and value */
	    isValueLeft  = evalExpr (exprNode->left, isFuzzyRun, aql, level+1);
	    isValueRight = evalExpr (exprNode->right, isFuzzyRun, aql, level+1);
	    
	    if (isValueLeft == NULLVAL || isValueRight == NULLVAL)
	      /* one of the expression operands is NULL, 
	       * so the result of the expression is NULL too. */
		exprNode->isValue = NULLVAL;
	    else if (isFuzzyRun &&
		     (isValueLeft == FUZZYVAL || isValueRight == FUZZYVAL))
	      /* one operand is fuzzy, so the result will be fuzzy too */
		exprNode->isValue = FUZZYVAL;
	    else
	      {
		/* both sides of the epxression have definitivly KNOWN values */
		/* the result depends on the operator now */

		if (isValueLeft == FUZZYVAL || isValueRight == FUZZYVAL)
		  aqlError (aql, 610, 0,"XXX expr eval with fuzzy vals XXX");

		exprNode->isValue = KNOWNVAL; /* default eval-result for most ops, only date-functions can fail now */

		if (exprNode->op == oPLUS ||
		    exprNode->op == oMINUS ||
		    exprNode->op == oTIMES ||
		    exprNode->op == oDIVIDE)
		  {
		    /* for now only numerictypes are allowed for the arithmetic operators */

		    if (exprNode->op == oDIVIDE && /* catch division by zero very early on */
			 ((exprNode->right->vtype == 'i' && exprNode->right->value.i == 0) ||
			  (exprNode->right->vtype == 'f' && exprNode->right->value.f == 0.0) ))
		      aqlError (aql, 895, exprNode->pos, "Arithmetic exception - division by zero");

		    if (exprNode->left->vtype == 'f' || exprNode->right->vtype == 'f')
		      exprNode->vtype = 'f'; /* arithmetic with floats results in floats */
		    else
		      exprNode->vtype = 'i';

		    /**  integer  op  integer **/
		    if (exprNode->left->vtype == 'i' && exprNode->right->vtype == 'i')
		      exprNode->value.i = evalIntArithmetic (aql, exprNode->op, exprNode->left->value.i, exprNode->right->value.i);

		    /** float  op  integer  **/
		    else if (exprNode->left->vtype == 'f' && exprNode->right->vtype == 'i')
		      exprNode->value.f = evalFloatArithmetic (aql, exprNode->op, exprNode->left->value.f, (float)exprNode->right->value.i);

		    /**   integer  op  float   **/
		    else if (exprNode->left->vtype == 'i' && exprNode->right->vtype == 'f')
		      exprNode->value.f = evalFloatArithmetic (aql, exprNode->op, (float)exprNode->left->value.i, exprNode->right->value.f);

		    /**  float  op  float   **/
		    else if (exprNode->left->vtype == 'f' && exprNode->right->vtype == 'f')
		      exprNode->value.f = evalFloatArithmetic (aql, exprNode->op, exprNode->left->value.f, exprNode->right->value.f);

		    else
		      aqlError (aql, 896, exprNode->pos
				, messprintf ("invalid expression types (%c). (%c) for arithmetic operator",
					      exprNode->left->vtype, exprNode->right->vtype)
				);
		  }
		else if (exprNode->op == oMOD)
		  {
		    /**  integer  mod  integer **/
		    /* C-implementation of x % y is bugged - modulo must be a periodic function, so
		     *    -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7
		     * modulo 3 gives
		     *     1  2  0  1  2  0  1  2  0  1  2  0  1
		     * (jean 23-feb-99) */
		    
		    if (exprNode->left->vtype == 'i' && exprNode->right->vtype == 'i')
		      {
			int tmp, x, m;
			
			x = exprNode->left->value.i, m =  exprNode->right->value.i;
			tmp = x - ((int) (x/m) * m) ;
			exprNode->value.i = (tmp + m) - ((int) (tmp+m)/m) * m ;
			exprNode->vtype = 'i';
		      }
		    else
		      aqlError (aql, 897, exprNode->pos
				, messprintf ("invalid expression types (%c). (%c) for modulo operator (integer only)",
					      exprNode->left->vtype, exprNode->right->vtype)
				);
		  }
		else
		  /* the remaining operators are date-functions */
		  {
		    BOOL isTimeDiffOK = FALSE;

		    if (exprNode->left->vtype != 't' || exprNode->right->vtype != 't')
		      aqlError (aql, 890,  exprNode->left->vtype != 't' ? exprNode->left->pos : exprNode->right->pos,
				"operand of date-function is not a DateType");

		    exprNode->vtype = 'i'; /* result is always Integer */

		    switch (exprNode->op)
		      {
		      case oYEARDIFF:
			isTimeDiffOK = timeDiffYears (exprNode->left->value.t, exprNode->right->value.t, &exprNode->value.i); break;
			
		      case oMONTHDIFF:
			isTimeDiffOK = timeDiffMonths (exprNode->left->value.t, exprNode->right->value.t, &exprNode->value.i); break;
		
		      case oWEEKDIFF:
			aqlError (aql, 798, exprNode->pos, "Sorry : Unimplemented date-function weekdiff ()"); break;

		      case oDAYDIFF:
			isTimeDiffOK = timeDiffDays (exprNode->left->value.t, exprNode->right->value.t, &exprNode->value.i); break;

		      case oHOURDIFF:
			isTimeDiffOK = timeDiffHours (exprNode->left->value.t, exprNode->right->value.t, &exprNode->value.i); break;

		      case oMINDIFF:
			isTimeDiffOK = timeDiffMins (exprNode->left->value.t, exprNode->right->value.t, &exprNode->value.i); break;
			
		      case oSECDIFF:
			isTimeDiffOK = timeDiffSecs (exprNode->left->value.t, exprNode->right->value.t, &exprNode->value.i); break;

		      default:
			aqlError (aql, 630, 0,
				  "invalid expression operator type in evalExpr () - EXPR_OP");
		      } /* end-switch exprNode->op */
		    
		    if (isTimeDiffOK == TRUE)
		      exprNode->isValue = KNOWNVAL;
		    else
		      exprNode->isValue = NULLVAL;
		  } /* end else date-functions */
	      }	/* end else both values KNOWN */
	  } /* end else (expr's with two operands) */
      }
      break; /* EXPR_OP */


    case nEXPR_TABLE_FUNC_FIELD:
      /*
	 the FIELD node specifies a specific column to work on,
	 but the number is in user-format, i.e. the leftmost column is No 1.
	 We use the variable tableFuncColumn to specify the column in the
	 subsequent code, but subtract 1 to get to C-numbers starting at 0.
      */
      tableFuncColumn = exprNode->number - 1;

      /* NOTE: fall-through with tableFuncColumn other than the pre-initialised ZERO */

    case nEXPR_TABLE_FUNC:
      {
	/* use a function on a table in column 'tableFuncColumn' (init to ZERO or set just above) */

	/* node->left is either a TABLE_SFW expression or a TABLE_VAR */
	evalTable (exprNode->left, isFuzzyRun, aql, level+1);

	if (exprNode->left->isValue == KNOWNVAL)
	  {
	    int    rowCount, numElements;
	    TABLE *theTable = exprNode->left->value.T;

	    if (!theTable)
	      aqlError (aql, 620, 0, "evalExpr - EXPR_TABLE_FUNC ->left->value.T == 0x0");
	    
	    exprNode->isValue = KNOWNVAL;

	    /* the table has a proper value, is non-empty and guaranteed to look at only only 1 column,
	       so lets find out the expression return value depending on the function specified */
	    switch (exprNode->op)
	      {
	      case oCOUNT:
		{
		  int rowCount;
		  
		  /* count-expression value returned is an integer */
		  exprNode->vtype = 'i';
		  
		  exprNode->value.i = 0;

		  /* count all elements that are non-empty */
		  for (rowCount = 0; rowCount < tabMax (theTable, tableFuncColumn); rowCount++)
		    {
		      if (tabEmpty (theTable, rowCount, tableFuncColumn))
                  continue;
		      
		      exprNode->value.i += 1;
		    }
		}
		break; /* COUNT */

	      case oMIN:
	      case oMAX:
	      case oSUM:
	      case oAVG:
		{
		  
		  if (theTable->type[tableFuncColumn] == 'i') /* perform function over integer values */
		    {
		      int intCurrent, intSum, intMin, intMax;

		      rowCount = 0;
		      while (tabEmpty (theTable, rowCount, tableFuncColumn)) /* find first nonEmpty row */
			rowCount++;
		      intMin = tabInt (theTable, rowCount, tableFuncColumn), /* init min with first value */
		      intMax = tabInt (theTable, rowCount, tableFuncColumn); /* init max with first value */
		      intSum = 0;

		      for (rowCount = 0, numElements = 0; rowCount < tabMax (theTable, tableFuncColumn); rowCount++)
			{
			  if (tabEmpty (theTable, rowCount, tableFuncColumn))
			    continue;

			  numElements++;

			  intCurrent = tabInt (theTable, rowCount, tableFuncColumn);
			  intSum += intCurrent;
			  
			  if (intCurrent < intMin) 
			    intMin = intCurrent;
			  
			  if (intCurrent > intMax)
			    intMax = intCurrent;
			}

		      exprNode->vtype = 'i'; /* default for all apart from AVG */

		      /* finished trawling through list items,
		       * now decide which characteristic of the list to return */
		      switch (exprNode->op)
			{
			case oSUM: exprNode->value.i = intSum; break;
			  
			case oMIN: exprNode->value.i = intMin; break;
			  
			case oMAX: exprNode->value.i = intMax; break;
			  
			case oAVG:
			  exprNode->vtype = 'f';
			  exprNode->value.f = (float) (intSum / numElements); break;
			default:;
			}
		    }
		  else if (theTable->type[tableFuncColumn] == 'f') /* perform function over float values */
		    {
		      float floatCurrent, floatSum, floatMin, floatMax;

		      rowCount = 0;
		      while (tabEmpty (theTable, rowCount, tableFuncColumn)) /* find first nonEmpty row */
			rowCount++;
		      floatMin = tabFloat (theTable, rowCount, tableFuncColumn), /* init min with first value */
		      floatMax = tabFloat (theTable, rowCount, tableFuncColumn), /* init max with first value */
		      floatSum = 0;

		      for (rowCount = 0, numElements = 0; rowCount < tabMax (theTable, tableFuncColumn); rowCount++)
			{
			  if (tabEmpty (theTable, rowCount, tableFuncColumn))
			    continue;

			  numElements++;

			  floatCurrent = tabFloat (theTable, rowCount, tableFuncColumn);
			  floatSum += floatCurrent;
			  
			  if (floatCurrent < floatMin) 
			    floatMin = floatCurrent;
			  
			  if (floatCurrent > floatMax)
			    floatMax = floatCurrent;
			}
		      
		      exprNode->vtype = 'f'; /* all results have  value-type 'float' */
		      
		      /* finished trawling through list items,
		       * now decide which characteristic of the list to return */
		      switch (exprNode->op)
			{
			case oSUM: exprNode->value.f = floatSum; break;
			  
			case oMIN: exprNode->value.f = floatMin; break;
			  
			case oMAX: exprNode->value.f = floatMax; break;
			  
			case oAVG: exprNode->value.f = (float) (floatSum / numElements);  break;

			default:;
			}
		    }
		  else if (theTable->type[tableFuncColumn] == 't')
		    /* perform function over date value, but only min & max allowed */
		    {
		      mytime_t dateCurrent,
			dateMin = tabDate (theTable, 0, tableFuncColumn), /* init min with first value */
			dateMax = tabDate (theTable, 0, tableFuncColumn); /* init max with first value */
		      
		      if (! (exprNode->op == oMIN || exprNode->op == oMAX))
			/* generated by e.g. 'select avg (select m->released from m in class movie)' */
			aqlError (aql, 887, exprNode->pos, /* point to the function name */
				  "Only min () or max () function can operate over DateType values!");
		      
		      for (rowCount = 0; rowCount < tabMax (theTable, tableFuncColumn); rowCount++)
			{
			  dateCurrent = tabDate (theTable, rowCount, tableFuncColumn);
			  
			  if (dateCurrent < dateMin) 
			    dateMin = dateCurrent;
			  
			  if (dateCurrent > dateMax)
			    dateMax = dateCurrent;
			}
		      
		      /* result is also a DateType */
		      exprNode->vtype = 't';
		      
		      if (exprNode->op == oMIN)
			exprNode->value.t = dateMin;
		      else if (exprNode->op == oMAX)
			exprNode->value.t = dateMax;
		    }
		  else if (theTable->type[tableFuncColumn] == 's')
		    {
		      char *stringCurrent,
			*stringMin = tabString (theTable, 0, tableFuncColumn), /* init min with first value */
			*stringMax = tabString (theTable, 0, tableFuncColumn); /* init max with first value */
		      
		      if (! (exprNode->op == oMIN || exprNode->op == oMAX))
			/* generated by e.g. 'select avg (select p->Full_name from p in class person)' */
			aqlError (aql, 888, exprNode->pos, /* point to the function name */
				  "Only min () or max () function can operate over Text values!");
		      
		      for (rowCount = 0; rowCount < tabMax (theTable, tableFuncColumn); rowCount++)
			{
			  stringCurrent = tabString (theTable, rowCount, tableFuncColumn);
			  
			  if (lexstrcmp (stringCurrent, stringMin) < 0)
			    stringMin = stringCurrent;
			  
			  if (lexstrcmp (stringCurrent, stringMin) > 0)
			    stringMax = stringCurrent;
			}
		      
		      /* result is also of type Text */
		      exprNode->vtype = 's';
		      
		      if (exprNode->op == oMIN)
			exprNode->value.s = stringMin;
		      else if (exprNode->op == oMAX)
			exprNode->value.s = stringMax;
		    }
		  else
		    aqlError (aql, 889, exprNode->left->pos,
			      "Table selection returns wrong type of value.\n"
			      "Functions min, max only work over columns of Int, Float, Text or DateType values,\n"
			      "Functions sum, avg work over numeric values only.");
		}
		break;		/* end MIN/MAX/SUM/AVG */


		/* NOTE, that we wouldn't need to evaluate the 
		   whole table for FIRST and LAST, needs optimising */
	      case oFIRST:		
	      case oLAST:
		{
		  int first = -1, last = 0, rowCount = 0;
		  
		  /* count all elements that are non-empty */
		  for (rowCount = 0; rowCount < tabMax (theTable, tableFuncColumn); rowCount++)
		    {
		      if (tabEmpty (theTable, rowCount, tableFuncColumn))
			continue;

		      if (first == -1) first = rowCount;
		      last = rowCount;
		    }

		  if (exprNode->op == oFIRST)
		    rowCount = first;
		  else
		    rowCount = last;

		  if (first > -1) /* nonEmpty fields were found */
		    {
		      exprNode->vtype = theTable->type[tableFuncColumn];
		      
		      switch (exprNode->vtype)
			{
			case 'k':
			case 'g':
			  exprNode->value.k = tabKey (theTable, rowCount, tableFuncColumn); break;
			  
			case 'i':
			  exprNode->value.i = tabInt (theTable, rowCount, tableFuncColumn); break;
			  
			case 'f':
			  exprNode->value.f = tabFloat (theTable, rowCount, tableFuncColumn); break;
			  
			case 't':
			  exprNode->value.t = tabDate (theTable, rowCount, tableFuncColumn); break;
			  
			case 'b':
			  exprNode->value.b = tabBool (theTable, rowCount, tableFuncColumn); break;
			  
			case 's':
			  exprNode->value.s = tabString (theTable, rowCount, tableFuncColumn); break;
			}
		    }
		  else
		    exprNode->isValue = NULLVAL;
		}
		break;


	      default:
		aqlError (aql, 620, 0, "evalExpr - invalid tableFuncNode->operator");
	      }	/* end-switch exprNode->op */
	  }
	else
	  /* evalTable failed or still fuzzy */
	  {
	    exprNode->isValue = exprNode->left->isValue;
	  }
      }
      break; /* EXPR_TABLE_FUNC */


      /* whether the locator is KNOWN (the result is never NULL, just TRUE or FALSE) */
    case nBOOL_EXISTS:
    {
        exprNode->vtype = 'b';

        if (aql->evalDebugLevel >= 4)
	  {
            indent (aql, level);
            aceOutf (aql->debug_out, "testing exists (exprNode->left=%lx && exprNode->left->loc=%lx)\n",
                         exprNode->left, exprNode->left->loc);
	  }
	
        if (exprNode->left && exprNode->left->loc)
	  {
            AqlFuzzyType child_value = exprNode->left->loc->isValue;
            /* KNOWNVAL means it exists, NULLVAL means it doesn't exist */
            exprNode->isValue = (( (child_value == KNOWNVAL) || (child_value == NULLVAL)) ?
                                 KNOWNVAL :
                                 FUZZYVAL); /* Note: may be still fuzzy */
            exprNode->value.b = (child_value == KNOWNVAL);
            if (aql->evalDebugLevel >= 4) {
	      indent (aql, level);
	      aceOutf (aql->debug_out, "child_value = %s\n",
			   child_value == KNOWNVAL ? "known":
			   child_value == NULLVAL ? "null" :
			   child_value == FUZZYVAL ? "fuzzy" : "other");
            }
	  }
	else
	  aqlError (aql, 620, 0, "BOOL_EXISTS node should have ->left and left->loc");
        }
      break; /* BOOL_EXISTS */



    case nBOOL_NOT:		/* give TRUE if the boolean evaluates FALSE */
      {
	exprNode->vtype = 'b';
	exprNode->isValue = evalExpr (exprNode->left, isFuzzyRun, aql, level+1);

	if (exprNode->isValue == NULLVAL)
	  exprNode->value.b = TRUE; /* assume 'NOT null' is TRUE */
	else
	  exprNode->value.b = ! (valueAsBool (aql, exprNode->left));

      }
      break; /* BOOL_NOT */



    case nBOOL_COMPARISON:
      {
	/* evaluate the expression on either side of the operator first */
	AqlFuzzyType isValueLeft  = evalExpr (exprNode->left, isFuzzyRun, aql, level+1);
	AqlFuzzyType isValueRight = evalExpr (exprNode->right, isFuzzyRun, aql, level+1);

	exprNode->vtype = 'b';

	if (isValueLeft == NULLVAL || isValueRight == NULLVAL ||
	    isValueLeft == FUZZYVAL || isValueRight == FUZZYVAL)
	  /* a comparison involving NULL or FUZZY values */
	  {
	    if (isValueLeft == NULLVAL || isValueRight == NULLVAL)
	      exprNode->isValue = NULLVAL; /* comparison with null values is NULL */
	    else
	      exprNode->isValue = FUZZYVAL; /* a comparison involving a fuzzy value has a FUZZY result */
	    
	    /* Although this node is NULL or FUZZY in an expression context, we may still
	     *   need its proper boolean value (in a where-clause context) */

	    if (exprNode->left->vtype == 'b' || exprNode->right->vtype == 'b')
	      /* if the comarison involves booleans we can make a better decision */
	      {
		if (! (exprNode->op == oEQ || exprNode->op == oNE))
		  aqlError (aql, 885, 0, "Comparison involving boolean values only allows = and != operators");
		
		if (exprNode->left->vtype == 'b' && exprNode->right->vtype == 'b')
		  /* both sides are boolean */
		  exprNode->value.b = boolComparison (aql, exprNode->op, exprNode->left->value.b, exprNode->right->value.b);
		else if (exprNode->left->vtype == 'b')
		  /* right is non-boolean */
		  exprNode->value.b = boolComparison (aql, exprNode->op,
						      exprNode->left->value.b,
						      isValueRight == KNOWNVAL ? TRUE : FALSE); /* non-booleans are TRUE if non-null */
		else if (exprNode->right->vtype == 'b')
		  /* left is non-boolean */
		  exprNode->value.b = boolComparison (aql, exprNode->op,
						      isValueLeft == KNOWNVAL ? TRUE : FALSE, /* non-booleans are TRUE if non-null */
						      exprNode->right->value.b);  
	      }
	    else
	      {
		/* The comparison involves non-booleans */
		if (exprNode->isValue == NULLVAL)
		  exprNode->value.b = FALSE;
		else
		  exprNode->value.b = TRUE; /* comparing fuzzy values is assumed, 
					     * although the fuzzy value may turn out different next time round */
	      }
	  }
	else
	  /* no NULL or FUZZY values from here on */
	  {
	    char lType = exprNode->left->vtype, rType = exprNode->right->vtype ;
	    exprNode->isValue = KNOWNVAL;


	    /*** comparing strings and keys/tags (by their string representation) ****/
	    
	    if ((lType == 'k' || lType == 'g') && rType == 's')	/* comparing key to string */
	      {
		if (lType == 'k' && pickCaseSensitive (exprNode->left->value.k))
		  exprNode->value.b = stringComparison (aql, exprNode->op, name (exprNode->left->value.k), exprNode->right->value.s, TRUE);
		else
		  exprNode->value.b = stringComparison (aql, exprNode->op, name (exprNode->left->value.k), exprNode->right->value.s, FALSE);
	      }

	    else if (lType == 's' && (rType == 'k' || rType == 'g')) /* comparing string to key */
	      {
		char *lString, *rString;

		if (exprNode->op == oLIKE)
		  {
		    /* we're comparing a string template to an object or tag name, so we flip the argument, 
		     * because the template is expected to be the second string */
		    lString = name (exprNode->right->value.k);
		    rString = exprNode->left->value.s;
		  }
		else
		  {
		    lString = exprNode->left->value.s;
		    rString = name (exprNode->right->value.k);
		  }
		
		if (rType == 'k'  && pickCaseSensitive (exprNode->right->value.k))
		  exprNode->value.b = stringComparison (aql, exprNode->op, lString, rString, TRUE);
		else
		  exprNode->value.b = stringComparison (aql, exprNode->op, lString, rString, FALSE);
	      }
	    
	    else if (lType == 's' && rType == 's') /* comparing two strings */
	      exprNode->value.b = stringComparison (aql, exprNode->op, exprNode->left->value.s, exprNode->right->value.s, FALSE);

	    /* past this point, no more string comparison is undertaken, 
	     * so the LIKE operator becomes invalid */
	    else if (exprNode->op == oLIKE)
	      aqlError (aql, 883, exprNode->pos, "The like operator can only be used on Text types");
	    
	    /* comparing two keys */
	    else if ((lType == 'k' || lType == 'g') && (rType == 'k' || rType == 'g'))
	      {
		if (! (exprNode->op == oEQ || exprNode->op == oNE))
		  /* only == and != allowed */
		  aqlError (aql, 886, exprNode->pos, "Comparison of two object values only allows = and != operators");
		
		exprNode->value.b = keyComparison (aql, exprNode->op, exprNode->left->value.k, exprNode->right->value.k);
	      }
	    
	    /*** comparing two numbers ****/
	    
	    else if (lType == 'i' && rType == 'i') /* both integers */
	      exprNode->value.b = integerComparison (aql, exprNode->op, exprNode->left->value.i, exprNode->right->value.i);

	    else if (lType == 'f' && rType == 'f') /* both floats */
	      exprNode->value.b = floatComparison (aql, exprNode->op, exprNode->left->value.f, exprNode->right->value.f);
	    
	    else if (lType == 'i' && rType == 'f') /* typecast the integer to a float for comparison */
	      exprNode->value.b = floatComparison (aql, exprNode->op, (float)exprNode->left->value.i, exprNode->right->value.f);
	    
	    else if (lType == 'f' && rType == 'i') /* typecast the integer to a float for comparison */
	      exprNode->value.b = floatComparison (aql, exprNode->op, exprNode->left->value.f, (float)exprNode->right->value.i);

	    /*** comparing two dates ***/

	    else if (lType == 't' && rType == 't')
	      {
		switch (exprNode->op)
		  {
		  case oEQ:		/* = */
		    exprNode->value.b = timeComparison (0, exprNode->left->value.t, exprNode->right->value.t);
		    break;
			
		  case oLT:		/* < */
		    exprNode->value.b = timeComparison (-1, exprNode->left->value.t, exprNode->right->value.t);
		    break;
		    
		  case oGT:		/* > */
		    exprNode->value.b = timeComparison (1, exprNode->left->value.t, exprNode->right->value.t);
		    break;
		    
		  case oNE:		/* != */
		    exprNode->value.b = (!timeComparison (0, exprNode->left->value.t, exprNode->right->value.t));
		    break;
		    
		  case oLE:		/* <= */
		    exprNode->value.b = ((timeComparison (0,  exprNode->left->value.t, exprNode->right->value.t)) || 
					 (timeComparison (-1, exprNode->left->value.t, exprNode->right->value.t)));
		    break;
		    
		  case oGE:		/* >= */
		    exprNode->value.b = ((timeComparison (0, exprNode->left->value.t, exprNode->right->value.t)) || 
					 (timeComparison (1, exprNode->left->value.t, exprNode->right->value.t)));
		    break;
		    
		  default:
		    aqlError (aql, 620, 0, "evalBool () - invalid date-comparison operator");
		  }
	      }

	    /*** comparison involving booleans ***/
	    else if (lType == 'b' || rType == 'b')
	      {
		if (! (exprNode->op == oEQ || exprNode->op == oNE))
		  aqlError (aql, 885, 0, "Comparison involving boolean values only allows = and != operators");
		
		if (lType == 'b' && rType == 'b')
				/* both booleans */
		  exprNode->value.b = boolComparison (aql, exprNode->op, exprNode->left->value.b, exprNode->right->value.b);
		else if (lType == 'b')
				/* left is boolean */
		  exprNode->value.b = boolComparison (aql, exprNode->op, exprNode->left->value.b, valueAsBool (aql, exprNode->right));  
		else 
				/* right is boolean */
		  exprNode->value.b = boolComparison (aql, exprNode->op, valueAsBool (aql, exprNode->left), exprNode->right->value.b);  
	      }
	    else
	      /* all legal combinations of value types for boolean comparisons have been dealt with now */
	      aqlError (aql, 884, exprNode->pos, "illegal value types for comparison operator");
	  }
      }
      break; /* BOOL_COMPARISON */



    case nBOOL_OP:		/* and, or & xor */
      {
	AqlFuzzyType isValueLeft;
	AqlFuzzyType isValueRight ;
	
	exprNode->vtype = 'b';

	isValueLeft = evalExpr (exprNode->left, isFuzzyRun, aql, level+1);

	if (isValueLeft == NULLVAL)
	  {
	    exprNode->isValue = NULLVAL;
	    exprNode->value.b = FALSE;
	  }
	else
	  /* the left value is known or fuzzy (we may not have to look at the right value) */
	  {
	    switch (exprNode->op)
	      {
	      case oAND:
		/* both have to be TRUE */
		/* evaluate lazily - 
		 * only need to look at the right expr if the left was TRUE */
		
		if (isValueLeft == KNOWNVAL && valueAsBool (aql, exprNode->left) == FALSE)
		  {
		    exprNode->isValue = KNOWNVAL;
		    exprNode->value.b = FALSE;	/* left is known to be false, result is FALSE */
		  }
		else
		  /* left value is known to be true or it's fuzzy, so we need to look at the right node */
		  {
		    isValueRight = evalExpr (exprNode->right, isFuzzyRun, aql, level+1);
		    
		    if (isValueRight == NULLVAL)
		      {
			exprNode->isValue = NULLVAL;
			exprNode->value.b = FALSE;
		      }
		    else if (isValueRight == KNOWNVAL && valueAsBool (aql, exprNode->right) == FALSE)
		      {
			exprNode->isValue = KNOWNVAL;
			exprNode->value.b = FALSE;	/* right is known to be false, result is FALSE */
		      }
		    else if ((isValueLeft == KNOWNVAL && valueAsBool (aql, exprNode->left) == TRUE) &&
			      (isValueRight == KNOWNVAL && valueAsBool (aql, exprNode->right) == TRUE) )
		      {
			exprNode->isValue = KNOWNVAL;
			exprNode->value.b = TRUE; /* both are known to be TRUE, logical and condition fullfilled */
		      }
		    else
		      /* neither side is NULL or KNOWN, so both are fuzzy
		       * and we'll say that this case is true in a boolean context */
		      {
			exprNode->isValue = FUZZYVAL;
			exprNode->value.b = (valueAsBool (aql, exprNode->left) && valueAsBool (aql, exprNode->right)); /* best guess */
		      }
		  }
		break;
		
	      case oOR:
		/* at least one is TRUE */
		if (isValueLeft == KNOWNVAL && valueAsBool (aql, exprNode->left) == TRUE)
		  {
		    exprNode->isValue = KNOWNVAL;
		    exprNode->value.b = TRUE; /* left is known to be true - that's enough */
		  }
		else
		  /* left value is known to be false or fuzzy, so we need to look at the right node */
		  {
 		    isValueRight = evalExpr (exprNode->right, isFuzzyRun, aql, level+1);
		    
		    if (isValueRight == NULLVAL)
		      {
			exprNode->isValue = NULLVAL;
			exprNode->value.b = FALSE; /* assume 'FALSE or NULL' is FALSE */
		      }
		    else if (isValueRight == KNOWNVAL && valueAsBool (aql, exprNode->right) == TRUE)
		      /* right is KNOWN to be true, that's fine to satisfy OR condition*/
		      {
			exprNode->isValue = KNOWNVAL;
			exprNode->value.b = TRUE;
		      }
		    else if (isValueLeft == KNOWNVAL && isValueRight == KNOWNVAL)
		      /* either side is either known to be false */
		      {
			exprNode->isValue = KNOWNVAL;
			exprNode->value.b = FALSE;
		      }
		    else
		      /* either side must be fuzzy */
		      {
			exprNode->isValue = FUZZYVAL;
			exprNode->value.b = (valueAsBool (aql, exprNode->left) || valueAsBool (aql, exprNode->right)); /* best guess */
		      }
		  }
		break;
		
		
	      case oXOR:
		/* only one TRUE, but not both */
		{
		  BOOL left, right;

		  /* we always need to look at both branches */
		  isValueRight = evalExpr (exprNode->right, isFuzzyRun, aql, level+1);

		  if (isValueRight == NULLVAL)
		    {
		      exprNode->isValue = NULLVAL;
		      exprNode->value.b = FALSE;
		    }
		  else
		    {
		      left = valueAsBool (aql, exprNode->left);
		      right = valueAsBool (aql, exprNode->right);

		      if (isValueLeft == KNOWNVAL && isValueRight == KNOWNVAL)
			exprNode->isValue = KNOWNVAL;
		      else
			/* one or both of the values if fuzzy, but none is null */
			exprNode->isValue = FUZZYVAL;

		      /* result is exclusive-OR */
		      exprNode->value.b = (left && !right) || (!left && right);
		    }
		}
		break;

	      default:
		aqlError (aql, 620, 0, "evalBool () - invalid operator in BOOL_OP-node");
	      } /* end-switch exprNode-op */
	  } /* end else non-NULL */
      }
      break; /* BOOL_OP */
      
    default:
      aqlError (aql, 620, 0, "evalExpr - invalid exprNode->type");
      /* end-switch exprNode->type */
    }

  if (!isFuzzyRun && exprNode->isValue == FUZZYVAL)
    aqlError (aql, 610, 0,"evalExpr () - produced FUZZYVAL in the non-fuzzy re-run");

  if (aql->evalDebugLevel >= 4)
    {
      indent (aql, level);
      aceOutf (aql->debug_out, "evalExpr result, isValue = %s:",
                   exprNode->isValue == NULLVAL ? "null" :
		   exprNode->isValue == KNOWNVAL ? "known" :
		   exprNode->isValue == FUZZYVAL ? "fuzzy" :
		   "other");
      aqlNodeValueOut (aql->debug_out, exprNode);
      aceOutf (aql->debug_out, "\n");
    }
  
  return exprNode->isValue;
} /* evalExpr */



/* returns TRUE is the operator yields TRUE for the two strings,
   otherwise false,

   The extra isCaseSensitive parameter is used when testing for equality
   between strings and also for the template matching with the LIKE
   operator. Currently it can only be TRUE when comparing object names
   from a class to be configured case-sensitive in options.wrm.

   For <greater-than> and <less-than> we use lexstrcmp,
   so numbers in strings are sorted properly, e.g. unc-7 before unc-13
*/
static BOOL stringComparison (AQL aql, AqlOpType boolOp, char *leftString, char *rightString, BOOL isCaseSensitive)
{
  int (*my_strcmp) (const char *s1, const char *s2);

  if (isCaseSensitive)
    my_strcmp = strcmp;
  else
    my_strcmp = strcasecmp;


  switch (boolOp)
    {
    case oEQ:		/* = */
      return ((my_strcmp) (leftString, rightString) == 0);

    case oLT:		/* < */
      return (lexstrcmp (leftString, rightString) < 0);

    case oGT:		/* > */
      return (lexstrcmp (leftString, rightString) > 0);

    case oNE:		/* != */
      return (! ((my_strcmp) (leftString, rightString) == 0));

    case oLE:		/* <= */
      if (( (my_strcmp) (leftString, rightString) == 0) ||
	  (lexstrcmp (leftString, rightString) < 0))
	return TRUE;
      else 
	return FALSE;

    case oGE:		/* >= */
      if (( (my_strcmp) (leftString, rightString) == 0) ||
	  (lexstrcmp (leftString, rightString) > 0))
	return TRUE;
      else 
	return FALSE;

    case oLIKE:
      if (isCaseSensitive)
	return (pickMatchCaseSensitive (leftString, rightString, TRUE) > 0);
      else
	return (queryRegExpMatch (leftString, rightString, FALSE) > 0);

    default:
      aqlError (aql, 610, 0, "stringComparison () - invalid comparison operator");
    }

  return FALSE;			/* can't get here, but keeps compiler happy */
} /* stringComparison */


/* returns TRUE if the operator yields TRUE for the two integers, 
   otherwise FALSE 
*/
static BOOL integerComparison (AQL aql, AqlOpType boolOp, int leftInt, int rightInt)
{
  switch (boolOp)
    {
    case oEQ:
      return (leftInt == rightInt);

    case oLT:
      return (leftInt < rightInt);

    case oGT:
      return (leftInt > rightInt);

    case oNE:
      return (leftInt != rightInt);

    case oLE:
      return (leftInt <= rightInt);

    case oGE:
      return (leftInt >= rightInt);

    default:
      aqlError (aql, 610, 0, "integerComparison () - invalid comparison operator");
    }

  return FALSE;			/* can't get here, but keeps compiler happy */
} /* integerComparison */


/* returns TRUE if the operator yields TRUE for the two floats, 
   otherwise FALSE 
*/
static BOOL floatComparison (AQL aql, AqlOpType boolOp, float leftFloat, float rightFloat)
{
  switch (boolOp)
    {
    case oEQ:
      return (leftFloat == rightFloat);

    case oLT:
      return (leftFloat < rightFloat);

    case oGT:
      return (leftFloat > rightFloat);

    case oNE:
      return (leftFloat != rightFloat);

    case oLE:
      return (leftFloat <= rightFloat);

    case oGE:
      return (leftFloat >= rightFloat);

    default:
      aqlError (aql, 610, 0, "floatComparison () - invalid comparison operator");
    }

  return FALSE;			/* can't get here, but keeps compiler happy */
} /* floatComparison */


/* returns TRUE if both bools are the same, otherwise FALSE */
static BOOL boolComparison (AQL aql, AqlOpType boolOp, BOOL leftBool, BOOL rightBool)
{
  switch (boolOp)
    {
    case oEQ:
      return (leftBool == rightBool);

    case oNE:
      return (leftBool != rightBool);

    default:
      aqlError (aql, 610, 0, "boolComparison () - invalid comparison op - "
		"should have been caught by aqlError 885");
    }

  return FALSE;			/* can't get here, but keeps compiler happy */
} /* boolComparison */



static BOOL keyComparison (AQL aql, AqlOpType boolOp, KEY leftKey, KEY rightKey)
{
  switch (boolOp)
    {
    case oEQ:
      return (lexAliasOf (leftKey) == lexAliasOf (rightKey));

    case oNE:
      return (lexAliasOf (leftKey) != lexAliasOf (rightKey));

    default:
      aqlError (aql, 610, 0, "keyComparison () - invalid comparison op - "
		"should have been caught by aqlError 886");
    }

  return FALSE;			/* can't get here, but keeps compiler happy */
} /* keyComparison */


/* returns the value of a non-NULL node in boolean context */
static BOOL valueAsBool (AQL aql, AqlNode *inNode)
{
  BOOL value;

  if (inNode->isValue == NULLVAL)
    aqlError (aql, 610, 0, "valueAsBool () - inNode is NULL");

  if (inNode->isValue == KNOWNVAL && inNode->vtype == 'b')
    value = inNode->value.b;	/* known boolean */
  else if (inNode->isValue == KNOWNVAL && inNode->vtype == 'i')
    value = inNode->value.i ? TRUE : FALSE; /* known integer */
  else
    value = TRUE;		/* value is non-NULL, which is enough
				 * for it to count as TRUE */

  return value;
} /* valueAsBool */


/***********************************************************************/


/* returns the result of the expression depending on the operator of type wh/aql_.h:AqlOp */
/* this function is only used when both sides of the expression are integers */
static int evalIntArithmetic (AQL aql, AqlOpType exprOp, int leftInt, int rightInt)
{
  switch (exprOp)
    {
    case oPLUS:
      return leftInt + rightInt;
    case oMINUS:
      return leftInt - rightInt;
    case oTIMES:
      return leftInt * rightInt;
    case oDIVIDE:
      return leftInt / rightInt;
    default:
      aqlError (aql, 610, 0, "evalIntArithmetic () - invalid exprOp");
    }
  return FALSE;			/* can't get here, but keeps compiler happy */
} /* evalIntArithmetic */



/* returns the result of the expression depending on the operator */
/* this function is also used, when one of the value is an integer,
   it is typecast to a float before it is passed to this function. */
static float evalFloatArithmetic (AQL aql, AqlOpType exprOp, float leftFloat, float rightFloat)
{
  switch (exprOp)
    {
    case oPLUS:
      return leftFloat + rightFloat;
    case oMINUS:
      return leftFloat - rightFloat;
    case oTIMES:
      return leftFloat * rightFloat;
    case oDIVIDE:
      return leftFloat / rightFloat;
    default:
      aqlError (aql, 620, 0, "evalFloatArithmetic () - invalid exprOp");
    }

  return FALSE;			/* can't get here, but keeps compiler happy */
} /* evalFloatArithmetic */

static void indent (AQL aql, int level)
{
  int i;

  for (i = 0; i < level; i++)
    {
      aceOut (aql->debug_out, "|   ");
    }

  return ;
}

/**************************** eof **************************************/
